# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import logging
import traceback
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Counter, Tuple, List, Dict
import yaml
import numpy as np
import time
from skimage.filters.edges import sobel
from pwem.emlib.image.image_readers import ImageReadersRegistry
from pyworkflow import BETA
from pyworkflow.protocol import STEPS_PARALLEL, BooleanParam, ProtStreamingBase, LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, GE, LE, EnumParam
from pyworkflow.object import Set, Pointer, String
from pyworkflow.utils import cyanStr, Message, redStr, yellowStr
from pwem.protocols import EMProtocol
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

logger = logging.getLogger(__name__)

# Form variables
IN_TS_SET = 'inTsSet'
MIN_TILT = 'minTilt'
MAX_TILT = 'maxTilt'
MAX_SX = 'maxShiftX'
MAX_SY = 'maxShiftY'
MIN_DOSE = 'minDose'
MAX_DOSE = 'maxDose'
QUALITY_FILTER = 'qualityFilterFactor'
DARK_SENSITIVITY = 'darkSensitivity'
DARK_LOCAL_RATION_TH = 'darkLocalRatioTh'
DARK_WINDOW = 'darkWindow'
DARK_COS_FLOOR = 'darkCosFloor'
DARK_CROP_FRACTION = 'darkCropFraction'
EXTREME_CONTRAST_TH = 'extremeContrastTh'
MIN_ENTROPY = 'minEntropy'
MIN_EDGE_ENERGY = 'minEdgeEnergy'
MIN_VIEWS = 'minViews'
DO_RESTACK = 'doReStack'


# Tilt-series annotation keys
@dataclass
class TiMetadata:
    by_max_shift: bool = False
    by_tilt_angle: bool = False
    by_dose: bool = False
    by_dark: bool = False
    by_extreme_contrast_split: bool = False
    by_low_info_or_flat_image: bool = False
    by_edge_energy: bool = False


class TiLabel(Enum):
    BY_MAX_SHIFT = 'by_max_shift'
    BY_TILT_ANGLE = 'by_tilt_angle'
    BY_DOSE = 'by_dose'
    BY_DARK = 'by_dark'
    BY_EXTREME_CONTRAST_SPLIT = 'by_extreme_contrast_split'
    BY_LOW_INFO = 'by_low_info_or_flat_image'
    BY_EDGE_ENERGY = 'by_edge_energy'  # Very unfocused or lack of sample


class QualityFilterModes(Enum):
    disabled = 0
    conservative = 1
    balanced = 2
    aggressive = 3
    custom = 4


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries()
    failedTiltSeries = SetOfTiltSeries()


class ProtExclViewFilter(EMProtocol, ProtStreamingBase):
    """
    This protocol allows to filter a set of aligned tilt series according to a set of parameter as they are:
     * Maximum allowed shift after the tilt series alignment
     * Range of accumulated dose
     * Range of tilt angle
     * Minimum number of views
    """
    _label = 'exclude views filter'
    _devStatus = BETA
    _possibleOutputs = outputObjects
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.itemTsIdReadList = []
        self.failedItems = []
        self.sRate = -1
        self.removedTsIds = String('')
        self.tiLabelDict = dict()
        self.darkSensitivityVal = None
        self.darkLocalRatioThVal = None
        self.darkWindowVal = None
        self.darkCosFloorVal = None
        self.darkCropFractionVal = None
        self.extremeContrastThVal = None
        self.minEntropyVal = None
        self.minEdgeEnergyVal = None

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)

        form.addParam(IN_TS_SET,
                      PointerParam,
                      important=True,
                      label="Tilt series",
                      pointerClass='SetOfTiltSeries',
                      help='Select several sets of tilt-series where to evaluate the consensus in their alignment. '
                           'Output set will bring the information from the first selected set.')

        group = form.addGroup('Filter criteria')
        lineShift = group.addLine('Filter by max shift (%)',
                                  help='This is the minimum/maximum shift allowed along the X or Y direction. '
                                       'A value of 0.1 means that a 10% of the dimensions of the tilt image '
                                       'is allowed.')
        lineShift.addParam(MAX_SX, FloatParam, default=0.0, validators=[GE(0), LE(1)], label="X")
        lineShift.addParam(MAX_SY, FloatParam, default=0.0, validators=[GE(0), LE(1)], label="Y   ")

        lineTilt = group.addLine('Filter by tilt angle (deg)',
                                 help='This is the minimum/maximum tilt angle allowed.'
                                      'Only the tilt images with tilt angles in the range min<tilt<max are'
                                      'allowed.')
        lineTilt.addParam(MIN_TILT, FloatParam, default=-30.0, label="Min")
        lineTilt.addParam(MAX_TILT, FloatParam, default=30.0, label="Max")

        lineDose = group.addLine('Filter by dose (e/A^2)',
                                 help='This is the minimum/maximum dose per tilt image.'
                                      'Only the tilt images with accumulated dose in this range will be kept.')
        lineDose.addParam(MIN_DOSE, FloatParam, default=0.0, label="Min")
        lineDose.addParam(MAX_DOSE, FloatParam, default=70.0, label="Max")

        lineDark = group.addLine('Filter by image quality',
                                 help='Disable low-quality tilt images (dark/blank frames, partial obstructions, '
                                      'flat/low-texture or structureless views) using a strictness preset '
                                      '(Conservative / Balanced / Aggressive / Custom).\n\n'
                                      '*Conservative:*\n '
                                      'Removes only clear, catastrophic failures (e.g., strong blackouts '
                                      'or obvious obstructions). Recommended when you want to keep as many views as '
                                      'possible, especially at high tilts, and rely on downstream steps to handle '
                                      'minor quality variations.\n\n'
                                      '*Balanced:*\n'
                                      'A practical default for most datasets. Targets typical punctual failures (dark '
                                      'frames, heavy partial obstructions) while ignoring normal high-tilt darkening '
                                      'and modest contrast changes.\n\n'
                                      '*Aggressive:*\n'
                                      'Stricter filtering intended to remove borderline low-quality views (milder '
                                      'dark dips, stronger local contrast anomalies, lower texture). Recommended when you '
                                      'see frequent intermittent acquisition problems or when you prefer a cleaner '
                                      'but smaller set of views\n\n'
                                      '*Custom:*\n'
                                      'Uses the preset as a starting point, then applies user-defined thresholds. '
                                      'Choose this when you want precise control over the sensitivity of the dark '
                                      'detector and the thresholds for extreme contrast, low information (entropy), '
                                      'and edge-energy.')

        lineDark.addParam(QUALITY_FILTER, EnumParam,
                          display=EnumParam.DISPLAY_COMBO,
                          choices=[item.name for item in QualityFilterModes],
                          default=QualityFilterModes.balanced.value,
                          label='Quality preset')

        customQualityCond = f'{QUALITY_FILTER} == {QualityFilterModes.custom.value}'
        group.addParam(DARK_SENSITIVITY, FloatParam,
                       label='Dark sensitivity',
                       default=2.0,
                       validators=[GE(0.5), LE(6)],
                       condition=customQualityCond,
                       help='Controls how extreme a brightness drop must be to classify an image as a dark outlier.'
                            '\n\t- Lower values = more aggressive (more images flagged).'
                            '\n\t- Higher values = more conservative (only strong blackouts).'
                            '\n\t- Typical range: 1.5 – 3.0 (valid range: 0.5 – 6.0).')
        group.addParam(DARK_LOCAL_RATION_TH, FloatParam,
                       label='Dark local ratio threshold',
                       default=0.65,
                       validators=[GE(0.30), LE(0.90)],
                       expertLevel=LEVEL_ADVANCED,
                       help='Requires the image to be sufficiently darker than its neighboring tilts.'
                            '\n\t- Higher values = more aggressive (flags smaller local drops).'
                            '\n\t- Lower values = more conservative (only strong, punctual drops).'
                            '\n\t- Typical range: 0.6 – 0.7 (valid range: 0.3 – 0.9).')
        group.addParam(DARK_WINDOW, IntParam,
                       label='Dark neighborhood tilts',
                       default=2,
                       validators=[GE(1), LE(5)],
                       expertLevel=LEVEL_ADVANCED,
                       help='Number of neighboring tilts used to compare local brightness.'
                            '\n\t- Small values detect very local drops.'
                            '\n\t- Larger values smooth the comparison but may miss punctual failures.'
                            '\n\t- Typical range: 2 – 3 (valid range: 1 – 5).')
        group.addParam(DARK_COS_FLOOR, FloatParam,
                       label='Cosine floor (high tilt stabilization)',
                       default=0.20,
                       validators=[GE(0.05), LE(0.50)],
                       expertLevel=LEVEL_ADVANCED,
                       help='Prevents instability of the brightness model at extreme tilt angles.'
                            '\n\t- Lower values increase sensitivity at high tilts but may be less stable.'
                            '\n\t- Higher values stabilize the model but reduce tilt dependence.'
                            '\n\t- Typical range: 0.1 – 0.3 (valid range: 0.05 – 0.5).')
        group.addParam(DARK_CROP_FRACTION, FloatParam,
                       label='Central crop fraction',
                       default=0.70,
                       validators=[GE(0.40), LE(1.00)],
                       expertLevel=LEVEL_ADVANCED,
                       help='Fraction of the image used to compute brightness (center crop).'
                            '\n\t- Lower values reduce edge artifacts (e.g., grid bars) but may be noisier.'
                            '\n\t- Higher values use more of the image but may include unwanted edges.'
                            '\n\t- Typical range: 0.6 – 0.8 (valid range: 0.4 – 1.0).')
        group.addParam(EXTREME_CONTRAST_TH, FloatParam,
                       label='Extreme contrast threshold',
                       default=1.0,
                       validators=[GE(0.5), LE(3.0)],
                       condition=customQualityCond,
                       help='Flags images with unusually uneven local contrast (e.g., grid bars, shadows).'
                            '\n\t- Lower values = more aggressive (detects mild obstructions).'
                            '\n\t- Higher values = more conservative (only strong artifacts).'
                            '\n\t- Typical range: 0.8 – 1.3 (valid range: 0.5 – 3.0).')
        group.addParam(MIN_ENTROPY, FloatParam,
                       label='Low information / flat image threshold',
                       default=3.0,
                       validators=[GE(1.0), LE(6.0)],
                       condition=customQualityCond,
                       help='Minimum entropy. Flags images with low structural content (flat or low-texture images).'
                            '\n\t- Higher values = more aggressive (requires more texture). '
                            '\n\t- Lower values = more conservative (only very flat images).'
                            '\n\t- Typical range: 2.5 – 3.5 (valid range: 1.0 – 6.0).')
        group.addParam(MIN_EDGE_ENERGY, FloatParam,
                       label='Edge energy threshold',
                       default=0.0,
                       validators=[GE(0.0), LE(0.05)],
                       condition=customQualityCond,
                       help='Minimum edge energy. Flags images with very low edge content (blurred or no sample '
                            'signal).'
                            '\n\t- Higher values = more aggressive (detects mild blur).'
                            '\n\t- Lower values = more conservative (only near-zero edge content).'
                            '\n\t- Typical range: 0.005 – 0.02 (valid range: 0.0 – 0.05).')

        form.addParam(MIN_VIEWS,
                      IntParam,
                      default=30,
                      label="Minimum number of tilts allowed",
                      help='Minimum number of views to include a tilt series.')

        form.addParam(DO_RESTACK, BooleanParam,
                      default=False,
                      label='Re-stack the output tilt-series?',
                      help='If set to No, the output tilt-series will be filtered at metadata level '
                           'instead of generating a new re-stacked file for each tilt-series.')

        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions ---------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        self._initialize()
        inTsSet = self._getInTsSet()
        self.sRate = inTsSet.getSamplingRate()
        self.readingOutput()

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            # In the if statement below, Counter is used because in the tsId comparison the order doesn't matter
            # but duplicates do. With a direct comparison, the closing step may not be inserted because of the order:
            # ['ts_a', 'ts_b'] != ['ts_b', 'ts_a'], but they are the same with Counter.
            if not inTsSet.isStreamOpen() and Counter(self.itemTsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.itemTsIdReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                excId = self._insertFunctionStep(self.excludeViewFilteringStep,
                                                 ts,
                                                 prerequisites=[],
                                                 needsGPU=False)
                closeSetStepDeps.append(excId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.itemTsIdReadList.append(tsId)

            time.sleep(10)
            if inTsSet.isStreamOpen():
                with self._lock:
                    inTsSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions ----------------------------
    def _initialize(self):
        qualityFilterMode = self.getAttribValue(QUALITY_FILTER)
        if qualityFilterMode == QualityFilterModes.balanced.value:
            self.darkSensitivityVal = 2.0
            self.darkLocalRatioThVal = 0.65
            self.darkWindowVal = 2
            self.darkCosFloorVal = 0.20
            self.darkCropFractionVal = 0.70
            self.extremeContrastThVal = 1.0
            self.minEntropyVal = 3.0
            self.minEdgeEnergyVal = 0.0
        elif qualityFilterMode == QualityFilterModes.conservative.value:
            self.darkSensitivityVal = 3.0
            self.darkLocalRatioThVal = 0.55
            self.darkWindowVal = 2
            self.darkCosFloorVal = 0.20
            self.darkCropFractionVal = 0.70
            self.extremeContrastThVal = 1.3
            self.minEntropyVal = 2.5
            self.minEdgeEnergyVal = 0.0
        elif qualityFilterMode == QualityFilterModes.aggressive.value:
            self.darkSensitivityVal = 1.5
            self.darkLocalRatioThVal = 0.75
            self.darkWindowVal = 2
            self.darkCosFloorVal = 0.20
            self.darkCropFractionVal = 0.70
            self.extremeContrastThVal = 0.8
            self.minEntropyVal = 3.5
            self.minEdgeEnergyVal = 0.01
        elif qualityFilterMode == QualityFilterModes.custom.value:
            self.darkSensitivityVal = self.getAttribValue(DARK_SENSITIVITY)
            self.darkLocalRatioThVal = self.getAttribValue(DARK_LOCAL_RATION_TH)
            self.darkWindowVal = self.getAttribValue(DARK_WINDOW)
            self.darkCosFloorVal = self.getAttribValue(DARK_COS_FLOOR)
            self.darkCropFractionVal = self.getAttribValue(DARK_CROP_FRACTION)
            self.extremeContrastThVal = self.getAttribValue(EXTREME_CONTRAST_TH)
            self.minEntropyVal = self.getAttribValue(MIN_ENTROPY)
            self.minEdgeEnergyVal = self.getAttribValue(MIN_EDGE_ENERGY)

    def excludeViewFilteringStep(self, ts: TiltSeries):
        try:
            self._registerOutput(ts)
        except Exception as e:
            logger.error(
                redStr(f'tsId = {ts.getTsId()} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries):
        angleMin = 999.
        angleMax = -999.
        accumDose = 0.
        initialDose = 999.
        doReStack = self.getAttribValue(DO_RESTACK)
        sxThreshold, syThreshold = self._getMaxShiftThresholds(ts)

        imgStack = ImageReadersRegistry.open(ts.getFirstItem().getFileName())

        # Compute image quality analysis once for the whole stack
        darkDict = {}
        if self.getAttributeValue(QUALITY_FILTER) != QualityFilterModes.disabled.value:
            darkDict = self._compute_dark_flags(ts, imgStack, self.darkSensitivityVal)

        finalNoImgs = 0
        tiList = []
        # Tilt-images
        for i, ti in enumerate(ts.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD)):
            tiltAngle = ti.getTiltAngle()
            self._genTiDict(ti)
            # Filter by tilt angle
            self._filterByTiltAngle(ti)
            # Filter by max shift
            if ts.hasAlignment():
                self._filterByMaxShifts(ti, sxThreshold, syThreshold)
            # Filter by dose
            self._filterByDose(ti)
            # Filter by quality
            tiData = imgStack.getImage(i)
            tqd = TiltImageQualityDetector(tiData)
            metricsDict = tqd.analyze_image(tiltAngle)
            if darkDict:
                metricsDict.update(darkDict.get(int(ti.getAcquisitionOrder()), {}))
            self._filterByImgQuality(ti, metricsDict)

            if ti.isEnabled():
                angleMin = min(tiltAngle, angleMin)
                angleMax = max(tiltAngle, angleMax)
                accumDose = max(ti.getAcquisition().getAccumDose(), accumDose)
                initialDose = min(ti.getAcquisition().getDoseInitial(), initialDose)
                newTi = ti.clone()
                tiList.append(newTi)
                finalNoImgs += 1
            elif not doReStack:
                newTi = ti.clone()
                tiList.append(newTi)

        # Now we know the excluded views. Re-stack if requested
        if doReStack:
            presentAcqOrders = set([ti.getAcquisitionOrder() for ti in tiList if ti.isEnabled()])
            inFileName = ts.getFirstItem().getFileName()
            reStackedFn = self._getExtraPath(f'{ts.getTsId()}.mrcs')
            ts.reStack(inFileName, reStackedFn, presentAcqOrders)

        with self._lock:
            minNoViewsAllowed = self.getAttribValue(MIN_VIEWS)
            # Successful TS
            if finalNoImgs >= minNoViewsAllowed:
                # Set of tilt-series
                outTsSet = self.getOutputSetOfTS()
                # Tilt-series
                outTs = TiltSeries()
                outTs.copyInfo(ts)
                outTsSet.append(outTs)

                if doReStack:
                    self._populateRestackedTs(outTs, tiList, reStackedFn, angleMin, angleMax, accumDose, initialDose)
                else:
                    self._populateFinalTs(outTs, tiList)

            # Failed TS
            else:
                # Set of tilt-series
                outTsSet = self.getOutputSetOfTS(failedTs=True)
                # Tilt-series
                outTs = ts.clone()
                outTsSet.append(outTs)
                outTs.copyItems(ts)

                tsId = ts.getTsId()
                logger.info(yellowStr(f'tsId = {tsId} was removed because the number '
                                      f'of tilt-images after filtering [{finalNoImgs}] '
                                      f'is lower than the minimum specified [{minNoViewsAllowed}].'))
                self._updateRemovedTsIds(tsId)

            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)

            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()

    def closeOutputSetsStep(self):
        self._closeOutputSet()
        outputName = self._possibleOutputs.tiltSeries.name
        output = getattr(self, outputName, [])
        if not output or (output and len(output) == 0):
            raise Exception(f'No output {outputName} was generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')
        # Generate a report with the results
        self._genYamlReport()

    # --------------------------- UTILS functions ----------------------------
    def getAttribValue(self, attribName: str):
        return getattr(self, attribName).get()

    def readingOutput(self) -> None:
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None) or []
        for item in outTsSet:
            self.itemTsIdReadList.append(item.getTsId())
        if outTsSet:
            self.info(cyanStr(f'TsIds processed: {self.itemTsIdReadList}'))
        else:
            self.info(cyanStr('No tilt-series have been processed yet'))

    def _getInTsSetPointer(self) -> Pointer:
        return getattr(self, IN_TS_SET)

    def _getInTsSet(self) -> SetOfTiltSeries:
        return self._getInTsSetPointer().get()

    def getOutputSetOfTS(self, failedTs: bool = False) -> SetOfTiltSeries:
        if failedTs:
            attribName = self._possibleOutputs.failedTiltSeries.name
            suffix = 'failed'
        else:
            attribName = self._possibleOutputs.tiltSeries.name
            suffix = ''

        outputSet = getattr(self, attribName, None)
        if isinstance(outputSet, SetOfTiltSeries):
            outputSet.enableAppend()
        else:
            outputSet = SetOfTiltSeries.create(self._getPath(),
                                               template='tiltseries',
                                               suffix=suffix)
            outputSet.copyInfo(self._getInTsSet())
            outputSet.setStreamState(Set.STREAM_OPEN)
            # # Write set properties, otherwise it may expose the set (sqlite) without properties.
            # outputSet.write()
            # Define outputs and relations
            self._defineOutputs(**{attribName: outputSet})
            self._defineSourceRelation(self._getInTsSetPointer(), outputSet)
        return outputSet

    def _getMaxShiftThresholds(self, ts: TiltSeries) -> Tuple[int, int]:
        xdim, ydim, _ = ts.getFirstItem().getDim()
        xdimThreshold = self.getAttribValue(MAX_SX) * xdim
        ydimThreshold = self.getAttribValue(MAX_SY) * ydim
        return xdimThreshold, ydimThreshold

    @staticmethod
    def _getTiId(ti: TiltImage) -> str:
        return f'{ti.getAcquisitionOrder()}@{ti.getTsId()}'

    def _genTiDict(self, ti: TiltImage) -> None:
        self.tiLabelDict[self._getTiId(ti)] = TiMetadata()

    def _updateTiDict(self, ti: TiltImage, label: str) -> None:
        setattr(self.tiLabelDict[self._getTiId(ti)], label, True)

    def _filterByTiltAngle(self, ti: TiltImage) -> None:
        tiltAngle = ti.getTiltAngle()
        if tiltAngle > self.getAttribValue(MAX_TILT) or tiltAngle < self.getAttribValue(MIN_TILT):
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_TILT_ANGLE.value)

    def _filterByMaxShifts(self, ti: TiltImage, sxThreshold: int, syThreshold: int) -> None:
        tm = ti.getTransform().getMatrix()
        sx = tm[0, 2]
        sy = tm[1, 2]
        if abs(sx) > sxThreshold or abs(sy) > syThreshold:
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_MAX_SHIFT.value)

    def _filterByDose(self, ti: TiltImage) -> None:
        dose = ti.getAcquisition().getAccumDose()
        if dose < self.getAttribValue(MIN_DOSE) or dose > self.getAttribValue(MAX_DOSE):
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_DOSE.value)

    def _filterByImgQuality(self, ti: TiltImage, tiMetricsDict: Dict) -> None:
        # Dark outlier (punctual blackout / shutter failure). Controlled by 'darkSensitivity'.
        if tiMetricsDict.get('is_dark', 0) and float(tiMetricsDict.get('is_dark')) > 0:
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_DARK.value)

        if tiMetricsDict["extreme_contrast_score"] > self.extremeContrastThVal:
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_EXTREME_CONTRAST_SPLIT.value)

        if tiMetricsDict["entropy"] < self.minEntropyVal:
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_LOW_INFO.value)

        if tiMetricsDict["edge_energy"] <= self.minEdgeEnergyVal:
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_EDGE_ENERGY.value)

    def _updateRemovedTsIds(self, tsId: str) -> None:
        updatedMsg = self.removedTsIds.get() + f' {tsId}'
        self.removedTsIds.set(updatedMsg)
        self._store(self.removedTsIds)

    @staticmethod
    def _populateRestackedTs(
            outTs: TiltSeries,
            tiList: List[TiltImage],
            reStackedFn: str,
            angleMin: float,
            angleMax: float,
            accumDose: float,
            initialDose: float) -> None:
        # Update the acquisition minAngle and maxAngle values of the tilt-series
        acq = outTs.getAcquisition()
        acq.setAngleMin(angleMin)
        acq.setAngleMax(angleMax)
        acq.setAccumDose(accumDose)
        acq.setDoseInitial(initialDose)
        outTs.setAcquisition(acq)
        # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
        # specific accum and initial dose values
        for tiOut in tiList:
            if tiOut.isEnabled():
                tiOut.setFileName(reStackedFn)
                tiAcq = tiOut.getAcquisition()
                tiAcq.setAngleMin(angleMin)
                tiAcq.setAngleMax(angleMax)
                tiOut.setAcquisition(tiAcq)
                outTs.append(tiOut)
        outTs.setAnglesCount(len(outTs))

    @staticmethod
    def _populateFinalTs(outTs: TiltSeries, tiList: List[TiltImage]) -> None:
        for tiOut in tiList:
            outTs.append(tiOut)

    def _genYamlReport(self):
        fileName = self._getExtraPath('exclude_tilt_series.yaml')
        try:
            with open(fileName, 'w', encoding='utf-8') as file:
                data_to_save = {k: asdict(v) for k, v in self.tiLabelDict.items()}
                with open(fileName, 'w', encoding='utf-8') as f:
                    yaml.dump(data_to_save, f, default_flow_style=False, sort_keys=False)
        except Exception as e:
            print(f"Error generating the yaml report {fileName} with the exception -> {e}")

    def closeOutputsForStreaming(self):
        # Close explicitly the outputs (for streaming)
        for output in self._possibleOutputs:
            output = getattr(self, output.name, None)
            if isinstance(output, Set):
                output.close()

    # --------------------------- INFO functions ---------------------------------------
    def _summary(self):
        summary = []
        if self.isFinished() and self.removedTsIds.get():
            summary.append(f'Some tilt-series were removed: *{self.removedTsIds.get()}*')
        return summary

    # -------------------------- IMAGE QUALITY: DARK OUTLIERS --------------------------
    def _robust_brightness(self,
                           image_data: np.ndarray,
                           p_low: float = 5.0,
                           p_high: float = 95.0,
                           eps: float = 1e-8) -> Dict[str, float]:
        """Compute a robust, scale-invariant brightness descriptor.

        No assumption of stack normalization. Uses a central crop to reduce border / grid-bar influence.
        """
        image_data = np.squeeze(image_data)
        h, w = image_data.shape
        cf = float(np.clip(self.darkCropFractionVal, 0.2, 1.0))
        dh = int(h * cf)
        dw = int(w * cf)
        y0 = max((h - dh) // 2, 0)
        x0 = max((w - dw) // 2, 0)
        crop = image_data[y0:y0 + dh, x0:x0 + dw]

        p05 = float(np.percentile(crop, p_low))
        p50 = float(np.median(crop))
        p95 = float(np.percentile(crop, p_high))
        scale = (p95 - p05)
        brightness = (p50 - p05) / (scale + eps)
        return {'p05': p05, 'p50': p50, 'p95': p95, 'brightness': float(brightness)}

    @staticmethod
    def _robust_linear_fit(x: np.ndarray,
                           y: np.ndarray,
                           iters: int = 6,
                           huber_c: float = 1.345,
                           eps: float = 1e-8) -> Tuple[float, float]:
        """Robust linear fit y ~ a + b*x using IRLS (Iteratively Reweighted Least Squares) with Huber
        weights to prevent outliers influence."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        b, a = np.polyfit(x, y, 1)  # slope, intercept
        for _ in range(max(1, iters)):
            resid = y - (a + b * x)
            med = np.median(resid)
            mad = np.median(np.abs(resid - med))
            sigma = 1.4826 * mad + eps
            u = resid / (huber_c * sigma)
            w = 1.0 / np.maximum(1.0, np.abs(u))
            X = np.vstack([np.ones_like(x), x]).T
            XtW = X.T * w
            beta = np.linalg.lstsq(XtW @ X, XtW @ y, rcond=None)[0]
            a, b = float(beta[0]), float(beta[1])
        return a, b

    def _compute_dark_flags(self,
                            ts: TiltSeries,
                            imgStack,
                            eps: float = 1e-8) -> Dict[int, Dict[str, float]]:
        """Compute per-image dark outlier scores for a tilt-series.

        Assumption: the stack is sorted by tilt angle from min to max.
        Strategy: robust global model log(brightness) vs log(cos(|tilt|)) + local neighborhood check.

        Main idea: An image is considered dark if it is much darker than expected for its angle
        (a global outlier), and it is also much darker than its neighbors (a local outlier).

        Returns a dict keyed by acquisitionOrder with fields: is_dark, dark_z, local_ratio, brightness, p05/p50/p95.
        """
        window = self.darkWindowVal
        ti_by_angle = list(ts.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD))
        n = len(ti_by_angle)
        angles = np.asarray([float(ti.getTiltAngle()) for ti in ti_by_angle], dtype=float)

        # Robust per‑image brightness (without normalization) ----------------------------------------------------

        brightness = np.zeros(n, dtype=float)
        p05 = np.zeros(n, dtype=float)
        p50 = np.zeros(n, dtype=float)
        p95 = np.zeros(n, dtype=float)

        for k, ti in enumerate(ti_by_angle):
            img = imgStack.getImage(k)
            bm = self._robust_brightness(img, eps=eps)
            brightness[k] = max(float(bm['brightness']), eps)
            p05[k], p50[k], p95[k] = bm['p05'], bm['p50'], bm['p95']

        # Global ‘expected’ model versus angle (physical/geometric trend) ---------------------------------------
        cosv = np.cos(np.radians(np.abs(angles)))
        cosv = np.maximum(cosv, self.darkCosFloorVal)
        x = np.log(cosv + eps)
        y = np.log(brightness + eps)
        a, b = self._robust_linear_fit(x, y)
        resid = y - (a + b * x)

        # Global outlier score based on a robust z‑score --------------------------------------------------------
        med_r = np.median(resid)
        mad_r = np.median(np.abs(resid - med_r))
        sigma_r = 1.4826 * mad_r + eps
        z = (resid - med_r) / sigma_r

        # Local check (neighboring-window analysis) -------------------------------------------------------------
        local_ratio = np.ones(n, dtype=float)
        for pos in range(n):
            lo = max(0, pos - window)
            hi = min(n, pos + window + 1)
            neigh = [j for j in range(lo, hi) if j != pos]
            if not neigh:
                continue
            ref = float(np.median(brightness[neigh]))
            local_ratio[pos] = brightness[pos] / (ref + eps)

        # Final decision based on global + local analysis -----------------------------------------------------
        is_dark = (z < -float(self.darkSensitivityVal)) & (local_ratio < self.darkLocalRatioThVal)
        out = {}
        for k, ti in enumerate(ti_by_angle):
            out[int(ti.getAcquisitionOrder())] = {
                'is_dark': float(is_dark[k]),
                'dark_z': float(z[k]),
                'local_ratio': float(local_ratio[k]),
                'brightness': float(brightness[k]),
                'p05': float(p05[k]),
                'p50': float(p50[k]),
                'p95': float(p95[k]),
                'exp_log_brightness': float(a + b * x[k])
            }
        return out


class TiltImageQualityDetector:
    """Per-image quality descriptors for cryo-ET tilt images.

    Darkness detection needs context from the whole tilt-series and is computed at the protocol level,
    then appended to the per-image metrics.
    """

    def __init__(self, image_data: np.ndarray) -> None:
        self.raw_image = np.squeeze(image_data)
        img_min = float(np.min(self.raw_image))
        img_max = float(np.max(self.raw_image))
        self.norm_image = (self.raw_image - img_min) / (img_max - img_min + 1e-8)

    def get_extreme_contrast_metrics(self, grid_size: Tuple[int, int] = (6, 6)) -> Dict[str, float]:
        """Detect extreme partial obstructions via patch-contrast variability."""
        h, w = self.raw_image.shape
        ph, pw = max(h // grid_size[0], 1), max(w // grid_size[1], 1)

        patch_contrast_ranges = []
        for i in range(grid_size[0]):
            for j in range(grid_size[1]):
                p = self.raw_image[i * ph:(i + 1) * ph, j * pw:(j + 1) * pw]
                if p.size == 0:
                    continue
                p10, p90 = np.percentile(p, [10, 90])
                patch_contrast_ranges.append(float(p90 - p10))

        if not patch_contrast_ranges:
            return {'mean_patch_contrast': 0.0, 'std_patch_contrast': 0.0, 'extreme_contrast_score': 0.0}

        mean_contrast = float(np.mean(patch_contrast_ranges))
        std_contrast = float(np.std(patch_contrast_ranges))
        extreme_contrast_score = std_contrast / (mean_contrast + 1e-8)
        return {'mean_patch_contrast': mean_contrast,
                'std_patch_contrast': std_contrast,
                'extreme_contrast_score': extreme_contrast_score}

    def get_entropy(self) -> float:
        """Shannon entropy for texture analysis."""
        img_8bit = (self.norm_image * 255).astype(np.uint8)
        hist, _ = np.histogram(img_8bit, bins=256, range=(0, 256))
        prob = hist / (hist.sum() + 1e-8)
        prob = prob[prob > 0]
        return float(-np.sum(prob * np.log2(prob)))

    def get_edge_energy(self) -> float:
        """Mean edge energy using Sobel filter."""
        edges = sobel(self.norm_image)
        return float(np.mean(edges))

    def analyze_image(self, tilt_angle: float) -> Dict[str, float]:
        """Return per-image metrics (protocol may append dark outlier fields)."""
        extreme_contrast = self.get_extreme_contrast_metrics()
        # bright = self.get_robust_brightness_metrics()
        return {
            'entropy': self.get_entropy(),
            'edge_energy': self.get_edge_energy(),
            'extreme_contrast_score': extreme_contrast['extreme_contrast_score'],
            'tilt_angle': float(tilt_angle)
        }
