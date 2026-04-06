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
from typing import Counter, Tuple, List, Optional
import numpy as np
import yaml
import time
from pwem.emlib.image.image_readers import ImageReadersRegistry
from pyworkflow import BETA
from pyworkflow.protocol import STEPS_PARALLEL, BooleanParam, ProtStreamingBase
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, GE, LE
from pyworkflow.object import Set, Pointer, String
from pyworkflow.utils import cyanStr, Message, redStr, yellowStr
from pwem.protocols import EMProtocol
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

logger = logging.getLogger(__name__)

EXCL_VIEWS_SUFFIX = '_exclViews'

# Form variables
IN_TS_SET = 'inTsSet'
MIN_TILT = 'mintilt'
MAX_TILT = 'maxtilt'
MAX_SX = 'maxShiftX'
MAX_SY = 'maxShiftY'
MIN_DOSE = 'minDose'
MAX_DOSE = 'maxDose'
DARK_TH = 'darkSensitivity'
MIN_VIEWS = 'minViews'
DO_RESTACK = 'doReStack'

# Tilt-series annotation keys
@dataclass
class TiMetadata:
    by_max_shift: bool = False
    by_tilt_angle: bool = False
    by_dose: bool = False
    by_dark: bool = False


class TiLabel(Enum):
    BY_MAX_SHIFT = 'by_max_shift'
    BY_TILT_ANGLE = 'by_tilt_angle'
    BY_DOSE = 'by_dose'
    BY_DARK = 'by_dark'


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
        lineShift.addParam(MAX_SX, FloatParam, default=0.0, validators=[GE(0), LE(1)], label="X   ")
        lineShift.addParam(MAX_SY, FloatParam, default=0.0, validators=[GE(0), LE(1)], label="Y    ")

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

        lineDark = group.addLine('Filter by dark sensitivity',
                                 help='Values lower than 0 means that no dark filter will be applied. Behavior:\n\n'
                                      '- *High Sensitivity: Low Factor --> 1.0 - 1.5.* It might flag images '
                                      'that are just slightly darker than the average, such as high-tilt images '
                                      'where the ice is naturally thicker.\n\n'
                                      '- *Balanced: around 2.0.* It ignores the natural darkening of high '
                                      'tilts but will catch "heavy" shadows or partial grid bars.\n\n'
                                      '- *Low Sensitivity: High Factor --> greater than 3.0.* It will only flag '
                                      '"catastrophic" failures, like a solid copper grid bar completely blocking the '
                                      'electron beam (total blackouts).')

        lineDark.addParam(DARK_TH, FloatParam,
                          default=2.0,
                          label='Dark sensitivity factor', )

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
                excId = self._insertFunctionStep(self.excludeViewFilteringStep, ts,
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
        sxThreshold, syThreshold = self._getMaxShiftThresholds(ts)
        darkImgIndices = self._getDarkImgIndices(ts)
        with self._lock:
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS()
            # Tilt-series
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            outTsSet.append(outTs)
            finalNoImgs = 0
            tiList = []
            # Tilt-images
            for ti in ts:
                newTi = TiltImage()
                newTi.copyInfo(ti)
                self._genTiDict(ti)
                # Filter by tilt angle
                self._filterByTiltAngle(ti)
                # Filter by max shift
                if ts.hasAlignment():
                    self._filterByMaxShifts(ti, sxThreshold, syThreshold)
                # Filter by dose
                self._filterByDose(ti)
                # Filter dark images
                self._filterByDarkImgs(ti, darkImgIndices)

                tiltAngle = ti.getTiltAngle()
                angleMin = min(tiltAngle, angleMin)
                angleMax = max(tiltAngle, angleMax)
                accumDose = max(ti.getAcquisition().getAccumDose(), accumDose)
                initialDose = min(ti.getAcquisition().getDoseInitial(), initialDose)

                tiList.append(ti)
                finalNoImgs += 1

            minNoViewsAllowed = self.getAttribValue(MIN_VIEWS)
            if finalNoImgs >= minNoViewsAllowed:
                if self.getAttribValue(DO_RESTACK):
                    self._populateRestackedTs(outTs, tiList, angleMin, angleMax, accumDose, initialDose)
                else:
                    self._populateFinalTs(outTs, tiList)

                outTs.write()
                outTsSet.update(outTs)
                outTsSet.write()
                self._store(outTsSet)
            else:
                tsId = ts.getTsId()
                logger.info(yellowStr(f'tsId = {tsId} was removed because the number '
                                      f'of tilt-images after filtering [{finalNoImgs}] '
                                      f'is lower than the minimum specified [{minNoViewsAllowed}].'))
                self._updateRemovedTsIds(tsId)
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

    def getOutputSetOfTS(self) -> SetOfTiltSeries:
        attrName = self._possibleOutputs.tiltSeries.name
        outputSet = getattr(self, attrName, None)
        if isinstance(outputSet, SetOfTiltSeries):
            outputSet.enableAppend()
        else:
            outputSet = SetOfTiltSeries.create(self._getPath(),
                                               template='tiltseries',
                                               suffix=EXCL_VIEWS_SUFFIX)
            outputSet.copyInfo(self._getInTsSet())
            outputSet.setStreamState(Set.STREAM_OPEN)
            # Write set properties, otherwise it may expose the set (sqlite) without properties.
            outputSet.write()
            # Define outputs and relations
            self._defineOutputs(**{attrName: outputSet})
            self._defineSourceRelation(self._getInTsSetPointer(), outputSet)
        return outputSet

    def _getMaxShiftThresholds(self, ts: TiltSeries) -> Tuple[int, int]:
        xdim, ydim, _ = ts.getFirstItem().getDim()
        xdimThreshold = self.getAttribValue(MAX_SX) * xdim
        ydimThreshold = self.getAttribValue(MAX_SY) * ydim
        return xdimThreshold, ydimThreshold

    def _getDarkImgIndices(self, ts: TiltSeries) -> Optional[List[int]]:
        darkSensitivity = self.getAttribValue(DARK_TH)
        if darkSensitivity > 0.0:
            imgStack = ImageReadersRegistry.open(ts.getFirstItem().getFileName())
            # Statistic indicators
            medians = np.fromiter((np.median(img) for img in imgStack), dtype=np.float64)
            q1, q3 = np.percentile(medians, [25.0, 75.0])
            iqr = q3 - q1
            lowLimit = q1 - (darkSensitivity * iqr)
            # Dark images indices
            darkImgsIndices = np.where(medians < lowLimit)[0]
            return darkImgsIndices.tolist()
        return None

    @staticmethod
    def _getTiId(ti: TiltImage) -> str:
        return f'{ti.getAcquisitionOrder()}@{ti.getTsId()}'

    def _genTiDict(self, ti:TiltImage) -> None:
        self.tiLabelDict[self._getTiId(ti)] = TiMetadata()

    def _updateTiDict(self, ti:TiltImage, label: str) -> None:
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

    def _filterByDarkImgs(self, ti: TiltImage, darkImgIndices: Optional[List[int]]) -> None:
        if not isinstance(darkImgIndices, list):
            return
        if ti.getIndex() in darkImgIndices:
            ti.setEnabled(False)
            self._updateTiDict(ti, TiLabel.BY_DARK.value)

    def _updateRemovedTsIds(self, tsId: str) -> None:
        updatedMsg = self.removedTsIds.get() + f' {tsId}'
        self.removedTsIds.set(updatedMsg)
        self._store(self.removedTsIds)

    @staticmethod
    def _populateRestackedTs(
            outTs: TiltSeries,
            tiList: List[TiltImage],
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

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if self.isFinished() and self.removedTsIds.get():
            summary.append(f'Some tilt-series were removed: *{self.removedTsIds.get()}*')
        return summary
