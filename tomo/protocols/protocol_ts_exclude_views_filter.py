# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
from enum import Enum
import time
from typing import Union, List, Counter

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam
from pyworkflow.object import Set, Pointer
from pyworkflow.utils import Message, makePath, cyanStr, redStr
from pwem.protocols import EMProtocol
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries
from tomo.protocols import ProtTomoBase

logger = logging.getLogger(__name__)

#Form variables
IN_TS_SET = 'inTsSet'

OUTPUT_TS_FAILED_NAME = "FailedTiltSeries"


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtExclViewFilter(EMProtocol, ProtTomoBase):
    """
    This protocol allows to filter a set of aligned tilt series according to a set of parameter as they are:\n
     * Maximum allowed shift after the tilt series alignment\n
     * Range of accumulated dose\n
     * Range of tilt angle
     * Minimum number of views
    """
    _label = 'exclude views filter'
    _devStatus = BETA
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.itemTsIdReadList = []
        self.failedItems = []
        self.sRate = -1

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam(IN_TS_SET,
                      PointerParam,
                      important=True,
                      label="Tilt series",
                      pointerClass='SetOfTiltSeries',
                      help='Select several sets of tilt-series where to evaluate the consensus in their alignment. '
                           'Output set will bring the information from the first selected set.')

        lineShift = form.addLine('Filter by shift (px)',
                                help='This is the minimum/maximum shift allowed along the X or Y direction.'
                                     ' A value of 0.1 means that a 10% of the dimensions of the tilt image '
                                     'is allowed.')
        lineShift.addParam('maxShiftX', FloatParam, default=0.1, label="Max shift in X")
        lineShift.addParam('maxShiftY', FloatParam, default=0.1, label="Max shift in Y")

        lineTilt = form.addLine('Filter by tilt angle (deg)',
                                help='This is the minimum/maximum tilt angle allowed.'
                                     'Only the tilt images with tilt angles in the range min<tilt<max are'
                                     'allowed.')
        lineTilt.addParam('mintilt', FloatParam, default=-30.0, label="Min")
        lineTilt.addParam('maxtilt', FloatParam, default=30.0,  label="Max")

        lineDose = form.addLine('Filter by dose (e/A^2)',
                                help='This is the minimum/maximum dose per tilt image.'
                                     'Only the tilt images with accumulated dose in this range will be kept.')
        lineDose.addParam('minDose', FloatParam, default=0.0, label="Min")
        lineDose.addParam('maxDose', FloatParam, default=70.0,  label="Max")

        form.addParam('minViews',
                      IntParam,
                      default = 30,
                      label="Min number of views",
                      help='Minimum number of views to include a tilt series.')

    # -------------------------- INSERT steps functions ---------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        inTsSet = self._getInTsSet()
        self.sRate = self._getInTsSet().getSamplingRate()
        self.readingOutput()

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            # In the if statement below, Counter is used because in the tsId comparison the order doesnt matter
            # but duplicates do. With a direct comparison, the closing step may not be inserted because of the order:
            # ['ts_a', 'ts_b'] != ['ts_b', 'ts_a'], but they are the same with Counter.
            if not inTsSet.isStreamOpen() and Counter(self.itemTsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self._closeOutputSet,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.itemTsIdReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                excId = self._insertFunctionStep(self.excludeViewFilteringStep, tsId,
                                                     needsGPU=True)
                closeSetStepDeps.append(excId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.itemTsIdReadList.append(tsId)

            time.sleep(10)
            if inTsSet.isStreamOpen():
                with self._lock:
                    inTsSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions ----------------------------
    def excludeViewFilteringStep(self):
        inputData = self._getInTsSet()

        suffix = '_exclViews'
        newSetTs = self._createSet(SetOfTiltSeries, 'tiltseries%s.sqlite', suffix)
        newSetTs.copyInfo(inputData)
        newSetTs.copyAttributes(inputData)

        xdim, ydim, ndim = inputData.getFirstItem().getDim()
        xdimThreshold = self.maxShiftX.get()*xdim
        ydimThreshold = self.maxShiftY.get()*ydim
        minTiltThreshold = self.mintilt.get()
        maxTiltThreshold = self.maxtilt.get()
        minDoseThreshold = self.minDose.get()
        maxDoseThreshold = self.maxDose.get()
        #tsi = ts.getFirstEnabledITem()

        for ts in inputData:

            newTs = TiltSeries(tsId=ts.getTsId())
            newTs.copyInfo(ts)
            newSetTs.append(newTs)
            tomoAcq = ts.getAcquisition()
            tomoAcq.setAngleMax(tomoAcq.getAngleMax())
            tomoAcq.setAngleMax(tomoAcq.getAngleMin())
            newTs.setAcquisition(tomoAcq)
            counter = 0
            for ti in ts:
                newTi = TiltImage()
                newTi.copyInfo(ti)
                tiltAngle = ti.getTiltAngle()

                dose = ti.getAcquisition().getAccumDose()

                if tiltAngle>maxTiltThreshold or tiltAngle<minTiltThreshold:
                    newTi.setEnabled(False)
                if ts.hasAlignment():
                    tm = ti.getTransform().getMatrix()
                    sx = tm[0, 2]
                    sy = tm[1, 2]
                    if sx>xdimThreshold or sy>ydimThreshold:
                        newTi.setEnabled(False)

                if dose<minDoseThreshold or dose>maxDoseThreshold:
                    newTi.setEnabled(False)

                newTs.append(newTi)
                counter = counter + 1

            if counter>=self.minViews.get():
                newTs.setDim(ts.getDim())
                newTs.write()
                newSetTs.update(ts)

        newSetTs.write()
        self._store()

        self._defineOutputs(**{self._possibleOutputs.tiltSeries.name: newSetTs})
        self._defineSourceRelation(inputData, newSetTs)

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def readingOutput(self) -> None:
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None)
        if outTsSet:
            for item in outTsSet:
                self.itemTsIdReadList.append(item.getTsId())
            self.info(cyanStr(f'TsIds processed: {self.itemTsIdReadList}'))
        else:
            self.info(cyanStr('No tilt-series have been processed yet'))

    def _getInTsSet(self, returnPointer: bool = False) -> Union[SetOfTiltSeries, Pointer]:
        inTsPointer = getattr(self, IN_TS_SET)
        return inTsPointer if returnPointer else inTsPointer.get()
    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary


