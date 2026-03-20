# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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
import time
import traceback
import typing
from collections import Counter
from enum import Enum
from typing import List, Union
import numpy as np
from pwem.objects import Transform
from pyworkflow import BETA
from pwem.protocols import EMProtocol
from pyworkflow.object import Pointer, Set
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase, PointerParam
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtStreamingBase):
    """
    Assign the transformation matrices from an input set of tilt-series to a target one.
    """

    _label = 'Tilt-series assign alignment'
    _devStatus = BETA
    _possibleOutputs = outputObjects
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsIdsReadFrom = []
        self.tsIdsReadTo = []
        self.sRateRatio = None

    @classmethod
    def worksInStreaming(cls):
        return True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam('getTMSetOfTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series from which transformation matrices will be obtained.',
                      label='Tilt-series from which to take the alignment')

        form.addParam('setTMSetOfTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices will be assigned.',
                      label='Tilt-series to assign the alignment to')

        form.addParallelSection(threads=3, mpi=0)

    # -------------------------- INSERT steps functions ---------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None)
        inTsSetFrom = self.getInTsSetFrom()
        self.readingOutput(outTsSet)
        inTsSetTo = self.getInTsSetTo()
        self.readingOutput(outTsSet, tsSetFrom=False)
        self.sRateRatio = inTsSetTo.getSamplingRate() / inTsSetFrom.getSamplingRate()

        while True:
            with self._lock:
                inTsIdsFrom = set(inTsSetFrom.getTSIds())
                inTsIdsTo = set(inTsSetTo.getTSIds())
                presentTsIds = inTsIdsFrom & inTsIdsTo

            if ((not inTsSetFrom.isStreamOpen() and Counter(self.tsIdsReadFrom) == Counter(presentTsIds)) and
                    (not inTsSetTo.isStreamOpen() and Counter(self.tsIdsReadTo) == Counter(presentTsIds))):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIdsFrom = inTsIdsFrom - set(self.tsIdsReadFrom)
            nonProcessedTsIdsTo = inTsIdsTo - set(self.tsIdsReadTo)
            tsFrom2ProcessDict = {tsId: ts.clone() for ts in inTsSetFrom.iterItems()
                                  if (tsId := ts.getTsId()) in nonProcessedTsIdsFrom  # Only not processed tsIds (from)
                                  and ts.getSize() > 0}  # Avoid processing empty TS
            tsTo2ProcessDict = {tsId: ts.clone() for ts in inTsSetTo.iterItems()
                                if (tsId := ts.getTsId()) in nonProcessedTsIdsTo  # Only not processed tsIds (to)
                                and ts.getSize() > 0}  # Avoid processing empty CTFs

            for tsId, tsFrom in tsFrom2ProcessDict.items():
                tsTo = tsTo2ProcessDict.get(tsId, None)
                if not tsTo:
                    logger.info(yellowStr(f'tsId = {tsId} - no corresponding tsTo to tsFrom was found...'))
                    continue
                pId = self._insertFunctionStep(self.assignTrMatStep,
                                               tsId,
                                               tsFrom,
                                               tsTo,
                                               prerequisites=[],
                                               needsGPU=False)
                closeSetStepDeps.append(pId)
                logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                self.tsIdsReadFrom.append(tsId)
                self.tsIdsReadTo.append(tsId)

            self.refreshStreaming(inTsSetFrom)
            self.refreshStreaming(inTsSetTo)

    def refreshStreaming(self, inSet: SetOfTiltSeries) -> None:
        # Refresh status for the streaming
        time.sleep(10)
        if inSet.isStreamOpen():
            with self._lock:
                inSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions ----------------------------
    def assignTrMatStep(self, tsId: str, tsFrom: TiltSeries, tsTo: TiltSeries):
        logger.info(cyanStr(f"tsId = {tsId} - assigning alignment..."))
        try:
            outTsSet = self.getOutTsSet()
            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(tsTo)
            # The tilt axis angle may have been re-assigned, so it must be updated
            # to keep the coherence with the values of the transformation matrix assigned
            fromTsTAx = tsFrom.getAcquisition().getTiltAxisAngle()
            newTs.getAcquisition().setTiltAxisAngle(fromTsTAx)
            outTsSet.append(newTs)

            # Manage the possible previously excluded views or previous ts re-stacking
            matchingAcqOrders = self._getCommonAcqOrderInTsPair(tsFrom, tsTo)
            fromTsAcqDict = {ti.getAcquisitionOrder(): ti.clone() for ti in tsFrom}

            for tiTo in tsTo.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD):
                newTi = self._processTiltImage(tiTo, fromTsAcqDict, matchingAcqOrders)
                newTs.append(newTi)

            newTs.setDim(tsTo.getDim())
            newTs.write()
            outTsSet.update(newTs)
            outTsSet.write()
            self._store()

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> failed: {e}'))
            logger.error(traceback.format_exc())

    def _processTiltImage(self, tiTo, fromTsAcqDict, matchingAcqOrders):
        acqOrder = tiTo.getAcquisitionOrder()

        if acqOrder in matchingAcqOrders:
            tiFrom = fromTsAcqDict[acqOrder]
            newTi = TiltImage()
            newTi.copyInfo(tiFrom)
            newTi.setFileName(tiTo.getFileName())
            newTi.setAcquisition(tiTo.getAcquisition())

            # The tilt axis angle may have been re-assigned or even refined at tilt-image
            # level (and updated consequently in the tilt axis angle field in the metadata),
            # so it must be updated to keep the coherence with the values of the transformation
            # matrix assigned
            fromTiTAx = tiFrom.getAcquisition().getTiltAxisAngle()
            newTi.getAcquisition().setTiltAxisAngle(fromTiTAx)
            newTi.setTiltAngle(tiFrom.getTiltAngle())
            self.updateTiTrMatrix(newTi)
            return newTi

        # Case of disabled views
        newTi = tiTo.clone()
        t = Transform()
        t.setMatrix(np.identity(3))
        newTi.setTransform(t)
        newTi.setEnabled(False)
        return newTi

    def closeOutputSetsStep(self):
        self._closeOutputSet()
        attribName = self._possibleOutputs.tiltSeries.name
        output = getattr(self, attribName, None)
        if not output or (output and len(output) == 0):
            raise Exception(f'No output/s {attribName} were generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')

    # --------------------------- UTILS functions ----------------------------
    def getInTsSetFrom(self, asPointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        return self.getTMSetOfTiltSeries if asPointer else self.getTMSetOfTiltSeries.get()

    def getInTsSetTo(self, asPointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        return self.setTMSetOfTiltSeries if asPointer else self.setTMSetOfTiltSeries.get()

    def readingOutput(self,
                      outSet: SetOfTiltSeries,
                      tsSetFrom: bool = True) -> None:
        if outSet:
            if tsSetFrom:
                tsIdList = self.tsIdsReadFrom
                inObjStr = 'tsFrom'
            else:
                tsIdList = self.tsIdsReadTo
                inObjStr = 'tsTo'
            for item in outSet:
                tsIdList.append(item.getTsId())
            self.info(cyanStr(f'{inObjStr}: items processed {tsIdList}'))
        else:
            self.info(cyanStr('No items have been processed yet'))

    @staticmethod
    def _getCommonAcqOrderInTsPair(ts1: TiltSeries, ts2: TiltSeries) -> typing.Set[int]:
        tsAcqOrderSet1 = {ti.getAcquisitionOrder() for ti in ts1}
        tsAcqOrderSet2 = {ti.getAcquisitionOrder() for ti in ts2}
        return tsAcqOrderSet1 & tsAcqOrderSet2

    def getOutTsSet(self):
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None)
        if outTsSet:
            outTsSet.enableAppend()
        else:
            outTsSet = SetOfTiltSeries.create(self._getPath(),
                                              template='tiltseries',
                                              suffix='assignedTransform')
            fromTsSet = self.getInTsSetFrom()
            toTsSet = self.getInTsSetTo()
            outTsSet.copyInfo(toTsSet)
            outTsSet.setDim(toTsSet.getDim())
            # The tilt axis angle may have been re-assigned, so it must be updated to
            # keep the coherence with the values of the transformation matrix assigned
            fromTsSetTAx = fromTsSet.getAcquisition().getTiltAxisAngle()
            outTsSet.getAcquisition().setTiltAxisAngle(fromTsSetTAx)
            outTsSet.setStreamState(Set.STREAM_OPEN)
            # Write set properties, otherwise it may expose the set (sqlite) without properties.
            outTsSet.write()

            self._defineOutputs(**{self._possibleOutputs.tiltSeries.name: outTsSet})
            self._defineSourceRelation(self.getInTsSetFrom(asPointer=True), outTsSet)
            self._defineSourceRelation(self.getInTsSetTo(asPointer=True), outTsSet)
        return outTsSet

    @staticmethod
    def _getTsSize(ts: TiltSeries) -> int:
        stackSize = ts.getSize()
        metadataSize = len([enabled for ti in ts.iterItems() if (enabled := ti.isEnabled())])
        return min(stackSize, metadataSize)

    def updateTiTrMatrix(self, ti: TiltImage) -> None:
        """ Scale the transform matrix shifts. """
        transform = ti.getTransform()
        matrix = transform.getMatrix()
        matrix[0][2] /= self.sRateRatio
        matrix[1][2] /= self.sRateRatio
        transform.setMatrix(matrix)
        ti.setTransform(transform)

    # --------------------------- INFO functions ----------------------------
    def _validate(self) -> List[str]:
        validateMsgs = []
        fromTsSet = self.getInTsSetFrom()
        # The "from" TS set is expected to have alignment
        for ts in fromTsSet.iterItems():
            if not ts.hasAlignment():
                validateMsgs.append("Tilt-series %s from the input set do not have a "
                                    "transformation matrix assigned." % ts.getTsId())
                break
        return validateMsgs

    def _summary(self):
        summary = []
        outputTSName = self._possibleOutputs.tiltSeries.name
        if hasattr(self, outputTSName):
            outTsSet = getattr(self, outputTSName)
            summary.append(f"\nTransformation matrices assigned: {outTsSet.getSize()}\n")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

