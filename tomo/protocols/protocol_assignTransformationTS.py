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
import sqlite3
import traceback
import typing
from enum import Enum
from os.path import exists
from typing import List, Union, Dict
import numpy as np
from pwem import getExecStatusDir, appendStreamItem
from pwem.objects import Transform
from pyworkflow import BETA
from pwem.protocols import EMProtocol
from pyworkflow.object import Pointer, Set
from pyworkflow.protocol import STEPS_PARALLEL, PointerParam
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from tomo.protocols.protocol_base_streaming_tomo import ProtocolBaseStreamingTomo
from tomo.utils import getTsIdsIntersection, getTsIdsDicts, writeTsSidecar

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtocolBaseStreamingTomo):
    """
    Assign the transformation matrices from an input set of tilt-series to a target one.
    """

    _label = 'Tilt-series assign alignment'
    _devStatus = BETA
    _possibleOutputs = outputObjects
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsFromDict = None
        self.tsToDict = None
        self.tsIdsReadFrom = []
        self.tsIdsReadTo = []
        self.sRateRatio = None
        self.failedItems = []

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
    def _insertAllSteps(self) -> None:
        inTsSetFrom = self.getInTsSetFrom()
        inTsSetTo = self.getInTsSetTo()
        self.sRateRatio = inTsSetTo.getSamplingRate() / inTsSetFrom.getSamplingRate()
        if inTsSetFrom.isStreamOpen() or inTsSetTo.isStreamOpen():
            self._insertFunctionStep(self.stepsGeneratorStep,
                                     prerequisites=[],
                                     needsGPU=False)
        else:
            self._insertNonStreamingSteps()

    # Streaming Hooks ############################
    def _getStreamingInputSets(self):
        return [self.getInTsSetFrom(), self.getInTsSetTo()]

    def _getStreamingOutputNames(self) -> str:
        return self._possibleOutputs.tiltSeries.name

    def _discoverReadyWork(self, tsIds, inputSets):
        # Rebuild the ready TSfrom and TSTo from their OWN producers' sidecars
        # (no live-DB read) and join by tsId. A tsId whose CTF is not yet
        # materialisable is skipped and retried next cycle.
        tsFromDict = self.getInTsSetFrom().fetchNewItems(tsIds)
        tsToDict = self.getInTsSetTo().fetchNewItems(tsIds)
        work = {}
        for tsId, tsFrom in tsFromDict.items():
            tsTo = tsToDict.get(tsId)
            if tsTo is None:
                logger.info(yellowStr(f'tsId = {tsId} - no corresponding tsTo found yet, retrying...'))
                continue
            work[tsId] = (tsFrom, tsTo)
        return work
    # End of streaming hooks #####################

    def _insertNonStreamingSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId, tsFrom in self.tsFromDict.items():
            tsTo = self.tsToDict[tsId]
            self._insertCommonSteps(tsFrom, tsTo, closeSetStepDeps=closeSetStepDeps)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    def _insertCommonSteps(self, *stepsInputs, closeSetStepDeps: List[int]) -> None:
        tsFrom, tsTo = stepsInputs
        pId = self._insertFunctionStep(self.assignTrMatStep,
                                       tsFrom,
                                       tsTo,
                                       prerequisites=[],
                                       needsGPU=False)
        closeSetStepDeps.append(pId)

    # --------------------------- STEPS functions ----------------------------
    def _initialize(self):
        inTsSetFrom = self.getInTsSetFrom()
        inTsSetTo = self.getInTsSetTo()
        commonTsIds = getTsIdsIntersection(inTsSetFrom, inTsSetTo)
        self.tsFromDict, self.tsToDict = getTsIdsDicts(inTsSetFrom, inTsSetTo, present_ts_ids=commonTsIds)

    def assignTrMatStep(self, tsFrom: TiltSeries, tsTo: TiltSeries):
        tsId = tsFrom.getTsId()
        logger.info(cyanStr(f"tsId = {tsId} - assigning alignment..."))
        try:
            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(tsTo)

            # The tilt axis angle may have been re-assigned, so it must be updated
            # to keep the coherence with the values of the transformation matrix assigned
            fromTsTAx = tsFrom.getAcquisition().getTiltAxisAngle()
            newTs.getAcquisition().setTiltAxisAngle(fromTsTAx)
            newTs.setDim(tsTo.getDim())
            newTs.setAlignment2D()

            # Manage the possible previously excluded views or previous ts re-stacking
            matchingAcqOrders = self._getCommonAcqOrderInTsPair(tsFrom, tsTo)
            fromTsAcqDict = {ti.getAcquisitionOrder(): ti.clone() for ti in tsFrom}

            newTiList = []
            for tiTo in tsTo.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD):
                newTi = self._processTiltImage(tiTo, fromTsAcqDict, matchingAcqOrders)
                newTiList.append(newTi)

            self._registerOutput(newTs, newTiList)

            # Streaming only: publish the per-TS metadata sidecar (built from the
            # in-memory ts/tiltImages, no DB read) and the journal id. Without this,
            # a downstream streaming consumer (e.g. ProtImodCtfCorrection) sees this
            # producer's status dir -- created by the inherited stepsGeneratorStep --
            # and reads readiness from its (otherwise empty) stream journal, so it
            # never discovers any ready tsId and never generates per-item steps.
            execStatusDir = getExecStatusDir(self)
            if exists(execStatusDir):
                writeTsSidecar(execStatusDir, newTs, newTiList)
                appendStreamItem(self, tsId)
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> failed: {e}'))
            logger.error(traceback.format_exc())
            self.failedItems.append(tsId)

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, newTs: TiltSeries, tiList: List[TiltImage]) -> None:
        with self._lock:
            try:
                # Set of tilt-series
                outTsSet = self.getOutTsSet()
                # Tilt-series
                outTsSet.append(newTs)
                # Tilt-images
                for newTi in tiList:
                    newTs.append(newTi)
                # Data persistance
                newTs.write()
                outTsSet.update(newTs)
                outTsSet.write()
                self._store(outTsSet)
            except sqlite3.OperationalError as e:
                # Release the write lock and reset the in-memory append state so
                # the @retry_on_sqlite_lock retry is a clean, non-hogging redo
                # (covers the later commits -- newTs.write/outTsSet.write -- not
                # just the append phase) and never trips the duplicate-tsId guard.
                self._releaseOutputWriteLock(outTsSet, newTs.getTsId())
                raise e

    def _processTiltImage(self,
                          tiTo: TiltImage,
                          fromTsAcqDict: Dict[int, TiltImage],
                          matchingAcqOrders: typing.Set[int]):
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

    @staticmethod
    def _getCommonAcqOrderInTsPair(ts1: TiltSeries, ts2: TiltSeries) -> typing.Set[int]:
        tsAcqOrderSet1 = {ti.getAcquisitionOrder() for ti in ts1}
        tsAcqOrderSet2 = {ti.getAcquisitionOrder() for ti in ts2}
        return tsAcqOrderSet1 & tsAcqOrderSet2

    def getOutTsSet(self) -> SetOfTiltSeries:
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
