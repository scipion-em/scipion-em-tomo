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
from typing import Tuple, Set, List

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import String
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from tomo.protocols import ProtTomoBase

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtTomoBase):
    """
    Assign the transformation matrices from an input set of tilt-series to a target one.
    """

    _label = 'Tilt-series assign alignment'
    _devStatus = BETA
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fromTsDict = None
        self.toTsDict = None
        self.nonMatchingTsIdsMsg = String()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam('getTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series from which transformation matrices will be obtained.',
                      label='Tilt-series from which to take the alignment')

        form.addParam('setTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices will be assigned.',
                      label='Tilt-series to assign the alignment to')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        commonTsIds = self._initialize()
        for tsId in commonTsIds:
            self._insertFunctionStep(self.assignTrMat, tsId,
                                     needsGPU=False)
        self._insertFunctionStep(self.closeOutputSetStep,
                                 needsGPU=False)

    # --------------------------- STEPS functions ----------------------------
    def _initialize(self) -> Set[str]:
        fromTsSet = self.getTMSetOfTiltSeries.get()
        toTsSet = self.setTMSetOfTiltSeries.get()
        commonTsIds, nonCommonTsIds = self._matchTsIds(fromTsSet, toTsSet)
        self.fromTsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in fromTsSet if ts.getTsId() in commonTsIds}
        self.toTsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in toTsSet if ts.getTsId() in commonTsIds}
        if nonCommonTsIds:
            msg = f'Non-matching tsIds: *{nonCommonTsIds}*'
            logger.info(cyanStr(msg))
            self.nonMatchingTsIdsMsg.set(msg)
        return commonTsIds

    def assignTrMat(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - assigning alignment..."))
        fromTs = self.fromTsDict[tsId]
        toTs = self.toTsDict[tsId]
        outTsSet = self.getOutTsSet()

        fromTsSize = fromTs.getSize()
        toTsSize = toTs.getSize()
        if fromTsSize != toTsSize:
            logger.info(cyanStr(f"tsId = {tsId} - The number of tilt-images in the source [{fromTsSize}] "
                                f"and target [{toTsSize}] sets is different."))

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(toTs)
        # The tilt axis angle may have been re-assigned, so it must be updated to keep the coherence with the
        # values of the transformation matrix assigned
        fromTsTAx = fromTs.getAcquisition().getTiltAxisAngle()
        newTs.getAcquisition().setTiltAxisAngle(fromTsTAx)
        outTsSet.append(newTs)

        # Manage the possible previously excluded views or previous ts re-stacking
        presentAcqOrdersFrom = fromTs.getTsPresentAcqOrders()
        presentAcqOrdersTo = toTs.getTsPresentAcqOrders()
        matchingAcqOrders = presentAcqOrdersFrom & presentAcqOrdersTo
        fromTsAcqDir = {ti.getAcquisitonOrder(): ti.clone() for ti in fromTs
                        if ti.getAcquisitonOrder() in matchingAcqOrders}
        toTsAcqDir = {ti.getAcquisitonOrder(): ti.clone() for ti in toTs
                      if ti.getAcquisitonOrder() in matchingAcqOrders}

        for acqOrder in matchingAcqOrders:
            tiFrom = fromTsAcqDir[acqOrder]
            tiTo = toTsAcqDir[acqOrder]
            newTi = TiltImage()
            newTi.copyInfo(tiTo, copyId=True)
            newTi.setLocation(tiTo.getLocation())

            # The tilt axis angle may have been re-assigned or even refined at tilt-image level (and updated
            # consequently in the tilt axis angle field in the metadata), so it must be updated to keep the
            # coherence with the values of the transformation matrix assigned
            fromTiTAx = tiFrom.getAcquisition().getTiltAxisAngle()
            newTi.getAcquisition().setTiltAxisAngle(fromTiTAx)
            newTi.setTiltAngle(tiFrom.getTiltAngle())
            newTransform = self.updateTM(tiFrom.getTransform())
            newTi.setTransform(newTransform)
            newTs.append(newTi)

        newTs.setDim(toTs.getDim())
        newTs.write()
        outTsSet.update(newTs)
        outTsSet.write()
        self._store()

    def closeOutputSetStep(self):
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None)
        if not outTsSet:
            raise Exception('No transformation matrix assignment was carried out. Please '
                            'check the Output Log > run.stdout and run.stderr')
        self._closeOutputSet()

    # --------------------------- UTILS functions ----------------------------
    def getOutTsSet(self):
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None)
        if outTsSet:
            outTsSet.enableAppend()
        else:
            outTsSet = SetOfTiltSeries.create(self._getPath(),
                                              template='tiltseries',
                                              suffix='AssignedTransform')
            toTsSet = self.setTMSetOfTiltSeries.get()
            fromTsSet = self.getTMSetOfTiltSeries.get()
            outTsSet.copyInfo(toTsSet)
            outTsSet.setDim(toTsSet.getDim())
            # The tilt axis angle may have been re-assigned, so it must be updated to keep the coherence with the
            # values of the transformation matrix assigned
            fromTsSetTAx = fromTsSet.getAcquisition().getTiltAxisAngle()
            outTsSet.getAcquisition().setTiltAxisAngle(fromTsSetTAx)

            self._defineOutputs(**{self._possibleOutputs.tiltSeries.name: outTsSet})
            self._defineSourceRelation(self.getTMSetOfTiltSeries, outTsSet)
            self._defineSourceRelation(self.setTMSetOfTiltSeries, outTsSet)
        return outTsSet

    def getSamplingRatio(self):
        return self.setTMSetOfTiltSeries.get().getSamplingRate() / self.getTMSetOfTiltSeries.get().getSamplingRate()

    def updateTM(self, transform):
        """ Scale the transform matrix shifts. """
        matrix = transform.getMatrix()

        sr = self.getSamplingRatio()

        matrix[0][2] /= sr
        matrix[1][2] /= sr

        transform.setMatrix(matrix)

        return transform

    @staticmethod
    def _matchTsIds(fromTsSet: SetOfTiltSeries, toTsSet: SetOfTiltSeries) -> Tuple[Set[str], Set[str]]:
        fromTsIds = fromTsSet.getTSIds()
        toTsIds = toTsSet.getTSIds()
        setCastedFromTsIds = set(fromTsIds)
        setCastedToTsIds = set(toTsIds)
        commonTsIds = setCastedFromTsIds & setCastedToTsIds  # Intersection, common elements
        nonCommonTsIds = setCastedFromTsIds ^ setCastedToTsIds  # Symmetric difference, non-common elements
        return commonTsIds, nonCommonTsIds

    # --------------------------- INFO functions ----------------------------
    def _validate(self) -> List[str]:
        validateMsgs = []
        # The from TS set is expected to have aligment
        for tsGetTM in self.getTMSetOfTiltSeries.get():
            if not tsGetTM.hasAlignment():
                validateMsgs.append("Tilt-series %s from the input set do not have a "
                                    "transformation matrix assigned." % tsGetTM.getTsId())
                break
        # Check the tsId matching between the from and to TS sets
        fromTsSet = self.getTMSetOfTiltSeries.get()
        toTsSet = self.setTMSetOfTiltSeries.get()
        commonTsIds, _ = self._matchTsIds(fromTsSet, toTsSet)
        if not commonTsIds:
            validateMsgs.append('There are no matching tsIds between the sets of tilt-series introduced.')
        return validateMsgs

    def _summary(self):
        summary = []
        nonMatchingTsIdsMsg = self.nonMatchingTsIdsMsg.get()
        outputTSName = self._possibleOutputs.tiltSeries.name
        if nonMatchingTsIdsMsg:
            summary.append(nonMatchingTsIdsMsg)
        if hasattr(self, outputTSName):
            outTsSet = getattr(self, outputTSName)
            summary.append(f"Input tilt-series:"
                           f"\n\t- Get Transform: {self.getTMSetOfTiltSeries.get().getSize()}"
                           f"\n\t- Set Transform: {self.setTMSetOfTiltSeries.get().getSize()}"
                           f"\nTransformation matrices assigned: {outTsSet.getSize()}\n")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name)
        if outTsSet:
            methods.append("The transformation matrix has been assigned to %d tilt-series from the input set.\n"
                           % (outTsSet.getSize()))
        else:
            methods.append("Outputs are not ready yet.")
        return methods
