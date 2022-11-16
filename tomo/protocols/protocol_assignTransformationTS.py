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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol

import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtTomoBase):
    """
    Assign the transformation matrices from an input set of tilt-series to a target one.
    """

    _label = 'Tilt-series assign alignment'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('getTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series from which transformation matrices will be obtained.',
                      label='Tilt-series WITH alignment')

        form.addParam('setTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices will be assigned.',
                      label='Tilt-series WITHOUT alignment')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.getTMSetOfTiltSeries.get():
            self._insertFunctionStep('assignTransformationMatricesStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def assignTransformationMatricesStep(self, tsObjId):
        outputAssignedTransformSetOfTiltSeries = self.getOutputAssignedTransformSetOfTiltSeries()

        setTMTS = self.setTMSetOfTiltSeries.get()[tsObjId]
        getTMTS = self.getTMSetOfTiltSeries.get()[tsObjId]
        tsId = getTMTS.getTsId()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(setTMTS)
        outputAssignedTransformSetOfTiltSeries.append(newTs)

        for tiltImageGetTM, tiltImageSetTM in zip(getTMTS, setTMTS):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImageSetTM, copyId=True)
            newTi.setLocation(tiltImageSetTM.getLocation())
            newTransform = self.updateTM(tiltImageGetTM.getTransform())
            newTi.setTransform(newTransform)
            newTs.append(newTi)

        newTs.setDim(setTMTS.getDim())
        newTs.write()

        outputAssignedTransformSetOfTiltSeries.update(newTs)
        outputAssignedTransformSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputAssignedTransformSetOfTiltSeries(self):
        if not hasattr(self, "outputAssignedTransformSetOfTiltSeries"):
            outputAssignedTransformSetOfTiltSeries = self._createSetOfTiltSeries(suffix='AssignedTransform')
            outputAssignedTransformSetOfTiltSeries.copyInfo(self.setTMSetOfTiltSeries.get())
            outputAssignedTransformSetOfTiltSeries.setDim(self.setTMSetOfTiltSeries.get().getDim())

            self._defineOutputs(outputAssignedTransformSetOfTiltSeries=outputAssignedTransformSetOfTiltSeries)
            self._defineSourceRelation(self.getTMSetOfTiltSeries, outputAssignedTransformSetOfTiltSeries)
        return self.outputAssignedTransformSetOfTiltSeries

    def getSamplingRatio(self):
        return self.setTMSetOfTiltSeries.get().getSamplingRate() / self.getTMSetOfTiltSeries.get().getSamplingRate()

    def updateTM(self, transform):
        """ Scale the transform matrix shifts. """
        matrix = transform.getMatrix()

        matrix[0][2] /= self.getSamplingRatio()
        matrix[1][2] /= self.getSamplingRatio()

        transform.setMatrix(matrix)

        return transform

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for tsGetTM, tsSetTM in zip(self.getTMSetOfTiltSeries.get(), self.setTMSetOfTiltSeries.get()):
            if not tsGetTM.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set do not have a "
                                    "transformation matrix assigned.")

            if tsGetTM.getSize() != tsSetTM.getSize():
                validateMsgs.append("Every input tilt-series and its target "
                                    "must have the same number of elements.")

        if self.getTMSetOfTiltSeries.get().getSize() != self.setTMSetOfTiltSeries.get().getSize():
            validateMsgs.append("Both input sets must have the same "
                                "number of tilt-series.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            summary.append("Input tilt-series: %d\nTransformation matrices assigned: %d\n"
                           % (self.getTMSetOfTiltSeries.get().getSize(),
                              self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            methods.append("The transformation matrix has been assigned to %d tilt-series from the input set.\n"
                           % (self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            methods.append("Outputs are not ready yet.")
        return methods
