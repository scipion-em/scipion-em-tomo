# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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

import os
import numpy as np
import itertools
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from pwem.emlib.image import ImageHandler


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtTomoBase):
    """
    Assign the transformation matrices from an input set of tilt-series to a target one.
    """

    _label = 'Tilt-series assign alignment'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series from which transformation matrices will be obtained.',
                      label='Input set of tilt-series')

        form.addParam('assignTransformSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices will be assigned.',
                      label='Assign transformation set of tilt-series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('assignTransformationMatricesStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def assignTransformationMatricesStep(self, tsObjId):
        outputAssignedTransformSetOfTiltSeries = self.getOutputAssignedTransformSetOfTiltSeries()

        ts = self.assignTransformSetOfTiltSeries.get()[tsObjId]
        tsTM = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputAssignedTransformSetOfTiltSeries.append(newTs)

        for tiltImage, tiltImageTM in zip(ts, tsTM):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(tiltImage.getLocation())
            newTi.setTransform(tiltImageTM.getTransform())
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))
        newTs.write()

        outputAssignedTransformSetOfTiltSeries.update(newTs)
        outputAssignedTransformSetOfTiltSeries.updateDim()
        outputAssignedTransformSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputAssignedTransformSetOfTiltSeries(self):
        if not hasattr(self, "outputAssignedTransformSetOfTiltSeries"):
            outputAssignedTransformSetOfTiltSeries = self._createSetOfTiltSeries(suffix='AssignedTransform')
            outputAssignedTransformSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputAssignedTransformSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

            self._defineOutputs(outputAssignedTransformSetOfTiltSeries=outputAssignedTransformSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputAssignedTransformSetOfTiltSeries)
        return self.outputAssignedTransformSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for ts in self.inputSetOfTiltSeries.get():
            if not ts.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set of tilt-series does not have a "
                                    "transformation matrix assigned.")

            if ts.getSize() != self.assignTransformSetOfTiltSeries.get()[ts.getObjId()].getSize():
                validateMsgs.append("Some tilt-series from the input set of tilt-series and its target in the assign "
                                    "transfomration set of tilt-series size's do not match. Every input tilt-series "
                                    "and its target must have the same number of elements")

        if self.inputSetOfTiltSeries.get().getSize() != self.assignTransformSetOfTiltSeries.get().getSize():
            validateMsgs.append("Both input sets of tilt-series size's do not match. Both sets must have the same "
                                "number of elements.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices assigned: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            methods.append("The transformation matrix has been assigned to %d Tilt-series from the input set.\n"
                           % (self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
