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

import numpy as np
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
    _possibleOutputs = {"assignedTransformSetOfTiltSeries": tomoObj.SetOfTiltSeries}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.outputAssignedTransformSetOfTiltSeries = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
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

        form.addParam('combineAlignment',
                      params.BooleanParam,
                      default=0,
                      label='Combine alignments',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='If this option is selected and the the tilt-series to assign the alignment has a previous '
                           'alignment, both precious and new alignment are combined.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        if self.combineAlignment.get():
            self._insertFunctionStep(self.combineTransformationMatricesStep)
        else:
            self._insertFunctionStep(self.assignTransformationMatricesStep)

    # --------------------------- STEPS functions ----------------------------
    def combineTransformationMatricesStep(self):
        self.info("Combining alignments mode started...")

        for getTS in self.getTMSetOfTiltSeries.get():
            tsId = getTS.getTsId()

            try:
                setTS = self.setTMSetOfTiltSeries.get()[{'_tsId': tsId}]

                self.getOutputAssignedTransformSetOfTiltSeries()

                if getTS.getSize() != setTS.getSize():
                    self.info(
                        "Number of tilt-images in source and target set differ for tilt-series %s. Ignoring ts "
                        "(will not appear in output set)." % tsId)

                else:
                    newTs = tomoObj.TiltSeries(tsId=tsId)
                    newTs.copyInfo(setTS)
                    self.outputAssignedTransformSetOfTiltSeries.append(newTs)

                    for tiltImageGetTM, tiltImageSetTM in zip(getTS, setTS):
                        newTi = tomoObj.TiltImage()
                        newTi.copyInfo(tiltImageSetTM, copyId=True)
                        newTi.setLocation(tiltImageSetTM.getLocation())

                        if tiltImageSetTM.hasTransform():
                            previousTransform = tiltImageSetTM.getTransform().getMatrix()
                            newTransform = self.updateTM(tiltImageGetTM.getTransform())
                            previousTransformArray = np.array(previousTransform)
                            newTransformArray = np.array(newTransform.getMatrix())

                            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
                            newTransform.setMatrix(outputTransformMatrix)

                        else:
                            newTransform = self.updateTM(tiltImageGetTM.getTransform())

                        newTi.setTransform(newTransform)
                        newTs.append(newTi)

                    newTs.setDim(setTS.getDim())
                    newTs.write()

                    self.outputAssignedTransformSetOfTiltSeries.update(newTs)
                    self.outputAssignedTransformSetOfTiltSeries.write()
                    self._store()

            except:
                self.info("Tilt-series %s not found in target set. Ignoring ts (will not appear in output set)." % tsId)

    def assignTransformationMatricesStep(self):
        self.info("Assigning alignments mode started...")

        for getTS in self.getTMSetOfTiltSeries.get():
            tsId = getTS.getTsId()

            try:
                setTS = self.setTMSetOfTiltSeries.get()[{'_tsId': tsId}]

                self.getOutputAssignedTransformSetOfTiltSeries()

                if getTS.getSize() != setTS.getSize():
                    self.info(
                        "Number of tilt-images in source and target set differ for tilt-series %s. Ignoring ts "
                        "(will not appear in output set)." % tsId)

                else:
                    newTs = tomoObj.TiltSeries(tsId=tsId)
                    newTs.copyInfo(setTS)
                    self.outputAssignedTransformSetOfTiltSeries.append(newTs)

                    for tiltImageGetTM, tiltImageSetTM in zip(getTS, setTS):
                        newTi = tomoObj.TiltImage()
                        newTi.copyInfo(tiltImageSetTM, copyId=True)
                        newTi.setLocation(tiltImageSetTM.getLocation())

                        newTransform = self.updateTM(tiltImageGetTM.getTransform())
                        newTi.setTransform(newTransform)

                        newTs.append(newTi)

                    newTs.setDim(setTS.getDim())
                    newTs.write()

                    self.outputAssignedTransformSetOfTiltSeries.update(newTs)
                    self.outputAssignedTransformSetOfTiltSeries.write()
                    self._store()

            except:
                self.info("Tilt-series %s not found in target set. Ignoring ts (will not appear in output set)." % tsId)

    # --------------------------- UTILS functions ----------------------------
    def getOutputAssignedTransformSetOfTiltSeries(self):
        if self.outputAssignedTransformSetOfTiltSeries:
            self.outputAssignedTransformSetOfTiltSeries.enableAppend()
        else:
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

        sr = self.getSamplingRatio()

        matrix[0][2] /= sr
        matrix[1][2] /= sr

        transform.setMatrix(matrix)

        return transform

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for tsGetTM in self.getTMSetOfTiltSeries.get():
            if not tsGetTM.hasAlignment():
                validateMsgs.append("Tilt-series %s from the input set do not have a "
                                    "transformation matrix assigned." % tsGetTM.getTsId())
                break

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
