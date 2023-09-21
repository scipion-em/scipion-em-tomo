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
    _possibleOutputs = {"assignedTransformSetOfTiltSeries": tomoObj.SetOfTiltSeries, }

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
                      label='Tilt-series WITH alignment')

        form.addParam('setTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series on which transformation matrices will be assigned.',
                      label='Tilt-series WITHOUT alignment')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('assignTransformationMatricesStep')

    # --------------------------- STEPS functions ----------------------------
    def assignTransformationMatricesStep(self):
        self.getOutputAssignedTransformSetOfTiltSeries()

        getTMTSdict = self.getTSDict(self.getTMSetOfTiltSeries.get())
        setTMTSdict = self.getTSDict(self.setTMSetOfTiltSeries.get())

        for key in getTMTSdict:
            if key in setTMTSdict.keys():
                print(key)
                print(getTMTSdict[key])
                print(setTMTSdict[key])
                print(getTMTSdict[key].getSize())
                print(setTMTSdict[key].getSize())

                if getTMTSdict[key].getSize() != setTMTSdict[key].getSize():
                    self.info("Number of tilt-images in source and target set differ for tilt-series %s. Ignoring ts "
                              "(will not appear in output set)." % key)

                else:
                    newTs = tomoObj.TiltSeries(tsId=key)
                    newTs.copyInfo(setTMTSdict[key])
                    self.outputAssignedTransformSetOfTiltSeries.append(newTs)

                    for tiltImageGetTM, tiltImageSetTM in zip(getTMTSdict[key], setTMTSdict[key]):
                        newTi = tomoObj.TiltImage()
                        newTi.copyInfo(tiltImageSetTM, copyId=True)
                        newTi.setLocation(tiltImageSetTM.getLocation())
                        newTransform = self.updateTM(tiltImageGetTM.getTransform())
                        newTi.setTransform(newTransform)
                        newTs.append(newTi)

                    newTs.setDim(setTMTSdict[key].getDim())
                    newTs.write()

                    outputAssignedTransformSetOfTiltSeries.update(newTs)
                    outputAssignedTransformSetOfTiltSeries.write()
                    self._store()

            else:
                self.info("No matching tilt-series in target set for tilt-series %s. Ignoring ts (will not appear in "
                          "output set)." % key)

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def getTSDict(tsSet):
        tsDict = {}

        for ts in tsSet:
            tsId = ts.getTsId()
            tsDict[tsId] = ts

        return tsDict

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

        for tsGetTM in self.getTMSetOfTiltSeries.get():
            if not tsGetTM.getFirstItem().hasTransform():
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
