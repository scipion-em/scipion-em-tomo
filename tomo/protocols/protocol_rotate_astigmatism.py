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
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
import math


class ProtRotateAstigmatism(EMProtocol, ProtTomoBase):
    """
    Rotate the astigmatism of a set of ctf tilt-series estimation given a set of transformation matrices coming from
    a set of tilt series.
    """

    _label = 'Astigmatism rotation'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('getTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series from which transformation matrices will be obtained.',
                      label='Set of tilt-series from which get transform')

        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="input tilt-series CTF estimation",
                      important=True,
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation whose astigmatism estimation will be rotated.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.getTMSetOfTiltSeries.get():
            self._insertFunctionStep(self.rotateAstimatism, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def rotateAstimatism(self, tsObjId):
        getTMTS = self.getTMSetOfTiltSeries.get()[tsObjId]
        tsId = getTMTS.getTsId()

        match = False

        for inputCtfTomoSeries in self.inputSetOfCtfTomoSeries.get():
            if tsId == inputCtfTomoSeries.getTsId():
                match = True

                self.getOutputSetOfCTFTomoSeries()

                newCTFTomoSeries = tomoObj.CTFTomoSeries()
                newCTFTomoSeries.copyInfo(inputCtfTomoSeries)
                newCTFTomoSeries.setTiltSeries(getTMTS)
                newCTFTomoSeries.setTsId(tsId)
                newCTFTomoSeries.setObjId(tsObjId)

                # Check IMOD specific fields
                if hasattr(inputCtfTomoSeries, '_IMODDefocusFileFlag'):
                    newCTFTomoSeries.setIMODDefocusFileFlag(inputCtfTomoSeries.getIMODDefocusFileFlag())

                if hasattr(inputCtfTomoSeries, '_estimationsInRange'):
                    newCTFTomoSeries.setNumberOfEstimationsInRange(inputCtfTomoSeries.getNumberOfEstimationsInRange())

                self.outputSetOfCTFTomoSeries.append(newCTFTomoSeries)

                for index, (tiltImageGetTM, inputCtfTomo) in enumerate(zip(getTMTS, inputCtfTomoSeries)):
                    newCTFTomo = tomoObj.CTFTomo()
                    newCTFTomo.copyInfo(inputCtfTomo)

                    rotationAngle = self.calculateRotationAngleFromTM(tiltImageGetTM)

                    if newCTFTomo.hasAstigmatismInfoAsList():
                        defocusAngleList = pwobj.CsvList(pType=float)

                        for angle in inputCtfTomo.getDefocusAngleList().split(','):
                            defocusAngleList.append(pwobj.Float(round(float(angle) + rotationAngle, 2)))

                        newCTFTomo.setDefocusAngleList(defocusAngleList)

                        newCTFTomo.completeInfoFromList()

                    else:
                        newCTFTomo.setDefocusAngle(pwobj.Float(inputCtfTomo.getDefocusAngle() + rotationAngle))

                        newCTFTomo.standardize()

                    newCTFTomoSeries.append(newCTFTomo)

                newCTFTomoSeries.setIsDefocusUDeviationInRange(inputCtfTomoSeries.getIsDefocusUDeviationInRange())
                newCTFTomoSeries.setIsDefocusVDeviationInRange(inputCtfTomoSeries.getIsDefocusVDeviationInRange())

                if not (newCTFTomoSeries.getIsDefocusUDeviationInRange() and
                        newCTFTomoSeries.getIsDefocusVDeviationInRange()):
                    newCTFTomoSeries.setEnabled(False)

                newCTFTomoSeries.write(properties=False)

                self.outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
                self.outputSetOfCTFTomoSeries.write()

                self._store()

        if not match:
            raise Exception("There is no matching CtfTomoSeries for a tilt-series %s"
                            % tsId)

    def closeOutputSetsStep(self):
        self.getOutputSetOfCTFTomoSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def calculateRotationAngleFromTM(ti):
        """ This method calculates que tilt image rotation angle from its associated transformation matrix."""
        tm = ti.getTransform().getMatrix()
        cosRotationAngle = tm[0][0]
        sinRotationAngle = tm[1][0]
        rotationAngle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

        return rotationAngle

    def getOutputSetOfCTFTomoSeries(self):
        if hasattr(self, "outputSetOfCTFTomoSeries"):
            self.outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = tomoObj.SetOfCTFTomoSeries.create(self._getPath(),
                                                                         template='CTFmodels%s.sqlite')
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self.getTMSetOfTiltSeries.get())
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCTFTomoSeries=outputSetOfCTFTomoSeries)
        return self.outputSetOfCTFTomoSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for tsGetTM, inputCtfTomoSeries in zip(self.getTMSetOfTiltSeries.get(), self.inputSetOfCtfTomoSeries.get()):
            if not tsGetTM.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set of tilt-series does not have a "
                                    "transformation matrix assigned.")

            if tsGetTM.getSize() != inputCtfTomoSeries.getSize():
                validateMsgs.append("Some tilt-series from the input set of tilt-series and its target in the assign "
                                    "transformation set of tilt-series size's do not match. Every input tilt-series "
                                    "and its target must have the same number of elements")

        if self.getTMSetOfTiltSeries.get().getSize() != self.inputSetOfCtfTomoSeries.get().getSize():
            validateMsgs.append("Input sets of tilt-series and input set of CTF estimation differ in size. Both sets "
                                "must have the same number of elements.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            summary.append("Input pairs of Tilt-Series and CTF estimations: %d.\n"
                           "CTF estimations astigmatically rotated: %d.\n"
                           % (self.getTMSetOfTiltSeries.get().getSize(),
                              self.outputSetOfCTFTomoSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            methods.append("The astigmatism estimation has been rotated for %d CTF estimations.\n"
                           % (self.outputSetOfCTFTomoSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
