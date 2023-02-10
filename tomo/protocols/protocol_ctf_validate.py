# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow import BETA
import pyworkflow.object as pwobj

from tomo.objects import SetOfCTFTomoSeries

logger = logging.getLogger(__name__)

OUTPUT_BAD_CTF_SERIE = "badCTFTomoSeries"
OUTPUT_GOOD_CTF_SERIE = "goodCTFTomoSeries"


class ProtCTFTomoSeriesValidate(EMProtocol):
    """
     Validate a set of CTF tomo series and separate into two sets (good and
     bad tomo series )
    """
    _label = 'ctf validate'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('inputCtfTomoSeries', params.PointerParam, important=True,
                      pointerClass='SetOfCTFTomoSeries',
                      label='Input ctf tomo series')

        form.addSection(label='CTF validation')
        form.addParam('validationType', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Global", "Per tilt"],
                      label="Validation type",
                      default=0,
                      help="Global mode: the series with at least "
                           "one image that does not satisfy the criteria are "
                           "rejected \n"
                           "Per tilt mode: the series containing a certain "
                           "number of images that does not satisfy the "
                           "criteria are rejected")
        form.addParam('numberImages', params.IntParam,
                      label='Number of images to rejected',
                      condition='validationType==1', default=None,
                      allowsNull=True,
                      help="Number of images taking into account to rejected a "
                           "ctf series")
        self.addCriteriaParams(form)

    def addCriteriaParams(self, form):
        form.addParam('defocusCriteria', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Yes", "No"],
                      label="Defocus tolerance",
                      default=1, help="Validate the defocus deviation taking "
                                      "into account a threshold(tolerance) "
                                      "respect to a defocus expected value.")
        form.addParam('defocusValue', params.FloatParam,
                      label='Expected value (Å)',
                      condition='defocusCriteria==0', default=None,
                      allowsNull=True,
                      help="Defocus expected value in Å")
        form.addParam('defocusTolerance', params.FloatParam,
                      label='Tolerance value (Å)',
                      condition='defocusCriteria==0', default=None,
                      allowsNull=True,
                      help="Defocus tolerance value in Å")

        form.addParam('astigmatismCriteria', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Yes", "No"],
                      label="Astigmatism",
                      default=1, help="Validate the astigmatism taking into "
                                      "account a tolerance value.")
        form.addParam('astigmatismTolerance', params.FloatParam,
                      label='Tolerance value',
                      condition='astigmatismCriteria==0', default=1.1,
                      help="Astigmatism tolerance value")

        form.addParam('resolutionCriteria', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Yes", "No"],
                      label="Resolution",
                      default=1, help="Validate the resolution taking into "
                                      "account a expected resolution.")
        form.addParam('resolutionTolerance', params.FloatParam,
                      label='Expected value',
                      allowsNull=True,
                      condition='resolutionCriteria==0', default=None,
                      help="Expected resolution value")

    def _insertAllSteps(self):
        self._insertFunctionStep(self.ctfValidateStep)
        self._insertFunctionStep(self.createOutputStep)

    # ----------------------------STEPS --------------------------------------
    def ctfValidateStep(self):
        """
        Validate all ctf tomo series and separate into two sets(good and bad
        following the selected criteria)
        """

        self.goodCTFTomoSeries = None
        self.badCTFTomoSeries = None
        self.validateDict = {}
        imagesToReject = 1  # Case of global validation type
        ctfSeries = self.inputCtfTomoSeries.get()

        if self.validationType.get() == 1:
            imagesToReject = self.numberImages.get()

        self.printTable(imagesToReject)
        printStr = ""
        for ctfSerie in ctfSeries:
            newCTFTomoSeries = ctfSerie.clone()
            ctfEstItems = []
            faileDefocusCriteria = 0
            failedAstigmatismCriteria = 0
            failedResolutionCriteria = 0
            for item in ctfSerie.iterItems():
                ctfEstItem = item.clone()
                ctfEstItems.append(ctfEstItem)

            # Defocus angles criteria
            if self.defocusCriteria.get() == 0:
                faileDefocusCriteria = self._validateCtfDefocusDeviation(ctfEstItems,
                                                                         self.defocusValue.get(),
                                                                         self.defocusTolerance.get())

            # Astigmatism criteria
            if self.astigmatismCriteria.get() == 0:
                failedAstigmatismCriteria = self._validateAstigmatism(ctfEstItems,
                                                                      self.astigmatismTolerance.get())

            # Resolution criteria
            if self.resolutionCriteria.get() == 0:
                # Some estimation methods don't calculate the resolution
                if ctfEstItems[0].getResolution() is not None:
                    failedResolutionCriteria = self._validateResolution(ctfEstItems,
                                                                        self.resolutionTolerance.get())

            validate = self.validateCtf(ctfEstItems,
                                        imagesToReject=imagesToReject)

            result = "Pass" if validate else "Rejected"
            printStr += "{0:25} {1:30} {2:30} {3:30} {4:^50}".format(ctfSerie.getTsId(),
                                                                     faileDefocusCriteria,
                                                                     failedAstigmatismCriteria,
                                                                     failedResolutionCriteria,
                                                                     result)
            printStr += "\n"

            if validate:
                newCTFTomoSeries.setEnabled(True)
                output = self.getOutputSetOfCTFTomoSeries(OUTPUT_GOOD_CTF_SERIE)
                output.append(newCTFTomoSeries)
            else:
                newCTFTomoSeries.setEnabled(False)
                output = self.getOutputSetOfCTFTomoSeries(OUTPUT_BAD_CTF_SERIE)
                output.append(newCTFTomoSeries)

            for ctfItem in ctfEstItems:
                newCTFTomoSeries.append(ctfItem)

        logger.info(printStr)

    def createOutputStep(self):
        if self.goodCTFTomoSeries is not None:
            self.goodCTFTomoSeries.setStreamState(pwobj.Set.STREAM_CLOSED)
            self.goodCTFTomoSeries.write()
            self._store()
        if self.badCTFTomoSeries is not None:
            self.badCTFTomoSeries.setStreamState(pwobj.Set.STREAM_CLOSED)
            self.badCTFTomoSeries.write()
            self._store()

    # --------------------------UTILS methods ----------------------------------
    def _validateCtfDefocusDeviation(self, ctfEstItems, defocusValue, tolerance):
        """
        Validate the set of ctf tomo series taking into account de defocus angle
        deviation
        """
        failedCTF = self.validateDefocusUDeviation(ctfEstItems, defocusValue,
                                                   defocusUTolerance=tolerance)
        failedCTF += self.validateDefocusVDeviation(ctfEstItems, defocusValue,
                                                    defocusVTolerance=tolerance)
        return failedCTF

    def _validateAstigmatism(self, ctfEstItems, astigmatismTolerance):
        """
         Validate the set of ctf tomo series taking into account de astigmatism
         threshold:
         the _objEnable property is set as True if the astigmatism is in range or False
         in other case, in each ctfTomoItems
         """
        failedCTF = 0
        for ctfTomo in ctfEstItems:
            isAstigmatismInRange = True if ctfTomo.getDefocusAngle() < astigmatismTolerance else False
            if not isAstigmatismInRange:
                failedCTF += 1
                ctfTomo.setEnabled(isAstigmatismInRange)
        return failedCTF

    def _validateResolution(self, ctfEstItems, resolutionTolerance):
        failedCTF = 0
        for ctfTomo in ctfEstItems:
            isAstigmatismInRange = True if ctfTomo.getResolution() < resolutionTolerance else False
            if not isAstigmatismInRange:
                failedCTF += 1
                ctfTomo.setEnabled(isAstigmatismInRange)
        return failedCTF

    @staticmethod
    def validateCtf(ctfItemsList, imagesToReject=1):
        """
         Validate the ctftomo serie.
         Return the number of images marked as disable
        """
        count = 0
        for ctfItem in ctfItemsList:
            if not ctfItem.isEnabled():
                count += 1
            if count == imagesToReject:
                return False
        return True

    def validateDefocusVDeviation(self, ctfTomoItems, defocusVValue, defocusVTolerance=20):
        """
          Set _objEnable property as True if the deviation is in range or False
          in other case, in each ctfTomoItems
        """
        failedCTF = 0
        for ctfTomo in ctfTomoItems:
            defocusVdeviation = abs(ctfTomo.getDefocusV() - defocusVValue)
            isDefocusVDeviationInRange = True if defocusVdeviation < defocusVTolerance else False
            if not isDefocusVDeviationInRange and ctfTomo.isEnabled():
                failedCTF += 1
                ctfTomo.setEnabled(isDefocusVDeviationInRange)
        return failedCTF

    def validateDefocusUDeviation(self, ctfTomoItems, defocusUValue, defocusUTolerance=20):
        """
        Set _objEnable property as True if the deviation is in range or False
        in other case, in each ctfTomoItems
        """
        failedCTF = 0
        for ctfTomo in ctfTomoItems:
            defocusUdeviation = abs(ctfTomo.getDefocusU() - defocusUValue)
            isDefocusUDeviationInRange = True if defocusUdeviation < defocusUTolerance else False
            if not isDefocusUDeviationInRange and ctfTomo.isEnabled():
                failedCTF += 1
                ctfTomo.setEnabled(isDefocusUDeviationInRange)
        return failedCTF

    def _getSetOfTiltSeries(self, pointer=False):
        return self.inputCtfTomoSeries.get().getSetOfTiltSeries(pointer=pointer)

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 prefix=outputSetName)
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self._getSetOfTiltSeries())
            outputSetOfCTFTomoSeries.setStreamState(pwobj.Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})

        return outputSetOfCTFTomoSeries

    def allowsDelete(self, obj):
        return True

    def _validate(self):
        return []

    def _summary(self):
        summary = []
        if hasattr(self, 'goodCTFTomoSeries'):
            summary.append("Number of good ctf series: %d." % (self.goodCTFTomoSeries.getSize()))
        if hasattr(self, 'badCTFTomoSeries'):
            summary.append("Number of bad ctf series: %d." % (self.badCTFTomoSeries.getSize()))

        return summary

    def printTable(self, imagesToReject):
        printStr = "\n-------------------------------------------------------\n"
        printStr += "Number of failed images to reject: %s\n" % imagesToReject
        printStr += '--------------------------------------------------------\n'
        printStr += "{0:35} {1:25} {2:25} {3:25} {4:20}".format('CTFTomoSerie',
                                                                'Defocus tolerance',
                                                                'Astigmatism',
                                                                'Resolution',
                                                                'Validate').ljust(5, ' ')
        printStr += "\n"
        printStr += '------------------------' * 7

        logger.info(printStr)
