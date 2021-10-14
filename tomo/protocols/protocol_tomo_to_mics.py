# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
import enum
from os.path import basename

from pwem.emlib.image import ImageHandler
from pwem.objects import SetOfMicrographs, Micrograph
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam
from pwem.protocols import EMProtocol
from pyworkflow.utils import replaceExt

import tomo.objects as tomoObjs

class ProtTomoToMicsOutput(enum.Enum):
    outputMicrographs = SetOfMicrographs()

class ProtTomoToMics(EMProtocol):
    """ Turns tomograms into set of micrographs to apply SPA picking methods."""
    _label = 'tomograms to micrographs'
    _devStatus = BETA
    _possibleOutputs = ProtTomoToMicsOutput

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Select the tomograms to be turned into micrographs')
        form.addParam('slicesGap', IntParam, label="Slices gap", default=10,
                      help='Number of slices to skip when turning tomogram slices into micrographs.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        input = self.input.get()

        output = SetOfMicrographs.create(self._getPath())
        # output.copyInfo(input)
        output.setSamplingRate(input.getSamplingRate())

        # For each tomogram
        for tomo in input:
            self.appendMicsFromTomogram(output, tomo)

        self._defineOutputs(**{ProtTomoToMicsOutput.outputMicrographs.name: output})
        self._defineSourceRelation(self.input, output)

    # --------------------------- UTILS functions --------------------------------------------
    @staticmethod
    def tomoSliceToMicName(tomo, slice):
        return "S%s_%s" % (slice, basename(tomo.getFileName()))

    def appendMicsFromTomogram(self, output, tomo):

        self.info("Creating micrographs for %s" % tomo.getFileName())
        # Load the tomogram
        ih = ImageHandler()
        img = ih.read(tomo.getFileName())
        data = img.getData()

        # For each slice
        for index in range(0, len(data),self.slicesGap.get()):
            self.debug("Creating micrograph for slice %s" % index)
            outputMicName = self._getExtraPath(self.tomoSliceToMicName(tomo, index))
            outputMicName = replaceExt(outputMicName, "mrc")
            slice = data[index]
            micImg = ImageHandler()
            micImg._img.setData(slice)
            micImg.write(micImg._img, outputMicName)

            # Create the micrograph metadata object
            newMic = Micrograph()
            newMic.setFileName(outputMicName)
            newMic.setMicName()
            newMic.setSamplingRate(tomo.getSamplingRate())

            # Append it
            output.append(newMic)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validateMsgs = []
        return validateMsgs
