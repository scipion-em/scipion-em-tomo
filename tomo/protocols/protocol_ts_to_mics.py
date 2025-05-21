# **************************************************************************
# *
# * Authors:     Federico de Isidro Gomez (fepe.ig@gmail.com) [1]
# *              Oier Lauzirika Zarrabeitia (oierlauzi@bizkaia.eu) [2]
# *
# * [1] Astex Therapeutics, UK
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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

from pwem.objects import SetOfMicrographs, Micrograph, Acquisition, String, Integer
from pwem.protocols import EMProtocol
from pwem import emlib
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam
from tomo.objects import SetOfTiltSeries


class ProtTsToMicsOutput(enum.Enum):
    outputMicrographs = SetOfMicrographs


class ProtTsToMics(EMProtocol):
    """ Turns tilt-series into set of micrographs to apply SPA picking methods."""
    _label = 'tilt-series to micrographs'
    _devStatus = BETA
    _possibleOutputs = ProtTsToMicsOutput

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam,
                      pointerClass=SetOfTiltSeries,
                      label='Tilt series',
                      important=True,
                      help='Select the tilt-series to be turned into micrographs')
        form.addParam('maxTiltAngle', FloatParam,
                      default=-1, label='Max tilt angle',
                      help='Maximum tilt angle of the tilt-series to be converted to micrographs. -1 to consider all images.')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        input = self.input.get()
        output = SetOfMicrographs.create(self._getPath())
        output.setSamplingRate(input.getSamplingRate())

        # Acquisition
        micAcq = Acquisition()
        tomoAcq = input.getAcquisition()
        micAcq.setMagnification(tomoAcq.getMagnification())
        micAcq.setSphericalAberration(tomoAcq.getSphericalAberration())
        micAcq.setVoltage(tomoAcq.getVoltage())
        micAcq.setAmplitudeContrast(tomoAcq.getAmplitudeContrast())
        micAcq.setDoseInitial(tomoAcq.getDoseInitial())
        micAcq.setDosePerFrame(tomoAcq.getDosePerFrame())

        output.setAcquisition(micAcq)

        ih = emlib.image.ImageHandler()
        for tiltSeries in input:
            tsId = tiltSeries.getTsId()
            for tiltImage in tiltSeries:
                if self.maxTiltAngle > 0 and abs(tiltImage.getTiltAngle()) < self.maxTiltAngle:
                    tiltId = tiltImage.getObjId()
                    micName = '%s_%03d' % (tsId, tiltId)
                    micrograph = Micrograph(location=self._getMicrographFilename(micName))
                    micrograph.setMicName(micName)
                    micrograph._tiltSeriesId = String(tsId)
                    micrograph._tiltImageId = Integer(tiltId)
                    
                    ih.convert(tiltImage, micrograph)
                    output.append(micrograph)
            
        self._defineOutputs(**{ProtTsToMicsOutput.outputMicrographs.name: output})
        self._defineSourceRelation(self.input, output)

    # --------------------------- UTILS functions ------------------------------
    def _getMicrographFilename(self, micName) -> str:
        return self._getExtraPath('%s.mrc' % micName)
