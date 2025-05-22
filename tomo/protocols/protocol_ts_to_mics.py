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

from pwem.objects import (SetOfMicrographs, Micrograph, Acquisition, 
                          String, Integer, Set)
from pwem.protocols import EMProtocol
from pwem import emlib
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage


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
        stepIds = []
        for ts in self.input.get():
            tsId = ts.getTsId()
            self._insertFunctionStep(self.createOutputStep, tsId, prerequisites=[])
        self._insertFunctionStep(self._closeOutputSet, prerequisites=stepIds)

    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self, tsId):
        with self._lock:
            inputTs = self._getInputTiltSeries(tsId)
            output = self._getOutputMicrographs()
            ih = emlib.image.ImageHandler()
            tiltImage: TiltImage = None
            for tiltImage in inputTs:
                if self.maxTiltAngle < 0 or abs(tiltImage.getTiltAngle()) < self.maxTiltAngle:
                    acqOrder = tiltImage.getAcquisitionOrder()
                    micName = '%s_%03d' % (tsId, acqOrder)
                    micrograph = Micrograph(location=self._getMicrographFilename(micName))
                    micrograph.setMicName(micName)
                    setattr(micrograph, TiltImage.TS_ID_FIELD, String(tsId))
                    setattr(micrograph, TiltImage.ACQ_ORDER_FIELD , Integer(acqOrder))
                    
                    ih.convert(tiltImage, micrograph)
                    output.append(micrograph)
            
            self._store()

    # --------------------------- UTILS functions ------------------------------
    def _getInputTiltSeries(self, tsId) -> TiltSeries:
        setOfTs: SetOfTiltSeries = self.input.get()
        return setOfTs.getTiltSeriesFromTsId(tsId)
    
    def _getMicrographFilename(self, micName) -> str:
        return self._getExtraPath('%s.mrc' % micName)

    def _getOutputMicrographs(self) -> SetOfMicrographs:
        result = getattr(self, ProtTsToMicsOutput.outputMicrographs.name, None)
        if result is None:
            setOfTs: SetOfTiltSeries = self.input.get()
            result = SetOfMicrographs.create(self._getPath())
            result.setStreamState(Set.STREAM_OPEN)
            result.setSamplingRate(setOfTs.getSamplingRate())

            # Acquisition
            micAcq = Acquisition()
            tomoAcq = setOfTs.getAcquisition()
            micAcq.setMagnification(tomoAcq.getMagnification())
            micAcq.setSphericalAberration(tomoAcq.getSphericalAberration())
            micAcq.setVoltage(tomoAcq.getVoltage())
            micAcq.setAmplitudeContrast(tomoAcq.getAmplitudeContrast())
            micAcq.setDoseInitial(tomoAcq.getDoseInitial())
            micAcq.setDosePerFrame(tomoAcq.getDosePerFrame())
            result.setAcquisition(micAcq)

            self._defineOutputs(**{ProtTsToMicsOutput.outputMicrographs.name: result})
            self._defineSourceRelation(self.input, result)
            
        return result