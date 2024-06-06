# **************************************************************************
# *
# * Authors:     Alberto García Mena (alberto.garcia@cnb.csic.es)
# *
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


from pwem.protocols.protocol_import.base import ProtImport
import pyworkflow as pw
from pyworkflow.protocol import params, STEPS_PARALLEL, ProtStreamingBase
import pyworkflow.protocol.constants as cons
from pyworkflow.object import Set
from tomo.convert.mdoc import MDoc
import pwem.objects as emobj
import pwem.objects as SetOfCTF
import tomo.objects as tomoObj
from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries


OUT_SCTF = "CTFSeries"


class ProtComposeCTF(ProtImport, ProtStreamingBase):
    """ Compose in streaming a set of tilt series based on a set of micrographs and mdoc files.
    Two time parameters are available for the streaming behaviour:
    Time to next tilt and time to next tilt series
    """
    _devStatus = pw.BETA
    _label = 'Compose CTF Series'
    _possibleOutputs = {OUT_SCTF: SetOfCTFTomoSeries}

    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('inputTS', params.PointerParam, allowsNull=False,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label="Input Tilt Series",
                      help='Select the SetOfTiltSeries to import')
        
        form.addParam('inputSetCTF', params.PointerParam, allowsNull=False,
                      pointerClass='SetOfCTF',
                      important=True,
                      label="Input CTF Series",
                      help='Select the SetOfCTF to import')
        
    def stepsGeneratorStep(self):
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        #while self.inputTS.isStreamOpen() or self.inputSetCTF.isStreamOpen():
        while True:
            self.info('True')
            self.TS = self.inputTS.get()
            self.CTFs = self.inputSetCTF.get()
            matchCTF = self._insertFunctionStep(self.matchCTF, prerequisites=[])
            outputs = self._insertFunctionStep(self.outputs, prerequisites=[matchCTF])
            if not self.TS.isStreamOpen() and not self.CTFs.isStreamOpen():
                self.info('Not more micrographs are expected, set closed')
                break
			    
    def matchCTF(self):
        self.info('Match')
        for TS in self.TS:
            for tilt in TS:
                micName = tilt.getMicName()
                for CTF in self.CTFs:
                    if micName == CTF.getMicrograph().getMicName():
                        self.info('match!: {}'.format(micName))
    
    def outputs(self):
        pass
    
    def _validate(self):
        pass


    def _summary(self):
        pass


