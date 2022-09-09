# **************************************************************************
# *
# * Authors:     Alberto García Mena (alberto.garcia@cnb.csic.es) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

from pwem.protocols.protocol_import.base import ProtImportFiles, ProtImport
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
import pyworkflow as pw
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params, Positive, STATUS_NEW, STEPS_PARALLEL
import pyworkflow.protocol.constants as cons
import time
import os
from glob import glob
import subprocess


class ProtImportTsStreaming(ProtImport):
    """ class for Tilt-Series import protocols in streaming.
    """
    _devStatus = pw.BETA
    mdocFilesList = []
    mdocFilesRead = []
    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL # Defining that the protocol contain parallel steps
        self.newSteps = []
        self.time4NextTilt_current = time.time()
        self.time4NextTS_current = time.time()
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Root directory of the tilt-series "
                           "(or movies) files.")

        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the path "
                           "should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addSection('Streaming')

        form.addParam('dataStreaming', params.BooleanParam, default=True,
                      label="Process data in streaming?",
                      help="Select this option if you want import data as it "
                           "is generated and process on the fly by next "
                           "protocols. In this case the protocol will "
                           "keep running to check new files and will "
                           "update the output Set, which can "
                           "be used right away by next steps.")

        form.addParam('time4NextTilt', params.IntParam, default=180,
                      condition='dataStreaming',
                      label="Time for next Tilt (secs)",
                      help="Delay until the next tilt is registered. After "
                           "timeout,\n if there is no new tilt, the tilt serie is considered as completed\n")
        form.addParam('time4NextTS', params.IntParam, default=1800,
                      condition='dataStreaming',
                      label="Time for next TiltSerie (secs)",
                      help="Interval of time (in seconds) after which, "
                           "if no new tilt serie is detected, the protocol will "
                           "end. When finished, the output Set will be "
                           "closed and no more data will be "
                           "added to it. \n"
                           "Note 1:  The default value is  high (30 min) to "
                           "avoid the protocol finishes during the acq of the "
                           "microscope. You can also stop it from right click "
                           "and press STOP_STREAMING.\n")
        form.addParallelSection(threads=3, mpi=1)



    def initializeParams(self):
        pass

    def _insertAllSteps(self):
        print('_insertAllSteps')
        self.newSteps.append(self._insertFunctionStep('initializeParams'))
        self.CloseStep_ID = self._insertFunctionStep('createOutputStep',
                                                     prerequisites=[],
                                                     wait=True)
        self.newSteps.append(self.CloseStep_ID)
        print('END')

        #wait to prevent finish the protocol

    def _stepsCheck(self):
        print('stepsCheck-------------------')
        currentTime = time.time()
        print(str(int(int(currentTime) - self.time4NextTS_current)) + ' segs')
        if int(currentTime - self.time4NextTS_current) > int(self.time4NextTS.get()):
            print('Timeout reached!!')
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)
            #self.CloseStep_ID.setStatus(cons.STATUS_NEW)
            self._steps[self.CloseStep_ID].setStatus(cons.STATUS_NEW)

        if self.findMdoc():
            newStepID = self._insertFunctionStep('readMdoc',prerequisites=[])
            self.newSteps.append(newStepID)
            self.updateSteps()

    def createOutputStep(self):
        pass

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overwritten in subclasses
        # (eg in Xmipp 'sortPSDStep')
        return 'createOutputStep'


    def findMdoc(self):
        fpath = self.filesPath.get()
        self.mdocFilesList = glob(os.path.join(fpath))
        if self.mdocFilesList == []: return False
        else:
            print(self.mdocFilesList)
            return True

    def readMdoc(self):
        print('readMdoc')

    def registerTilt(self):
        pass

    def sortTS(self):
        pass

    def stakTS(self):
        pass

    def _summary(self):
        summary = []

        return summary