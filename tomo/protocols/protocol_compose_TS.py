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
from .protocol_ts_import import ProtImportTsBase
import pyworkflow as pw
from pyworkflow.protocol import params, Positive, STATUS_NEW, STEPS_PARALLEL
import pyworkflow.protocol.constants as cons
from tomo.convert.mdoc import MDoc
import pyworkflow.utils as pwutils
import pwem.objects as emobj
from pwem.objects import Transform
import tomo.objects as tomoObj
from collections import OrderedDict
from tomo.protocols import ProtTomoBase


import time
import os
from glob import glob


class ProtComposeTS(ProtImport, ProtTomoBase):
    """ class for Tilt-Series import protocols in streaming.
    """
    _devStatus = pw.BETA
    _label = 'Compose Tilt Serie'
    mdocFilesList = []
    mdocFilesRead = []
    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL # Defining that the protocol contain parallel steps
        self.newSteps = []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported movies.')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory ot the tiltSerie files",
                      help="Root directory of the tilt-series. "
                           "Will be search the *.mdoc file for each Tilt Serie")



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

    def _initialize(self):
        self.listMdocsRead = []
        self.time4NextTS_current = time.time()
        self._counterTS = 0
        self.SOTS = self._createSetOfTiltSeries(suffix='')
        self.summary = []
    def _insertAllSteps(self):
        self._insertFunctionStep(self._initialize)
        self.CloseStep_ID = self._insertFunctionStep('createOutputStep',
                                                     prerequisites=[],
                                                     wait=True)
        self.newSteps.append(self.CloseStep_ID)

    def _stepsCheck(self):
        currentTime = time.time()
        #print('stepsCheck ' + str(int(int(currentTime) - self.time4NextTS_current)) + ' segs')
        listCurrent = self.findMdoc()
        listRemain = [x for x in listCurrent if x not in self.listMdocsRead]

        if int(currentTime - self.time4NextTS_current) > int(self.time4NextTS.get()):
            print('Timeout reached!!')
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

        elif listRemain != []:
            self.listMdocsRead = listCurrent
            self.time4NextTS_current = time.time()
            newStepID = self._insertFunctionStep('readMdoc', listRemain,
                                        prerequisites=[], wait=False)
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
        self.MDOC_DATA_SOURCE = glob(os.path.join(fpath, '*.mdoc'))
        self.MDOC_DATA_SOURCE.sort(key=os.path.getmtime)
        return self.MDOC_DATA_SOURCE

    def readMdoc(self, listRemains):
        for file2Read in listRemains:
            self.debug('mdoc file to read: {}'.format(file2Read))
            mdocObj = MDoc(file2Read)
            validationError = mdocObj.read()
            if validationError:
                print(validationError)
            else:
                fileOrderAngleList = []
                for tiltMetadata in mdocObj.getTiltsMetadata():
                    fileOrderAngleList.append((
                        tiltMetadata.getAngleMovieFile(),             # Filename
                        '{:03d}'.format(tiltMetadata.getAcqOrder()),  # Acquisition
                        tiltMetadata.getTiltAngle()))
                while time.time() - self.readDateFile(file2Read) < \
                        2 * self.time4NextTilt.get():
                    self.debug('waiting...')
                    time.sleep(self.time4NextTilt.get() / 2)
                if len(fileOrderAngleList) < 4:
                    print('Mdoc error. Less than 2 tilts')
                    break
                else:
                    self.matchTS(fileOrderAngleList)
                    tsObj = self.createTS(mdocObj)
                    self.manageSetOfTS(tsObj)

    def readDateFile(self, file):
        return os.path.getmtime(file)

    def matchTS(self, fileOrderAngleList):
        self._loadInputList()
        print('listOfMdocMovies: {}, listOfFolderMovies: {}'.format(
             len(fileOrderAngleList), len(self.listOfMovies)))
        listMdocFiles = [os.path.basename(fp[0]) for fp in fileOrderAngleList]
        for x, movie in enumerate(self.listOfMovies):
            if os.path.basename(movie.getFileName()) not in listMdocFiles:
                self.debug('deleting: {}'.format(self.listOfMovies[x].getFileName()))
                del self.listOfMovies[x]

    def _loadInputList(self):
        """ Load the input set of movies and create a list. """
        moviesFile = self.inputMovies.get().getFileName()
        self.debug("Loading input db: %s" % moviesFile)
        movieSet = emobj.SetOfMovies(filename=moviesFile)
        movieSet.loadAllProperties()
        self.listOfMovies = [m.clone() for m in movieSet]
        movieSet.close()
        self.debug("Closed db.")

    def createTS(self, mdocObj):
        fileOrderAngleList = []
        accumulatedDoseList = []
        incomingDoseList = []

        for tiltMetadata in mdocObj.getTiltsMetadata():
            fileOrderAngleList.append((
                tiltMetadata.getAngleMovieFile(),  # Filename
                '{:03d}'.format(tiltMetadata.getAcqOrder()),  # Acquisition
                tiltMetadata.getTiltAngle()))
            accumulatedDoseList.append(tiltMetadata.getAccumDose())
            incomingDoseList.append(tiltMetadata.getIncomingDose())

        fileOrderedAngleList = sorted(fileOrderAngleList, key=lambda angle: float(angle[2]))
        tsObj = tomoObj.TiltSeries(tsId='TS_' + str(self._counterTS))
        origin = Transform()  # tilt serie
        tsObj.setOrigin(origin)
        tsObj.setAnglesCount(len(fileOrderedAngleList))

        counter = 0
        for f, to, ta in fileOrderedAngleList:
            for movie in self.listOfMovies:
                if f in movie.getFileName():
                    ti = tomoObj.TiltImage()
                    ti.clone(movie)
                    ti = tomoObj.TiltImage(
                        location=f, _acqOrder=to, _tiltAngle=ta)
                    ti.setAcquisition(tsObj.getAcquisition().clone())
                    ti.setLocation(counter, movie.getFileName())
                    dosePerFrame = incomingDoseList[counter]
                    accumDose = accumulatedDoseList[counter]
                    counter += 1
                    ti.getAcquisition().setDosePerFrame(dosePerFrame)
                    # Accumulated dose in current ti
                    ti.getAcquisition().setAccumDose(accumDose)
                    tsObj.append(ti)

        self._counterTS += 1
        return tsObj

    def manageSetOfTS(self, tsObj):
        print('manageSetOfTS')
        self._defineOutputs(SOTS=self.SOTS)
        self.SOTS.enableAppend()
        newTs = tsObj.clone()
        newTs.copyInfo(tsObj)
        self.SOTS.append(newTs)
        self.SOTS.write()
        self._store(self.SOTS)

    def _validate(self):
        errors = [] if len(self.inputMovies.get()) > 1 else \
            ["More than one Input movie is needed to run."]
        return errors

    def _summary(self):
        self.summary.append('Finished')
        return self.summary