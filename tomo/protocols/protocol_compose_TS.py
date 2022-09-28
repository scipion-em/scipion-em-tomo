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

from pwem.protocols.protocol_import.base import ProtImport
import pyworkflow as pw
from pyworkflow.protocol import params, STEPS_PARALLEL
import pyworkflow.protocol.constants as cons
from tomo.convert.mdoc import MDoc
import pwem.objects as emobj
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from pwem.emlib.image import ImageHandler
import time
import os
from glob import glob


class ProtComposeTS(ProtImport, ProtTomoBase):
    """ Compose in streaming a set of tilt series based on a sets of micrographs and mdoc files.
    """
    _devStatus = pw.BETA
    _label = 'Compose Tilt Serie'
    mdocFilesList = []
    mdocFilesRead = []

    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL # Defining that the protocol contain parallel steps
        self.newSteps = []
        self.TiltSeries = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      important=True,
                      label="Input micrographs",
                      help='Select the SetOfMicrographs to import')

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
        self.ih = ImageHandler()

    def _insertAllSteps(self):
        self._insertFunctionStep(self._initialize)
        self.CloseStep_ID = self._insertFunctionStep('closeSet',
                                                     prerequisites=[],
                                                     wait=True)
        self.newSteps.append(self.CloseStep_ID)


    def _stepsCheck(self):
        currentTime = time.time()
        self.debug('stepsCheck ' +
            str(int(int(currentTime) - self.time4NextTS_current)) + ' segs')
        listCurrent = self.findMdoc()
        listRemain = [x for x in listCurrent if x not in self.listMdocsRead]

        if int(currentTime - self.time4NextTS_current) \
                > int(self.time4NextTS.get()):
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

    def closeSet(self):
        print('Hello')
        self.TiltSeries.setStreamState(self.TiltSeries.STREAM_CLOSED)

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
        return 'closeSet'

    def findMdoc(self):
        """
        :return: return a sorted by date list of all mdoc files in the path
        """
        """ return a sorted by date list of all mdoc files in the path """
        fpath = self.filesPath.get()
        self.MDOC_DATA_SOURCE = glob(os.path.join(fpath, '*.mdoc'))
        self.MDOC_DATA_SOURCE.sort(key=os.path.getmtime)
        return self.MDOC_DATA_SOURCE

    def readMdoc(self, listRemains):
        """
        Main function to launch the match with the set of micrographs and
        launch the create the SetOfTiltSeries and each TiltSerie
        :param listRemains: list of mdoc files in the path
        """
        for file2Read in listRemains:
            mdocObj = MDoc(file2Read)
            validationError = mdocObj.read()
            if validationError:
                self.debug(validationError)
            else:
                self.info('mdoc file to read: {}'.format(file2Read))
                fileOrderAngleList = []
                for tiltMetadata in mdocObj.getTiltsMetadata():
                    fileOrderAngleList.append((
                        tiltMetadata.getAngleMovieFile(),           # Filename
                        '{:03d}'.format(tiltMetadata.getAcqOrder()),# Acquisition
                        tiltMetadata.getTiltAngle()))

                while time.time() - self.readDateFile(file2Read) < \
                        2 * self.time4NextTilt.get():
                    self.debug('waiting...')
                    time.sleep(self.time4NextTilt.get() / 2)
                if len(fileOrderAngleList) < 4:
                    self.error('Mdoc error. Less than 4 tilts in the serie')
                    break
                else:
                    self.matchTS(fileOrderAngleList)
                    self.createTS(mdocObj)

    def readDateFile(self, file):
        return os.path.getmtime(file)

    def matchTS(self, fileOrderAngleList):
        """
        Edit the self.listOfMics with the ones in the mdoc file
        :param fileOrderAngleList: for each tilt:
                filename, acquisitionOrder, Angle
        """
        self._loadInputList()
        self.info('Tilts in the mdoc file: {}\n'
                  'Micrographs abailables: {}'.format(
             len(fileOrderAngleList), len(self.listOfMics)))
        listMdocFiles = [os.path.basename(fp[0]) for fp in fileOrderAngleList]
        listMicsMatched = []
        for x, mic in enumerate(self.listOfMics):
            if mic.getMicName() in listMdocFiles:
                listMicsMatched.append(mic)
        self.listOfMics = listMicsMatched
        self.info('Micrographs matched for the mdoc file: {}'.format(
            len(self.listOfMics)))

    def _loadInputList(self):
        """ Load the input set of mics and create a list. """
        micFile = self.inputMicrographs.get().getFileName()
        self.debug("Loading input db: %s" % micFile)
        mic_Set = emobj.SetOfMicrographs(filename=micFile)
        mic_Set.loadAllProperties()
        self.listOfMics = [m.clone() for m in mic_Set]
        mic_Set.close()

    def createTS(self, mdocObj):
        """
        Create the SetOfTiltSeries and each TiltSerie
        :param mdocObj: mdoc object to manage
        """
        if self.TiltSeries == None:
            SOTS = self._createSetOfTiltSeries(suffix='Set')
            SOTS.setStreamState(SOTS.STREAM_OPEN)
            SOTS.enableAppend()
            self._defineOutputs(TiltSeries=SOTS) #generates self.TiltSeries
            self._defineSourceRelation(self.inputMicrographs, SOTS)
            self._store(SOTS)
        else:
            SOTS = self.TiltSeries
            SOTS.setStreamState(SOTS.STREAM_OPEN)
            SOTS.enableAppend()
            self._store(SOTS)

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

        fileOrderedAngleList = sorted(fileOrderAngleList,
                                      key=lambda angle: float(angle[2]))

        tsObj = tomoObj.TiltSeries()#alom,ejor tengo k append todas al inicio
        tsObj.setTsId(mdocObj.getTsId())
        tsObj.getAcquisition().setTiltAxisAngle(mdocObj.getTiltAxisAngle())

        SOTS.append(tsObj)

        counterTi = 0
        for f, to, ta in fileOrderedAngleList:
            to_a = int(to) - 1
            try:
                for mic in self.listOfMics:
                    if tsObj.getSamplingRate() == None:
                        tsObj.setSamplingRate(mic.getSamplingRate())
                    if SOTS.getSamplingRate() == None:
                        SOTS.setSamplingRate(mic.getSamplingRate())
                    if os.path.basename(f) in mic.getMicName():
                        ti = tomoObj.TiltImage(
                            location=mic.getFileName(),
                            _acqOrder=to_a,
                            _tiltAngle=ta)
                        ti.setTsId(counterTi)
                        ti.setIndex(counterTi)
                        ti.setSamplingRate(mic.getSamplingRate())
                        ti.setAcquisition(tsObj.getAcquisition().clone())
                        ti.getAcquisition().setDosePerFrame(incomingDoseList[to_a])
                        ti.getAcquisition().setAccumDose(accumulatedDoseList[to_a])

                        # Create stack and append
                        self.addTiltImage(ti.getFileName(), tsObj, mic.getMicName(), ti.getTiltAngle(),
                                          ti.getAcquisitionOrder(), ti.getAcquisition(),
                                          ti.getTsId(), ti.getSamplingRate(),
                                          counterTi, counterTi)
                        #tsObj.append(ti)
                        counterTi += 1
            except Exception as e:
                self.info(e)

        SOTS.updateDim()
        tsObj.write(properties=False)
        SOTS.update(tsObj)
        SOTS.write()
        #tsObj.close()

        self._store(SOTS)
        #self._store()
        #tsObj.clear()
        #SOTS.close()

    def addTiltImage(self, tiFile, tsObject, suffix, ti_ang, ti_Ord, acq,
                     tsIde, samplingRate, objId, index):
        """
        :param tiFile: aligned tilt image file
        :param tsObject: Tilt Series to which the new Ti Image will be added
        :param suffix: used for the location.
        :param ti_ang:
        :param ti_Ord:
        :param acq:
        :param tsIde: Tilt Series Movies object
        :param samplingRate:  current Tilt Series sampling rate
        :param objId: location of the Tilt Image which will be added
        :param index: position of the slice in the generated slack
        """
        ti = tomoObj.TiltImage(tiltAngle=ti_ang, tsId=tsIde, acquisitionOrder=ti_Ord)
        ti.setSamplingRate(samplingRate)
        ti.setIndex(index)
        ti.setAcquisition(acq)
        newLocation = (self._getExtraPath(str(tsIde) + '_' + suffix + '.mrcs'))
        self.ih.convert(inputObj=tiFile, outputObj= newLocation)
        ti.setLocation(newLocation)
        tsObject.append(ti)
        #pw.utils.cleanPath(tiFile)



    def _validate(self):
        errors = [] if len(self.inputMicrographs.get()) > 1 else \
            ["More than one Input micrographs is needed to run."]
        return errors

    def _summary(self):
        summary = []
        text = 'Set of Tilt Serie composed'
        summary.append(text)
        return summary