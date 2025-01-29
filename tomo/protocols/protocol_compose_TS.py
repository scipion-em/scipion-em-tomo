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

import time
import os
from glob import glob
import re

from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from pwem.protocols.protocol_import.base import ProtImport
import pyworkflow as pw
from pyworkflow.protocol import params, STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.object import Set
from tomo.convert.mdoc import MDoc
import pwem.objects as emobj
import tomo.objects as tomoObj
from tomo.objects import SetOfTiltSeries
from pwem.objects.data import Transform
from tomo.protocols import ProtTomoBase
from pwem.emlib.image import ImageHandler

OUT_STS = "TiltSeries"


class ProtComposeTS(ProtImport, ProtTomoBase, ProtStreamingBase):
    """ Compose in streaming a set of tilt series based on a set of micrographs and mdoc files.
    One time parameters are available for the streaming behaviour: Time to next tilt
    """
    _devStatus = pw.BETA
    _label = 'Compose Tilt Series'
    _possibleOutputs = {OUT_STS: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL  # Defining that the protocol contain parallel steps

    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        self.MDOC_DATA_SOURCE = None
        self.TiltSeries = None
        self.waitingMdoc = True
        self.time4NextMic = 10
        self.time4NextTS_current = time.time()
        self.timeNextLoop = 30

    # -------------------------- DEFINES AND STEPS -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('inputMicrographs', params.PointerParam, allowsNull=False,
                      pointerClass='SetOfMicrographs',
                      important=True,
                      label="Input micrographs",
                      help='Select the SetOfMicrographs to import')

        form.addParam('filesPath', params.PathParam,
                      label="Path with the *.mdoc files for each tilt series",
                      help="Root directory of the tilt-series. "
                           "Use of * will work for multiple characters or ? for a single one. Also [ ] can specify ranges.")
        form.addParam('mdocPattern', params.PathParam,
                      label="Mdoc pattern",
                      default="*.mdoc",
                      help="Pattern that should match for mdoc files."
                           "Use of * will work for multiple characters or ? for a single one. Also [ ] can specify ranges.")
        form.addParam('excludedWords', params.StringParam,
                      label="Exclusion words",
                      default="",
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Space separated words that will be used to exclude mdoc files that could be listed with the above parameters.")

        form.addParam('isTomo5', params.BooleanParam, default=False,
                      label="Tomography 5 mdoc?",
                        help = "If these mdocs were generated by the Tomography 5 software, check this box to ensure that "
                        "the tilt axis angle is converted properly: -1 * TiltAxisAngle - 90")

        form.addParam('mdoc_bug_Correction', params.BooleanParam, default=False,
                      label="mdoc bug Correction",
                      help="Setting True, the mdoc generated by SerialEM "
                           "will be read considering the bug: filepath formatting and ensures tiltA is rounded and consistent.")

        form.addSection('Streaming')

        form.addParam('dataStreaming', params.BooleanParam, default=True,
                      label="Process data in streaming?",
                      help="Select this option if you want import data as it "
                           "is generated and process on the fly by next "
                           "protocols. In this case the protocol will "
                           "keep running to check new files and will "
                           "update the output Set, which can "
                           "be used right away by next steps.")

        form.addParam('time4NextTilt', params.StringParam, default="3m",
                      condition='dataStreaming',
                      label="Time for next Tilt",
                      help="Delay until the next tilt is "
                           "registered in the mdoc file. After "
                           "timeout, the mdoc file is not updated, the tilt series "
                           "is considered as proccessed."
                           "Minimum time recommended 20 secs (20s). For PACEtomo propose, please increase this time acording your acquisition."
                           "A correct format is an integer number in "
                           "seconds or the following syntax: {days}d {hours}h "
                           "{minutes}m {seconds}s separated by spaces "
                           "e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25")

        form.addParallelSection(threads=3, mpi=1)

    def _initialize(self):
        self.ih = ImageHandler()


    def isMdocBanned(self, mdoc):

        for bannedWord in self.excludedWords.getListFromValues(caster=str):
            if bannedWord in mdoc:
                self.info("mdoc %s contains the exclusion word %s. Skipping it." % (mdoc, bannedWord))
                return True
        return False

    def stepsGeneratorStep(self):
        """
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met
        call the self._insertFunctionStep method.
        """
        self._initialize()
        list_reading = []

        whileRunning = True
        while whileRunning:
            list_current = self.findMdocs()
            self.len_mics_input_1, _ = self._loadInputList()
            for mdocFile in list_current:
                # Exclusion
                if self.isMdocBanned(mdocFile):
                    continue
                if mdocFile not in list_reading:
                    self._insertFunctionStep('readMdoc', mdocFile, prerequisites=[], wait=False, needsGPU=False)
                    list_reading.append(mdocFile)

            if self.inputMicrographs.get().isStreamOpen() == False: #TODO maybe update the set
                self.info('The set of micrographs is closed')
                whileRunning = False

            time.sleep(self.timeNextLoop)

        self.debug('Happy EMProcessing')
        # self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        # self.closeSet()


    # def closeSet(self):
    #     self.info('Closing sets')
    #     for _, output_set in self.iterOutputAttributes():
    #         output_set.setStreamState(Set.STREAM_CLOSED)
	#
    #     self._store()


    # -------------------------- MAIN FUNCTIONS -----------------------
    def findMdocs(self):
        """
        :return: return a sorted by date list of all mdoc files in the path
        """
        fpath = self.filesPath.get()
        self.MDOC_DATA_SOURCE = glob(os.path.join(fpath, self.mdocPattern.get()))
        self.MDOC_DATA_SOURCE.sort(key=os.path.getmtime)
        return self.MDOC_DATA_SOURCE

    def readMdoc(self, file2read):
        """
        Main function to launch the match with the set of micrographs and
        launch the creation of the SetOfTiltSeries and each tilt series

        :param file2read: mdoc file in the path

        """
        self.info('Reading mdoc file: {}'.format(file2read))
        # checking time after last mdoc file update to consider it closed
        time4NextTilt = self.time4NextTilt.toSeconds()
        while time.time() - self.readDateFile(file2read) < time4NextTilt:
            self.debug('Waiting next tilt...)')
            time.sleep(time4NextTilt / 2)

        statusMdoc, mdoc_order_angle_list = self.readingMdocTiltInfo(file2read)
        self.info(f'mdoc file {os.path.basename(file2read)} with {len(mdoc_order_angle_list)} tilts considered closed')

        if statusMdoc:
            if len(mdoc_order_angle_list) < 3:
                self.info('Mdoc error. Less than 3 tilts on the series')
            elif self.matchTS(mdoc_order_angle_list, file2read):
                self.createTS(self.mdoc_obj)
                self.info("Tilt serie ({} tilts) composed from mdoc file: {}\n".
                          format(len(mdoc_order_angle_list), os.path.basename(file2read)))
        else:
            self.info('Mdoc file did not pass the format validation')


    def readingMdocTiltInfo(self, file2read):
        """
        :param file2read: mdoc file to read
        :return: Bool: if the validation of the mdoc goes good or bad
                 mdoc_order_angle_list: list with info for each tilt
                    file, acquisition order and tilt Angle
        """
        mdoc_order_angle_list = []
        self.mdoc_obj = MDoc(file2read)
        validation_error = self.mdoc_obj.read(ignoreFilesValidation=True)
        if validation_error:
            self.debug(validation_error)
            return False, mdoc_order_angle_list
        for tilt_metadata in self.mdoc_obj.getTiltsMetadata():
            filepath = tilt_metadata.getAngleMovieFile()
            tiltA = tilt_metadata.getTiltAngle()
            if self.mdoc_bug_Correction.get():
                filepath, tiltA = self.fixingMdocBug(filepath, tiltA)

            mdoc_order_angle_list.append((filepath,
                                          '{:03d}'.format(tilt_metadata.getAcqOrder()), tiltA))
        return True, mdoc_order_angle_list


    @staticmethod
    def fixingMdocBug(filepath, tiltA):
        idx = filepath.find(']_')
        filepath = filepath[:idx + 2] + filepath[idx + 2].upper() + filepath[idx + 3:]
        if float(tiltA) - round(float(tiltA), 0) != 0:
            filepath = filepath.replace(str(tiltA), str(round(float(tiltA))) + '.00')
            tiltA = str(round(float(tiltA))) + '.00'
        return filepath, tiltA

    @staticmethod
    def readDateFile(file):
        return os.path.getmtime(file)

    def matchTS(self, mdoc_order_angle_list, file2read):
        """
        Edit the self.listOfMics with the ones in the mdoc file

        :param mdoc_order_angle_list: for each tilt:
                filename, acquisitionOrder, Angle

        :param file2read: mdoc file to read

        """
        streameState = True
        while streameState:
            streameState = self.inputMicrographs.get().isStreamOpen()
            #Each self.timeNextLoop secs the self.listOfMics is updated (in stepsGeneratorStep)
            self.info(f'Tilts on the mdoc file: {len(mdoc_order_angle_list)}\n'
                      f'Micrographs available: {len(self.listOfMics)}')

            #MATCH
            list_mdoc_files = [os.path.splitext(os.path.basename(fp[0]))[0] for fp in mdoc_order_angle_list]
            list_mics_matched = []
            for x, mic in enumerate(self.listOfMics):
                if os.path.splitext(mic.getMicName())[0] in list_mdoc_files:
                    list_mics_matched.append(mic)

            if len(list_mics_matched) < len(mdoc_order_angle_list):
                    self.info(f"{len(self.listOfMics) - len(mdoc_order_angle_list)} micrographs are not abailable to compose the Tilt Serie"
                              "The Tilt serie will not be generated because not all mirographs are available."
                              f"The mdoc file {file2read} will not provide a TiltSerie")
                    time.sleep(self.timeNextLoop) #time until next check is run.
            else:
                self.info(f'Micrographs matched for the mdoc file: {len(list_mics_matched)}')
                return True

    def _loadInputList(self):
        """ Load the input set of mics and create a list. """
        mic_file = self.inputMicrographs.get().getFileName()
        self.debug("Loading input db: %s" % mic_file)
        mic_set = emobj.SetOfMicrographs(filename=mic_file)
        mic_set.loadAllProperties()
        self.listOfMics = [m.clone() for m in mic_set]
        mic_set.close()
        return len(self.listOfMics), self.inputMicrographs.get().isStreamOpen()

    def createTS(self, mdoc_obj):
        """
        Create the SetOfTiltSeries and each tilt series
        :param mdoc_obj: mdoc object to manage
        """

        self.info('Tilt series {} being composed...'.format(mdoc_obj.getTsId()))
        time.sleep(10)
        if self.TiltSeries is None:
            SOTS = self._createSetOfTiltSeries(suffix='')
            SOTS.setStreamState(SOTS.STREAM_OPEN)
            SOTS.enableAppend()
            self._defineOutputs(TiltSeries=SOTS)
            self._defineSourceRelation(self.inputMicrographs, SOTS)
            self._store(SOTS)
        else:
            SOTS = self.TiltSeries
            SOTS.setStreamState(SOTS.STREAM_OPEN)
            SOTS.enableAppend()
            self._store(SOTS)

        file_order_angle_list = []
        accumulated_dose_list = []
        incoming_dose_list = []
        for tilt_metadata in mdoc_obj.getTiltsMetadata():
            filepath = tilt_metadata.getAngleMovieFile()
            tiltAngle = tilt_metadata.getTiltAngle()
            if self.mdoc_bug_Correction.get():
                filepath, tiltAngle = self.fixingMdocBug(filepath, tiltAngle)

            file_order_angle_list.append((filepath,  # Filename
                                          '{:03d}'.format(tilt_metadata.getAcqOrder()),  # Acquisition
                                          tiltAngle))
            accumulated_dose_list.append(tilt_metadata.getAccumDose())
            incoming_dose_list.append(tilt_metadata.getIncomingDose())

        file_ordered_angle_list = sorted(file_order_angle_list,
                                         key=lambda angle: float(angle[2]))
        # Tilt series object
        ts_obj = tomoObj.TiltSeries()
        ts_obj.setTsId(mdoc_obj.getTsId())
        acq = ts_obj.getAcquisition()
        acq.setVoltage(mdoc_obj.getVoltage())
        acq.setMagnification(mdoc_obj.getMagnification())
        acq.setSphericalAberration(self.listOfMics[0].getAcquisition().getSphericalAberration())
        acq.setAmplitudeContrast(self.listOfMics[0].getAcquisition().getAmplitudeContrast())
        if self.isTomo5.get():
            ts_obj.getAcquisition().setTiltAxisAngle(-1 * mdoc_obj.getTiltAxisAngle() - 90)
        else:
            ts_obj.getAcquisition().setTiltAxisAngle(mdoc_obj.getTiltAxisAngle())

        origin = Transform()
        ts_obj.setOrigin(origin)
        SOTS.setAcquisition(acq)
        SOTS.append(ts_obj)

        self.settingTS(SOTS, ts_obj, file_ordered_angle_list, incoming_dose_list)

        SOTS.write()
        self._store(SOTS)

    def settingTS(self, SOTS, ts_obj, file_ordered_angle_list, incoming_dose_list):
        """
        Set all the info in each tilt and set the ts_obj information with all
        the tilts

        :param SOTS: Set of tilt series.
        :param ts_obj: Tilt series object to add tilts too.
        :param file_ordered_angle_list: list of files sorted by angle.
        :param incoming_dose_list: list of dose per tilt.
        :return:
        """

        ts_fn = self._getOutputTiltSeriesPath(ts_obj)
        counter_ti = 0

        TSAngleFile = self._getExtraPath("{}.rawtlt".format(ts_obj.getTsId()))
        TSAngleFile = open(TSAngleFile, "a")
        for n in file_ordered_angle_list:
            TSAngleFile.write('{}\n'.format(str(n[2])))
        TSAngleFile.close()
        sr = self.listOfMics[0].getSamplingRate()
        properties = {"sr": sr}
        newStack = ImageStack(properties=properties)
        ti = None
        for f, to, ta in file_ordered_angle_list:
            try:
                for mic in self.listOfMics:
                    if ts_obj.getSamplingRate() is None:
                        ts_obj.setSamplingRate(sr)
                    if SOTS.getSamplingRate() is None:
                        SOTS.setSamplingRate(sr)
                    if os.path.basename(f) in mic.getMicName():
                        ti = tomoObj.TiltImage()
                        ti.setTsId(ts_obj.getTsId())
                        new_location = (counter_ti, ts_fn)
                        ti.setLocation(new_location)
                        ti.setObjId(counter_ti + 1)
                        ti.setIndex(counter_ti + 1)
                        ti.setAcquisition(ts_obj.getAcquisition())
                        ti.setAcquisitionOrder(int(to))
                        ti.setTiltAngle(ta)
                        ti.setSamplingRate(sr)
                        ti.setAcquisition(ts_obj.getAcquisition().clone())
                        dosePerFrame = incoming_dose_list[int(to) - 1]
                        ti.getAcquisition().setDosePerFrame(dosePerFrame)
                        ti.getAcquisition().setAccumDose(int(to)*dosePerFrame)
                        newStack.append(ImageReadersRegistry.open(mic.getFileName()))
                        ts_obj.append(ti)
                        counter_ti += 1
            except Exception as e:
                self.error(e)
                return
        ImageReadersRegistry.write(newStack, ts_fn, isStack=True)
        ts_obj._setFirstDim(ti)

        SOTS.update(ts_obj)

    # -------------------------- AUXILIARY FUNCTIONS -----------------------
    def _getOutputTiltSeriesPath(self, ts, suffix=''):
        return self._getExtraPath('%s%s.mrcs' % (ts.getTsId(), suffix))

    def _getOutputTiltImagePaths(self, tilt_image):
        """ Return expected output path for correct movie and DW one.
        """
        base = self._getExtraPath(self._getTiltImageMRoot(tilt_image))
        return base + '.mrc', base + '_Out.mrc'


    @staticmethod
    def _getTiltImageMRoot(ti):
        return '%s_%02d' % (ti.getTsId(), ti.getObjId())

    def getTimeOutInSeconds(self, timeOut):
        timeOutFormatRegexList = {r'\d+s': 1, r'\d+m': 60, r'\d+h': 3600,
                                  r'\d+d': 72000}
        try:
            return int(timeOut)
        except Exception:
            seconds = 0
        for regex, secondsUnit in timeOutFormatRegexList.items():
            matchingTimes = re.findall(regex, timeOut)
            for matchTime in matchingTimes:
                seconds += int(matchTime[:-1]) * secondsUnit
        return seconds


    def _validate(self):
        pass


    def _summary(self):
        summary = []
        summary.append('Path with the *.mdoc files for each tilt serie:{}\n'.format(self.filesPath.get()))
        if not hasattr(self, 'TiltSeries'):
            summary.append("Output SetOfTiltSeries not ready yet.")
        else:
            try:
                summary.append("{} tilt series added".format(self.TiltSeries.getSize()))
            except Exception as e:
                print(e)
        return summary


