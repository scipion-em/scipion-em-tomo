# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
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

import os
import re
from glob import glob
import time
from datetime import timedelta, datetime
from collections import OrderedDict
from os.path import join
from pathlib import PureWindowsPath

import numpy as np
from sqlite3 import OperationalError

import pyworkflow as pw
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.objects import Acquisition
from pyworkflow.utils import removeExt, getParentFolder, removeBaseExt
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtImport

from tomo.convert import parseMdoc, getAnglesFromHeader, getAnglesFromMdoc, getAnglesFromTlt
from .protocol_base import ProtTomoBase


class ProtImportTsBase(ProtImport, ProtTomoBase):
    """ Base class for Tilt-Series and Tilt-SeriesMovies import protocols.
    """
    IMPORT_FROM_FILES = 0

    # How to handle the input files into the project
    IMPORT_COPY_FILES = 0
    IMPORT_LINK_ABS = 1
    IMPORT_LINK_REL = 2

    ANGLES_FROM_FILENAME = 'Filename'
    ANGLES_FROM_HEADER = 'Header'
    ANGLES_FROM_MDOC = 'Mdoc'
    ANGLES_FROM_TLT = 'Tlt'
    ANGLES_FROM_RANGE = 'Range'

    NOT_MDOC_GUI_COND = 'filesPattern is None or (filesPattern is not None and "mdoc" not in filesPattern)'
    MDOC_DATA_SOURCE = False

    acquisitions = None
    sRates = None
    accumDoses = None
    meanDosesPerFrame = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Root directory of the tilt-series (or movies).")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      help="Pattern of the tilt series\n\n"
                           "The pattern can contain standard wildcards such as *, ?, etc.\n\n"
                           "It should also contains the following special tags:\n"
                           "   {TS}: tilt series identifier, which can be any UNIQUE part of the path. This must be an "
                           "alpha-numeric sequence (avoid symbols as -) that can not start with a number.\n"
                           "   {TO}: acquisition order, an integer value (important for dose).\n"
                           "   {TA}: tilt angle, a positive or negative float value.\n\n"
                           "Examples:\n"
                           "To import a set of image stacks (tilt-series or tilt-series movies) as: \n"
                           "TiltSeries_a_001_0.0.mrc\n"
                           "TiltSeries_a_002_3.0.mrc\n"
                           "TiltSeries_a_003_-3.0.mrc\n"
                           "...\n"
                           "TiltSeries_b_001_0.0.mrc\n"
                           "TiltSeries_b_002_3.0.mrc\n"
                           "TiltSeries_b_003_-3.0.mrc\n"
                           "...\n"
                           "The pattern TiltSeries_{TS}_{TO}_{TA}.mrc will identify:\n"
                           "{TS} as a, b, ...\n"
                           "{TO} as 001, 002, 003, ...\n"
                           "{TA} as 0.0, 3.0, -3.0, ...\n")
        self._defineAngleParam(form)
        form.addParam('importAction', params.EnumParam,
                      default=self.IMPORT_LINK_REL,
                      choices=['Copy files',
                               'Absolute symlink',
                               'Relative symlink'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Import action on files",
                      help="This parameters determine how the project will deal "
                           "with imported files. It can be: \n"
                           "*Copy files*: Input files will be copied into your "
                           "project. (this will duplicate the raw data)."
                           "*Absolute symlink*: Create symbolic links to the "
                           "absolute path of the files."
                           "*Relative symlink*: Create symbolic links as "
                           "relative path from the protocol run folder. ")

        self._defineAcquisitionParams(form, isTiltSeries=True)

        form.addSection('Streaming')

        form.addParam('dataStreaming', params.BooleanParam, default=False,
                      label="Process data in streaming?",
                      help="Select this option if you want import data as it is "
                           "generated and process on the fly by next protocols. "
                           "In this case the protocol will keep running to check "
                           "new files and will update the output Set, which can "
                           "be used right away by next steps.")

        form.addParam('timeout', params.IntParam, default=43200,
                      condition='dataStreaming',
                      label="Timeout (secs)",
                      help="Interval of time (in seconds) after which, if no new file "
                           "is detected, the protocol will end. When finished, "
                           "the output Set will be closed and no more data will be "
                           "added to it. \n"
                           "Note 1:  The default value is  high (12 hours) to avoid "
                           "the protocol finishes during the acq of the "
                           "microscope. You can also stop it from right click and press "
                           "STOP_STREAMING.\n"
                           "Note 2: If you're using individual frames when importing "
                           "movies, the timeout won't be refreshed until a whole "
                           "movie is stacked.")

        form.addParam('fileTimeout', params.IntParam, default=30,
                      condition='dataStreaming',
                      label="File timeout (secs)",
                      help="Interval of time (in seconds) after which, if a file has "
                           "not changed, we consider it as a new file. \n")

        self._defineBlacklistParams(form)

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        pass

    def _defineAcquisitionParams(self, form, isTiltSeries):
        """ Define acq parameters, it can be overridden
        by subclasses to change what parameters to include.
        """
        acqGroupCondition = True  # Tilt series movies case
        if isTiltSeries:
            acqGroupCondition = self.NOT_MDOC_GUI_COND

        group = form.addGroup('Acquisition info',
                              condition=acqGroupCondition)
        group.addParam('voltage', params.FloatParam, default=300,
                       label=Message.LABEL_VOLTAGE,
                       condition=self.NOT_MDOC_GUI_COND,
                       help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', params.FloatParam, default=2.7,
                       label=Message.LABEL_SPH_ABERRATION,
                       expertLevel=params.LEVEL_ADVANCED,
                       help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
                       label=Message.LABEL_AMPLITUDE,
                       expertLevel=params.LEVEL_ADVANCED,
                       help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', params.IntParam, default=50000,
                       label=Message.LABEL_MAGNI_RATE,
                       condition=self.NOT_MDOC_GUI_COND,
                       help=Message.TEXT_MAGNI_RATE)
        group.addParam('samplingRate', params.FloatParam, default=1.0,
                       important=True,
                       label=Message.LABEL_SAMP_RATE,
                       condition=self.NOT_MDOC_GUI_COND,
                       help=Message.TEXT_SAMP_RATE)

        return group

    def _defineBlacklistParams(self, form):
        """ Override to add options related to blacklist info.
        """
        pass

    # -------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('importStep')

    # -------------------------- STEPS functions -------------------------------
    def importStep(self):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        doseList = []
        counter = 0

        if not self.MDOC_DATA_SOURCE:
            self.info("Using glob pattern: '%s'" % self._globPattern)
            self.info("Using regex pattern: '%s'" % self._regexPattern)

        outputSet = getattr(self, self._outputName, None)

        if outputSet is None:
            createSetFunc = getattr(self, self._createOutputName)
            outputSet = createSetFunc()
        elif outputSet.getSize() > 0:
            outputSet.loadAllProperties()

            for ts in outputSet:
                self._existingTs.add(ts.getTsId())

        self._fillAcquisitionInfo(outputSet)
        tsClass = outputSet.ITEM_TYPE
        tiClass = tsClass.ITEM_TYPE

        finished = False
        lastDetectedChange = datetime.now()

        # Ignore the timeout variables if we are not really in streaming mode
        if self.dataStreaming:
            timeout = timedelta(seconds=self.timeout.get())
            fileTimeout = timedelta(seconds=self.fileTimeout.get())
        else:
            timeout = timedelta(seconds=5)
            fileTimeout = timedelta(seconds=5)

        while not finished:
            time.sleep(3)  # wait 3 seconds before check for new files
            someNew = False  # Check if some new TS has been found
            someAdded = False  # Check if some new were added
            # incompleteTs = False  # Check if there are incomplete TS

            matchingFiles = self.getMatchingFiles(fileTimeOut=fileTimeout)

            if self._existingTs:
                outputSet.enableAppend()

            for ts, tiltSeriesList in matchingFiles.items():
                someNew = True

                tsObj = tsClass(tsId=ts)
                # we need this to set mapper before adding any item
                outputSet.append(tsObj)

                if self.MDOC_DATA_SOURCE:
                    doseList = self.accumDoses[ts]
                    counter = 0

                # Add tilt images to the tiltSeries
                for f, to, ta in tiltSeriesList:
                    try:
                        ti = tiClass(location=f,
                                     acquisitionOrder=to,
                                     tiltAngle=ta)
                        if self.MDOC_DATA_SOURCE:
                            ti.setAcquisition(tsObj.getAcquisition())
                            ti.getAcquisition().setDosePerFrame(doseList[counter])  # Accumulated dose in current ti
                            counter += 1
                        tsObj.append(ti)
                    except OperationalError as e:

                        raise Exception("%s is an invalid for the {TS} field, it must be an alpha-numeric sequence "
                                        "(avoid symbols as -) that can not start with a number." % ts)

                if self.MDOC_DATA_SOURCE:
                    # Tilt series object dose per frame has been updated each time the tilt image dose per frame has
                    # been updated before, so the mean value is used to be the reference in the acquisition of the
                    # whole tilt series movie
                    tsObj.getAcquisition().setDosePerFrame(self.meanDosesPerFrame[ts])

                outputSet.update(tsObj)  # update items and size info
                self._existingTs.add(ts)
                someAdded = True

            if someAdded:
                self.debug('Updating output...')
                outputSet.updateDim()
                self._updateOutputSet(self._outputName, outputSet,
                                      state=outputSet.STREAM_OPEN)
                self.debug('Update Done.')

            self.debug('Checking if finished...someNew: %s' % someNew)

            now = datetime.now()

            if not someNew:
                # If there are no new detected files, we should check the
                # inactivity time elapsed (from last event to now) and
                # if it is greater than the defined timeout, we conclude
                # the import and close the output set
                # Another option is to check if the protocol have some
                # special stop condition, this can be used to manually stop
                # some protocols such as import movies
                finished = (now - lastDetectedChange > timeout
                            or self.streamingHasFinished())
                self.debug("Checking if finished:")
                self.debug("   Now - Last Change: %s"
                           % pwutils.prettyDelta(now - lastDetectedChange))
            else:
                # If we have detected some files, we should update
                # the timestamp of the last event
                lastDetectedChange = now
                finished = not self.isInStreaming()

            self.debug("Finished: %s" % finished)

        # Close the output set
        self._updateOutputSet(self._outputName, outputSet,
                              state=outputSet.STREAM_CLOSED)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        self._initialize()
        try:
            matching = self.getMatchingFiles()

        except Exception as e:
            errorStr = str(e)
            return [errorStr]
            # if 'Missing angles file: ' in errorStr:
            #     return [errorStr]
            # else:
            #     raise e

        if not matching:
            return ["There are no files matching the pattern %s"
                    % self._globPattern]

        self._firstMatch = list(matching.items())[0]
        self._tiltAngleList = self._getSortedAngles(self._firstMatch[1])

        return self._validateAngles()

    def _validateAngles(self):
        """ Function to be implemented in subclass to validate
        the angles range.
        """
        return []

    # -------------------------- BASE methods to be overridden -----------------
    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages form such as: xmipp3, eman2, relion...etc.
        """
        return ['files']

    # -------------------------- UTILS functions -------------------------------
    def _initialize(self):
        """ Initialize some internal variables such as:
        - patterns: Expand the pattern using environ vars or username
            and also replacing special character # by digit matching.
        - outputs: Output variable names
        """
        self.MDOC_DATA_SOURCE = 'mdoc' in self.filesPattern.get()
        if not self.MDOC_DATA_SOURCE:
            path = self.filesPath.get('').strip()
            pattern = self.filesPattern.get('').strip()
            self._pattern = os.path.join(path, pattern) if pattern else path

            def _replace(p, ts, to, ta):
                p = p.replace('{TS}', ts)
                p = p.replace('{TO}', to)
                p = p.replace('{TA}', ta)
                return p

            self._regexPattern = _replace(self._pattern.replace('*', '(.*)'),
                                          '(?P<TS>.*)',
                                          '(?P<TO>\d+)',
                                          '(?P<TA>[+-]?\d+(\.\d+)?)')
            self._regex = re.compile(self._regexPattern)
            self._globPattern = _replace(self._pattern, '*', '*', '*')

        # Set output names depending on the import type (either movies or images)
        self._outputName = 'outputTiltSeries'
        self._createOutputName = '_createSetOfTiltSeries'

        # Keep track of which existing tilt-series has already been found
        self._existingTs = set()

    def _anglesInPattern(self):
        """ This function should be called after a call to _initialize"""
        return '{TA}' in self._pattern and '{TO}' in self._pattern

    def _getMatchingFilesFromMdoc(self):
        """If the list of files provided by the user is a list of mdoc files, then the tilt series movies
        are built from them, following the considerations listed below:
            - For each mdoc file, it and the corresponding movie files must be in the same directory.
            - The tilt series id will be the base name of the mdoc file, by default, so the mdocs must have different
              base name. If another name is desired, the user can introduce the name structure (see advanced parameter)
            """
        fpath = self.filesPath.get()
        mdocList = glob(join(fpath, self.filesPattern.get()))
        if not mdocList:
            raise Exception('No mdoc files were found in the introduced path:\n%s' % fpath)

        def _validateMdocInfoRead():
            msg = ['Data not found in file\n%s:\n' % mdoc]
            if not acquisition.getVoltage():
                msg.append('Voltage')
            if not acquisition.getMagnification():
                msg.append('Magnification')
            if not sRate:
                msg.append('PixelSpacing')
            if None in tiltAngleList:
                indices = [str(i) for i, v in enumerate(tiltAngleList) if v is None]
                msg.append('TiltAngle: %s' % ' '.join(indices))
            if len(msg) > 1:
                validationErrorMsgList.append(msg)

            return validationErrorMsgList

        matchingFiles = OrderedDict()
        self.acquisitions = OrderedDict()
        self.sRates = OrderedDict()
        self.accumDoses = OrderedDict()
        self.meanDosesPerFrame = OrderedDict()
        validationErrorMsgList = []
        for mdoc in mdocList:
            accumulatedDose = 0
            zValueList = parseMdoc(mdoc)
            z0 = zValueList[0]
            # Acquisition and pixel size are read from the first slice (acquisition angle)
            acquisition = self._genTsAcquisitionFromMdoc(z0)
            sRate = z0.get('PixelSpacing', self.samplingRate.get())  # Get the pixel size

            # fileList = []
            fileOrderAngleList = []
            accumulatedDoseList = []
            tiltAngleList = []

            counter = 0
            for zSlice in zValueList:
                fileName = self._getAngleMovieFileName(zSlice['SubFramePath'])
                tiltAngle = zSlice.get('TiltAngle', None)
                tiltAngleList.append(tiltAngle)
                fileOrderAngleList.append((
                    join(getParentFolder(mdoc), fileName),  # Filename
                    '{:03d}'.format(counter),               # Acquisition order
                    tiltAngle,                              # Tilt angle
                ))
                accumulatedDose = self._getDoseFromMdoc(zSlice, accumulatedDose)
                accumulatedDoseList.append(accumulatedDose)
                # fileList.append(removeExt(fileName))
                counter += 1
            tsId = removeBaseExt(mdoc)  # TODO: possible parameter (advanced) used by the user to indicate the structure
            # self._getTsIdFromMdocData(fileList)
            matchingFiles[tsId] = fileOrderAngleList
            self.acquisitions[tsId] = acquisition
            self.sRates[tsId] = sRate
            self.accumDoses[tsId] = accumulatedDoseList
            self.meanDosesPerFrame[tsId] = accumulatedDoseList[-1]/counter

            # Check Mdoc info read
            validationErrorMsgList = _validateMdocInfoRead()

        if validationErrorMsgList:
            raise Exception(' '.join(validationErrorMsgList))

        return matchingFiles

    @staticmethod
    def _getDoseFromMdoc(zSlice, accumulatedDose):
        """It calculates the accumulated dose on the frames represented by zSlice, and add it to the
        previous accumulated dose"""
        EXPOSURE_DOSE = 'ExposureDose'  # Dose on specimen during camera exposure in electrons/sq. Angstrom
        FRAME_DOSES_AND_NUMBERS = 'FrameDosesAndNumbers'  # Dose per frame in electrons per square Angstrom followed
        # by number of frames at that dose
        DOSE_RATE = 'DoseRate'  # Dose rate to the camera, in electrons per unbinned pixel per second
        EXPOSURE_TIME = 'ExposureTime'  # Image exposure time
        MIN_MAX_MEAN = 'MinMaxMean'  # Minimum, maximum, and mean value for this image
        PIXEL_SIZE = 'PixelSpacing'  # Pixel spacing in Angstroms for individual image
        COUNTS_PER_ELECTRON = 'CountsPerElectron'

        def _keysInDict(listOfKeys):
            return all([key in zSlice.keys() for key in listOfKeys])

        # Different ways of calculating the dose, ordered by priority considering the possible variability between
        # different mdoc files
        newDose = 0

        # Directly from field ExposureDose
        if EXPOSURE_DOSE in zSlice:
            expDoseVal = zSlice[EXPOSURE_DOSE]
            if expDoseVal:
                newDose = float(expDoseVal)

        # Directly from field FrameDosesAndNumbers
        if not newDose and FRAME_DOSES_AND_NUMBERS in zSlice:
            frameDoseAndNums = zSlice[FRAME_DOSES_AND_NUMBERS]
            if frameDoseAndNums:
                newDose = float(list(frameDoseAndNums.strip())[-1])  # Get the mean from a string like '0 6'

        # Calculated from fields DoseRate and ExposureTime
        if not newDose and _keysInDict([DOSE_RATE, EXPOSURE_TIME]):
            doseRate = zSlice[DOSE_RATE]
            expTime = zSlice[EXPOSURE_TIME]
            if doseRate and expTime:
                newDose = float(doseRate) * float(expTime)

        # Calculated from fields MinMaxMean, PixelSpacing and CountsPerElectron
        if not newDose and _keysInDict([MIN_MAX_MEAN, PIXEL_SIZE, COUNTS_PER_ELECTRON]):
            minMaxMean = zSlice[MIN_MAX_MEAN]
            pixelSize = zSlice[PIXEL_SIZE]
            counts = zSlice[COUNTS_PER_ELECTRON]
            if all([minMaxMean, pixelSize, counts]):
                meanVal = list(minMaxMean.strip())[-1]  # Get the mean from a string like '-42 2441 51.7968'
                newDose = (float(meanVal)/int(counts)) / float(pixelSize)**2

        return newDose + accumulatedDose

    @staticmethod
    def _getAngleMovieFileName(subFramePath):
        # PureWindowsPath pathlib is ised to make possible deal with different path separators, like \\
        return PureWindowsPath(subFramePath).parts[-1]

    # @staticmethod
    # def _getTsIdFromMdocData(baseNameList):
    #     """Find the most common substring in the list of subFrame base names to get the tsId"""
    #     def _findStem(arr):
    #         """function to find the stem (longest common substring) from the string array"""
    #         # Determine size of the array
    #         n = len(arr)
    #         # Take first word from array as reference
    #         s = arr[0]
    #         l = len(s)
    #         res = ""
    #         for i in range(l):
    #             for j in range(i + 1, l + 1):
    #                 # generating all possible substrings
    #                 # of our reference string arr[0] i.e s
    #                 stem = s[i:j]
    #                 k = 1
    #                 for k in range(1, n):
    #                     # Check if the generated stem is
    #                     # common to all words
    #                     if stem not in arr[k]:
    #                         break
    #                 # If current substring is present in
    #                 # all strings and its length is greater
    #                 # than current result
    #                 if k + 1 == n and len(res) < len(stem):
    #                     res = stem
    #
    #         return res
    #
    #     d = _findStem(baseNameList).replace('_0', '')
    #     return d

    def _genTsAcquisitionFromMdoc(self, zSlice):
        acq = Acquisition()
        acq.setVoltage(zSlice.get('Voltage', self.voltage.get()))
        acq.setSphericalAberration(self.sphericalAberration.get())
        acq.setAmplitudeContrast(self.amplitudeContrast.get())
        acq.setMagnification(zSlice.get('Magnification', self.magnification.get()))
        if hasattr(self, 'doseInitial'):  # This field is only present in the form for TsM import
            acq.setDoseInitial(self.doseInitial.get())

        return acq

    def getMatchingFiles(self, fileTimeOut=None):
        """ Return an ordered dict with TiltSeries found in the files as key
        and a list of all tilt images of that series as value.
        """
        if self.MDOC_DATA_SOURCE:
            return self._getMatchingFilesFromMdoc()
        else:
            filePaths = glob(self._globPattern)
            filePaths.sort(key=lambda fn: os.path.getmtime(fn))
            fileTimeOut = fileTimeOut or timedelta(seconds=5)

            matchingFiles = OrderedDict()

            def _getTsId(match):
                """ Retrieve the TiltSerie ID from the matching object.
                We need to have tsId that starts with character, so
                let's add a prefix if it is not the case
                """
                tsId = match.group('TS')
                return tsId

            def _addOne(fileList, f, m):
                """ Add one file matching to the list. """
                fileList.append((f, int(m.group('TO')), float(m.group('TA'))))

            def _addMany(fileList, f, m):
                """ Add many 'files' (when angles in header or mdoc) to the list. """
                anglesFrom = self.getEnumText('anglesFrom')

                if anglesFrom == self.ANGLES_FROM_HEADER:
                    angles = getAnglesFromHeader(f)
                elif anglesFrom == self.ANGLES_FROM_MDOC:
                    mdocFn = os.path.splitext(f)[0] + '.mdoc'
                    if not os.path.exists(mdocFn):
                        raise Exception("Missing angles file: %s" % mdocFn)
                    angles = getAnglesFromMdoc(mdocFn)
                elif anglesFrom == self.ANGLES_FROM_TLT:
                    tltFn = os.path.splitext(f)[0] + '.tlt'
                    if not os.path.exists(tltFn):
                        raise Exception("Missing angles file: %s" % tltFn)
                    angles = getAnglesFromTlt(tltFn)
                elif anglesFrom == self.ANGLES_FROM_RANGE:
                    angles = self._getTiltAngleRange()
                else:
                    raise Exception('Invalid angles option: %s' % anglesFrom)

                for i, a in enumerate(angles):
                    fileList.append(((i+1, f), i+1, a))

            addFunc = _addOne if self._anglesInPattern() else _addMany

            # Handle special case of just one TiltSeries, to avoid
            # the user the need to specify {TS}
            if len(filePaths) == 1 and not self.isInStreaming():
                f = filePaths[0]
                ts = pwutils.removeBaseExt(f)  # Base name without extension
                matchingFiles[ts] = []
                _addMany(matchingFiles[ts], f, None)
            else:
                for f in filePaths:
                    if self.fileModified(f, fileTimeOut):
                        continue
                    m = self._regex.match(f)
                    if m is not None:
                        ts = _getTsId(m)
                        # Only report files of new tilt-series
                        if ts not in self._existingTs:
                            if ts not in matchingFiles:
                                matchingFiles[ts] = []
                            addFunc(matchingFiles[ts], f, m)

            return matchingFiles

    def getCopyOrLink(self):
        # Set a function to copyFile or createLink
        # depending in the user selected option
        if self.copyFiles:
            return pw.utils.copyFile
        else:
            return pw.utils.createAbsLink

    def fileModified(self, fileName, fileTimeout):
        """ Check if the fileName modification time is less
        than a given timeout.
        Params:
            fileName: input filename that will be checked.
            fileTimeout: timeout
        """
        self.debug('Checking file: %s' % fileName)
        mTime = datetime.fromtimestamp(os.path.getmtime(fileName))
        delta = datetime.now() - mTime
        self.debug('   Modification time: %s' % pw.utils.prettyTime(mTime))
        self.debug('   Delta: %s' % pw.utils.prettyDelta(delta))

        return delta < fileTimeout

    def isBlacklisted(self, fileName):
        """ Overwrite in subclasses """
        return False

    @classmethod
    def worksInStreaming(cls):
        # Import protocols always work in streaming
        return True

    def _fillAcquisitionInfo(self, inputTs):
        inputTs.setSamplingRate(self.samplingRate.get())
        acq = inputTs.getAcquisition()
        acq.setVoltage(self.voltage.get())
        acq.setSphericalAberration(self.sphericalAberration.get())
        acq.setAmplitudeContrast(self.amplitudeContrast.get())
        acq.setMagnification(self.magnification.get())

    def _getTiltAngleRange(self):
        """ Return the list with all expected tilt angles. """
        offset = 1
        if self.minAngle.get() > self.maxAngle.get():
            offset = -1 * offset
        return np.arange(self.minAngle.get(),
                         self.maxAngle.get() + offset,  # also include last angle
                         self.stepAngle.get())

    @staticmethod
    def _getSortedAngles(tiltSeriesList):
        """ Return the sorted angles from a given tiltSeriesList. """
        return sorted(item[2] for item in tiltSeriesList)

    def _sameTiltAngleRange(self, tiltAngleRange, tiltSeriesList):
        # allow some tolerance when comparing tilt-angles
        return np.allclose(tiltAngleRange,
                           self._getSortedAngles(tiltSeriesList),
                           atol=0.1)

    # --------------- Streaming special functions -----------------------
    def _getStopStreamingFilename(self):
        return self._getExtraPath("STOP_STREAMING.TXT")

    def getActions(self):
        """ This method will allow that the 'Stop import' action to appears
        in the GUI when the user right-click in the protocol import box.
        It will allow a user to manually stop the streaming.
        """
        # Only allow to stop if running and in streaming mode
        if self.dataStreaming and self.isRunning():
            return [('STOP STREAMING', self.stopImport)]
        else:
            return []

    def stopImport(self):
        """ Since the actual protocol that is running is in a different
        process that the one that this method will be invoked from the GUI,
        we will use a simple mechanism to place an special file to stop
        the streaming.
        """
        # Just place an special file into the run folder
        f = open(self._getStopStreamingFilename(), 'w')
        f.close()

    def streamingHasFinished(self):
        return (not self.isInStreaming() or
                os.path.exists(self._getStopStreamingFilename()))

    def isInStreaming(self):
        return self.dataStreaming.get()


class ProtImportTs(ProtImportTsBase):
    _label = 'import tilt-series'

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        group = form.addGroup('Tilt info',
                              condition=self.NOT_MDOC_GUI_COND)

        group.addParam('anglesFrom', params.EnumParam,
                       default=0,
                       choices=[self.ANGLES_FROM_RANGE,
                                self.ANGLES_FROM_HEADER,
                                # self.ANGLES_FROM_MDOC,
                                self.ANGLES_FROM_TLT],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Import angles from',
                       help='Select how the tilt angles will be imported. '
                            'It can be defined by range: Min, Max, Step '
                            'or from image header, or from complementary'
                            'mdoc or tlt files (should have the same filename '
                            'but with the .mdoc or .tlt extension).')

        line = group.addLine('Tilt angular range',
                             condition='anglesFrom==0',  # ANGLES_FROM_RANGE
                             help="Specify the tilting angular range. "
                                  "Depending on the collection schema, the "
                                  "order of the acquisition does not need to "
                                  "be the same order of the angular range. ")
        line.addParam('minAngle', params.FloatParam, default=-60, label='min')
        line.addParam('maxAngle', params.FloatParam, default=60, label='max')
        line.addParam('stepAngle', params.FloatParam, default=3, label='step')

    def _validateAngles(self):
        if not self.MDOC_DATA_SOURCE:
            ts, tiltSeriesList = self._firstMatch
            i, fileName = tiltSeriesList[0][0]
            x, y, z, n = ImageHandler().getDimensions(fileName)
            nImages = max(z, n)  # Just handle ambiguity with mrc format
            nAngles = len(self._tiltAngleList)
            if nAngles != nImages:
                return ['Tilt-series %s stack has different number of images (%d) '
                        'than the expected number of tilt angles (%d). '
                        % (fileName, nImages, nAngles)]
        return []


class ProtImportTsMovies(ProtImportTsBase):
    _label = 'import tilt-series movies'

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        group = form.addGroup('Tilt info',
                              condition=self.NOT_MDOC_GUI_COND)

        group.addParam('anglesFrom', params.EnumParam,
                       default=0,
                       choices=[self.ANGLES_FROM_FILENAME],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Import angles from',
                       help='Angles will be parsed from the filename pattern.'
                            'The special token {TA} should be specified as part'
                            'of the pattern, that will be used to match the '
                            'value of the angle for each TiltSeriesMovie.')

    def _defineAcquisitionParams(self, form, isTiltSeries):
        """ Add movie specific options to the acquisition section. """
        group = ProtImportTsBase._defineAcquisitionParams(self, form, isTiltSeries=False)

        line = group.addLine('Dose (e/A^2)',
                             condition=self.NOT_MDOC_GUI_COND,
                             help="Initial accumulated dose (usually 0) and "
                                  "dose per frame. ")
        line.addParam('doseInitial', params.FloatParam, default=0,
                      label='Initial')
        line.addParam('dosePerFrame', params.FloatParam, default=None,
                      allowsNull=True,
                      label='Per frame')
        group.addParam('gainFile', params.FileParam,
                       label='Gain image',
                       help='A gain reference related to a set of movies '
                            'for gain correction')
        group.addParam('darkFile', params.FileParam,
                       label='Dark image',
                       help='A dark image related to a set of movies')
        return group

    def _fillAcquisitionInfo(self, inputTs):
        inputTs.setGain(self.gainFile.get())
        inputTs.setDark(self.darkFile.get())

        if self.MDOC_DATA_SOURCE:
            # TODO Acquisition is historically expected as one per set of tilt series movies, so the first one is hte one used, at least for now
            firstTsId = list(self.acquisitions.keys())[0]
            acq = self.acquisitions[firstTsId]
            sRate = self.sRates[firstTsId]
            inputTs.setAcquisition(acq)
            inputTs.setSamplingRate(sRate)
        else:
            ProtImportTsBase._fillAcquisitionInfo(self, inputTs)
            acq = inputTs.getAcquisition()
            acq.setDoseInitial(self.doseInitial.get())
            acq.setDosePerFrame(self.dosePerFrame.get())

    def _initialize(self):
        ProtImportTsBase._initialize(self)
        self._outputName += 'M'
        self._createOutputName += 'M'

    def _validateAngles(self):
        """ Function to be implemented in subclass to validate
        the angles range.
        """
        if not self.MDOC_DATA_SOURCE and self.getEnumText('anglesFrom') == self.ANGLES_FROM_FILENAME:
            if not self._anglesInPattern():
                return ['When importing movies, {TA} and {TO} should be in the '
                        'files pattern.']
        return []
