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
import numpy as np

import pyworkflow as pw
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message

import tomo.convert
from .protocol_base import ProtTomoBase


class ProtImportTiltSeries(pwem.ProtImport, ProtTomoBase):
    """ Base class for other Import protocols.
    All imports protocols will have:
    1) Several options to import from (_getImportOptions function)
    2) First option will always be "from files". (for this option
      files with a given pattern will be retrieved  and the ### will
      be used to mark an ID part from the filename.
      - For each file a function to process it will be called
        (_importFile(fileName, fileId))
    """
    IMPORT_FROM_FILES = 0

    # How to handle the input files into the project
    IMPORT_COPY_FILES = 0
    IMPORT_LINK_ABS = 1
    IMPORT_LINK_REL = 2

    IMPORT_TYPE_MICS = 0
    IMPORT_TYPE_MOVS = 1

    ANGLES_FROM_FILES = 0
    ANGLES_FROM_HEADER = 1
    ANGLES_FROM_MDOC = 2

    _label = 'import tilt-series'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('importType', params.EnumParam, default=1,
                      choices=['Tilt series', 'Tilt series Movies'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Select type of input images',
                      help='If you import tilt-series movies, then'
                           'you need to provide also dose information'
                           'and maybe gain/dark images.')
        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Root directory of the tilt-series (or movies).")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      help="Pattern of the tilt series\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc.\n\n"
                           "It should also contains the following special tags:"
                           "   {TS}: tilt series identifier "
                           "         (can be any UNIQUE part of the path).\n"
                           "   {TO}: acq order"
                           "         (an integer value, important for dose).\n"
                           "   {TA}: tilt angle"
                           "         (positive or negative float value).\n\n"
                           "Examples:\n"
                           "")
        form.addParam('anglesFrom', params.EnumParam,
                      default=self.ANGLES_FROM_FILES,
                      choices=['Files', 'Header', 'Mdoc'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Import angles from',
                      help='By default, the angles should be provided as part '
                           'of the files pattern, using the {TA} token. '
                           'Additionally, in the case of importing tilt-series,'
                           'angles can be read from the header or from mdoc '
                           'files.')
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

        group = form.addGroup('Tilt info')
        line = group.addLine('Tilt angular range',
                             help="Specify the tilting angular range. "
                                  "Depending on the collection schema, the "
                                  "order of the acquisition does not need to "
                                  "be the same order of the angular range. ")
        line.addParam('minAngle', params.FloatParam, default=-60, label='min')
        line.addParam('maxAngle', params.FloatParam, default=60, label='max')
        line.addParam('stepAngle', params.FloatParam, default=3, label='step')

        self._defineAcquisitionParams(form)

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

    def _defineAcquisitionParams(self, form):
        """ Define acq parameters, it can be overridden
        by subclasses to change what parameters to include.
        """
        group = form.addGroup('Acquisition info')
        group.addParam('voltage', params.FloatParam, default=300,
                       label=Message.LABEL_VOLTAGE,
                       help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', params.FloatParam, default=2.7,
                       label=Message.LABEL_SPH_ABERRATION,
                       help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
                       label=Message.LABEL_AMPLITUDE,
                       help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', params.IntParam, default=50000,
                       label=Message.LABEL_MAGNI_RATE,
                       help=Message.TEXT_MAGNI_RATE)
        group.addParam('samplingRate', params.FloatParam,  default=1.0,
                       important=True,
                       label=Message.LABEL_SAMP_RATE,
                       help=Message.TEXT_SAMP_RATE)

        moviesCond = 'importType==%d' % self.IMPORT_TYPE_MOVS
        line = group.addLine('Dose (e/A^2)',
                             condition=moviesCond,
                             help="Initial accumulated dose (usually 0) and "
                                  "dose per frame. ")
        line.addParam('doseInitial', params.FloatParam, default=0,
                      label='Initial')
        line.addParam('dosePerFrame', params.FloatParam, default=None,
                      allowsNull=True,
                      label='Per frame')
        group.addParam('gainFile', params.FileParam,
                       condition=moviesCond,
                       label='Gain image',
                       help='A gain reference related to a set of movies '
                            'for gain correction')
        group.addParam('darkFile', params.FileParam,
                       condition=moviesCond,
                       label='Dark image',
                       help='A dark image related to a set of movies')
        return group

    def _defineBlacklistParams(self, form):
        """ Override to add options related to blacklist info.
        """
        pass

    # -------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('importStep', self._pattern,
                                 self.voltage.get(),
                                 self.sphericalAberration.get(),
                                 self.amplitudeContrast.get(),
                                 self.magnification.get())

    # -------------------------- STEPS functions -------------------------------
    def importStep(self, pattern, voltage, sphericalAberration,
                         amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self._initialize()
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
        tiltAngleRange = self._getTiltAngleRange()
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
            incompleteTs = False  # Check if there are incomplete TS

            matchingFiles = self.getMatchingFiles(fileTimeOut=fileTimeout)

            if self._existingTs:
                outputSet.enableAppend()

            for ts, tiltSeriesList in matchingFiles.iteritems():
                someNew = True

                if len(tiltAngleRange) != len(tiltSeriesList):
                    if incompleteTs:
                        raise Exception("More than one tilt-series seems "
                                        "incomplete regarding the expected "
                                        "number of tilt images.")
                    incompleteTs = True
                    continue

                if not self._sameTiltAngleRange(tiltAngleRange, tiltSeriesList):
                    raise Exception("Invalid tilt angle range for tilt series: "
                                    "%s" % ts)

                tsObj = tsClass(tsId=ts)
                # we need this to set mapper before adding any item
                outputSet.append(tsObj)
                print("After append...")
                pwutils.prettyDict(tsObj.getObjDict(includeClass=True))
                # Add tilt images to the tiltSeries
                for f, to, ta in tiltSeriesList:
                    tsObj.append(tiClass(location=f,
                                         acquisitionOrder=to,
                                         tiltAngle=ta))

                outputSet.update(tsObj)  # update items and size info
                print("After update....")
                pwutils.prettyDict(tsObj.getObjDict(includeClass=True))
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
                self.debug("Finished: %s" % finished)
            else:
                # If we have detected some files, we should update
                # the timestamp of the last event
                lastDetectedChange = now

        # Close the output set
        self._updateOutputSet(self._outputName, outputSet,
                              state=outputSet.STREAM_CLOSED)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        self._initialize()
        matching = self.getMatchingFiles()

        if self.importType == self.IMPORT_TYPE_MOVS:
            if self.anglesFrom != self.ANGLES_FROM_FILES:
                errors.append('Angle information could only be taken from '
                              'the files pattern when importing tilt-series'
                              'movies.')
            elif not self._anglesInPattern():
                errors.append('When importing movies, {TA} and {TO} should be '
                              'in the files pattern.')
        else:
            if self.anglesFrom == self.ANGLES_FROM_MDOC:
                errors.append('Importing angles from mdoc is not yet implemented.')

        if not matching:
            errors.append("There are no files matching the pattern %s"
                          % self._globPattern)

        tiltAngleRange = self._getTiltAngleRange()
        n = len(tiltAngleRange)
        ts = matching.keys()[0]
        tiltSeriesList = matching.values()[0]

        if (len(tiltSeriesList) != n or
            not self._sameTiltAngleRange(tiltAngleRange, tiltSeriesList)):
            errors.append('Tilt-series %s differs from expected tilt '
                          'angle range. ')

        return errors

    # -------------------------- BASE methods to be overridden -----------------
    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages form such as: xmipp3, eman2, relion...etc.
        """
        return ['files']

    # -------------------------- UTILS functions -------------------------------
    def doImportMovies(self):
        """ Return True if we are importing TiltSeries movies. """
        return self.importType == self.IMPORT_TYPE_MOVS

    def _initialize(self):
        """ Initialize some internal variables such as:
        - patterns: Expand the pattern using environ vars or username
            and also replacing special character # by digit matching.
        - outputs: Output variable names
        """
        self._pattern = os.path.join(self.filesPath.get('').strip(),
                                     self.filesPattern.get('').strip())

        def _replace(p, ts, to, ta):
            p = p.replace('{TS}', ts)
            p = p.replace('{TO}', to)
            p = p.replace('{TA}', ta)
            return p

        self._regexPattern = _replace(self._pattern.replace('*', '(.*)'),
                                      '(?P<TS>.*)', '(?P<TO>\d+)',
                                      '(?P<TA>[+-]?\d+(\.\d+)?)')
        self._regex = re.compile(self._regexPattern)
        self._globPattern = _replace(self._pattern, '*', '*', '*')

        # Set output names depending on the import type (either movies or images)
        self._outputName = 'outputTiltSeries'
        self._createOutputName = '_createSetOfTiltSeries'

        if self.doImportMovies():
            self._outputName += 'M'
            self._createOutputName += 'M'

        # Keep track of which existing tilt-series has already been found
        self._existingTs = set()

    def _anglesInPattern(self):
        """ This function should be called after a call to _initialize"""
        return '{TA}' in self._pattern and '{TO}' in self._pattern

    def getMatchingFiles(self, fileTimeOut=None):
        """ Return an ordered dict with TiltSeries found in the files as key
        and a list of all tilt images of that series as value.
        """
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
            return 'TS_%s' % tsId if tsId[0].isdigit() else tsId

        def _addOne(fileList, f, m):
            """ Return one file matching. """
            fileList.append((f, int(m.group('TO')), float(m.group('TA'))))

        def _addMany(fileList, f, m):
            """ Return many 'files' (when angles in header or mdoc). """
            _, _, _, n = pwem.ImageHandler().getDimensions(f)
            # FIXME, read proper angles
            angles = tomo.convert.getAnglesFromHeader(f)
            for i, a in enumerate(angles):
                fileList.append(((i+1, f), i+1, a))

        addFunc = _addOne if self._anglesInPattern() else _addMany

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

    def worksInStreaming(self):
        # Import protocols always work in streaming
        return True

    def _fillAcquisitionInfo(self, inputTs):
        inputTs.setSamplingRate(self.samplingRate.get())
        acq = inputTs.getAcquisition()
        acq.setVoltage(self.voltage.get())
        acq.setSphericalAberration(self.sphericalAberration.get())
        acq.setAmplitudeContrast(self.amplitudeContrast.get())
        acq.setMagnification(self.magnification.get())

        if self.doImportMovies():
            inputTs.setGain(self.gainFile.get())
            inputTs.setDark(self.darkFile.get())
            acq.setDoseInitial(self.doseInitial.get())
            acq.setDosePerFrame(self.dosePerFrame.get())

    def _getTiltAngleRange(self):
        """ Return the list with all expected tilt angles. """
        return np.arange(self.minAngle.get(),
                         self.maxAngle.get() + 1,  # also include max angle
                         self.stepAngle.get())

    def _sameTiltAngleRange(self, tiltAngleRange, tiltSeriesList):
        return np.allclose(tiltAngleRange,
                           sorted(item[2] for item in tiltSeriesList))

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
        return os.path.exists(self._getStopStreamingFilename())
