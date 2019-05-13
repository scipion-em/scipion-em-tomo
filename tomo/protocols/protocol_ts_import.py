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
from datetime import datetime
from collections import OrderedDict

import pyworkflow as pw
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params
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
                      help="By default ...")

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
        """ Define acq parameters, it can be overriden
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
        self.loadPatterns()
        self._insertFunctionStep('importStep', self._pattern,
                                 self.voltage.get(),
                                 self.sphericalAberration.get(),
                                 self.amplitudeContrast.get(),
                                 self.magnification.get())

    # -------------------------- STEPS functions -------------------------------
    def doImportMovies(self):
        """ Return True if we are importing TiltSeries movies. """
        return self.importType == self.IMPORT_TYPE_MOVS

    def _createOutputSet(self, suffix=''):
        if self.doImportMovies():
            self._outputSuffix = 'M'
            return self._createSetOfTiltSeriesM()
        else:
            self._outputSuffix = ''
            return self._createSetOfTiltSeries()

    def importStep(self, pattern, voltage, sphericalAberration,
                         amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.loadPatterns()
        self.info("Using glob pattern: '%s'" % self._globPattern)
        self.info("Using regex pattern: '%s'" % self._regexPattern)

        outputSet = self._createOutputSet()
        self._fillAcquisitionInfo(outputSet)
        tsClass = outputSet.ITEM_TYPE
        tiClass = tsClass.ITEM_TYPE

        self.info("Files: ")
        tiltSeriesDict = OrderedDict()

        for f, ts, to, ta in self.getMatchingFiles():
            if ts not in tiltSeriesDict:
                tiltSeriesDict[ts] = []

            tiltSeriesList = tiltSeriesDict[ts]
            order = len(tiltSeriesList) + 1
            tiltSeriesList.append(tiClass(location=f,
                                          acqOrder=order,
                                          tiltAngle=ta))

        for ts, tiltSeriesList in tiltSeriesDict.iteritems():
            tsObj = tsClass(tsId=ts)
            # we need this to set mapper before adding any item
            outputSet.append(tsObj)
            # Add tilt images to the tiltSeries
            for tim in tiltSeriesList:
                tsObj.append(tim)

            outputSet.update(tsObj)  # update items and size info
            tsObj.printAll()

        outputSet.updateDim()
        outputName = 'outputTiltSeries%s' % self._outputSuffix
        self._defineOutputs(**{outputName: outputSet})

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
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

        return errors

    # -------------------------- BASE methods to be overridden -----------------
    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages form such as: xmipp3, eman2, relion...etc.
        """
        return ['files']

    # -------------------------- UTILS functions ------------------------------
    def loadPatterns(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
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

    def _anglesInPattern(self):
        """ This function should be called after a call to loadPatterns"""
        return '{TA}' in self._pattern and '{TO}' in self._pattern

    def getMatchingFiles(self):
        """ Return a sorted list with the paths of files that
        matched the pattern.
        """
        self.loadPatterns()

        filePaths = glob(self._globPattern)
        filePaths.sort()

        def _addOne(fileList, f, m):
            """ Return one file matching. """
            fileList.append(
                (f, m.group('TS'), int(m.group('TO')), float(m.group('TA'))))

        def _addMany(fileList, f, m):
            """ Return many 'files' (when angles in header or mdoc). """
            _, _, _, n = pwem.ImageHandler().getDimensions(f)
            # FIXME, read proper angles
            angles = tomo.convert.getAnglesFromHeader(f)
            for i, a in enumerate(angles):
                fileList.append(((i+1, f), m.group('TS'), i+1, a))

        addFunc = _addOne if self._anglesInPattern() else _addMany

        matchingFiles = []
        for f in filePaths:
            m = self._regex.match(f)
            if m is not None:
                addFunc(matchingFiles, f, m)

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
            fileTimeout: timeout """
        self.debug('Checking file: %s' % fileName)
        mTime = datetime.fromtimestamp(os.path.getmtime(fileName))
        delta = datetime.now() - mTime
        self.debug('   Modification time: %s' % pw.utils.prettyTime(mTime))
        self.debug('   Delta: %s' % pw.utils.prettyDelta(delta))

        return delta < fileTimeout

    def isBlacklisted(self, fileName):
        """ Overwrite in subclasses """
        return False

    def iterFiles(self):
        """ Iterate through the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()

        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'"
                                    % (fileName, self._pattern))

                fileId = int(match.group(1))

            else:
                fileId = None

            yield fileName, fileId

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
