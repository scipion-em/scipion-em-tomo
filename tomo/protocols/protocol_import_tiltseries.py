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

from tomo.objects import TiltSeriesM, TiltImageM


class ProtImportTiltSeries(pwem.ProtImport):
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

    _label = 'import tilt-series'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

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
                           "   {TO}: acquisition order"
                           "         (an integer value, important for dose).\n"
                           "   {TA}: tilt angle"
                           "         (positive or negative float value).\n\n"
                           "Examples:\n"
                           "")
        form.addParam('importAction', params.EnumParam,
                      choices=['Copy files',
                               'Absolute symlink',
                               'Relative symlink'],
                      default=self.IMPORT_LINK_REL,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Import action on files",
                      help="By default ...")

        self._defineImportParams(form)

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
                           "the protocol finishes during the acquisition of the "
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

    def _defineImportParams(self, form):
        """ Override to add options related to the different types
        of import that are allowed by each protocol.
        """
        pass

    def _defineAcquisitionParams(self, form):
        """ Define acquisition parameters, it can be overriden
        by subclasses to change what parameters to include.
        """
        group = form.addGroup('Acquisition info')
        group.addParam('haveDataBeenPhaseFlipped', params.BooleanParam,
                       default=False,
                       label='Have data been phase-flipped?',
                       help='Set this to Yes if the images have been ctf-phase '
                            'corrected.')
        group.addParam('acquisitionWizard', params.LabelParam, important=True,
                       label='Use the wizard button to import acquisition.',
                       help='Depending on the import Format, the wizard\n'
                            'will try to import the acquisition values.\n'
                            'If not found, required ones should be provided.')
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
        return group

    def _defineBlacklistParams(self, form):
        """ Override to add options related to blacklist info.
        """
        pass

    # -------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):

        # Only the import movies has property 'inputIndividualFrames'
        # so let's query in a non-intrusive manner
        # inputIndividualFrames = getattr(self, 'inputIndividualFrames', False)

        # if self.dataStreaming or inputIndividualFrames:
        #     funcName = 'importImagesStreamStep'
        # else:
        #     funcName = 'importStep'

        self._insertFunctionStep('importStep', self.getPattern(),
                                 self.voltage.get(),
                                 self.sphericalAberration.get(),
                                 self.amplitudeContrast.get(),
                                 self.magnification.get())

    # -------------------------- STEPS functions -------------------------------

    def __createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pw.utils.cleanPath(setFn)
        from pyworkflow.mapper.sqlite_db import SqliteDb
        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj

    def _createSetOfTiltSeries(self, suffix=''):
        return self.__createSet(SetOfTiltSeriesM,
                                'tiltseries_set%s.sqlite', suffix)

    def importStep(self, pattern, voltage, sphericalAberration,
                         amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.getPattern()
        self.info("Using glob pattern: '%s'" % self._globPattern)
        self.info("Using regex pattern: '%s'" % self._regexPattern)

        self.info("Files: ")
        tiltSeriesDict = OrderedDict()
        for f, ts, to, ta in self.getMatchingFiles():
            if ts not in tiltSeriesDict:
                tiltSeriesDict[ts] = TiltSeriesM(filename=':memory:')

            tiltSeriesM = tiltSeriesDict[ts]
            tiltSeriesM.append(TiltImageM(location=f,
                                          acquisitionOrder=to,
                                          tiltAngle=ta))

        for ts, tiltSeriesM in tiltSeriesDict.iteritems():
            self.info("Tilt series: %s" % ts)
            for tim in tiltSeriesM.iterItems():#orderBy='_tiltAngle', direction='DESC'):
                self.info("   acq order: %02d, tilt angle: %0.2f"
                          % (tim.getAcquisitionOrder(), tim.getTiltAngle()))

        return

        createSetFunc = getattr(self, '_create' + self._outputClassName)
        imgSet = createSetFunc()
        imgSet.setIsPhaseFlipped(self.haveDataBeenPhaseFlipped.get())
        acquisition = imgSet.getAcquisition()

        self.fillAcquisition(acquisition)

        # Call a function that should be implemented by each subclass
        self.setSamplingRate(imgSet)

        outFiles = [imgSet.getFileName()]
        imgh = ImageHandler()
        img = imgSet.ITEM_TYPE()
        img.setAcquisition(acquisition)
        n = 1
        copyOrLink = self.getCopyOrLink()
        alreadyWarned = False  # Use this flag to warn only once

        for i, (fileName, fileId) in enumerate(self.iterFiles()):
            if self.isBlacklisted(fileName):
                continue
            uniqueFn = self._getUniqueFileName(fileName)
            dst = self._getExtraPath(uniqueFn)
            if ' ' in dst:
                if not alreadyWarned:
                    self.warning('Warning: your file names have white spaces!')
                    self.warning('Removing white spaces from copies/symlinks.')
                    alreadyWarned = True
                dst = dst.replace(' ', '')
            copyOrLink(fileName, dst)
            # Handle special case of Imagic images, copying also .img or .hed
            self.handleImgHed(copyOrLink, fileName, dst)

            if self._checkStacks:
                _, _, _, n = imgh.getDimensions(dst)

            if n > 1:
                for index in range(1, n + 1):
                    img.cleanObjId()
                    img.setMicId(fileId)
                    img.setFileName(dst)
                    img.setIndex(index)
                    self._addImageToSet(img, imgSet)
            else:
                img.setObjId(fileId)
                img.setFileName(dst)
                # Fill the micName if img is either Micrograph or Movie
                uniqueFn = uniqueFn.replace(' ', '')
                self._fillMicName(img, uniqueFn)
                self._addImageToSet(img, imgSet)

            outFiles.append(dst)

            sys.stdout.write("\rImported %d/%d\n" % (i + 1, self.numberOfFiles))
            sys.stdout.flush()

        print("\n")

        args = {}
        outputSet = self._getOutputName()
        args[outputSet] = imgSet
        self._defineOutputs(**args)

        return outFiles

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        self.getPattern()
        matching = self.getMatchingFiles()

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
    def getPattern(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
        """
        pattern = os.path.join(self.filesPath.get('').strip(),
                               self.filesPattern.get('').strip())

        def _replace(p, ts, to, ta):
            p = p.replace('{TS}', ts)
            p = p.replace('{TO}', to)
            p = p.replace('{TA}', ta)
            return p

        self._regexPattern = _replace(pattern.replace('*', '(.*)'),
                                      '(?P<TS>.*)', '(?P<TO>\d+)',
                                      '(?P<TA>[+-]?\d+(\.\d+)?)')
        self._regex = re.compile(self._regexPattern)
        self._globPattern = _replace(pattern,
                                     '*', '*', '*')

    def getMatchingFiles(self):
        """ Return a sorted list with the paths of files that
        matched the pattern.
        """
        self.getPattern()

        filePaths = glob(self._globPattern)
        filePaths.sort()

        matchingFiles = []
        for f in filePaths:
            m = self._regex.match(f)
            if m is not None:
                f += ' (matched)'
                matchingFiles.append((f, m.group('TS'),
                                      int(m.group('TO')),
                                      float(m.group('TA'))))

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
                                    % (fileName, self.getPattern()))

                fileId = int(match.group(1))

            else:
                fileId = None

            yield fileName, fileId

    def worksInStreaming(self):
        # Import protocols always work in streaming
        return True