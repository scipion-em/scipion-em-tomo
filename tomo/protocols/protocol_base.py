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
import logging
import re
from glob import glob
from os.path import join, getmtime

import pyworkflow as pw
from pyworkflow.protocol.params import (PointerParam, EnumParam, PathParam,
                                        FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.mapper.sqlite_db import SqliteDb
from pyworkflow.utils.properties import Message
from pwem.protocols import ProtImport, EMProtocol, ProtImportFiles

import tomo.objects
from tomo.constants import TS_LABEL
from tomo.convert.mdoc import normalizeTSId

logger = logging.getLogger(__name__)


class ProtTomoBase:
    def _createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be deleted. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pw.utils.cleanPath(setFn)

        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj

    def _createSetOfTiltSeriesM(self, suffix='') -> tomo.objects.SetOfTiltSeriesM:
        return self._createSet(tomo.objects.SetOfTiltSeriesM,
                               'tiltseriesM%s.sqlite', suffix)

    def _createSetOfTiltSeries(self, suffix='') -> tomo.objects.SetOfTiltSeries:
        self._ouputSuffix = ''
        return self._createSet(tomo.objects.SetOfTiltSeries,
                               'tiltseries%s.sqlite', suffix)

    def _createSetOfCoordinates3D(self, volSet, suffix='') -> tomo.objects.SetOfCoordinates3D:
        coord3DSet = self._createSet(tomo.objects.SetOfCoordinates3D,
                                     'coordinates%s.sqlite', suffix,
                                     indexes=['_volId'])
        coord3DSet.setPrecedents(volSet)
        return coord3DSet

    def _createSetOfTomograms(self, suffix='') -> tomo.objects.SetOfTomograms:
        return self._createSet(tomo.objects.SetOfTomograms,
                               'tomograms%s.sqlite', suffix)

    def _createSetOfSubTomograms(self, suffix='') -> tomo.objects.SetOfSubTomograms:
        return self._createSet(tomo.objects.SetOfSubTomograms,
                               'subtomograms%s.sqlite', suffix)

    def _createSetOfAverageSubTomograms(self, suffix='') -> tomo.objects.SetOfAverageSubTomograms:
        return self._createSet(tomo.objects.SetOfAverageSubTomograms,
                               'avgSubtomograms%s.sqlite', suffix)

    def _createSetOfClassesSubTomograms(self, subTomograms, suffix='') -> tomo.objects.SetOfClassesSubTomograms:
        classes = self._createSet(tomo.objects.SetOfClassesSubTomograms,
                                  'subtomogramClasses%s.sqlite', suffix)
        classes.setImages(subTomograms)

        return classes

    def _createSetOfLandmarkModels(self, suffix='') -> tomo.objects.SetOfLandmarkModels:
        return self._createSet(tomo.objects.SetOfLandmarkModels, 'setOfLandmarks%s.sqlite', suffix)

    def _createSetOfMeshes(self, volSet, suffix='') -> tomo.objects.SetOfMeshes:
        meshSet = self._createSet(tomo.objects.SetOfMeshes,
                                  'meshes%s.sqlite', suffix)
        meshSet.setPrecedents(volSet)
        return meshSet

    def _getOutputSuffix(self, cls):
        """ Get the name to be used for a new output.
        For example: output3DCoordinates7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(cls):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1  # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter + 1) if maxCounter > 0 else ''  # empty if not output


class ProtTomoPicking(ProtImport, ProtTomoBase):
    OUTPUT_PREFIX = 'output3DCoordinates'

    """ Base class for Tomogram boxing protocols. """

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputTomograms', PointerParam, label="Input Tomograms", important=True,
                      pointerClass='SetOfTomograms',
                      help='Select the Tomogram to be used during picking.')

    def _summary(self):
        summary = []
        if self.isFinished() and self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:*\n%s" % (key, output.getSummary()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary


class ProtTomoImportFiles(ProtImportFiles, ProtTomoBase):
    """Base protocol to import tomography files. There are two modes that should be implemented to import files,
    based on the value of the filesPattern form parameter:

    1. If empty (a single file is imported) or using the classic wildcard patterns.
    2. If the pattern contains the label {TS}, to represent the part of the name desired to be considered as the tsId.

    How to implement an import protocol in Scipion tomo:

    * The method initializeParsing must be called before the steps are generated, directly or as part of an
    initialization method. This method manages all the functionality required to deal with the {TS} pattern. If this
    label is not present, it does nothing, and the files matching must be carried out with the method iterFiles.

    * When processing the data, an if statement must be implemented, asking for a pattern {TS} introduced or not (it's
    the value of the protocol attribute self.regEx, which is filled properly in the execution of the method
    initializeParsing. Then:

        > if self.regEx: the iterator is obtained by calling the method getMatchingFilesFromRegEx
        > else: the iterator is obtained via the method iterFiles.

    * With the iterators defined, the rest of the code to generate the corresponding scipion objects is common to both
    modes.

    A good example is the protocol ProtImportTomograms.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.regEx = None
        self.regExPattern = None
        self.globPattern = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        self._defineImportParams(form)
        self._defineAcquisitionParams(form)

    @staticmethod
    def addExclusionWordsParam(form):
        form.addParam('exclusionWords', StringParam,
                      label='Exclusion words:',
                      expertLevel=LEVEL_ADVANCED,
                      help="List of words separated by a space that "
                           "the path should not have.")

    def _defineImportParams(self, form):
        """ Override to add options related to the different types
        of import that are allowed by each protocol.
        """
        importChoices = self._getImportChoices()

        form.addSection(label='Import')
        if len(importChoices) > 1:  # not only from files
            form.addParam('importFrom', EnumParam,
                          choices=importChoices, default=self._getDefaultChoice(),
                          label='Import from',
                          help='Select the type of import.')
        else:
            form.addHidden('importFrom', EnumParam,
                           choices=importChoices, default=self.IMPORT_FROM_FILES,
                           label='Import from',
                           help='Select the type of import.')

        form.addParam('filesPath', PathParam,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/Tomograms/data/day??_tomograms/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/Tomograms/data/day*_tomograms/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/Tomograms/data/day#_tomograms/\n"
                           "'#' represents one digit that will be used as "
                           "tomogram ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', StringParam,
                      label='Pattern',
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addLine('* NOTE: About files import in Scipion Tomo\n(click on the help icon to the right)',
                     help='There are two modes that should be implemented to import files, based on the value of '
                          'the Pattern parameter:\n'
                          '\n1. *Empty* (a single file is imported if File directory parameter is the full path of a '
                          'single file) or *using the classic wildcard patterns*, as in Scipion for SPA, and explain '
                          'in the help of the parameter Pattern.\n'
                          '\n*2. The pattern contains the label {TS}*, used to represent the part of the name desired '
                          'to be considered as the tsId (Unique identifier of a tilt-series and its derived products, '
                          'such as tomograms, ctfs, alignment, etc).\n\n'
                          '*Examples*:\n\n'
                          '*Example 1*: import only one tomogram named tomoAbc.mrc located if the directory '
                          '/home/exampleDir.\n'
                          '\n\tOption 1a: set Files directory to /home/exampleDir/tomoAbc.mrc and leave Pattern empty.\n'
                          '\n\tOption 1b: set set Files directory to /home/exampleDir and Pattern to tomoAbc.mrc.\n'
                          '\n\tOption 2: set set Files directory to /home/exampleDir and Pattern to the substring '
                          'desired to be the tsId:\n\n'
                          'If the tsId is desired to be tomoAbc, the pattern should be {TS}.mrc, while if the tsId we '
                          'want, for example to match other objects in the project labelled with tsId based on what '
                          'follows the word "tomo", it should be tomo{TS}.mrc. In the first case the tsId would be the '
                          'whole basename, "tomoAbc", while in the second it would be only "Abc".\n\n'
                          '*Example 2*: import alignment files (.xf) generated with IMOD, named TS_01_binned.xf, '
                          'TS_02_binned.xf, TS_03_binned.xf, located in the directory /home/aliExample/binnedTS.\n\n'
                          'Let us consider that we have already the tilt-series imported in Scipion and we know that '
                          'their tsIds are, respectively, TS_01, TS_02, TS_03. To import the alignment files we have '
                          'obtaiend using IMOD before having installed Scipion, the\n'
                          '\n\tFiles directory should be /home/aliExample/binnedTS, and the'
                          '\n\tPattern {TS}_binned.xf.\n\n'
                          'Thus, the transformation matrices would be imported with tsIds TS_01, TS_02, TS_03, that '
                          'matches to the existing tsIds of the tilt-series we want to assign the imported alignments '
                          'to.')
        form.addParam('copyFiles', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Copy files?",
                      help="By default the files are not copied into the "
                           "project to avoid data duplication and to save "
                           "disk space. Instead of copying, symbolic links are "
                           "created pointing to original files. This approach "
                           "has the drawback that if the project is moved to "
                           "another computer, the links need to be restored.")

    def _defineAcquisitionParams(self, form):
        """ Override to add options related to acquisition info.
        """
        form.addParam('samplingRate', FloatParam,
                      label=Message.LABEL_SAMP_RATE)

    # --------------------------- UTILS functions -----------------------------
    def initializeParsing(self):
        pattern = self.filesPattern.get()
        if pattern:
            if TS_LABEL in pattern:
                logger.info('Importing using a pattern.')
                path = self.filesPath.get().strip()
                pattern = pattern.strip()
                pattern = join(path, pattern)
                regExPattern = pattern.replace(TS_LABEL, r'(?P<TS>.*)')  # regex pattern for TS
                self.regEx = re.compile(regExPattern)
                self.regExPattern = regExPattern
                globPattern = pattern.replace(TS_LABEL, '*')
                # Glob module does not handle well the brackets (it does not list them)
                self.globPattern = globPattern.replace('[', '*').replace(']', '*')
            else:
                logger.info(f'Direct import. Pattern {TS_LABEL} not introduced.')

    def _excludeByWords(self, files):
        exclusionWords = self.exclusionWords.get()

        if exclusionWords is None:
            return files

        exclusionWordList = exclusionWords.split()

        allowedFiles = []

        for file in files:
            if any(bannedWord in file for bannedWord in exclusionWordList):
                logger.info("%s excluded. Contains any of %s" %
                            (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles

    def getMatchingFilesFromRegEx(self):
        filePaths = glob(self.globPattern)
        filePaths = self._excludeByWords(filePaths)
        filePaths.sort(key=lambda fn: getmtime(fn))
        matchingFilesDict = dict()
        for f in filePaths:
            matchRes = self.regEx.match(f)
            if matchRes is not None:
                tsId = matchRes.group('TS')  # Return the complete matched subgroup
                logger.info("Raw tilt series id is %s." % tsId)
                tsId = normalizeTSId(tsId)
                logger.info("Normalized tilt series id is %s." % tsId)
                matchingFilesDict[tsId] = f
        return matchingFilesDict


    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        pass


class ProtTomoSubtomogramAveraging(EMProtocol, ProtTomoBase):
    """ Base class for subtomogram averaging protocols. """
    pass


class ProtTomoImportAcquisition:
    MANUAL_IMPORT = 0
    FROM_FILE_IMPORT = 1

    def _defineParams(self, form):

        """ Override to add options related to acquisition info.
        """

        importAcquisitionChoices = ['Manual', 'From file']

        form.addSection(label='Acquisition Info')

        form.addParam('importAcquisitionFrom', EnumParam,
                      choices=importAcquisitionChoices,
                      default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('acquisitionData', PathParam,
                      label="Acquisition parameters file",
                      help="File with the acquisition parameters for each "
                           "tomogram or subtomogram to import. File must be in plain format."
                           " The file must contain a row per file to be imported "
                           "and have the following parameters in order: \n"
                           "\n"
                           "'File_name AcquisitionAngleMin AcquisitionAngleMax Step TiltAxisAngle' \n"
                           "\n"
                           "An example would be:\n"
                           "subtomo1.em -40 40 3 85\n"
                           "subtomo2.em -45 50 2 85\n",
                      condition="importAcquisitionFrom == %d" % self.FROM_FILE_IMPORT)

        form.addParam('acquisitionAngleMax', FloatParam,
                      default=60,
                      label='Acquisition angle max',
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      help='Enter the positive limit of the acquisition angle')

        form.addParam('acquisitionAngleMin', FloatParam,
                      default=-60,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Acquisition angle min',
                      help='Enter the negative limit of the acquisition angle')

        form.addParam('step', FloatParam,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Step',
                      help='Enter the step size for the import')

        form.addParam('tiltAxisAngle', FloatParam,
                      label='Tilt axis angle (deg.)',
                      allowsNull=True,
                      help="The rotation angle is the angle from the vertical "
                           "to the axis of tilting, where counterclockwise is "
                           "positive.\n See "
                           "https://bio3d.colorado.edu/imod/doc/tomoguide.html#UnknownAxisAngle")

        form.addParam('voltage', FloatParam, default=300,
                      label=Message.LABEL_VOLTAGE,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      help=Message.TEXT_VOLTAGE)

        form.addParam('sphericalAberration', FloatParam, default=2.7,
                      label=Message.LABEL_SPH_ABERRATION,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      help=Message.TEXT_SPH_ABERRATION)

        form.addParam('amplitudeContrast', FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      help=Message.TEXT_AMPLITUDE)

    def _parseAcquisitionData(self):
        if self.importAcquisitionFrom.get() == self.MANUAL_IMPORT:
            self.acquisitionParameters = {
                'angleMin': self.acquisitionAngleMin.get(),
                'angleMax': self.acquisitionAngleMax.get(),
                'step': self.step.get(),
                'tiltAxisAngle': self.tiltAxisAngle.get(),
                'voltage': self.voltage.get(),
                'sphericalAberration': self.sphericalAberration.get(),
                'amplitudeContrast': self.amplitudeContrast.get(),
            }
        else:
            params = open(self.acquisitionData.get(), "r")
            self.acquisitionParameters = {}
            for line in params.readlines():
                param = line.split()
                try:
                    self.acquisitionParameters.update({param[0]: {
                        'angleMin': float(param[1]),
                        'angleMax': float(param[2]),
                        'step': int(param[3]),
                        'tiltAxisAngle': float(param[4]),
                    }})
                except Exception as e:
                    print('Wrong acquisition data file format', e)

    def _extractAcquisitionParameters(self, fileName):
        if self.importAcquisitionFrom.get() == self.FROM_FILE_IMPORT:
            onlyName = fileName.split('/')[-1]
            acquisitionParams = self.acquisitionParameters[onlyName]
        else:
            acquisitionParams = self.acquisitionParameters

        return tomo.objects.TomoAcquisition(**acquisitionParams)

    def _summary(self, summary, setOfObject):
        for obj in setOfObject:
            if obj.hasAcquisition():
                summary.append(u"File: %s" % obj.getFileName())
                summary.append(u"Acquisition angle max: *%0.2f*" % obj.getAcquisition().getAngleMax())

                summary.append(u"Acquisition angle min: *%0.2f*" % obj.getAcquisition().getAngleMin())
                if obj.getAcquisition().getStep():
                    summary.append(u"Step: *%d*" % obj.getAcquisition().getStep())
                if obj.getAcquisition().getTiltAxisAngle():
                    summary.append(u"Tilt axis angle: *%0.2f*" % obj.getAcquisition().getTiltAxisAngle())
