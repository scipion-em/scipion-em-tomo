# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
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
from os.path import abspath, basename, join
from pwem.convert.headers import Ccp4Header
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.utils.path import createAbsLink, removeBaseExt, getExt
import pyworkflow.protocol.params as params
from .protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from ..convert.mdoc import normalizeTSId
from ..objects import Tomogram, SetOfTomograms

logger = logging.getLogger(__name__)
OUTPUT_NAME = 'Tomograms'


class ProtImportTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfTomograms'
    _label = 'import tomograms'
    _possibleOutputs = {OUTPUT_NAME: SetOfTomograms}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.Tomograms = None
        self.ih = None


    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)
        ProtTomoImportFiles.addExclusionWordsParam(form)

        ProtTomoImportAcquisition._defineParams(self, form)

        form.addSection('Origin Info')
        form.addParam('setOrigCoord', params.BooleanParam,
                      condition='importFrom == IMPORT_FROM_FILES',
                      label="Set origin of coordinates",
                      help="Option YES:\nA new volume will be created with "
                           "the "
                           "given ORIGIN of coordinates. This ORIGIN will be "
                           "set in the map file header.\nThe ORIGIN of "
                           "coordinates will be placed at the center of the "
                           "whole volume if you select n(x)/2, n(y)/2, "
                           "n(z)/2 as "
                           "x, y, z coordinates (n(x), n(y), n(z) are the "
                           "dimensions of the whole volume). However, "
                           "selecting "
                           "0, 0, 0 as x, y, z coordinates, the volume will be "
                           "placed at the upper right-hand corner.\n\n"
                           "Option NO:\nThe ORIGIN of coordinates will be "
                           "placed at the center of the whole volume ("
                           "coordinates n(x)/2, n(y)/2, n(z)/2 by default). "
                           "This "
                           "ORIGIN will NOT be set in the map file header.\n\n"
                           "WARNING: In case you want to process "
                           "the volume with programs requiring a specific "
                           "symmetry regarding the origin of coordinates, "
                           "for example the protocol extract unit "
                           "cell, check carefully that the coordinates of the "
                           "origin preserve the symmetry of the whole volume. "
                           "This is particularly relevant for loading "
                           "fragments/subunits of the whole volume.\n",
                      default=False)

        form.addBooleanParam('fromMrcHeader', label='From mrc header',
                             help='Use origin information in mrc headers of the tomograms.',
                             default=False, condition='setOrigCoord', )

        form.addLine('Manual offset',
                     help="A wizard will suggest you possible "
                          "coordinates for the ORIGIN. In MRC volume "
                          "files, the ORIGIN coordinates will be "
                          "obtained from the file header.\n "
                          "In case you prefer set your own ORIGIN "
                          "coordinates, write them here. You have to "
                          "provide the map center coordinates in "
                          "Angstroms (pixels x sampling).\n",
                     condition='setOrigCoord and not fromMrcHeader')
        # line.addParam would produce a nicer looking form
        # but them the wizard icon is drawn outside the visible
        # window. Until this bug is fixed form is a better option
        form.addParam('x', params.FloatParam, condition='setOrigCoord and not fromMrcHeader',
                      label="x", help="offset along x axis (Angstroms)")
        form.addParam('y', params.FloatParam, condition='setOrigCoord and not fromMrcHeader',
                      label="y", help="offset along y axis (Angstroms)")
        form.addParam('z', params.FloatParam, condition='setOrigCoord and not fromMrcHeader',
                      label="z", help="offset along z axis (Angstroms)")

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.importTomogramsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.ih = ImageHandler()
        self.initializeParsing()

    def importTomogramsStep(self):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        samplingRate = self.samplingRate.get()

        # Create a Volume template object
        tomo = Tomogram()
        tomo.setSamplingRate(samplingRate)
        tomoSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        tomoSet.setSamplingRate(samplingRate)

        self._parseAcquisitionData()
        if self.importAcquisitionFrom.get() != self.FROM_FILE_IMPORT:
            tomoSet.setAcquisition(self._extractAcquisitionParameters(None))

        if self.regEx:
            logger.info("Using regex pattern: '%s'" % self.regExPattern)
            logger.info("Generated glob pattern: '%s'" % self.globPattern)
            for tsId, fileName in self.getMatchingFilesFromRegEx().items():
                self.addTomoToSet(fileName, tsId, tomo, tomoSet)
        else:
            inPattern = self.filesPattern.get()
            pattern = inPattern.strip() if inPattern else ''
            logger.info("Using direct pattern: '%s'" % join(self.filesPath.get().strip(), pattern))
            filePaths = [fileName[0] for fileName in self.iterFiles()]
            fileList = self._excludeByWords(filePaths)
            for fileName in fileList:
                tsId = normalizeTSId(removeBaseExt(fileName))
                self.addTomoToSet(fileName, tsId, tomo, tomoSet)

        self._defineOutputs(**{OUTPUT_NAME: tomoSet})

    # --------------------------- UTILS functions ------------------------------
    def _getOrigCoord(self):
        return -1. * self.x.get(), -1. * self.y.get(), -1. * self.z.get()

    def setDefaultOrigin(self, fileName, origin):
        samplingRate = self.samplingRate.get()
        x, y, z, n = self.ih.getDimensions(fileName)
        origin.setShifts(x / -2. * samplingRate,
                         y / -2. * samplingRate,
                         z / -2. * samplingRate)

    def getTomoNewFileName(self, tsId, ext):
        return self._getExtraPath(f'{tsId}{ext}')

    def addTomoToSet(self, fileName: str, tsId: str, tomoObj: Tomogram, tomoSet: SetOfTomograms) -> None:
        origin = Transform()
        if self.setOrigCoord.get():
            if self.fromMrcHeader.get():
                if Ccp4Header.isCompatible(fileName):
                    ccp4Header = Ccp4Header(fileName, readHeader=True)
                    origin.setShiftsTuple(ccp4Header.getOrigin())
                else:
                    logger.info("File %s not compatible with mrc format. Setting default origin: geometrical center "
                                "of it." % fileName)
                    self.setDefaultOrigin(fileName, origin)
            else:
                origin.setShiftsTuple(self._getOrigCoord())
        else:
            self.setDefaultOrigin(fileName, origin)

        tomoObj.setOrigin(origin)
        tomoObj.setTsId(tsId)
        newFileName = self.getTomoNewFileName(tsId, getExt(fileName))
        createAbsLink(abspath(fileName), abspath(newFileName))
        tomoObj.setAcquisition(self._extractAcquisitionParameters(fileName))
        tomoObj.cleanObjId()
        tomoObj.setFileName(newFileName)
        tomoSet.append(tomoObj)
        tomoSet.update(tomoObj)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return self.Tomograms is not None

    def _getTomMessage(self):
        return "Tomograms %s" % self.getObjectTag(OUTPUT_NAME)

    def _summary(self):
        summary = []
        try:
            if self._hasOutput():
                summary.append("%s imported from:\n%s"
                               % (self._getTomMessage(), self.getPattern()))

                if self.samplingRate.get():
                    summary.append(u"Sampling rate: *%0.2f* (Å/px)" % self.samplingRate.get())

                ProtTomoImportAcquisition._summary(self, summary, self.Tomograms)

                x, y, z = self.Tomograms.getFirstItem().getShiftsFromOrigin()
                summary.append(u"Tomograms Origin (x,y,z):\n"
                               u"    x: *%0.2f* (Å/px)\n"
                               u"    y: *%0.2f* (Å/px)\n"
                               u"    z: *%0.2f* (Å/px)" % (x, y, z))

        except Exception as e:
            print(e)

        return summary

    def _methods(self):
        methods = []
        if self._hasOutput():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getTomMessage(), self.samplingRate.get()), )
        return methods

    def _getVolumeFileName(self, fileName, extension=None):
        if extension is not None:
            baseFileName = "import_" + str(basename(fileName)).split(".")[0] + ".%s" % extension
        else:
            baseFileName = "import_" + str(basename(fileName)).split(":")[0]

        return self._getExtraPath(baseFileName)

    def _validate(self):
        errors = []
        self._initialize()
        if self.regEx:
            matchingFileDict = self.getMatchingFilesFromRegEx()
            if not matchingFileDict:
                errors.append('No files matching the pattern %s were found.' % self.globPattern)
        else:
            try:
                next(self.iterFiles())
            except StopIteration:
                errors.append('No files matching the pattern %s were found.' % self.getPattern())
        return errors
