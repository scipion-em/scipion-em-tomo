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
    """
    Import Tomograms (ProtImportTomograms) — User Manual

    Overview

    The Import Tomograms protocol imports reconstructed tomographic volumes into
    Scipion and converts them into a native `SetOfTomograms`.

    Its main purpose is to register tomographic volumes together with their spatial
    origin, sampling information, and acquisition metadata so they can be used in
    downstream cryo-electron tomography workflows.

    For a biological user, this protocol is typically the starting point of a
    tomographic analysis pipeline, since all subsequent particle picking,
    subtomogram extraction, visualization, and averaging depend on a correct
    tomogram spatial reference.

    Inputs and General Workflow

    The protocol requires one or more tomogram files as input.

    These files are imported from a directory or from a filename pattern.

    During execution, the protocol scans the input path, identifies matching files,
    and creates a `SetOfTomograms` in the Scipion project.

    Each imported tomogram receives:

    - a tilt-series identifier (`tsId`),
    - a sampling rate,
    - an origin definition,
    - acquisition metadata.

    The protocol does not duplicate the original tomogram data. Instead, it creates
    internal links to the source files.

    This makes the import lightweight while preserving direct access to the original
    reconstructed volumes.

    File Matching Strategy

    The protocol supports two file discovery modes.

    Direct pattern matching:
    - files are imported directly from the specified path and filename pattern.

    Regular-expression matching:
    - filenames are parsed through a regular expression and converted into
      normalized tilt-series identifiers.

    This flexibility is particularly useful in tomography facilities where
    reconstructed tomograms may come from different reconstruction pipelines and
    follow different naming conventions.

    Sampling Rate and Spatial Interpretation

    Every imported tomogram is assigned a sampling rate.

    This sampling rate defines the voxel size in Angstroms per pixel and becomes
    the spatial scale used by all downstream tomography protocols.

    From a biological perspective, correct sampling rate assignment is essential
    because all measurements, particle coordinates, and subtomogram extraction
    boxes depend directly on it.

    Acquisition Metadata

    The protocol can also import tomography acquisition metadata.

    Depending on the user configuration, acquisition information may be obtained:

    - from the input files, or
    - from manually defined acquisition parameters.

    When acquisition parameters are assigned at import time, the resulting tomogram
    set stores experimental metadata such as microscope conditions in a form that
    remains available for downstream processing.

    This is especially important in cryo-ET workflows where accurate acquisition
    metadata may later influence reconstruction interpretation or subtomogram
    analysis.

    Origin of Coordinates

    One of the most biologically important features of the protocol is the handling
    of tomogram origin coordinates.

    The spatial origin determines how all particle coordinates and extracted
    subtomograms are interpreted inside the tomographic volume.

    The protocol provides three possible origin strategies.

    Default geometric center:
    - if no manual origin is requested, the origin is placed at the geometric
      center of the tomogram.

    MRC header origin:
    - if enabled, the protocol reads origin information directly from the MRC
      header whenever the file is compatible.

    Manual origin:
    - the user may explicitly provide X, Y, and Z shifts in Angstroms.

    This flexibility is biologically important because different reconstruction
    software packages may define the tomogram origin differently.

    A wrong origin may not affect visualization immediately, but it can strongly
    impact downstream coordinate interpretation, particle extraction, symmetry
    operations, and subtomogram averaging.

    Origin Handling in Practice

    When origin information is requested from the MRC header:

    - the protocol reads the origin directly from the file metadata.

    If the file is not compatible with MRC origin metadata:

    - the protocol automatically falls back to the geometric center.

    When a manual origin is introduced:

    - the specified coordinates are converted into internal shifts.

    If no explicit origin is provided:

    - the protocol places the origin at the center of the tomogram volume.

    This ensures that every imported tomogram always receives a valid spatial
    reference frame.

    Tomogram Registration

    For each imported tomogram, the protocol performs several registration steps.

    It assigns:

    - the tomogram origin,
    - the tilt-series identifier,
    - acquisition metadata,
    - the internal file reference.

    The tomogram is then appended to the output set.

    This produces a Scipion-native tomogram object that can immediately be used by
    downstream tomography protocols.

    Output Generation

    The protocol generates a `SetOfTomograms`.

    The output set contains:

    - all successfully imported tomograms,
    - a common sampling rate,
    - optional acquisition metadata.

    Each tomogram preserves its own identity and spatial reference.

    This organization is especially important in experiments involving multiple
    tomograms, where each reconstructed volume must remain individually traceable.

    Validation of Input Consistency

    Before execution, the protocol checks whether matching tomogram files can be
    found.

    Validation depends on the selected matching mode.

    For direct pattern import:
    - at least one matching file must exist.

    For regular-expression import:
    - at least one filename must satisfy the regular-expression matching rule.

    If no files are found, the protocol stops before execution.

    This prevents the creation of empty tomogram sets.

    Practical Recommendations

    In routine cryo-electron tomography workflows, this protocol should be used
    with particular attention to spatial consistency.

    A few practical considerations are especially important:

    - Verify the sampling rate carefully before import.
    - Use MRC-header origin only when the reconstruction software is known to
      write correct origin metadata.
    - For symmetry-sensitive downstream analyses, verify that the origin preserves
      the expected geometric symmetry.
    - When importing tomograms from multiple reconstruction pipelines, confirm that
      tilt-series identifiers remain consistent.

    In most biological applications, correct origin assignment is often more
    important than the import itself.

    Final Perspective

    The Import Tomograms protocol is much more than a file-loading utility.

    It defines the spatial and experimental reference frame on which the entire
    tomographic workflow will depend.

    In modern cryo-ET pipelines, correct tomogram import is a foundational step
    because all downstream coordinate interpretation, particle extraction, and
    biological analysis inherit the spatial assumptions established here.
    """
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
