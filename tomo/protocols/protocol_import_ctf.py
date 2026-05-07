# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from enum import Enum, IntEnum
import logging
from os.path import join

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pwem import fnMatching
from pyworkflow.plugin import Domain
from pyworkflow.utils import removeBaseExt
from ..convert.mdoc import normalizeTSId, deNormalizeTSId

from ..objects import SetOfCTFTomoSeries, TiltSeries, CTFTomoSeries
from .protocol_base import ProtTomoImportFiles

logger = logging.getLogger(__name__)


class ImportChoice(IntEnum):
    CTFFIND = 0
    IMOD = 1
    GCTF = 2
    ARETOMO = 3

class outputs(Enum):
    CTFs = SetOfCTFTomoSeries

class ProtImportTsCTF(ProtTomoImportFiles):
    """
    Import Tomography CTFs (ProtImportTsCTF) — User Manual

    Overview

    The Import Tomography CTFs protocol imports contrast transfer function (CTF)
    estimations for tilt series from external tomography software into Scipion.

    Its main purpose is to associate externally estimated defocus and CTF parameters
    with previously imported tilt series so they can be used in downstream cryo-electron
    tomography workflows.

    For a biological user, this protocol provides the link between raw tilt-series
    acquisition and CTF-aware tomographic processing, where accurate defocus information
    becomes essential for reliable reconstruction, subtomogram extraction, and
    subsequent biological interpretation.

    Inputs and General Workflow

    The protocol requires two principal inputs:

    - A `SetOfTiltSeries`, representing the tomographic tilt series already available
      in the project.
    - A collection of external CTF estimation files.

    During execution, the protocol attempts to associate each tilt series with its
    corresponding defocus or CTF file.

    Once a match is found, the selected external parser reads the defocus information
    and creates a corresponding `CTFTomoSeries` object linked to that tilt series.

    This workflow preserves the identity of each tilt series and ensures that CTF
    information remains associated with the correct tomographic dataset.

    Supported Import Formats

    The protocol supports several external tomography packages.

    Available import modes include:

    - CTFFIND
    - IMOD
    - GCTF
    - AreTomo

    However, GCTF import is currently not supported for execution, even though it
    appears as an available option.

    This design reflects the fact that tomography users often estimate CTF parameters
    using different software environments before integrating results into Scipion.

    File Matching Strategy

    A key part of the protocol is the association between tilt series and defocus files.

    Two matching modes are supported:

    Direct filename matching:
    - files are matched against tilt-series identifiers using normalized names.

    Regular-expression matching:
    - files are matched using a user-defined pattern that extracts tilt-series
      identifiers.

    This flexibility is especially useful in facility-scale tomography pipelines where
    naming conventions may differ across acquisition or preprocessing systems.

    If no file can be matched to a tilt series, the protocol reports the missing
    association.

    CTF Import and Parsing

    For every successfully matched tilt series, the protocol creates a new
    `CTFTomoSeries`.

    The output object inherits the identity and metadata of the original tilt series.

    Then, the selected import parser reads the external defocus file and populates the
    CTF information for that tilt series.

    From a biological perspective, this step transfers optical information estimated
    outside Scipion into the native tomographic framework, enabling downstream
    correction-aware processing.

    IMOD-Specific Handling

    When importing from IMOD, the protocol creates additional IMOD-specific metadata.

    These internal fields are initialized before parsing and updated later during
    import.

    This specialized handling reflects the fact that different tomography software
    packages may encode defocus information differently.

    Output Generation

    The protocol generates a `SetOfCTFTomoSeries` as output.

    For each successfully imported tilt series, the output contains:

    - the original tilt-series identity,
    - the associated CTF estimations,
    - package-specific metadata when relevant.

    The output set maintains a direct relation to the input tilt-series set.

    This ensures compatibility with downstream tomography protocols requiring CTF-aware
    processing.

    Quality Control and Automatic Disabling

    After importing CTF information, the protocol performs a quality consistency check.

    It evaluates whether the imported defocus estimations fall within acceptable
    deviation ranges.

    If the imported CTF values fail this internal consistency check:

    - the corresponding `CTFTomoSeries` is automatically disabled.

    This is biologically important because highly inconsistent defocus estimations may
    indicate failed fitting, unreliable tilt images, or acquisition artifacts.

    Validation of Input Consistency

    Before execution, the protocol validates that input files can actually be found.

    Validation checks include:

    - whether files matching the requested pattern exist;
    - whether regular-expression matching produces valid candidates;
    - whether the selected import mode is currently supported.

    If no matching files are found, the protocol stops before execution.

    This prevents accidental generation of incomplete or meaningless CTF metadata.

    Output Interpretation

    After completion, the summary reports:

    - the number of imported CTF tomography series;
    - which input tilt series could not be matched to external CTF files.

    This feedback is particularly useful in large tomography datasets where some tilt
    series may be missing corresponding CTF estimations.

    Practical Recommendations

    In routine cryo-electron tomography workflows, this protocol is most useful when
    CTF estimation has been performed outside Scipion.

    A few practical considerations are important:

    - Keep tilt-series identifiers and defocus filenames as consistent as possible.
    - Use regular-expression matching when external software introduces complex naming
      conventions.
    - Carefully inspect unmatched tilt series after import.
    - Pay attention to automatically disabled CTF series, since they may indicate
      problematic defocus estimation.

    In most biological applications, correct association between tilt series and CTF
    files is more important than the specific import format.

    Final Perspective

    The Import Tomography CTFs protocol is a key interoperability tool between external
    CTF estimation software and Scipion tomography workflows.

    Rather than simply reading defocus files, it preserves tilt-series identity,
    validates optical consistency, and integrates CTF metadata into the native
    tomographic framework.

    In modern cryo-ET pipelines, this protocol provides an essential preparatory step
    before reconstruction, subtomogram analysis, and biologically meaningful
    interpretation.
    """
    _possibleOutputs = outputs
    _label = 'import tomo CTFs'

    def _getImportChoices(self):
        """ Return a list of possible choices from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.)
        """
        return [ImportChoice.CTFFIND.name,
                ImportChoice.IMOD.name,
                ImportChoice.GCTF.name,
                ImportChoice.ARETOMO.name]

    def _getDefaultChoice(self):
        return 0

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        self._defineImportParams(form)
        self.addExclusionWordsParam(form)
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      help='Select the corresponding tilt-series for '
                           'which you want to update the CTF parameters.')

    # --------------------------- INSERT functions ----------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.importCTFStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.initializeParsing()

    def importCTFStep(self):
        ci = self.getImportClass()
        atLeastOneMatch = False

        defocusFileDict = {}
        if self.regEx:
            logger.info("Using regex pattern: '%s'" % self.regExPattern)
            logger.info("Generated glob pattern: '%s'" % self.globPattern)
            tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inputSetOfTiltSeries.get()}
            defocusFilesDict = self.getMatchingFilesFromRegEx()
            for tsId, ts in tsDict.items():
                defocusFile = defocusFilesDict.get(tsId, None)
                if defocusFile:
                    atLeastOneMatch = True
                    self.genCtfSet(ts, defocusFile, ci)
                else:
                    logger.warning(f'tsId = {tsId}: No defocus file was found.')

        else:
            logger.info("Using direct pattern: '%s'" % join(self.filesPath.get().strip(),
                                                            self.filesPattern.get().strip()))
            defocusFiles = [fn for fn in self.iterFiles()]
            for defocusFn in defocusFiles:
                defocusFileDict[removeBaseExt(defocusFn)] = defocusFn

            for ts in self._getInputTs():
                tsId = ts.getTsId()
                _, defocusFn = fnMatching(deNormalizeTSId(tsId), defocusFileDict, objType='Tilt-series')
                if defocusFn is not None:
                    atLeastOneMatch = True
                    self.genCtfSet(ts, defocusFn, ci)
                else:
                    logger.warning(f'tsId = {tsId}: No defocus file was found.')

        if not atLeastOneMatch:
            raise Exception('There are no files that match any of the Tilt series')

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        outSet = self._getOutputSet()
        if outSet:
            summary.append("Imported CTF tomo series: *%i*" % (self._getOutputSet().getSize()))
            inTsIdSet = set(self._getInputTs().getUniqueValues(TiltSeries.TS_ID_FIELD))
            outTsIdSet = set(outSet.getUniqueValues(TiltSeries.TS_ID_FIELD))
            nonCommonOutTsIds = inTsIdSet.difference(outTsIdSet)
            if nonCommonOutTsIds:
                summary.append(f'Some input TS did not match the imported files:\n'
                               f'\t*Non-matching tsIds = {nonCommonOutTsIds}*')
        else:
            summary.append("Output CTFs not ready yet.")
        return summary

    def _validate(self):
        errorMsg = []
        self._initialize()
        if self.regEx:
            matchingFileDict = self.getMatchingFilesFromRegEx()
            if not matchingFileDict:
                errorMsg.append('No files matching the pattern %s were found.' % self.globPattern)
        else:
            matchingFiles = self.getMatchFiles()
            if not matchingFiles:
                errorMsg.append('Unable to find the files provided:\n\n'
                                '\t-filePath = %s\n'
                                '\t-pattern = %s\n' % (self.filesPath.get(),
                                                       self.filesPattern.get()))

        if self.importFrom.get() == ImportChoice.GCTF.value:
            errorMsg.append("Import from GCTF is not supported yet.")

        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _getInputTs(self, pointer=False):
        return (self.inputSetOfTiltSeries.get() if not pointer
                else self.inputSetOfTiltSeries)

    def _createOutputSet(self):
        outputCtfs = SetOfCTFTomoSeries.create(self._getPath(),
                                               template='CTFmodels%s.sqlite')
        outputCtfs.setSetOfTiltSeries(self._getInputTs(pointer=True))
        self._defineOutputs(**{outputs.CTFs.name: outputCtfs})
        return outputCtfs

    def _getOutputSet(self):
        return getattr(self, self._possibleOutputs.CTFs.name, None)

    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        importFrom = self.importFrom.get()

        if importFrom == ImportChoice.IMOD.value:
            importFunc = Domain.importFromPlugin('imod.convert',
                                                 'ImodCtfParser',
                                                 doRaise=True)
        elif importFrom == ImportChoice.CTFFIND.value:
            importFunc = Domain.importFromPlugin('cistem.convert',
                                                 'GrigorieffLabImportCTF',
                                                 doRaise=True)
        elif importFrom == ImportChoice.GCTF.value:
            importFunc = Domain.importFromPlugin('gctf.convert',
                                                 'GctfImportCTF',
                                                 doRaise=True)
        elif importFrom == ImportChoice.ARETOMO.value:
            importFunc = Domain.importFromPlugin('aretomo.convert',
                                                 'AretomoCtfParser',
                                                 doRaise=True)
        else:
            importFunc = None

        return importFunc(self) or None

    def _importFromImod(self):
        return self.importFrom.get() == ImportChoice.IMOD.value

    def iterFiles(self):
        """ Iterate through the files matched with the pattern. """
        filePaths = self.getMatchFiles()
        filePaths = self._excludeByWords(filePaths)

        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise ValueError("File '%s' doesn't match the pattern '%s'"
                                     % (fileName, self.getPattern()))

            yield fileName

    def genCtfSet(self, ts, defocusFn, ci):
        outputCtfs = self._getOutputSet()
        if outputCtfs is None:
            outputCtfs = self._createOutputSet()

        newCTFTomoSeries = CTFTomoSeries()
        newCTFTomoSeries.copyInfo(ts)
        newCTFTomoSeries.setTiltSeries(ts)
        newCTFTomoSeries.setTsId(ts.getTsId())

        if self._importFromImod():
            # Create IMOD-specific attrs that will be updated later
            newCTFTomoSeries.setIMODDefocusFileFlag(None)
            newCTFTomoSeries.setNumberOfEstimationsInRange(None)

        outputCtfs.append(newCTFTomoSeries)

        ci.parseTSDefocusFile(ts, defocusFn, newCTFTomoSeries)

        if not (newCTFTomoSeries.getIsDefocusUDeviationInRange() and
                newCTFTomoSeries.getIsDefocusVDeviationInRange()):
            newCTFTomoSeries.setEnabled(False)

        outputCtfs.update(newCTFTomoSeries)
        outputCtfs.write()
        self._store()