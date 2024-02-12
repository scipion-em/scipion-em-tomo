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

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.plugin import Domain
from pyworkflow.utils import removeBaseExt
from ..convert.mdoc import normalizeTSId

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
    """Common protocol to import CTF estimation of a tilt-series. """
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
        super()._defineImportParams(form)
        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that "
                           "the path should not have",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      help='Select the corresponding tilt-series for '
                           'which you want to update the CTF parameters.')

    # --------------------------- INSERT functions ----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.importCTFStep)

    # --------------------------- STEPS functions -----------------------------
    def importCTFStep(self):
        ci = self.getImportClass()

        outputCtfs = self._getOutputSet()
        if outputCtfs is None:
            outputCtfs = self._createOutputSet()

        defocusFiles = [fn for fn in self.iterFiles()]
        for ts in self._getInputTs():
            tsId = ts.getTsId()
            tsObjId = ts.getObjId()

            for defocusFn in defocusFiles:
                if tsId == normalizeTSId(removeBaseExt(defocusFn).replace('_ctf', '').replace('_avrot', '')):
                    logger.info("Parsing file: " + defocusFn)

                    newCTFTomoSeries = CTFTomoSeries()
                    newCTFTomoSeries.copyInfo(ts)
                    newCTFTomoSeries.setTiltSeries(ts)
                    newCTFTomoSeries.setObjId(tsObjId)
                    newCTFTomoSeries.setTsId(tsId)

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
                    defocusFiles.remove(defocusFn)
                    break

        outputCtfs.write()
        self._store()

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
        matchingFiles = self.getMatchFiles()
        if matchingFiles:
            tsIdList = self._getInputTs().getUniqueValues(TiltSeries.TS_ID_FIELD)
            defocusBNames = [normalizeTSId(pwutils.removeBaseExt(defocusFn).replace('_ctf', '').replace('_avrot', ''))
                             for defocusFn in self.iterFiles()]
            matchResults = list(set(tsIdList) & set(defocusBNames))
            if not matchResults:
                errorMsg.append('No matching files found.\n'
                                'CTF filenames are expected to include tsId e.g. tsId_ctf_avrot.txt. The suffixes '
                                '"_ctf" or "_avrot" are not mandatory.\n'
                                'The tsIds detected in the tilt series introduced are:\n'
                                '%s\n'
                                'The defocus files base names detected are (excluding the suffixes "_ctf" and '
                                '"_avrot"):\n'
                                '%s' % (tsIdList, defocusBNames))
        else:
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
            importFunc = Domain.importFromPlugin('imod.protocols',
                                                 'ProtImodBase',
                                                 doRaise=True)
            return importFunc()
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

    def _excludeByWords(self, files):
        exclusionWords = self.exclusionWords.get()

        if exclusionWords is None:
            return files

        exclusionWordList = exclusionWords.split()
        allowedFiles = []

        for file in files:
            if any(bannedWord in file for bannedWord in exclusionWordList):
                print("%s excluded. Contains any of %s" % (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles
