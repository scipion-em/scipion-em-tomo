# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
import logging
from enum import Enum
from os.path import abspath, join
from pwem.emlib.image import ImageHandler
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import yellowStr
from pyworkflow.utils.path import removeBaseExt, getExt, createAbsLink
from .protocol_base import ProtTomoImportFiles
from ..constants import ERR_NO_TOMOMASKS_GEN
from ..convert.mdoc import normalizeTSId
from ..objects import TomoMask, SetOfTomoMasks


logger = logging.getLogger(__name__)


class importTomoMasksOutputs(Enum):
    tomomasks = SetOfTomoMasks


class ProtImportTomomasks(ProtTomoImportFiles):
    """Protocol to import a set of tomomasks (segmentations) to the project"""

    _label = 'import tomomasks (segmentations)'
    _possibleOutputs = importTomoMasksOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.warnMsg = String('')
        self.matchingTomoMaskDict = None  # keys = filenames of matching tomomasks, values = matching Tomogram

    def _defineParams(self, form):
        ProtTomoImportFiles()._defineImportParams(form)
        ProtTomoImportFiles.addExclusionWordsParam(form)
        form.addParam('inputTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Select the tomograms to be assigned to the input tomo masks.')

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.importStep)

    def _initialize(self):
        self.ih = ImageHandler()
        self.initializeParsing()

    def importStep(self):
        inTomos = self.inputTomos.get()
        samplingRate = inTomos.getSamplingRate()
        tomoMask = TomoMask()
        tomoMask.setSamplingRate(samplingRate)

        if self.regEx:
            logger.info("Using regex pattern: '%s'" % self.regExPattern)
            logger.info("Generated glob pattern: '%s'" % self.globPattern)
            filesDict = self.getMatchingFilesFromRegEx()
        else:
            inPattern = self.filesPattern.get()
            pattern = inPattern.strip() if inPattern else ''
            logger.info("Using direct pattern: '%s'" % join(self.filesPath.get().strip(), pattern))
            filePaths = [fileName[0] for fileName in self.iterFiles()]
            fileList = self._excludeByWords(filePaths)
            filesDict = {}
            for fileName in fileList:
                tsId = normalizeTSId(self.getTomoMaskName(fileName))
                filesDict[tsId] = fileName

        tomoDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos}
        tomoMaskSet = self._genOutputSetOfTomoMasks(filesDict, tomoDict)

        warnMsg = self.warnMsg
        if warnMsg.get():
            self._store(warnMsg)
            print(yellowStr('WARNING!') + '\n' + warnMsg.get())
        if len(tomoMaskSet) == 0:
            raise Exception(ERR_NO_TOMOMASKS_GEN)

        self._defineOutputs(**{self._possibleOutputs.tomomasks.name: tomoMaskSet})
        self._defineSourceRelation(self.inputTomos, tomoMaskSet)

    # --------------------------- INFO functions ------------------------------
    def _getTomMessage(self):
        return "Tomomasks %s" % self.getObjectTag(self._possibleOutputs.tomomasks.name)

    def _summary(self):
        try:
            summary = []
            if self.isFinished():
                summary.append("%s imported from:\n%s" % (self._getTomMessage(), self.getPattern()))

                if self.warnMsg:
                    summary.append('Some tomograms or tomomasks were excluded. Check the log for more details.')

            return summary

        except Exception as e:
            print(e)

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

    # --------------------------- UTILS functions ------------------------------
    @staticmethod
    def _tomoHasValidTsId(tomo):
        return True if getattr(tomo, '_tsId', None) else False

    def _genOutputSetOfTomoMasks(self, filesDict, tomoDict):
        tomoMasksNonMatchingDims = []
        nonMatchingDimsList = []
        tomoMaskSet = SetOfTomoMasks.create(self._getPath(), template='tomomasks%s.sqlite')
        inTomoSet = self.inputTomos.get()
        sRate = inTomoSet.getSamplingRate()
        tomoMaskSet.setSamplingRate(sRate)
        counter = 1
        for tsId in filesDict.keys():
            tomoMaskFile = filesDict[tsId]
            tomo = tomoDict.get(tsId, None)
            if tomo:
                x, y, z, _ = self.ih.getDimensions(tomoMaskFile)
                xt, yt, zt, _ = self.ih.getDimensions(tomo.getFileName())
                if (xt, yt, zt) == (x, y, z):
                    tomoMask = TomoMask()
                    tomoMask.setTsId(tsId)
                    tomoMask.setSamplingRate(sRate)
                    newFileName = self.getTomoMaskNewFileName(tsId, getExt(tomoMaskFile))
                    createAbsLink(abspath(tomoMaskFile), abspath(newFileName))
                    tomoMask.setFileName(tomoMaskFile)
                    tomoMask.setVolName(tomo.getFileName())
                    tomoMaskSet.append(tomoMask)
                    counter += 1
                else:
                    tomoMasksNonMatchingDims.append(self.getTomoMaskName(tomoMaskFile))
                    nonMatchingDimsList.append((x, y, z))

        if tomoMasksNonMatchingDims:
            nNonMatchingDimsMasks = len(tomoMasksNonMatchingDims)
            msgNonMatchingDimsMasks = yellowStr('[%i] tomomasks have different dimensions than the ones from the '
                                                'introduced set of tomograms (x, y, z) = (%i, %i, %i):' %
                                                (nNonMatchingDimsMasks, xt, yt, zt))

            for nonMatchingDimsMask, nonMatchingDims in zip(tomoMasksNonMatchingDims, nonMatchingDimsList):
                msgNonMatchingDimsMasks += '\n\t- %s (x, y, z) = (%i, %i, %i)\n' % \
                                           (nonMatchingDimsMask,
                                            nonMatchingDims[0],
                                            nonMatchingDims[1],
                                            nonMatchingDims[2])

            self.warnMsg.set(self.warnMsg.get() + msgNonMatchingDimsMasks + '\n\n')

        return tomoMaskSet

    def addTomoMaskToSet(self, fileName: str, tsId: str, tomoMask: TomoMask, tomoMaskSet: SetOfTomoMasks) -> None:
        tomoMask.setTsId(tsId)
        newFileName = self.getTomoMaskNewFileName(tsId, getExt(fileName))
        createAbsLink(abspath(fileName), abspath(newFileName))
        tomoMask.cleanObjId()
        tomoMask.setFileName(newFileName)
        tomoMaskSet.append(tomoMask)
        tomoMaskSet.update(tomoMask)

    def getTomoMaskNewFileName(self, tsId: str, ext: str) -> str:
        return self._getExtraPath(f'{tsId}{ext}')

    @staticmethod
    def getTomoMaskName(maskFileName):
        # These suffixes are added by the memb. annotator and membrain, respectively
        return removeBaseExt(maskFileName).replace('_materials', '').replace('_segmented', '')

