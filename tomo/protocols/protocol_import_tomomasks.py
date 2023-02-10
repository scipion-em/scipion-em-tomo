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
from pwem.emlib.image import ImageHandler
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import yellowStr
from pyworkflow.utils.path import removeBaseExt
from .protocol_base import ProtTomoImportFiles
from ..constants import ERR_NO_TOMOMASKS_GEN, ERR_NON_MATCHING_TOMOS
from ..objects import TomoMask, SetOfTomoMasks


class ProtImportTomomasks(ProtTomoImportFiles):
    """Protocol to import a set of tomomasks (segmentations) to the project"""

    _outputClassName = 'SetOfTomoMasks'
    _label = 'import tomomasks (segmentations)'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.warnMsg = ''
        self.matchingTomoMaskDict = None  # keys = filenames of matching tomomasks, values = matching Tomogram

    def _defineParams(self, form):
        ProtTomoImportFiles()._defineImportParams(form)
        form.addParam('inputTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Select the tomograms to be assigned to the input tomo masks.')

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.importStep)

    def importStep(self):
        self._checkMatchingFiles()
        tomoMaskSet = self._genOutputSetOfTomoMasks()
        if self.warnMsg:
            print(yellowStr('WARNING!') + '\n' + self.warnMsg)
        if not tomoMaskSet:
            raise Exception(ERR_NO_TOMOMASKS_GEN)

        self._defineOutputs(outputTomoMasks=tomoMaskSet)
        self._defineSourceRelation(self.inputTomos.get(), tomoMaskSet)

    # --------------------------- INFO functions ------------------------------
    def _getTomMessage(self):
        return "Tomomasks %s" % self.getObjectTag('outputTomoMasks')

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
        try:
            next(self.iterFiles())
        except StopIteration:
            errors.append('No files matching the pattern %s were found.' % self.getPattern())
        return errors

    # --------------------------- UTILS functions ------------------------------

    def _checkMatchingFiles(self):
        nonMatchingTomoMaskNames = []
        matchingTomoMaskDict = {}
        inTomoSet = self.inputTomos.get()
        tomoIds, tomoBaseNames = zip(*[(tomo.getTsId(), removeBaseExt(tomo.getFileName())) for tomo in inTomoSet])

        if self._tomoHasValidTsId(inTomoSet[1]):
            # Look for the tsId of each tomogram to be contained in the tomoFileNames
            list2check = tomoIds
        else:
            # The same, but considering the tomo basename instead of the tsId
            list2check = tomoBaseNames

        for file, _ in self.iterFiles():
            tomoMaskName = self.getTomoMaskName(file)
            matches = [True if tomo == tomoMaskName else False for tomo in list(list2check)]
            if any(matches):
                matchingTomoMaskDict[file] = inTomoSet[matches.index(True) + 1].clone()
            else:
                nonMatchingTomoMaskNames.append(tomoMaskName)

        # Check if there are non-matching tomograms
        matchingIds = matchingTomoMaskDict.keys()
        pattern = '\t- {}\n'
        if matchingIds:
            uniqueMatchingIds = set([self.getTomoMaskName(matchingId) for matchingId in matchingIds])
            set2check = set(list2check)
            nonMatchingTomogramIds = set2check ^ uniqueMatchingIds  # ^ is (a | b) - (a & b) inverse intersection
            if nonMatchingTomogramIds:
                nOfNonMatchingTomos = len(nonMatchingTomogramIds)
                headMsg = yellowStr('[%i] tomograms did not match to any of the tomomasks introduced:' % nOfNonMatchingTomos)
                msgNonMatchingTomos = ('%s\n%s' % (headMsg, pattern * nOfNonMatchingTomos)).format(*nonMatchingTomogramIds)
                self.warnMsg += msgNonMatchingTomos + '\n\n'
        else:
            raise Exception(ERR_NON_MATCHING_TOMOS)

        # The same for the non-matching tomomasks
        if nonMatchingTomoMaskNames:
            nOfNonMatchingTomomasks = len(nonMatchingTomoMaskNames)
            headMsg = yellowStr('[%i] tomomasks did not match to any of the tomograms introduced:' % nOfNonMatchingTomomasks)
            msgNonMatchingTomoMasks = ('%s\n%s' % (headMsg, pattern * nOfNonMatchingTomomasks)).format(*nonMatchingTomoMaskNames)
            self.warnMsg += msgNonMatchingTomoMasks + '\n\n'

        self.matchingTomoMaskDict = matchingTomoMaskDict

    @staticmethod
    def _tomoHasValidTsId(tomo):
        return True if getattr(tomo, '_tsId', None) else False

    def _genOutputSetOfTomoMasks(self):
        ih = ImageHandler()
        tomoMasksNonMatchingDims = []
        nonMatchingDimsList = []
        tomoMaskSet = SetOfTomoMasks.create(self._getPath(), template='tomomasks%s.sqlite', suffix='annotated')
        inTomoSet = self.inputTomos.get()
        xt, yt, zt = inTomoSet.getDimensions()
        sRate = inTomoSet.getSamplingRate()
        tomoMaskSet.setSamplingRate(sRate)
        counter = 1
        for tomomaskFile, tomoObj in self.matchingTomoMaskDict.items():
            x, y, z, _ = ih.getDimensions(tomomaskFile)
            if (xt, yt, zt) == (x, y, z):
                tomoMask = TomoMask()
                tomoMask.setSamplingRate(sRate)
                tomoMask.setLocation(counter, tomomaskFile)
                tomoMask.setVolName(tomoObj.getFileName())
                if self._tomoHasValidTsId(tomoObj):
                    tomoMask.setTsId(tomoObj.getTsId())
                tomoMaskSet.append(tomoMask)
                counter += 1
            else:
                tomoMasksNonMatchingDims.append(self.getTomoMaskName(tomomaskFile))
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

            self.warnMsg += msgNonMatchingDimsMasks + '\n\n'

        return tomoMaskSet

    @staticmethod
    def getTomoMaskName(maskFileName):
        return removeBaseExt(maskFileName).replace('_materials', '')  # This suffix is added by the memb. annotator
