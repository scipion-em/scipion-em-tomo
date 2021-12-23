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
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam
from pyworkflow.utils.path import removeBaseExt
from .protocol_base import ProtTomoImportFiles
from ..objects import TomoMask, SetOfTomoMasks


class ProtImportTomomasks(ProtTomoImportFiles):
    """Protocol to import a set of tomomasks (segmentations) to the project"""

    _outputClassName = 'SetOfTomoMasks'
    _label = 'import tomomasks (segmentations)'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.warnMsg = None
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
        self._defineOutputs(outputTomoMasks=tomoMaskSet)
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _getTomMessage(self):
        return "Tomomasks %s" % self.getObjectTag('outputTomoMasks')

    def _summary(self):
        try:
            summary = []
            if self.isFinished():
                summary.append("%s imported from:\n%s" % (self._getTomMessage(), self.getPattern()))
                summary.append(u"Sampling rate: *%0.2f* (Å/px)" % self.samplingRate.get())

                if self.warnMsg:
                    summary.append(self.warnMsg)

            return summary

        except Exception as e:
            print(e)

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getTomMessage(), self.samplingRate.get()), )
        return methods

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
        msgNonMatchingTomos = ''
        msgNonMatchingTomoMasks = ''
        inTomoSet = self.inputTomos.get()
        tomoIds, tomoBaseNames = zip(*[(tomo.getTsId(), removeBaseExt(tomo.getFileName())) for tomo in inTomoSet])

        def isMember(x, y):
            # x will be a list of ids and y the basename of the current tomomask. Thus, all the ids will be mapped to
            # check if one of them is contained in the current basename
            return x in y

        if self._tomoHasValidTsId(inTomoSet[1]):
            # Look for the tsId of each tomogram to be contained in the tomoFileNames
            list2check = tomoIds
        else:
            # The same, but considering the tomo basename instead of the tsId
            list2check = tomoBaseNames

        for file, _ in self.iterFiles():
            tomoMaskName = removeBaseExt(file).replace('_materials', '')
            matches = list(map(isMember, list2check, [tomoMaskName]))
            if any(matches):
                matchingTomoMaskDict[file] = inTomoSet[matches.index(True) + 1]
            else:
                nonMatchingTomoMaskNames.append(tomoMaskName)

        # Check if there are non-matching tomograms
        setMatchingIds = set(matchingTomoMaskDict.keys())
        pattern = '\t- {}\n'
        if setMatchingIds:
            set2check = set(list2check)
            nonMatchingTomogramIds = set2check ^ setMatchingIds  # ^ is (a | b) - (a & b) inverse intersection
            if nonMatchingTomogramIds:
                nOfNonMatchingTomos = len(nonMatchingTomogramIds)
                msgNonMatchingTomos = ('*[%i]* tomograms did not match to any of the tomomasks introduced:%s' %
                                       (nOfNonMatchingTomos, pattern * nOfNonMatchingTomos)).format(*nonMatchingTomogramIds)
        else:
            raise Exception('No matching tomograms were found with the input tomomasks. The match is carried out by '
                            'checking if the tsId of each tomogram is contained in the basename of the tomomaks. If '
                            'the tomograms do not have that attribute, then the base name is considered.')

        # The same for the non-matching tomomasks
        if nonMatchingTomoMaskNames:
            nOfNonMatchingTomomasks = len(nonMatchingTomoMaskNames)
            msgNonMatchingTomoMasks = ('*[%i]* tomomask did not match to any of the tomograms introduced:%s' %
                                   (nOfNonMatchingTomomasks, pattern * nOfNonMatchingTomomasks)).format(*nonMatchingTomoMaskNames)

        self.warnMsg = msgNonMatchingTomos + '\n\n' + msgNonMatchingTomoMasks
        self.matchingTomoMaskDict = matchingTomoMaskDict

    @staticmethod
    def _tomoHasValidTsId(tomo):
        return True if getattr(tomo, '_tsId', None) else False

    def _genOutputSetOfTomoMasks(self):
        tomoMaskSet = SetOfTomoMasks.create(self._getPath(), template='tomomasks%s.sqlite', suffix='annotated')
        inTomoSet = self.inputTomos.get()
        sRate = inTomoSet.getSamplingRate()
        tomoMaskSet.setSamplingRate(sRate)
        counter = 1
        for tomomaskFile, tomoObj in self.matchingTomoMaskDict.items():
            tomoMask = TomoMask()
            tomoMask.setSamplingRate(sRate)
            tomoMask.setLocation(counter, tomomaskFile)
            tomoMask.setVolName(tomoObj.getFileName())
            if self._tomoHasValidTsId(tomoObj):
                tomoMask.setTsId(tomoObj.getTsId())
            tomoMaskSet.append(tomoMask)
            counter += 1

        return tomoMaskSet
