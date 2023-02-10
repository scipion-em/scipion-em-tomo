# **************************************************************************
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from enum import Enum
from os import symlink
from os.path import basename, abspath
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from tomo.objects import TomoMask, SetOfTomoMasks
from tomo.utils import isMatchingByTsId


class assignT2TMaskOutputs(Enum):
    tomoMasks = SetOfTomoMasks


class ProtAssignTomo2TomoMask(EMProtocol):
    """ This protocol assign tomograms to tomomasks (segmentations)."""

    _devStatus = BETA
    _label = 'assign tomograms to tomo masks (segmentations)'
    _possibleOutputs = assignT2TMaskOutputs
    MATERIALS_SUFFIX = '_materials'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inTomoMasks', PointerParam,
                      pointerClass='SetOfTomoMasks',
                      label='Tomo masks (segmentations)',
                      help='Select the tomo masks desired to be referred to the introduced tomograms. The match '
                           'between both sets is carried out firstly by tsId and if not possible, then it will try '
                           'to do it by filename.')
        form.addParam('inputTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Select the tomograms to be assigned to the input tomo masks.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------

    def createOutputStep(self):
        inTomos = self.inputTomos.get()
        inTomoMasks = self.inTomoMasks.get()
        outputSetOfTomoMasks = inTomoMasks.create(self._getPath(), template='tomomasks%s.sqlite')
        outputSetOfTomoMasks.setSamplingRate(inTomos.getSamplingRate())

        if isMatchingByTsId(self.inTomoMasks.get(), self.inputTomos.get()):
            tomoTsIds = [tomo.getTsId() for tomo in inTomos]
            for inTomoMask in inTomoMasks:
                outTomoMask = self.setMatchingTomogram(tomoTsIds, inTomoMask, inTomos)
                outputSetOfTomoMasks.append(outTomoMask)
        else:
            # Membrane Annotator tool adds suffix _materials to the generated tomomasks
            tomoBaseNames = [basename(tomo.getFileName().replace(self.MATERIALS_SUFFIX, '')) for tomo in inTomos]
            for inTomoMask in inTomoMasks:
                outTomoMask = self.setMatchingTomogram(tomoBaseNames, inTomoMask, inTomos, isMatchingByTsId=False)
                outputSetOfTomoMasks.append(outTomoMask)

        self._defineOutputs(**{assignT2TMaskOutputs.tomoMasks.name: outputSetOfTomoMasks})
        self._defineSourceRelation(inTomos, outputSetOfTomoMasks)
        self._defineSourceRelation(inTomoMasks, outputSetOfTomoMasks)

    # --------------------------- INFO functions --------------------------------------------

    def _validate(self):
        errors = []
        inTomoMasks = self.inTomoMasks.get()
        inTomos = self.inputTomos.get()
        tol = 0.01
        if abs(inTomos.getSamplingRate() - inTomoMasks.getSamplingRate()) > tol:
            errors.append('The sampling rate of input sets of tomomasks [%.2f Å/pix] and tomograms [%.2f Å/pix] differ '
                          'in more than %.2f Å/pix.' % (inTomoMasks.getSamplingRate(), inTomos.getSamplingRate(), tol))
        else:
            # Check match by tsId
            tomoMaskTsIds = [tomoMask.getTsId() for tomoMask in inTomoMasks]
            tomoTsIds = [tomo.getTsId() for tomo in inTomos]
            numberMatchesByTsId = len(set(tomoTsIds) & set(tomoMaskTsIds))  # Length of the intersection of both lists

            # Check match by basename
            tomoMaskBaseNames = [tomoMask.getFileName().replace(self.MATERIALS_SUFFIX, '') for tomoMask in inTomoMasks]
            tomoBaseNames = [tomo.getFileName() for tomo in inTomos]
            numberMatchesByBaseName = len(set(tomoMaskBaseNames) & set(tomoBaseNames))
            # Length of the intersection of both lists

            if numberMatchesByTsId == 0 and numberMatchesByBaseName:
                errors.append('Unable to find any match between the introduced datasets after checking the tsIds and '
                              'the basename of all elements.')
        return errors

    def _warnings(self):
        warnings = []
        inTomoMasks = self.inTomoMasks.get()
        inTomos = self.inputTomos.get()
        if isMatchingByTsId(inTomoMasks, inTomos):
            notMatchingMsg = ''
            # Check match by tsId
            tomoTsIds = [tomo.getTsId() for tomo in inTomos]
            for tomoMask in inTomoMasks:
                if tomoMask.getTsId() not in tomoTsIds:
                    notMatchingMsg += '\n\t-%s' % tomoMask.getFileName()
            if notMatchingMsg:
                warnings.append('Not able to find a tsId-based match for the following tomomasks in '
                                'the introduced tomograms:' + notMatchingMsg)
        else:
            notMatchingMsg = ''
            # Check match by basename
            tomoBaseNames = [basename(tomo.getFileName().replace(self.MATERIALS_SUFFIX, '')) for tomo in inTomos]
            for tomoMask in inTomoMasks:
                tomoMaskName = tomoMask.getFileName()
                tomoMaskBaseName = basename(tomoMaskName.replace(self.MATERIALS_SUFFIX, ''))
                if tomoMaskBaseName not in tomoBaseNames:
                    notMatchingMsg += '\n\t-%s' % tomoMaskName
            if notMatchingMsg:
                warnings.append('Not able to find a basename-based match for the following tomomasks in '
                                'the introduced tomograms:' + notMatchingMsg)

        return warnings

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append("%s tomograms assigned to %s tomomasks." %
                           (self.getObjectTag('inputTomos'), self.getObjectTag('inTomoMasks')))
        return summary

# --------------------------- UTILS functions --------------------------------------------
    def setMatchingTomogram(self, idList, inTomoMask, inTomos, isMatchingByTsId=True):
        inFileName = inTomoMask.getFileName()
        outFileName = self._getExtraPath(basename(inFileName))
        symlink(abspath(inFileName), abspath(outFileName))
        outTomoMask = TomoMask()
        outTomoMask.setLocation(inFileName)
        outTomoMask.copyInfo(inTomoMask)
        outTomoMask.setFileName = outFileName
        if isMatchingByTsId:
            outTomoMask.setVolName(inTomos[idList.index(inTomoMask.getTsId()) + 1].getFileName())
        else:
            outTomoMask.setVolName(inTomos[idList.index(basename(inTomoMask.getFileName().replace(
                self.MATERIALS_SUFFIX, ''))) + 1].getFileName())
        # TODO: if the volume has not been matched at this point
        #  try to find out if the tsId is contained in the basename
        return outTomoMask
