# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils import removeExt
from pwem.protocols import EMProtocol


class ProtAssignTomo2Subtomo(EMProtocol):
    """ This protocol assign tomograms to subtomograms that have been imported before without tomograms.
    Subtomograms should contain the name of the original tomogram in their own file name.
    """
    _label = 'assign tomos to subtomos'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass='SetOfSubTomograms', label='Subtomograms',
                      help='Select the subtomograms that you want to update with original tomograms as precedents.'
                           'The subtomograms should contain the original tomogram name in their own filename.')
        form.addParam('inputTomos', PointerParam, pointerClass='SetOfTomograms', label='Tomograms',
                      help='Select the tomograms to be assigned to the subtomograms.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        self.tomoDict = {}
        inputTomos = self.inputTomos.get()
        for tomo in inputTomos:
            self.tomoDict[tomo.getBaseName()] = tomo.getObjId()
        inputSubtomos = self.inputSubtomos.get()
        outputSubtomos = inputSubtomos.create(self._getPath())
        outputSubtomos.copyInfo(inputSubtomos)
        outputSubtomos.copyItems(inputSubtomos, updateItemCallback=self._updateItem)
        self._defineOutputs(outputSubtomograms=outputSubtomos)
        self._defineSourceRelation(self.inputSubtomos, outputSubtomos)
        self._defineSourceRelation(self.inputTomos, outputSubtomos)

    # --------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        for tomoName in self.tomoDict:
            if removeExt(tomoName) in item.getFileName():
                item.setVolName(tomoName)
                item.setVolId(self.tomoDict.get(tomoName))
                break

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("%s tomograms assigned to %s subtomograms." %
                           (self.getObjectTag('inputTomos'), self.getObjectTag('inputSubtomos')))
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputSubtomograms'):
            methods.append("Output subtomograms not ready yet.")
        else:
            methods.append("%d %s tomograms assigned to %d %s subtomograms." %
                           (self.inputTomos.get().getSize(), self.getObjectTag('inputTomos'),
                            self.inputSubtomos.get().getSize(), self.getObjectTag('inputSubtomos')))
        return methods
