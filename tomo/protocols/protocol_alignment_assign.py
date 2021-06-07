# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *              David Herreros Calero (dherreros@cnb.csic.es) (Coordinates3D compatibility)
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
from pwem.protocols import EMProtocol

import tomo.objects as tomoObjs


class ProtAlignmentAssignSubtomo(EMProtocol):
    """ Assign the alignment stored in a set of Subtomograms/Coordinates3D
    to another set.
    Both sets should have same pixel size (A/px).
    The Subtomograms/Coordinates3D with the alignment can also be a subset of a bigger set.
    """
    _label = 'assign alignment'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass='SetOfSubTomograms, SetOfCoordinates3D',
                      label='Input',
                      help='Select the Subtomograms/Coordinates3D that you want to update the new alignment.')
        form.addParam('inputAlignment', PointerParam, pointerClass='SetOfSubTomograms, SetOfCoordinates3D',
                      label="Alignments",
                      help='Select the Subtomograms/Coordinates3D with alignment to be apply to the other object.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        input = self.input.get()
        inputAlignment = self.inputAlignment.get()
        output = input.create(self._getPath())
        output.copyInfo(input)
        if isinstance(output, tomoObjs.SetOfSubTomograms):
            output.setAlignment(inputAlignment.getAlignment())
        output.copyItems(input, updateItemCallback=self._updateItem)
        self._defineOutputs(outputAligned=output)
        self._defineSourceRelation(self.input, output)
        self._defineSourceRelation(self.inputAlignment, output)

    # --------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        # Add alignment info from corresponding item on inputAlignment
        inputAlignment = self.inputAlignment.get()
        alignedObj = inputAlignment[item.getObjId()]
        # If alignment is found for this object set the alignment info
        # on the output object, if not do not write that item
        if alignedObj is not None:
            alignment = self._getAlignment(alignedObj)
            self._setAlignment(item, alignment)
        else:
            item._appendItem = False

    def _getAlignment(self, object):
        if isinstance(object, tomoObjs.Coordinate3D):
            return object.getMatrix()
        elif isinstance(object, tomoObjs.SubTomogram):
            return object.getTransform()

    def _setAlignment(self, object, alignment):
        if isinstance(object, tomoObjs.Coordinate3D):
            object.setMatrix(alignment)
        elif isinstance(object, tomoObjs.SubTomogram):
            return object.setTransform(alignment)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputAligned'):
            summary.append("Output not ready yet.")
        else:
            summary.append("Assigned alignment to %s objects from a total of %s." % (
                self.outputAligned.getSize(), self.input.get().getSize()))
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputAligned'):
            methods.append("Output not ready yet.")
        else:
            methods.append("We assigned alignment to %s objects from %s and produced %s."
                           % (self.outputAligned.getSize(), self.getObjectTag('input'),
                              self.getObjectTag('outputAligned')))
        return methods

    def _validate(self):
        validateMsgs = []
        inputAlignment = self.inputAlignment.get().getFirstItem()
        input = self.input.get().getFirstItem()
        if type(input) != type(inputAlignment):
            validateMsgs.append('*Input* and *Alignments* parameters must belong to the same type '
                                '(SetOfCoordinates3D or SetOfSubtomorgams).')
        else:
            if not inputAlignment.hasTransform():
                validateMsgs.append('Please provide Subtomograms/Coordinates3D which '
                                    'have transformation matrix in "inputAlignment".')
        return validateMsgs
