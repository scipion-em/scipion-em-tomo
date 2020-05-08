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


from pyworkflow.protocol.params import PointerParam
from tomo.protocols import ProtTomoBase


class ProtAlignmentAssignSubtomo(ProtTomoBase):
    """ Assign the alignment calculated for a set of subtomograms to another set.
    Both sets should have same pixel size (A/px).
    The subtomograms with the alignment can also be a subset of a bigger set.
    """
    _label = 'assign alignment subtomo'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass='SetOfSubTomograms', label='Input subtomograms',
                      help='Select the subtomograms that you want to update the new alignment.')
        form.addParam('inputAlignment', PointerParam, pointerClass='SetOfSubTomograms',
                      label="Input alignments",
                      help='Select the subtomograms with alignment to be apply to the other particles.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        inputSubtomos = self.inputSubtomos.get()
        inputAlignment = self.inputAlignment.get()
        outputSubtomos = self._createSetOfSubTomograms()
        outputSubtomos.copyInfo(inputSubtomos)
        outputSubtomos.setAlignment(inputAlignment.getAlignment())
        outputSubtomos.copyItems(inputSubtomos, updateItemCallback=self._updateItem)
        self._defineOutputs(outputSubtomograms=outputSubtomos)
        self._defineSourceRelation(self.inputSubtomos, outputSubtomos)
        self._defineSourceRelation(self.inputAlignment, outputSubtomos)

    # --------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        # Add alignment info from corresponding item on inputAlignment
        inputAlignment = self.inputAlignment.get()
        alignedParticle = inputAlignment[item.getObjId()]
        # If alignment is found for this particle set the alignment info
        # on the output particle, if not do not write that item
        if alignedParticle is not None:
            alignment = alignedParticle.getTransform()
            item.setTransform(alignment)
        else:
            item._appendItem = False

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("Assigned alignment to %s subtomograms from a total of %s." % (
                self.outputSubtomograms.getSize(), self.inputSubtomos.get().getSize()))
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputSubtomograms'):
            methods.append("Output subtomograms not ready yet.")
        else:
            methods.append("We assigned alignment to %s subtomograms from %s and produced %s."
                           % (self.outputSubtomograms.getSize(), self.getObjectTag('inputSubtomos'),
                              self.getObjectTag('outputSubtomograms')))
        return methods

    def _validate(self):
        # check that input set of aligned particles do have 3D alignment
        errors = []
        inputAlignmentSet = self.inputAlignment.get()
        if not inputAlignmentSet.hasAlignment():
            errors.append("Input alignment set should contains some kind of alignment (2D, 3D or Projection).")
        else:
            # Just for consistency, check that the subtomograms really contains Transform object
            first = inputAlignmentSet.getFirstItem()
            alignment = first.getTransform()
            if alignment is None:
                errors.append('Inconsistency detected in *Input alignment* !!!')
                errors.append('It has alignment: _%s_, but the alignment is missing!!!' %
                              inputAlignmentSet.getAlignment())
        return errors
