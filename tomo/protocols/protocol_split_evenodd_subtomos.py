# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from tomo.objects import SetOfTomograms
from tomo.protocols import ProtTomoBase


class ProtSplitEvenOddTomoSet(EMProtocol, ProtTomoBase):
    """
    Splits a set of tomograms or subtomograms into two independent subsets
    according to the parity of each element identifier.

    AI Generated:

    Split Even/Odd Tomograms or Subtomograms (ProtSplitEvenOddTomoSet) — User Manual
        Overview

        The Split Even/Odd Tomograms or Subtomograms protocol divides an input
        set into two separate subsets based on the object identifier of each
        element. Items with even identifiers are placed in one output set,
        while items with odd identifiers are placed in another.

        In tomography workflows, this type of partition is commonly used when
        preparing independent half-sets for validation, testing reproducibility,
        or running parallel downstream analyses. The protocol does not modify
        the tomograms or subtomograms themselves; it only reorganizes the input
        data into two complementary groups.

        Inputs and General Workflow

        The protocol accepts as input either a set of tomograms or a set of
        subtomograms.

        During execution, the protocol first determines the nature of the input
        object. If the input corresponds to tomograms, two new sets of
        tomograms are created. If the input corresponds to subtomograms, two
        new sets of subtomograms are created instead.

        Both output sets inherit the metadata and general information from the
        original input set, ensuring that acquisition parameters, sampling
        information, and other relevant attributes remain consistent.

        Splitting Criterion

        The separation is based exclusively on the object identifier
        (`objId`) of each element.

        Elements whose identifier is divisible by two are assigned to the
        *even* output set. Elements whose identifier is not divisible by two
        are assigned to the *odd* output set.

        This criterion guarantees a deterministic partition of the input set.
        Running the protocol multiple times on the same input produces the same
        split.

        Biological Interpretation

        From a biological perspective, this protocol does not perform any form
        of structural classification, quality filtering, or data selection
        based on image content. The split is purely technical.

        For that reason, the resulting even and odd subsets should be
        interpreted simply as two independent partitions of the same dataset.
        They are particularly useful in workflows where one wishes to maintain
        statistical independence between subsets, for example during
        validation-oriented reconstruction procedures.

        Outputs and Their Interpretation

        After execution, the protocol generates two output sets:

        - **outputset_even**: contains all elements with even object identifiers.
        - **outputset_odd**: contains all elements with odd object identifiers.

        The outputs preserve the same data type as the input. Therefore, if the
        input is a set of tomograms, both outputs are tomogram sets. If the
        input is a set of subtomograms, both outputs are subtomogram sets.

        Source relations are also preserved, allowing downstream protocols to
        trace both subsets back to the original dataset.

        Practical Recommendations

        This protocol is most useful when a simple and reproducible split of
        the dataset is required.

        Since the partition depends entirely on object identifiers, users
        should keep in mind that the even/odd division does not guarantee
        biological balance between subsets. If the input dataset was assembled
        in a non-random order, one subset could accidentally contain a biased
        representation of the data.

        In routine practice, this protocol is best suited for technical
        half-set generation rather than for biologically meaningful sampling.

        Final Perspective

        The Split Even/Odd Tomograms or Subtomograms protocol provides a fast,
        deterministic, and lightweight way to divide tomography datasets into
        two complementary subsets.

        Although computationally simple, it can play an important practical
        role in validation workflows, independent testing pipelines, and
        downstream processing strategies where reproducible dataset partitioning
        is required.
    """
    _label = 'split even/odd tomos/subtomos'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfSubTomograms, SetOfTomograms',
                      label="Set to split",
                      help='Select the set of tomograms or subtomograms that you '
                           'want to split in even/odd sets.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        if isinstance(inputSet, SetOfTomograms):
            evenSet = self._createSetOfTomograms(suffix='_even')
            oddSet = self._createSetOfTomograms(suffix='_odd')
        else:
            evenSet = self._createSetOfSubTomograms(suffix='_even')
            oddSet = self._createSetOfSubTomograms(suffix='_odd')

        evenSet.copyInfo(inputSet)
        oddSet.copyInfo(inputSet)

        for element in inputSet:
            if element.getObjId() % 2 == 0:
                evenSet.append(element)
            else:
                oddSet.append(element)

        self._defineOutputs(outputset_even=evenSet)
        self._defineSourceRelation(inputSet, evenSet)
        self._defineOutputs(outputset_odd=oddSet)
        self._defineSourceRelation(inputSet, oddSet)

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        if not self.isFinished():
            return["Output sets not ready yet."]
        else:
            return["We have split the input set in even and odd sets."]
