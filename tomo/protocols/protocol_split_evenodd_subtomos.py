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
    """ Protocol to split set of tomograms or subtomograms in even/odd sets by element id.
    """

    """
    ProtSplitEvenOddTomoSet — Split Tomograms/Subtomograms into Even and Odd Sets

    This protocol separates a set of tomograms or subtomograms into two
    independent subsets according to the object identifier parity. Elements
    with even identifiers are assigned to one output set, while elements
    with odd identifiers are assigned to another. The protocol is mainly
    intended for dataset partitioning, validation workflows, and independent
    processing strategies in cryo-electron tomography pipelines.

    AI Generated:

    Split Even/Odd Tomograms or Subtomograms (ProtSplitEvenOddTomoSet)
    — User Manual

        Overview

        The Split Even/Odd Tomograms protocol divides an input dataset into
        two separate subsets based on the parity of the internal object
        identifiers. Tomograms or subtomograms with even identifiers are
        stored in one output set, while elements with odd identifiers are
        stored in another.

        In practical cryo-electron tomography workflows, this type of split
        is commonly used to generate independent datasets for validation,
        benchmarking, testing reproducibility, or parallel processing
        strategies. By separating the data into two groups, users can
        evaluate consistency between independent reconstructions or compare
        processing outcomes under different conditions.

        Inputs and General Workflow

        The protocol accepts either a SetOfTomograms or a
        SetOfSubTomograms as input. During execution, the protocol inspects
        each element in the dataset and evaluates its object identifier.
        Elements whose identifiers are divisible by two are assigned to the
        even subset, while the remaining elements are assigned to the odd
        subset.

        The protocol preserves the metadata and acquisition information from
        the original dataset, ensuring that both output subsets remain fully
        compatible with downstream tomography workflows inside Scipion.

        Dataset Partitioning Strategy

        The splitting strategy is deterministic and entirely based on object
        identifiers. This means that repeated executions on the same dataset
        will always generate identical even and odd subsets, which is
        particularly useful for reproducible benchmarking and validation
        experiments.

        Since the protocol does not alter the tomograms or subtomograms
        themselves, the resulting subsets maintain the original spatial,
        biological, and acquisition properties of the input data.

        Outputs and Interpretation

        After execution, the protocol produces two output datasets:
        an even set and an odd set. Both outputs preserve the structure,
        metadata, and object relationships of the original input collection.

        These subsets can be processed independently in downstream cryo-ET
        workflows, including subtomogram averaging, classification,
        reconstruction validation, or testing of independent refinement
        strategies.

        From a biological perspective, the protocol does not introduce any
        transformation or modification to the underlying data. Its role is
        purely organizational, facilitating controlled experimental designs
        and reproducible computational analyses.

        Practical Recommendations

        In routine cryo-electron tomography workflows, splitting datasets
        into even and odd subsets is useful when users need to validate the
        robustness of a processing pipeline or compare independent
        reconstructions. It can also help distribute computational load
        across multiple processing branches.

        Users should keep in mind that the split is based solely on object
        identifiers and not on biological content, acquisition conditions,
        or structural similarity. Therefore, if the dataset contains ordered
        acquisitions or strongly heterogeneous populations, additional care
        may be needed to ensure balanced biological representation between
        subsets.

        Final Perspective

        Although simple in implementation, dataset partitioning is an
        important organizational step in many cryo-ET workflows. By
        generating reproducible even and odd subsets, this protocol provides
        a straightforward mechanism for validation, independent analysis,
        and reproducible testing within Scipion tomography pipelines.

    """