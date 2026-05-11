# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * National Center of Biotechnology, CSIC, Spain
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
import logging
from enum import Enum
from typing import Union
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTiltSeries, SetOfTomograms, TiltSeries, TiltImage

logger = logging.getLogger(__name__)
IN_TOMO_SET = 'inTomoSet'
IN_TS_SET = 'inTsSet'


class OutputsTsFromTomos(Enum):
    tiltSeries = SetOfTiltSeries


class ProtTsFromTomos(EMProtocol):
    """This protocol gets the tilt-series that corresponds to a given set of tomograms. It is very useful for
    fiducial-less samples, when the quality of alignment results is difficult to be observed from
    the tilt-series, but the tomograms. In that case, the undesired data objects would be discarded
    at tomogram level, but further processing may be desired to be carried out with the corresponding
    tilt-series before getting the final tomograms."""

    """
    Retrieves the tilt-series associated with a selected set of tomograms,
    allowing users to recover the original tilt-series corresponding to
    validated or selected tomographic reconstructions.

    AI Generated:

    Tilt-Series From Tomograms (ProtTsFromTomos) — User Manual
        Overview

        The Tilt-Series From Tomograms protocol is designed to recover
        the subset of tilt-series associated with a given collection of
        tomograms. In cryo-electron tomography workflows, it is common
        to evaluate data quality only after tomographic reconstruction,
        since reconstruction artifacts, alignment problems, missing wedge
        effects, or low contrast may not be sufficiently visible directly
        at the tilt-series level. This protocol provides a convenient way
        to trace validated tomograms back to their original tilt-series,
        enabling additional processing or refinement steps to continue
        from the corresponding raw or aligned projection data.

        From a biological perspective, this protocol is especially useful
        in fiducial-less workflows, where the quality of tilt-series
        alignment can be difficult to assess before reconstruction.
        Researchers often inspect reconstructed tomograms to identify
        datasets with sufficient structural preservation, contrast, or
        alignment quality. Once unsuitable tomograms are discarded, this
        protocol allows the user to automatically recover only the
        tilt-series linked to the accepted tomograms.

        Inputs and General Workflow

        The protocol requires two inputs: a set of tomograms and a set
        of tilt-series. Both datasets are expected to share common tilt-
        series identifiers (tsIds), which are used internally to establish
        the correspondence between tomograms and their originating
        tilt-series.

        During execution, the protocol compares the identifiers present
        in both datasets and computes their intersection. Only tilt-series
        whose identifiers are present in both collections are preserved
        in the final output. This behavior ensures consistency between
        tomographic reconstructions and their associated projection data.

        In practical cryo-ET workflows, this allows users to manually or
        automatically curate tomograms first, and then continue processing
        only the corresponding tilt-series. Typical downstream applications
        include improved CTF estimation, subtomogram extraction, particle
        picking, denoising workflows, or reconstruction refinement.

        Identifier Matching and Dataset Consistency

        A central aspect of the protocol is the use of tsIds to establish
        relationships between tomograms and tilt-series. The protocol
        validates the overlap between both datasets and raises an exception
        if no common identifiers are found. This prevents accidental
        propagation of unrelated or incompatible datasets.

        When only partial overlap exists, the protocol reports the non-
        matching identifiers while still processing the common subset.
        This behavior is biologically useful because tomography datasets
        are often curated incrementally, and some tomograms or tilt-series
        may have been removed in previous quality-control stages.

        From a data-management perspective, preserving identifier
        consistency throughout the workflow is extremely important.
        Misaligned identifiers may lead to incorrect downstream analysis,
        particularly in automated pipelines where metadata relationships
        are heavily relied upon.

        Output Generation

        The protocol creates a new SetOfTiltSeries containing only the
        tilt-series associated with the selected tomograms. Each tilt-
        series is cloned together with its corresponding tilt images,
        preserving acquisition metadata, ordering, and structural
        information.

        The resulting output behaves as a standard Scipion tilt-series
        dataset and can therefore be directly connected to subsequent
        tomography protocols. Since the protocol only filters and copies
        metadata references rather than recomputing projections, execution
        is computationally lightweight and very fast even for large
        datasets.

        Biological and Practical Applications

        In biological cryo-ET studies, this protocol becomes especially
        valuable after tomogram inspection and curation. Researchers often
        reconstruct large numbers of tomograms and later retain only those
        showing sufficient quality for downstream structural analysis.
        Once the curated tomograms are identified, recovering the matching
        tilt-series allows the workflow to continue from projection-space
        data without manually tracking datasets.

        This approach is particularly important in fiducial-less alignment
        strategies, where reconstruction quality is frequently the most
        reliable criterion for evaluating dataset usability. By linking
        validated tomograms back to their tilt-series, users can refine
        alignment procedures, repeat reconstruction under different
        parameters, or apply advanced processing methods only to the best
        datasets.

        Practical Recommendations

        In routine workflows, users should ensure that both the tomogram
        and tilt-series datasets originate from compatible processing
        pipelines and preserve consistent tsIds. Maintaining clean and
        traceable metadata throughout the workflow greatly simplifies
        dataset management and avoids ambiguity during selection steps.

        It is generally advisable to perform tomogram curation carefully
        before applying this protocol. Since the output directly reflects
        the selected tomograms, any filtering or quality-control decision
        at the tomogram level will propagate to all subsequent processing
        stages involving the recovered tilt-series.

        Final Perspective

        The Tilt-Series From Tomograms protocol serves as a bridge between
        tomogram-level quality assessment and tilt-series-level processing.
        Although technically simple, it fulfills an important organizational
        and biological role in cryo-electron tomography workflows by
        preserving the connection between reconstructed volumes and their
        originating projection data. This capability is particularly useful
        in large-scale tomography projects, fiducial-less workflows, and
        iterative refinement strategies where dataset traceability is
        essential for reliable structural interpretation.
    """