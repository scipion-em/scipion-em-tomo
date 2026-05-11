# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
logger = logging.getLogger(__name__)

import math
import numpy as np

from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.protocol.params import MultiPointerParam, FloatParam, EnumParam
from pyworkflow.object import Set, Float

from pwem.protocols import EMProtocol

from tomo.objects import TiltSeries, TiltImage
from tomo.protocols import ProtTomoBase


class ProtConsensusAlignmentTS(EMProtocol, ProtTomoBase):
    """
    Perform a consensus of a set of alignments for the same tilt series. Returns the average alignment matrix of the
    consensus alignments and its standard deviation of shift and angle.
    """

    """
    Computes a consensus alignment from multiple independently aligned tilt-series
    datasets corresponding to the same tomographic acquisition. The protocol evaluates
    the agreement between alignment transformations and generates a consensus
    alignment together with statistical measurements describing alignment variability
    across datasets.

    AI Generated:

    Tilt-Series Consensus Alignment (ProtConsensusAlignmentTS) — User Manual
        Overview

        The Tilt-Series Consensus Alignment protocol evaluates the consistency of
        multiple alignment solutions computed for the same tilt-series dataset.
        Its main objective is to identify reproducible geometric transformations
        across independent alignments and to generate a robust consensus alignment
        that minimizes the influence of outliers or unstable registrations.

        In cryo-electron tomography workflows, tilt-series alignment is one of the
        most critical preprocessing stages because inaccuracies in image registration
        directly affect tomogram reconstruction quality and downstream structural
        interpretation. Different alignment strategies, fiducial selections,
        preprocessing conditions, or optimization procedures may produce slightly
        different alignment solutions for the same experimental dataset. This
        protocol provides a systematic way to compare those solutions and determine
        whether they converge toward a biologically and geometrically consistent
        alignment.

        From a practical perspective, the protocol acts as a validation and
        robustness-analysis framework for tilt-series alignment. Rather than relying
        on a single alignment estimation, users can compare several independently
        generated alignment sets and identify stable consensus transformations.

        Inputs and General Workflow

        The protocol requires multiple sets of aligned tilt-series as input. Each
        set must contain the same tilt-series identifiers and corresponding
        transformation matrices assigned to each tilt-image.

        During execution, the protocol extracts the transformation matrices from
        every input tilt-series and compares them pairwise. The comparison evaluates
        rotational and translational agreement between alignment solutions while
        accounting for differences in sampling rate between datasets.

        The protocol supports two consensus strategies: global consensus and local
        consensus. In global consensus mode, the complete tilt-series alignment is
        evaluated as a single entity, meaning that all tilt-images must collectively
        agree within the specified tolerances. In local consensus mode, each
        tilt-image is evaluated independently, allowing consensus to be achieved
        even when only subsets of images agree across the different alignments.

        Once consensus is established, the protocol computes an averaged
        transformation matrix together with angular and shift standard deviations
        that quantify alignment variability.

        Global and Local Consensus Strategies

        The choice between global and local consensus has important practical and
        biological implications. Global consensus is the strictest approach because
        it assumes that a valid alignment should remain internally consistent across
        the entire tilt-series. If one alignment deviates significantly from the
        others, the protocol recursively removes inconsistent solutions until a
        stable consensus is achieved or no valid agreement remains.

        This strategy is particularly useful for high-quality datasets where
        alignment consistency is expected across all projections. It is also well
        suited for automated pipelines that require globally stable geometric
        reconstructions before tomogram generation.

        In contrast, local consensus evaluates each tilt-image independently. This
        approach is more flexible and tolerant to localized alignment failures,
        which are relatively common in cryo-ET datasets affected by contamination,
        low contrast, specimen deformation, or missing fiducials at specific tilt
        angles.

        Biologically, local consensus can preserve usable alignment information
        even in partially degraded tilt-series, avoiding the complete rejection of
        datasets where only a subset of projections is problematic.

        Shift and Angle Tolerances

        Consensus evaluation relies on user-defined angular and translational
        tolerances. The shift tolerance defines the maximum accepted translational
        deviation between alignments and is expressed in Angstroms, while the angle
        tolerance specifies the maximum accepted rotational difference in degrees.

        These tolerances directly determine the strictness of the consensus process.
        Small tolerance values enforce highly reproducible alignments and are
        appropriate for high-resolution workflows where geometric precision is
        essential. However, excessively strict thresholds may reject biologically
        acceptable alignments, particularly in noisy or heterogeneous datasets.

        More permissive tolerances increase robustness and allow consensus to be
        achieved in challenging experimental conditions, although they may also
        incorporate larger geometric variability into the final averaged alignment.

        In practice, moderate tolerances often provide a good balance between
        alignment reliability and robustness to experimental noise.

        Consensus Averaging and Variability Estimation

        When consensus is achieved, the protocol computes an averaged transformation
        matrix from all accepted alignments. This averaging process generates a
        representative alignment solution that reduces the influence of isolated
        errors or unstable optimization results.

        In addition to the consensus alignment itself, the protocol estimates the
        standard deviation of rotational angles and translational shifts for each
        tilt-image. These statistics provide a quantitative measurement of alignment
        stability across the compared datasets.

        Low standard deviations indicate highly reproducible alignments and suggest
        strong geometric consistency between independent alignment procedures.
        Conversely, large deviations may reveal problematic projections, unstable
        fiducial tracking, low signal regions, or inconsistencies introduced during
        preprocessing.

        From a biological perspective, these variability measurements can help
        identify regions of the tilt-series that may compromise tomogram quality or
        downstream subtomogram averaging analyses.

        Outputs and Their Interpretation

        The protocol produces two possible output sets. The first contains
        tilt-series for which consensus alignment was successfully achieved. These
        outputs include the averaged transformation matrices together with the
        associated angular and shift variability statistics.

        The second output contains tilt-series for which no reliable consensus
        could be established. In these cases, the original tilt-images are preserved
        without consensus alignment information.

        This separation allows users to distinguish between geometrically reliable
        datasets and potentially problematic acquisitions that may require further
        inspection or reprocessing.

        Practical Recommendations

        In routine cryo-ET workflows, it is generally advisable to compare alignment
        solutions generated using different preprocessing conditions or fiducial
        tracking strategies. Consensus analysis can reveal whether the final
        alignment is stable or highly dependent on a particular processing choice.

        Global consensus is typically recommended for high-quality datasets with
        homogeneous alignment behavior, while local consensus is often more suitable
        for noisy cellular tomograms or datasets with localized alignment failures.

        Users should carefully inspect tilt-images showing large alignment standard
        deviations, since these projections frequently correspond to regions with
        poor contrast, contamination, ice deformation, or inaccurate fiducial
        detection.

        When consensus cannot be achieved, the dataset should not necessarily be
        considered unusable. Instead, the absence of consensus should be interpreted
        as an indication that the alignment requires additional validation,
        parameter optimization, or alternative preprocessing strategies.

        Final Perspective

        The Tilt-Series Consensus Alignment protocol provides a robust framework
        for validating and consolidating alignment solutions in cryo-electron
        tomography workflows. By combining multiple independent alignment estimates
        into statistically evaluated consensus transformations, the protocol helps
        improve alignment reliability while providing quantitative measurements of
        geometric stability.

        Beyond simple averaging, the protocol introduces an important quality-control
        layer into tomographic processing pipelines, enabling users to identify
        unstable projections, assess reproducibility, and increase confidence in
        downstream tomogram reconstruction and structural interpretation.
    """