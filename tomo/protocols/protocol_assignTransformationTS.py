# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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
import time
import traceback
import typing
from collections import Counter
from enum import Enum
from typing import List, Union
import numpy as np
from pwem.objects import Transform
from pyworkflow import BETA
from pwem.protocols import EMProtocol
from pyworkflow.object import Pointer, Set
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase, PointerParam
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtStreamingBase):
    """
    Assign the transformation matrices from an input set of tilt-series to a target one.
    """

    """
    Assigns transformation matrices from one set of tilt-series to another,
    allowing alignment information obtained in a previously processed dataset
    to be transferred and reused in a different tilt-series collection.

    AI Generated:

    Assign Transformation Matrix Tilt-Series (ProtAssignTransformationMatrixTiltSeries) — User Manual

        Overview

        The Assign Transformation Matrix Tilt-Series protocol transfers alignment
        information from one tilt-series dataset to another by assigning the
        transformation matrices of a reference set onto a target set. Its main
        purpose is to preserve geometrical alignment relationships between
        corresponding tilt-images while avoiding the need to recompute alignment
        parameters from scratch.

        In cryo-electron tomography workflows, this protocol is especially useful
        when working with reprocessed datasets, re-stacked tilt-series, corrected
        image collections, or alternative preprocessing pipelines where the image
        content changes but the acquisition geometry remains compatible. Instead
        of repeating computationally expensive alignment procedures, the protocol
        propagates previously validated transformations to a new dataset.

        Inputs and General Workflow

        The protocol requires two input sets of tilt-series. The first set acts
        as the alignment source and must already contain valid transformation
        matrices. The second set acts as the destination dataset, receiving the
        alignment information from the source set.

        During execution, the protocol continuously monitors both datasets in
        streaming mode and identifies matching tilt-series using their tilt-series
        identifiers. Once matching datasets are detected, the protocol creates a
        new output tilt-series where the transformation matrices from the source
        dataset are assigned to the corresponding tilt-images of the target
        dataset.

        The protocol has been designed for streaming environments, allowing
        tilt-series to be processed incrementally as soon as they become available.
        This behavior is particularly useful in automated tomography pipelines or
        facility-level workflows where data may still be actively generated during
        processing.

        Matching Tilt-Series and Acquisition Consistency

        A key aspect of the protocol is the matching of acquisition orders between
        source and target tilt-images. Transformation matrices are only assigned
        when both tilt-series contain images with compatible acquisition orders.
        This ensures that geometrical consistency is preserved during the transfer
        process.

        The protocol also handles cases where certain views were excluded during
        previous processing stages or when tilt-series were re-stacked. If a
        tilt-image exists in the target dataset but no compatible acquisition
        order is found in the source dataset, the image is marked as disabled and
        assigned an identity transformation matrix instead of a potentially
        incorrect alignment.

        From a biological and geometrical perspective, maintaining acquisition-order
        consistency is critical because transformation matrices directly describe
        the spatial relationship between projections. Incorrect mapping between
        views could propagate alignment errors into downstream tomographic
        reconstruction or subtomogram analysis.

        Streaming Execution and Parallel Processing

        The protocol operates using a streaming-oriented execution model. It
        continuously checks whether new tilt-series have appeared in the input
        datasets and schedules alignment-transfer tasks dynamically. This allows
        processing to begin before the complete datasets are fully available.

        Parallel execution is supported, enabling several tilt-series to be
        processed simultaneously. This significantly improves throughput in
        high-volume tomography facilities or automated acquisition pipelines.

        The protocol keeps track of processed tilt-series identifiers to avoid
        duplicate processing and automatically closes the output set once all
        compatible tilt-series have been transferred and both input streams are
        finalized.

        Transformation Matrix Adaptation

        An important feature of the protocol is the automatic adaptation of
        translational shifts according to the sampling-rate ratio between the
        source and target datasets. When voxel sizes differ between datasets,
        translational components of the transformation matrices are rescaled to
        preserve geometrical correctness.

        This adjustment is particularly important when the target tilt-series has
        undergone binning, resizing, or resampling operations. Without this
        correction, alignment shifts would become inconsistent with the new pixel
        size and could lead to reconstruction artifacts or spatial distortions.

        The protocol also updates tilt-axis angles whenever these values were
        modified or refined during previous alignment procedures. This guarantees
        coherence between metadata and the assigned transformation matrices.

        Outputs and Interpretation

        The protocol produces a new set of tilt-series containing the images from
        the destination dataset together with the transformation matrices imported
        from the reference dataset. The output preserves the metadata and geometry
        of the target tilt-series while incorporating the alignment information
        from the source dataset.

        Disabled or unmatched views remain explicitly flagged, ensuring that
        incomplete correspondences are handled safely and transparently. This is
        especially important in tomography workflows where missing projections or
        discarded views are relatively common.

        The resulting output can be directly used for tomographic reconstruction,
        subtomogram averaging, or downstream alignment-sensitive analyses without
        requiring a new alignment estimation step.

        Validation and Reliability

        Before execution, the protocol validates that the source tilt-series
        actually contain transformation matrices. If no alignment information is
        available, execution is halted to prevent invalid assignments.

        From a practical perspective, users should ensure that both tilt-series
        datasets originate from compatible acquisition schemes and preserve the
        same acquisition ordering. Major inconsistencies between datasets may lead
        to incomplete assignments or disabled projections.

        Final Perspective

        For cryo-electron tomography users, this protocol provides an efficient
        mechanism for reusing previously computed alignment information across
        multiple processing branches. By preserving geometrical consistency while
        avoiding redundant alignment computation, the protocol simplifies complex
        tomography workflows and accelerates iterative data processing strategies.

        In practical biological workflows, this capability becomes particularly
        valuable when comparing preprocessing methods, generating alternative
        reconstructions, or integrating streaming acquisition pipelines into
        automated Scipion environments.
    """