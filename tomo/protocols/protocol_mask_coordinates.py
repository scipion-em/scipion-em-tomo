# **************************************************************************
# *
# * Authors:     Oier Lauzirika Zarrabeitia (olauzirika@cnb.csic.es)
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
import numpy as np

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam

from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.objects import Set

from tomo.objects import (SetOfCoordinates3D, Coordinate3D, SetOfTomograms,
                          SetOfTomoMasks, TomoMask, Tomogram)
from tomo.protocols import ProtTomoBase
import tomo.constants as const

COORDINATES = 'Coordinates'

class ProtMaskCoordinates(EMProtocol, ProtTomoBase):
    """ To be filled by Oier HAHAHA
    """

    """
    Filters 3D coordinates using tomographic segmentation masks in order
    to retain only coordinates located within biologically relevant
    segmented regions. The protocol evaluates each coordinate against
    a tomomask and generates a new filtered set of coordinates based on
    the selected segmentation criteria.

    AI Generated:

    Mask 3D Coordinates (ProtMaskCoordinates) — User Manual
        Overview

        The Mask 3D Coordinates protocol filters a set of 3D coordinates
        according to tomographic segmentation masks. Its primary purpose
        is to remove coordinates located outside biologically relevant
        segmented regions while preserving those associated with the
        structures of interest. In cryo-electron tomography workflows,
        this operation is particularly useful for restricting particle
        picking results to specific cellular compartments, membranes,
        organelles, or macromolecular assemblies identified through
        segmentation procedures.

        In practical biological applications, automated particle picking
        frequently produces coordinates distributed across both relevant
        and irrelevant regions of the tomogram. By combining coordinate
        datasets with segmentation masks, the protocol allows users to
        focus subsequent analyses only on spatial regions supported by
        prior structural or biological knowledge. This improves dataset
        specificity and reduces the introduction of false positives into
        subtomogram averaging or classification workflows.

        Inputs and General Workflow

        The protocol requires two main inputs: a set of 3D coordinates
        and a set of tomographic segmentation masks. The coordinate set
        provides the spatial positions to be evaluated, while the
        segmentation masks define the regions considered biologically
        relevant. Each tomogram is processed independently by matching
        coordinates and segmentations through their corresponding
        tilt-series identifiers.

        During execution, the protocol creates a new output coordinate
        set linked to the original tomograms. For every tomogram, the
        associated segmentation mask is loaded into memory and converted
        into a binary mask representation. Each coordinate position is
        then evaluated against this mask in order to determine whether
        it belongs to a segmented region.

        Segmentation Labels and Region Selection

        An important feature of the protocol is the possibility of
        filtering coordinates using specific segmentation labels.
        Segmentation datasets often contain multiple annotated regions
        represented by different integer values. The protocol allows
        users either to select a particular label or to consider all
        non-zero segmented regions simultaneously.

        When the segmentation label parameter is negative, every
        non-background voxel is treated as a valid segmented region.
        This mode is useful for general filtering tasks in which all
        annotated structures should be retained. Alternatively, when a
        positive label value is provided, only coordinates falling
        within voxels matching that exact label are preserved. This
        enables highly selective biological analyses focused on specific
        organelles, membrane systems, or structural compartments.

        Coordinate Evaluation Strategy

        For each coordinate, the protocol retrieves its spatial position
        relative to the tomogram reference frame and converts the
        coordinates into voxel indices compatible with the segmentation
        mask dimensions. The coordinate is accepted only if the
        corresponding voxel position in the mask evaluates to true.

        This voxel-based filtering strategy ensures direct spatial
        consistency between segmentation data and particle coordinates.
        As a result, only coordinates physically located inside the
        segmented structures are propagated into the output dataset.

        Handling Missing Segmentations

        Biological datasets are not always fully segmented, and some
        tomograms may lack associated masks. The protocol therefore
        provides configurable behavior for handling these cases.

        When exclusion of unsegmented tomograms is enabled, coordinates
        belonging to tomograms without a corresponding segmentation are
        discarded entirely. This behavior is appropriate for highly
        curated workflows where only validated segmented regions should
        contribute to downstream analysis.

        Alternatively, users may disable this exclusion behavior. In
        that case, coordinates from tomograms lacking segmentations are
        copied directly into the output dataset without filtering. This
        option is useful in partially annotated datasets where retaining
        unsegmented tomograms remains biologically meaningful.

        Validation and Dataset Consistency

        Before execution, the protocol validates the compatibility
        between coordinate and segmentation datasets. The sampling rate
        of the segmentations must match the sampling rate associated
        with the tomograms used during coordinate picking. Likewise,
        tomogram and segmentation dimensions must be identical to
        guarantee spatial consistency during voxel-based masking.

        These validation steps are biologically important because even
        small mismatches in voxel size or dimensions may shift the
        coordinate positions relative to the segmentation mask, leading
        to incorrect filtering decisions.

        Outputs and Their Interpretation

        After processing, the protocol generates a new set of filtered
        3D coordinates linked to the original tomograms and input
        segmentations. The output preserves all metadata associated with
        the original coordinate set while restricting the coordinates to
        the segmented spatial regions selected during execution.

        Biologically, the resulting dataset represents a spatially
        refined subset of particles or annotations that are consistent
        with the segmentation information. This refinement often improves
        the quality of downstream subtomogram averaging, classification,
        or structural interpretation workflows.

        Practical Recommendations

        In routine cryo-electron tomography workflows, this protocol is
        particularly effective when segmentation masks represent stable
        or biologically meaningful structures such as membranes,
        organelles, cytoskeletal networks, or viral assemblies. Careful
        selection of segmentation labels is essential because incorrect
        label usage may unintentionally remove relevant coordinates or
        retain unwanted regions.

        Users should also verify that segmentation masks and coordinate
        datasets share the same voxel size and geometric dimensions
        before execution. Even when validation passes, visually checking
        several filtered coordinates in a tomographic viewer is strongly
        recommended to confirm biological correctness.

        Final Perspective

        Coordinate masking is not simply a technical filtering step but
        a biologically driven spatial refinement process. By integrating
        segmentation knowledge with particle coordinate datasets, the
        protocol enables more selective and biologically coherent
        tomographic analyses while reducing noise and irrelevant spatial
        information in downstream workflows.
    """
