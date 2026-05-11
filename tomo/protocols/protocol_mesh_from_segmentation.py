# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
# *              Oier Lauzirika Zarrabeita (oierlauzi@bizkaia.eu)
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
import logging
import traceback
from enum import Enum
import math
from os.path import basename

import numpy as np
import skimage.morphology
from pwem.emlib.image.image_readers import ImageReadersRegistry, MRCImageReader
from pwem.protocols import EMProtocol
from pwem.objects import Set, Integer
from pyworkflow import BETA
from pyworkflow.protocol import STEPS_PARALLEL, LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam, GE, LE)
from pyworkflow.utils import yellowStr, cyanStr, redStr, createLink
from tomo.objects import (SetOfMeshes, MeshPoint, SetOfTomoMasks, TomoMask,
                          SetOfTomograms, Tomogram)
import tomo.constants as const

logger = logging.getLogger(__name__)
MRC_EXT = '.mrc'


class OutputMeshesFromSegmentation(Enum):
    meshes = SetOfMeshes


class ProtMeshFromSegmentation(EMProtocol):
    """
    Creates meshes based on segmentations or voxels values (TomoMasks).
    """

    """
    Generates 3D mesh representations from tomographic segmentations or
    probabilistic voxel masks (TomoMasks). The protocol converts labeled
    regions or density-based masks into sparse point-based meshes that can
    be visualized, analyzed, or used in downstream structural workflows.

    AI Generated:

    Mesh From Segmentation (ProtMeshFromSegmentation) — User Manual

        Overview

        The Mesh From Segmentation protocol transforms tomographic masks
        into collections of 3D mesh points associated with tomograms. Its
        main purpose is to convert segmented biological structures into a
        geometric representation that can be visualized efficiently and used
        for spatial analysis within cryo-electron tomography workflows.

        In practical cryo-ET studies, biological objects such as membranes,
        filaments, vesicles, organelles, or macromolecular assemblies are
        often segmented before quantitative interpretation. This protocol
        allows those segmented regions to be converted into point-based
        meshes that preserve the spatial organization of the original
        structures while significantly reducing data complexity.

        The generated meshes are particularly useful for visualization
        environments, geometric analyses, particle contextualization, and
        downstream tools requiring sparse structural representations rather
        than dense voxel volumes.

        Inputs and General Workflow

        The protocol requires two principal inputs: a set of tomographic
        masks and a corresponding set of tomograms. Each tomographic mask
        must share a common tsId with its associated tomogram so both objects
        can be matched correctly during processing.

        During initialization, the protocol identifies the tomograms and
        masks that share common identifiers. Elements without matching tsIds
        are ignored, ensuring that only coherent tomogram-mask pairs are
        processed. This design is especially important in large cryo-ET
        projects where datasets may originate from multiple acquisition or
        segmentation pipelines.

        Before processing begins, input masks and tomograms are either linked
        directly or converted into MRC format when necessary. This guarantees
        compatibility with the internal image processing workflow and avoids
        unnecessary duplication when files are already stored in a supported
        format.

        Segmentation-Based and Smooth-Mask Processing

        The protocol supports two biologically distinct processing modes.

        In segmentation mode, the input tomomask is interpreted as a labeled
        segmentation map where each integer value represents a distinct
        structural region. The user may define a background label, which is
        excluded from processing. Every remaining label is independently
        transformed into a mesh representation.

        This mode is particularly suitable for semantic segmentations
        generated by manual annotation tools or machine-learning approaches.
        Biological users commonly employ it for separating membranes,
        cytoskeletal filaments, organelles, or compartment boundaries into
        independent spatial objects.

        In smooth-mask mode, the tomomask is interpreted as a continuous
        density or probability map instead of a discrete segmentation. In
        this case, the protocol selects voxels whose values fall within a
        user-defined intensity interval.

        This approach is especially useful when working with probabilistic
        segmentations, confidence maps, or soft masks derived from neural
        network predictions. Rather than relying on discrete labels, the mesh
        is generated from regions satisfying the specified density thresholds.

        Morphological Processing

        Before mesh generation, the protocol optionally applies morphological
        operations to refine the binary mask.

        Dilation can be applied using a spherical footprint with a user-defined
        radius. Biologically, dilation may help reconnect fragmented regions,
        compensate for segmentation discontinuities, or enlarge thin structures
        that would otherwise generate sparse or disconnected meshes.

        Skeletonization can also be enabled. This operation reduces structures
        to their topological backbone while preserving connectivity. In cryo-ET
        workflows, skeletonization is particularly useful for filamentous or
        tubular structures such as actin networks, microtubules, or membrane
        traces where the central geometry is more important than the full
        segmented volume.

        Since these operations are applied sequentially, their combined effect
        strongly influences the geometry of the resulting mesh. Biological users
        should therefore evaluate the processed masks visually to ensure that
        relevant structural features are preserved.

        Mesh Density and Point Sampling

        Once the binary mask has been finalized, the protocol extracts voxel
        coordinates corresponding to the segmented region. Instead of converting
        every voxel into a mesh point, a random subsampling strategy is applied
        according to the user-defined density percentage.

        This parameter determines how many voxels are retained as mesh points.
        Lower densities generate lighter meshes that are easier to visualize and
        manipulate interactively, while higher densities preserve more geometric
        detail at the cost of larger outputs.

        From a biological perspective, sparse meshes are often sufficient for
        representing large cellular structures or membrane networks. However,
        denser sampling may be preferable when fine spatial details or local
        curvature analyses are required.

        Generated mesh points preserve their original spatial coordinates and
        remain associated with the corresponding tomogram. Each processed region
        is additionally assigned a group identifier, allowing independent labels
        or segmented structures to remain distinguishable in downstream analyses.

        Parallel Execution and Robustness

        The protocol executes independently for each tomogram-mask pair using
        parallel processing steps. This design improves scalability and makes
        the protocol suitable for large cryo-electron tomography datasets.

        During execution, failed tomograms are tracked individually. Errors in
        one dataset therefore do not necessarily interrupt processing of the
        remaining tomograms. Conversion failures, incompatible files, or mesh
        generation problems are reported through the logging system for later
        inspection.

        Validation and Consistency Checks

        To ensure geometric consistency, the protocol validates that the
        introduced tomograms and tomomasks share the same sampling rate within
        a small tolerance threshold. This verification is biologically critical
        because mismatched voxel sizes would distort the spatial interpretation
        of the generated meshes.

        The protocol also verifies that valid tomogram-mask associations exist
        through shared tsIds. Non-matching identifiers are excluded from the
        workflow and reported to the user.

        Outputs and Their Interpretation

        The protocol produces a SetOfMeshes object containing mesh points
        associated with their original tomograms. The resulting meshes preserve
        the global spatial organization of the segmented structures while
        reducing volumetric complexity.

        Depending on the selected parameters, the generated meshes may represent
        compact volumetric regions, skeletal traces, membrane surfaces, or sparse
        probabilistic structures. These outputs can subsequently be used for
        visualization, quantitative spatial analysis, structural contextualization,
        or integration with downstream tomographic workflows.

        Practical Recommendations

        For standard semantic segmentations, the default segmentation mode is
        generally appropriate. Users should carefully verify the selected
        background label to avoid accidentally including solvent or empty regions
        in the final mesh.

        When processing probabilistic masks, selecting biologically meaningful
        threshold ranges is essential. Excessively permissive thresholds may
        generate noisy meshes, whereas overly restrictive values may fragment
        biologically continuous structures.

        Skeletonization is particularly recommended for elongated biological
        objects such as cytoskeletal networks or membrane traces. Conversely,
        for compact organelles or volumetric compartments, disabling
        skeletonization may better preserve morphology.

        The density parameter should be adjusted according to the intended use
        of the mesh. Visualization-oriented workflows generally benefit from
        lower densities, while geometric analyses may require denser point
        sampling.

        Final Perspective

        In modern cryo-electron tomography workflows, segmentation alone is
        often insufficient for advanced spatial interpretation. By transforming
        dense voxel masks into structured geometric representations, the Mesh
        From Segmentation protocol provides an efficient bridge between
        segmentation, visualization, and quantitative structural analysis.

        Careful selection of thresholds, morphological operations, and mesh
        density allows biological users to tailor the generated meshes to the
        specific structural properties of their tomographic datasets and the
        scientific questions under investigation.
    """