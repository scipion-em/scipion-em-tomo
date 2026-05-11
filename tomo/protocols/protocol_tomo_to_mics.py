# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
import enum
from os.path import basename
import mrcfile
import numpy as np
from pwem.objects import SetOfMicrographs, Micrograph, Acquisition
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, GT
from pwem.protocols import EMProtocol
from pyworkflow.utils import replaceExt, removeExt
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfCoordinates3D, Coordinate3D, SetOfTomograms


class ProtTomoToMicsOutput(enum.Enum):
    outputMicrographs = SetOfMicrographs()


class ProtTomoToMics(EMProtocol):
    """ Turns tomograms into set of micrographs to apply SPA picking methods."""

"""
Tomograms to Micrographs and 2D Coordinates to 3D Coordinates — User Manual

    Overview

    These two protocols are designed to work together as a complementary workflow
    that bridges cryo-electron tomography and Single Particle Analysis (SPA)
    methodologies. The first protocol converts tomograms into sets of 2D
    micrographs generated from selected tomographic slices, while the second
    protocol reconstructs the corresponding 3D coordinates from particles picked
    on those micrographs.

    Together, they provide a practical strategy for applying mature 2D particle
    picking tools to tomographic data while preserving the ability to recover
    biologically meaningful spatial information inside the original volume.

    In cryo-electron tomography, direct particle picking within tomograms is often
    difficult because of low signal-to-noise ratios, specimen thickness, and
    crowded intracellular environments. By transforming tomographic slices into
    micrograph-like images, the workflow allows standard SPA particle pickers
    to operate efficiently on tomographic datasets.

    Once particles are detected in 2D, the second protocol restores their
    volumetric context by converting the detected positions back into 3D
    tomographic coordinates. This enables downstream subtomogram extraction,
    averaging, classification, and spatial analysis workflows.

    Inputs and General Workflow

    The workflow begins with a set of tomograms. The Tomograms to Micrographs
    protocol extracts slices along the Z axis and generates 2D micrographs from
    them. The extraction starts from the central slice of each tomogram and
    progresses symmetrically toward upper and lower regions using a configurable
    slice gap.

    For every selected slice, neighboring slices can optionally be averaged
    together. This averaging operation improves contrast and suppresses noise,
    often making biological structures easier to detect automatically.

    Each generated micrograph preserves metadata from the original tomogram,
    including sampling rate and acquisition parameters. Most importantly, the
    micrograph naming convention stores information about the originating
    tomogram and the exact Z slice used during generation.

    After particle picking is performed on the generated micrographs, the
    2D Coordinates to 3D Coordinates protocol reconstructs the corresponding
    volumetric coordinates. The X and Y values are inherited directly from
    the picked particle positions, while the Z coordinate is recovered from
    the encoded slice identifier contained in the micrograph name.

    The resulting output is a SetOfCoordinates3D associated with the original
    tomograms.

    Slice Averaging and Signal Enhancement

    One of the most biologically relevant aspects of this workflow is the
    possibility of averaging adjacent tomographic slices before generating
    each micrograph.

    Cryo-electron tomograms are typically noisy, especially in cellular studies
    where sample thickness and radiation damage strongly reduce image contrast.
    Averaging neighboring slices improves visibility of macromolecular features
    and frequently leads to more robust particle detection.

    This strategy is particularly useful for ribosomes, membrane proteins,
    viral assemblies, and other complexes embedded in crowded intracellular
    environments.

    However, averaging too many slices may blur axial information and reduce
    the precision of the reconstructed Z coordinate. Biological users should
    therefore balance noise reduction against spatial specificity depending
    on the intended downstream analysis.

    Spatial Sampling Along the Z Axis

    The slice gap parameter determines how densely the tomogram is sampled.
    Small gaps generate a larger number of micrographs and provide more complete
    volumetric coverage, although they increase computational and storage costs.

    Larger gaps reduce redundancy and accelerate downstream processing but may
    skip biologically relevant regions if particles are sparsely distributed
    throughout the tomogram thickness.

    Since the reconstructed Z coordinate directly depends on the selected
    slices, the choice of slice gap also influences the final spatial accuracy
    of the workflow.

    Coordinate Reconstruction and Biological Interpretation

    The reconstruction step restores the three-dimensional spatial organization
    that was temporarily simplified during the 2D particle-picking stage.

    This is biologically important because many cryo-ET studies rely not only
    on identifying particles, but also on understanding their spatial
    distribution relative to membranes, organelles, cytoskeletal structures,
    or neighboring macromolecular complexes.

    By reconnecting 2D particle detections with their original tomographic
    positions, the workflow enables both structural and spatial biological
    analyses while maintaining compatibility with standard subtomogram
    processing pipelines.

    Outputs and Their Interpretation

    The workflow produces two main outputs. First, a SetOfMicrographs generated
    from tomographic slices, intended for particle picking and visualization.
    Second, a SetOfCoordinates3D containing the reconstructed volumetric
    particle positions.

    The generated micrographs should not be interpreted as independent
    experimental acquisitions, but rather as computational representations
    derived from tomographic data to facilitate particle detection.

    The reconstructed coordinates preserve the association with the original
    tomograms and remain fully compatible with subtomogram extraction and
    averaging workflows.

    Practical Recommendations

    In most biological applications, moderate slice averaging combined with
    intermediate slice gaps provides a good balance between contrast enhancement,
    computational efficiency, and spatial precision.

    Highly noisy tomograms may benefit from stronger averaging, whereas studies
    focused on accurate spatial localization should minimize both slice averaging
    and slice gaps.

    Users should visually inspect both the generated micrographs and the final
    reconstructed coordinates to ensure that biological structures remain well
    represented throughout the workflow.

    Final Perspective

    Together, the Tomograms to Micrographs and 2D Coordinates to 3D Coordinates
    protocols provide an effective hybrid strategy that combines the robustness
    of SPA-style particle picking with the volumetric richness of cryo-electron
    tomography.

    By temporarily transforming tomographic information into a 2D representation
    and subsequently restoring its original 3D spatial context, the workflow
    enables efficient particle detection while preserving the biological meaning
    of the tomographic environment.
"""