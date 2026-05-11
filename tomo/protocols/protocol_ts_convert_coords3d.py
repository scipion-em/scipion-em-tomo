# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.constants import CENTER_GRAVITY
from tomo.protocols import ProtTomoBase

METADATA_INPUT_COORDINATES = "fiducialCoordinates.xmd"


class ProtTsConvertCoordinates3d(EMProtocol, ProtTomoBase):
    """
    Scipion protocol to convert a set of tilt-series coordinates 3d to a set of coordinates 3d associated to a set of
    tomograms.
    """

    """
    Converts 3D coordinates associated with tilt-series into 3D coordinates
    linked to reconstructed tomograms. The protocol transfers spatial coordinate
    information from the tilt-series reference frame into the tomogram reference
    system while preserving positional and scoring metadata required for
    downstream tomographic analysis.

    AI Generated:

    Tilt-Series Convert Coordinates3D (ProtTsConvertCoordinates3d) — User Manual
        Overview

        The Tilt-Series Convert Coordinates3D protocol transforms a set of 3D
        coordinates associated with tilt-series into a standard set of 3D
        coordinates linked directly to reconstructed tomograms. Its primary goal
        is to reconnect coordinate information generated during tilt-series
        alignment or fiducial-tracking procedures with the final tomographic
        reconstruction space.

        In cryo-electron tomography workflows, several intermediate processing
        stages operate in the coordinate system of the tilt-series rather than in
        the reconstructed tomogram itself. Fiducial positions, alignment markers,
        or tracked particles are often stored relative to the tilt-series geometry
        before reconstruction is completed. However, downstream structural
        interpretation, subtomogram extraction, visualization, and spatial analysis
        require coordinates expressed directly within the tomographic reference
        frame.

        This protocol bridges that transition by converting tilt-series-based
        coordinates into tomogram-associated coordinates while preserving their
        spatial meaning and metadata relationships.

        Inputs and General Workflow

        The protocol requires two inputs: a set of 3D coordinates associated with
        tilt-series and a corresponding set of reconstructed tomograms.

        During execution, the protocol identifies the tomogram associated with
        each coordinate using the tilt-series identifier. A correspondence table
        between tilt-series identifiers and tomograms is internally generated to
        ensure that every coordinate is assigned to the correct reconstructed
        volume.

        Once the appropriate tomogram is identified, the protocol converts the
        coordinate positions into the tomogram spatial reference system using the
        tomogram sampling rate. The transformed coordinates are then stored as
        standard Coordinate3D objects associated with their respective tomograms.

        In addition to positional information, the protocol preserves scoring
        metadata from the original coordinates, allowing confidence or quality
        measurements generated during alignment procedures to remain available in
        downstream analyses.

        Coordinate Conversion and Spatial Interpretation

        Biologically, this protocol restores the direct relationship between
        detected spatial features and the reconstructed tomographic volume. While
        tilt-series coordinates are geometrically meaningful during alignment,
        tomogram-associated coordinates are necessary for interpreting the spatial
        organization of macromolecular complexes inside the reconstructed specimen.

        The conversion process rescales coordinate values according to the
        tomogram sampling rate, ensuring that the resulting coordinates are
        correctly expressed in tomographic voxel space. This step is especially
        important when reconstruction binning or sampling modifications have been
        applied during tomogram generation.

        Because the protocol preserves the original spatial associations, the
        resulting coordinates remain compatible with subtomogram extraction,
        particle averaging, fiducial visualization, and spatial distribution
        studies.

        Sampling Rate and Coordinate Accuracy

        The accuracy of the converted coordinates depends directly on the
        consistency between the tilt-series geometry and the reconstructed
        tomograms. If tomograms have been reconstructed using different binning
        factors or sampling rates, proper coordinate scaling becomes essential to
        maintain geometric consistency.

        The protocol automatically uses the sampling rate of the input tomograms
        during coordinate transformation. This ensures that coordinate positions
        remain synchronized with the physical dimensions of the reconstructed
        volumes.

        From a practical perspective, users should verify that the input
        tilt-series coordinates and tomograms originate from compatible alignment
        and reconstruction workflows. Inconsistent sampling conventions or
        mismatched tilt-series identifiers may produce biologically incorrect
        spatial mappings.

        Outputs and Their Interpretation

        The protocol produces a SetOfCoordinates3D directly associated with the
        input tomograms. Each coordinate contains transformed X, Y, and Z
        positions together with the tomographic reference and any inherited score
        information.

        These coordinates can subsequently be used for subtomogram extraction,
        fiducial inspection, structural averaging, or spatial organization
        analyses within Scipion tomography workflows.

        From a biological perspective, the resulting output represents a validated
        mapping between alignment-derived spatial information and the final
        reconstructed tomographic environment.

        Practical Recommendations

        In routine cryo-ET processing, users should ensure that the tilt-series
        coordinates were generated from the same acquisition and alignment context
        as the tomograms used for conversion. Maintaining consistent tilt-series
        identifiers throughout the workflow is essential for accurate coordinate
        assignment.

        It is also advisable to visually inspect a subset of converted coordinates
        within the reconstructed tomograms to confirm that fiducials or particles
        appear correctly positioned in three-dimensional space.

        When tomograms have been reconstructed using aggressive binning or
        rescaling, users should pay particular attention to coordinate precision,
        especially in workflows requiring accurate subtomogram localization or
        quantitative spatial analysis.

        Final Perspective

        The Tilt-Series Convert Coordinates3D protocol provides an essential
        bridge between tilt-series alignment procedures and tomogram-centered
        structural analysis workflows. By converting alignment-derived spatial
        coordinates into tomogram-associated coordinates, the protocol enables
        coherent downstream interpretation of three-dimensional biological
        structures while preserving geometric consistency across the tomography
        processing pipeline.
    """