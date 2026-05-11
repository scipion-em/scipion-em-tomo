# **************************************************************************
# *
# * Authors:     Oier Lauzirika Zarrabeitia (oierlauzi@bizkaia.eu) [1]
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
from pwem.protocols import EMProtocol

from tomo.objects import (SetOfTiltSeries, TiltImage,
                          SetOfCoordinates3D, Coordinate3D,
                          SetOfLandmarkModels, LandmarkModel )
from tomo.protocols import ProtTomoBase
import tomo.constants as constants
from tomo.utils import getObjFromRelation

import enum
import numpy as np

class OutputProjectCoordinates(enum.Enum):
    landmarkModels = SetOfLandmarkModels


class ProtProjectCoordinates(EMProtocol, ProtTomoBase):
    """
    Project 3D coordinates into a set of landmarks.
    """

    """
    Projects 3D coordinates onto a tilt-series to generate landmark models.
    The protocol converts spatial 3D coordinate information into corresponding
    2D landmark positions across all images of a tilt-series, enabling the
    generation of fiducial-like landmarks suitable for tomographic alignment,
    visualization, or tracking workflows.

    AI Generated:

    Project Coordinates (ProtProjectCoordinates) — User Manual
        Overview

        The Project Coordinates protocol projects a set of 3D coordinates onto
        one or more tilt-series in order to generate landmark models. Its main
        purpose is to transform volumetric spatial information into a set of
        consistent 2D landmark positions distributed across the angular views
        of a tomographic acquisition. In cryo-electron tomography workflows,
        this protocol is particularly useful for simulating fiducial markers,
        validating alignment strategies, generating synthetic tracking datasets,
        or linking structural coordinates with experimental tilt-series geometry.

        From a biological and tomographic perspective, the protocol establishes
        the geometric relationship between 3D particle locations and their
        expected appearance across a tilt-series. This allows researchers to
        study how spatial coordinates propagate through tomographic projections
        and how alignment transformations affect their observed positions.

        Inputs and General Workflow

        The protocol requires a set of 3D coordinates as its primary input.
        These coordinates usually correspond to particle positions, structural
        landmarks, or biologically relevant points identified inside tomograms.
        Optionally, the user may also provide a set of tilt-series on which the
        coordinates will be projected. If the tilt-series are not explicitly
        specified, the protocol attempts to deduce them automatically from the
        coordinate relationships stored in the dataset.

        During execution, the protocol iterates through each tilt-series and
        identifies all coordinates associated with that tomographic dataset.
        Every 3D coordinate is converted into homogeneous coordinates and scaled
        to match the sampling rate of the tilt-series. The protocol then applies
        a geometric projection model based on the tilt angle of each image in
        the series, generating the corresponding 2D landmark position for every
        projection view.

        The resulting projected positions are stored as landmark models linked
        to the original tilt-series. These landmarks can later be used in
        alignment, visualization, or geometric validation workflows.

        Geometric Projection Model

        The protocol relies on a projection matrix derived from the tilt angle
        associated with each tilt image. This matrix models the tomographic
        acquisition geometry by rotating the 3D coordinate system according to
        the tilt angle and projecting the resulting positions into the 2D image
        plane.

        Biologically, this operation reproduces how a particle or structural
        feature would appear during the experimental acquisition process. The
        protocol therefore provides a mathematically consistent bridge between
        volumetric coordinates and experimental tilt projections.

        If the tilt image already contains alignment transformations, the
        protocol additionally compensates for these transformations by applying
        the inverse transformation matrix to the projected coordinates. This
        ensures that landmark positions remain consistent with the aligned or
        transformed state of the tilt-series.

        Coordinate Scaling and Sampling Consistency

        An important aspect of the protocol is the handling of sampling rates.
        Coordinates and tilt-series may originate from datasets with different
        binning levels or voxel sizes. To preserve geometric consistency, the
        protocol rescales the coordinates according to the ratio between the
        coordinate sampling rate and the tilt-series sampling rate.

        From a practical perspective, this step is essential when combining
        datasets generated at different resolutions or after preprocessing
        operations such as binning or resampling. Without proper scaling,
        projected landmarks would appear shifted or geometrically inconsistent
        across the tilt-series.

        Landmark Model Generation

        For every tilt-series, the protocol creates an independent landmark
        model containing the projected 2D positions associated with all input
        coordinates. Each landmark is linked to its originating coordinate
        through a chain identifier, preserving the correspondence between the
        original 3D point and all of its 2D projections.

        This organization allows the resulting landmarks to behave similarly to
        fiducial trajectories or tracked particles across angular views. Such
        information can be valuable for alignment benchmarking, motion analysis,
        or synthetic dataset generation.

        Outputs and Their Interpretation

        After execution, the protocol produces a set of landmark models
        associated with the input tilt-series. Each landmark model contains the
        projected 2D coordinates corresponding to the original 3D positions
        across all tilt images.

        Biologically and computationally, these outputs represent the expected
        trajectories of spatial features during tomographic acquisition. They
        can therefore be interpreted as synthetic fiducials, geometric
        references, or coordinate-based tracking models suitable for downstream
        tomographic workflows.

        The generated landmark models preserve the relationship with both the
        original coordinates and the associated tilt-series, ensuring full
        traceability throughout the processing pipeline.

        Practical Recommendations

        In practical cryo-ET workflows, users should ensure that the coordinate
        set and tilt-series belong to the same tomographic reference frame.
        Significant inconsistencies in alignment, sampling rate, or coordinate
        origin may produce inaccurate landmark projections.

        The protocol is especially useful for testing alignment algorithms,
        validating geometric transformations, or generating controlled synthetic
        datasets for methodological development. When working with transformed
        tilt-series, preserving accurate transformation matrices is important to
        ensure geometrically correct landmark placement.

        Care should also be taken when interpreting projected coordinates near
        the tomogram boundaries, since projection effects at high tilt angles
        may move landmarks outside the visible field of view.

        Final Perspective

        For cryo-electron tomography users, the Project Coordinates protocol is
        fundamentally a geometric transformation tool that connects volumetric
        structural information with experimental tilt-series projections. By
        translating 3D coordinates into consistent 2D landmark trajectories,
        the protocol enables realistic modeling of tomographic acquisition
        geometry and facilitates a wide range of alignment, validation, and
        simulation workflows. Proper handling of coordinate scaling, tilt
        geometry, and transformation matrices is essential for obtaining
        biologically and geometrically meaningful results.
    """