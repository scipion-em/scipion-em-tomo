# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
import numpy as np
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from .protocol_base import ProtTomoPicking
import tomo.constants as const
from ..objects import Coordinate3D, SetOfCoordinates3D, SetOfSubTomograms, SetOfTomograms, SubTomogram


class Output3dCoordExtraction(enum.Enum):
    coordinates3d = SetOfCoordinates3D


class ProtTomoExtractCoords(ProtTomoPicking):
    """
    Extract the coordinates information from a set of subtomograms.

    This protocol is useful when we want to re-extract the subtomograms
    (maybe resulting from classification) with the
    original dimensions. It can be also handy to visualize the resulting
    subtomograms in their location on the tomograms.
    """


class ProtTomoExtractCoords(ProtTomoPicking):
    """
    ProtTomoExtractCoords — Extract 3D Coordinates Protocol

    Overview
    --------
    Extracts 3D coordinates from a set of subtomograms or pre-existing coordinates.
    This protocol allows re-extraction of subtomograms at their original dimensions
    and visualization of their positions within the tomograms.
    Common use cases include:
        - Re-extracting subtomograms after classification.
        - Associating subtomograms with updated tomograms.
        - Visualizing coordinate distributions on tomograms.

    Inputs and Workflow
    -------------------
    - Subtomograms or 3D Coordinates: Input objects from which coordinates will be extracted.
      Practical tips:
        * Input can be either SetOfSubTomograms or SetOfCoordinates3D.
        * Ensure coordinates are already associated with initial tomograms.
    - Tomograms: Target tomograms to associate extracted coordinates.
    - Box Size: Optional; defines the extraction box size. Defaults are taken from the input coordinates.

    The main workflow:
        1. Define parameters (input sets and optional box size).
        2. Extract coordinates from each item, scaling positions and shifts as necessary.
        3. Create output SetOfCoordinates3D with updated positions and transformations.

    Coordinate Extraction Details
    -----------------------------
    - Input subtomograms: Coordinates are scaled to match tomogram sampling, and subtomogram shifts are applied.
    - Input coordinates: Both coordinates and shifts are scaled according to tomogram sampling.
    - Two scale factors are used: one for coordinates, one for shifts.

    Outputs
    -------
    - SetOfCoordinates3D: Contains all extracted coordinates with positions and transformations adjusted.
    - Maintains links to input subtomograms and tomograms for provenance.

    Practical Recommendations
    -------------------------
    - Check input subtomograms to ensure coordinates exist.
    - Verify tomogram associations if using filename-based matching.
    - Carefully select box size when re-extracting subtomograms for consistency.

    Utilities
    ---------
    - Supports temporary output paths and file suffixes for intermediate results.
    - Provides methods to retrieve tomograms and coordinates from inputs.
    - Handles both subtomogram and 3D coordinate inputs seamlessly.

    Biological Perspective
    ---------------------
    - Coordinates extraction is essential for precise localization of subtomograms.
    - Enables downstream analyses such as averaging, visualization, or further tomographic processing.
    - Proper scaling ensures spatial consistency between subtomograms and their tomograms.

    """