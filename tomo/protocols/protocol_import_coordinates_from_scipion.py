# *
# * Authors:     Scipion Team
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
from enum import Enum
from os.path import exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol import FileParam, IntParam, PointerParam
from pyworkflow.utils import Message, removeBaseExt, yellowStr
from ..constants import SCIPION, ERR_COORDS_FROM_SQLITE_NO_MATCH
from .protocol_base import ProtTomoBase
from ..objects import SetOfCoordinates3D


class outputObjs(Enum):
    coordinates = SetOfCoordinates3D


class ProtImportCoordinates3DFromScipion(EMProtocol, ProtTomoBase):
    """Protocol to import a set of 3d coordinates from Scipion sqlite file"""

    """
    ProtImportCoordinates3DFromScipion — Import 3D Coordinates Protocol

    Overview
    --------
    Imports 3D coordinates from a Scipion SQLite file and maps them to a set of tomograms. 
    This allows previously extracted coordinates to be reused for subtomogram extraction, 
    visualization, or further processing while maintaining spatial consistency.

    Inputs and Workflow
    -------------------
    - SQLite File: Contains the 3D coordinates to import.
      Practical tips:
        * Ensure the file exists and is readable.
        * Verify that coordinates reference tomograms using tsId/tomoId or filename.
    - Input Tomograms: Tomograms to which the coordinates will be associated.
      Practical tips:
        * Matching is attempted via tsId/tomoId first, then by filename if needed.
        * Coordinates are assumed to have the same sampling rate as the tomograms.
    - Box Size: Integer defining the extraction box size in pixels (default 20).

    Workflow Steps
    --------------
    1. Load coordinates from the SQLite file into a SetOfCoordinates3D object.
    2. Set sampling rate and box size from input tomograms.
    3. Match coordinates to tomograms by tsId/tomoId or filename.
    4. Exclude and log coordinates that do not match any tomogram.
    5. Assign volume pointers to each coordinate for spatial reference.
    6. Define outputs and maintain provenance with input tomograms.

    Matching and Validation
    -----------------------
    - Coordinates are matched to tomograms via tsId/tomoId or filename.
    - Non-matching coordinates are excluded and reported in detailed logs.
    - Raises an error if no coordinates match any input tomogram.

    Outputs
    -------
    - SetOfCoordinates3D: Coordinates successfully imported and assigned to volumes.
    - Provenance relationships to input tomograms are preserved.

    Practical Recommendations
    -------------------------
    - Confirm the SQLite file exists and contains valid coordinate data.
    - Ensure all tomograms referenced by coordinates are provided.
    - Maintain consistent sampling rates between coordinates and tomograms.

    Biological Perspective
    ---------------------
    - Enables reuse of previously extracted coordinates for consistent analysis.
    - Maintains spatial context for accurate subtomogram extraction and downstream analysis.
    - Supports reproducibility and integrity in structural biology workflows.
    """