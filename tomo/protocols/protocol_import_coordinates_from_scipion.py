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
    Imports a set of 3D particle coordinates directly from a Scipion
    sqlite database and associates them with a corresponding set of
    tomograms. The protocol reconstructs the spatial relationship
    between coordinates and tomograms using tomo identifiers or
    filename matching to ensure compatibility with downstream
    tomography workflows.

    AI Generated:

    Import Coordinates 3D From Scipion (ProtImportCoordinates3DFromScipion) — User Manual
        Overview

        The Import Coordinates 3D From Scipion protocol is designed to
        recover and reuse previously generated 3D particle coordinates
        stored in a Scipion sqlite database. Its primary objective is
        to reconnect coordinate information with a new or existing set
        of tomograms while preserving the spatial consistency required
        for subtomogram extraction, averaging, and structural analysis.

        In cryo-electron tomography workflows, coordinates are frequently
        generated during earlier processing stages or in independent
        projects. This protocol allows those coordinates to be imported
        back into the current workflow without the need to regenerate
        particle positions manually. From a biological perspective, this
        facilitates data reuse, reproducibility, and integration between
        multiple tomography processing pipelines.

        The protocol is especially useful in collaborative environments
        or long-term projects where coordinates may have been generated
        in different Scipion sessions but still need to remain associated
        with updated tomograms or refined reconstructions.

        Inputs and General Workflow

        The protocol requires two main inputs: a Scipion sqlite file
        containing a SetOfCoordinates3D object and a set of tomograms
        that will serve as the spatial reference for those coordinates.

        During execution, the protocol loads the coordinate set directly
        from the sqlite database and creates a new output coordinate set
        associated with the provided tomograms. The imported coordinates
        are assumed to already be expressed in the same sampling rate as
        the introduced tomograms. Consequently, no geometric rescaling is
        performed during import.

        The protocol then attempts to establish a relationship between
        coordinates and tomograms. This association is primarily based on
        the tomoId or tsId attributes. If no direct identifier match is
        found, the protocol performs a secondary comparison using the base
        filenames of the tomograms.

        This dual matching strategy increases robustness when coordinates
        originate from different projects or when metadata conventions
        vary slightly across workflows.

        Coordinate-to-Tomogram Matching

        One of the most important aspects of this protocol is the matching
        mechanism between imported coordinates and tomograms. Every 3D
        coordinate must be associated with a valid tomogram reference in
        order to preserve its spatial meaning.

        The protocol first searches for exact tomoId matches between the
        coordinate metadata and the tomogram set. This is the most reliable
        association strategy because it depends on explicit tomography
        identifiers rather than filenames.

        If no tomoId match is found, the protocol attempts to identify the
        corresponding tomogram using filename similarity. This fallback
        mechanism is particularly useful when coordinates were exported or
        migrated between projects where metadata identifiers may have been
        modified or lost.

        From a biological perspective, correct coordinate association is
        essential because an incorrect tomogram assignment would place
        particles into the wrong cellular or structural context. Such
        mismatches could invalidate subtomogram extraction and all
        downstream analyses derived from those particles.

        Handling Non-Matching Coordinates

        Coordinates that cannot be associated with any tomogram are
        automatically excluded from the output set. The protocol generates
        detailed warning messages describing the excluded coordinates,
        including their identifiers and spatial positions.

        This behavior prevents invalid coordinate assignments while still
        allowing valid particles to be imported successfully. The detailed
        reporting also helps users diagnose metadata inconsistencies or
        naming problems between projects.

        If none of the imported coordinates can be matched to the provided
        tomograms, the protocol raises an exception and terminates execution.
        This validation step ensures that biologically meaningless coordinate
        sets are not propagated into subsequent tomography workflows.

        Sampling Rate and Spatial Interpretation

        Unlike protocols that import coordinates from external software,
        this protocol assumes that the imported coordinates already share
        the same sampling rate as the introduced tomograms. As a result,
        coordinate positions are transferred directly without scaling.

        Biologically, this assumption is valid when coordinates and
        tomograms originate from the same Scipion workflow or from
        datasets processed with identical voxel sizes. However, users
        should exercise caution when importing coordinates generated from
        binned or resampled tomograms, since mismatched sampling rates
        could lead to systematic localization errors.

        The protocol also assigns a user-defined box size to the output
        coordinates. This box size defines the extraction region that
        will later be used during subtomogram extraction or particle
        analysis workflows.

        Outputs and Their Interpretation

        After execution, the protocol produces a new SetOfCoordinates3D
        object containing all successfully matched coordinates associated
        with the corresponding tomograms.

        Each coordinate is linked to its tomogram volume reference,
        allowing downstream protocols to correctly interpret particle
        positions in three-dimensional space. The output coordinates
        preserve the original spatial information stored in the sqlite
        database while adapting the dataset to the currently selected
        tomogram set.

        In cases where only a subset of tomograms contains matching
        coordinates, the protocol effectively filters the imported
        coordinate dataset to retain only biologically relevant entries.

        Practical Recommendations

        In practical cryo-ET workflows, it is strongly recommended to
        maintain stable tomoId conventions across projects whenever
        possible. Identifier-based matching is considerably more robust
        than filename matching and reduces the risk of ambiguous
        associations.

        Users should also verify that imported tomograms correspond to
        the same voxel size and preprocessing stage used when the original
        coordinates were generated. Even though the protocol assumes
        compatible sampling rates, incorrect assumptions may compromise
        downstream subtomogram extraction quality.

        After import, visual inspection of several coordinates within the
        tomograms is highly recommended. Confirming that particles appear
        correctly localized provides an effective validation of the matching
        process and helps detect metadata inconsistencies early.

        When importing coordinates from archived or external Scipion
        projects, preserving original filenames and tomography identifiers
        greatly improves workflow reproducibility and minimizes the risk of
        coordinate exclusion.

        Final Perspective

        For tomography users, coordinate import from Scipion sqlite files
        is more than a simple database recovery operation. It represents a
        mechanism for preserving and reusing biologically meaningful spatial
        annotations across workflows and projects. Accurate tomogram matching,
        consistent metadata management, and careful validation of imported
        coordinates are essential to ensure reliable downstream structural
        interpretation in cryo-electron tomography analyses.
    """