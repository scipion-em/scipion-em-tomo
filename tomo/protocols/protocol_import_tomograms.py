# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
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
from os.path import abspath, basename, join
from pwem.convert.headers import Ccp4Header
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.utils.path import createAbsLink, removeBaseExt, getExt
import pyworkflow.protocol.params as params
from .protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from ..convert.mdoc import normalizeTSId
from ..objects import Tomogram, SetOfTomograms

logger = logging.getLogger(__name__)
OUTPUT_NAME = 'Tomograms'


class ProtImportTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """Protocol to import a set of tomograms to the project"""


    """
    Imports a set of subtomograms into a Scipion project together with
    their associated acquisition metadata and spatial origin information.
    The protocol supports volumetric formats commonly used in cryo-electron
    tomography workflows and prepares the imported subtomograms for downstream
    visualization, classification, averaging, or refinement procedures.

    AI Generated:

    Import Subtomograms (ProtImportSubTomograms) — User Manual

        Overview

        The Import Subtomograms protocol is designed to incorporate a collection
        of previously generated subtomograms into a Scipion tomography project.
        In cryo-electron tomography workflows, subtomograms usually correspond
        to small 3D regions extracted from larger tomograms and centered on
        biological particles, macromolecular complexes, membranes, or repeating
        cellular structures. These subtomograms often represent the starting
        point for subtomogram averaging, alignment, classification, or structural
        heterogeneity analysis.

        From a biological perspective, this protocol acts as a bridge between
        external subtomogram extraction procedures and downstream structural
        analysis inside Scipion. The protocol does not perform particle picking
        or extraction itself. Instead, it imports already generated subtomogram
        volumes while preserving important metadata such as voxel size,
        acquisition parameters, and spatial origin.

        The protocol is especially useful when subtomograms were generated using
        external tomography software packages or when data produced in previous
        workflows must be integrated into a unified Scipion processing pipeline.

        Inputs and General Workflow

        The protocol requires a collection of volumetric files corresponding to
        subtomograms. These files are identified using a filename pattern, which
        allows the protocol to iterate automatically through all matching entries.
        Each imported volume is registered internally as a SubTomogram object and
        incorporated into a SetOfSubTomograms container.

        During the import process, the protocol determines the dimensions of each
        volume and calculates a spatial origin centered within the box. This
        origin definition is particularly important in subtomogram averaging
        workflows because alignment algorithms assume that particles are roughly
        centered within the subtomogram reference frame.

        The voxel size, provided as the sampling rate, is propagated to the full
        dataset and stored for downstream processing. Correct sampling rate
        definition is biologically critical because it directly affects
        structural interpretation, resolution estimation, and compatibility with
        subsequent refinement procedures.

        File Formats and Volume Interpretation

        The protocol supports standard volumetric formats commonly used in cryo-EM
        and tomography workflows, including MRC and MAP files. Special handling
        is implemented for stack-like MRC files in which multiple subtomograms
        may be stored within a single container.

        When the imported file contains multiple volumes, the protocol separates
        them logically into individual subtomogram entries while preserving the
        correct indexing information. This behavior is particularly useful for
        large-scale tomography projects where particle extraction software exports
        many particles into consolidated stacks rather than independent files.

        MAP files are internally interpreted as MRC-compatible volumes to ensure
        interoperability across cryo-EM software ecosystems.

        Acquisition Metadata Handling

        One of the central features of this protocol is the integration of
        acquisition metadata inherited from tomography experiments. The protocol
        parses acquisition parameters and associates them with each imported
        subtomogram.

        In biological workflows, acquisition metadata can become highly relevant
        when comparing datasets collected under different microscope conditions,
        magnifications, or imaging strategies. Preserving this information helps
        maintain consistency throughout downstream refinement and averaging
        procedures.

        The protocol therefore not only imports image data but also preserves
        experimental context, which is essential for reproducible cryo-ET
        analysis.

        Spatial Origin and Coordinate System

        The protocol automatically assigns a centered spatial origin to every
        subtomogram based on the dimensions of the imported volume and the
        sampling rate. In practice, this means that the geometric center of the
        subtomogram becomes the reference point for subsequent alignment and
        averaging operations.

        This origin assignment is biologically important because many subtomogram
        alignment methods assume that the target macromolecule is approximately
        centered in the extraction box. Incorrect origins may lead to unstable
        alignments, poor averages, or inaccurate structural interpretation.

        The protocol computes the shifts in physical units using the voxel size,
        ensuring consistency between image geometry and spatial calibration.

        Integration into Cryo-ET Workflows

        Imported subtomograms are commonly used in subtomogram averaging pipelines
        aimed at increasing signal-to-noise ratio and recovering high-resolution
        structural information from noisy tomographic data.

        In practical biological applications, subtomograms may correspond to
        ribosomes, viral spikes, membrane channels, cytoskeletal assemblies, or
        large molecular complexes inside native cellular environments. Once
        imported, these particles can be aligned, classified, or averaged to
        identify structural states and conformational variability.

        Because the protocol preserves acquisition information and spatial
        consistency, it facilitates reliable integration with downstream Scipion
        tomography protocols.

        Coordinate Association

        The code structure includes placeholders for future or optional
        association between subtomograms and previously imported 3D coordinates.
        Although this functionality is currently disabled, the design indicates
        support for workflows where each subtomogram can be directly linked to
        its original particle coordinate inside the tomogram.

        Such associations are biologically valuable because they preserve the
        spatial context of particles within the cellular environment, enabling
        correlation between structural information and native localization.

        Outputs and Their Interpretation

        After execution, the protocol generates a SetOfSubTomograms object
        containing all successfully imported subtomograms together with their
        associated metadata, origins, and acquisition parameters.

        Each subtomogram maintains its file reference and geometric information,
        allowing it to be used immediately in visualization, averaging,
        classification, or refinement workflows.

        The resulting dataset becomes a standardized Scipion-compatible container
        suitable for large-scale cryo-electron tomography analysis.

        Validation and Data Integrity

        Before import begins, the protocol validates the existence of files
        matching the provided pattern. If no matching files are detected, the
        execution stops with an error message.

        This validation step prevents incomplete imports and ensures that users
        are aware of incorrect file patterns or missing datasets before starting
        computationally expensive downstream analyses.

        Practical Recommendations

        In practical cryo-ET workflows, users should verify that all imported
        subtomograms share a consistent voxel size and box dimensions before
        proceeding to averaging or classification. Mixing particles extracted at
        different sampling rates may lead to incorrect alignments or unreliable
        structural interpretation.

        It is also advisable to confirm that particles are approximately centered
        within the extraction box. Although the protocol assigns a centered
        origin automatically, strongly off-centered particles may still require
        re-extraction or additional preprocessing.

        For large datasets stored as MRC stacks, users should ensure that the
        indexing and dimensionality are correctly interpreted after import,
        particularly when subtomograms originate from external software packages.

        Final Perspective

        For cryo-electron tomography users, importing subtomograms is more than a
        simple data-loading operation. It represents the transition from raw
        extracted particles to biologically meaningful structural analysis.
        Proper definition of voxel size, acquisition metadata, and spatial origin
        is essential for obtaining reliable averages and interpretable structural
        results.

        By standardizing subtomogram datasets inside Scipion, this protocol
        provides the foundation for robust subtomogram averaging workflows and
        integrative structural studies within native cellular environments.
    """