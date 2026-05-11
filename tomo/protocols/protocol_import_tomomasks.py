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
import logging
from enum import Enum
from os.path import abspath, join
from pwem.emlib.image import ImageHandler
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import yellowStr
from pyworkflow.utils.path import removeBaseExt, getExt, createAbsLink
from .protocol_base import ProtTomoImportFiles
from ..constants import ERR_NO_TOMOMASKS_GEN
from ..convert.mdoc import normalizeTSId
from ..objects import TomoMask, SetOfTomoMasks


logger = logging.getLogger(__name__)


class importTomoMasksOutputs(Enum):
    tomomasks = SetOfTomoMasks


class ProtImportTomomasks(ProtTomoImportFiles):
    """Protocol to import a set of tomomasks (segmentations) to the project"""


    """
    Imports a set of tomographic segmentation masks into a Scipion project
    and associates them with their corresponding tomograms. The protocol
    validates dimensional consistency between tomograms and masks and
    prepares the segmentation volumes for downstream visualization,
    annotation, structural interpretation, and quantitative analysis.

    AI Generated:

    Import Tomomasks (ProtImportTomomasks) — User Manual

        Overview

        The Import Tomomasks protocol is designed to incorporate segmented
        tomographic volumes, also referred to as tomomasks, into a Scipion
        tomography workflow. In cryo-electron tomography, segmentation masks
        are commonly used to identify biologically meaningful regions inside
        tomograms, such as membranes, ribosomes, cytoskeletal filaments,
        organelles, viral particles, vesicles, or other macromolecular
        assemblies.

        From a biological perspective, segmentation masks provide a way to
        transform noisy tomographic reconstructions into interpretable
        structural information. These masks are frequently generated using
        manual annotation tools, machine learning segmentation software, or
        specialized membrane annotation pipelines. The protocol enables these
        segmented regions to be imported and formally associated with the
        original tomograms inside the Scipion environment.

        The protocol therefore acts as a bridge between segmentation workflows
        and downstream tomographic analysis, facilitating visualization,
        particle extraction, contextual interpretation, and structural studies
        within native cellular environments.

        Inputs and General Workflow

        The protocol requires two main inputs: a set of tomograms already
        present in the project and a collection of segmentation mask files.
        The masks are identified using either direct filename patterns or
        regular expression matching strategies.

        During execution, the protocol scans the provided directory structure,
        identifies compatible mask files, and attempts to associate each mask
        with its corresponding tomogram through normalized tilt-series
        identifiers. This matching strategy is particularly important in large
        tomography projects where multiple segmentations may coexist for
        different tomograms or experimental conditions.

        Once the associations are established, the protocol generates a
        SetOfTomoMasks object that stores the imported segmentation masks
        together with their tomogram references and metadata.

        Segmentation Masks in Biological Context

        In practical cryo-ET workflows, tomomasks often represent biologically
        relevant structures isolated from the surrounding cellular environment.
        For example, segmentation masks may define membrane boundaries,
        cytoplasmic compartments, viral envelopes, protein assemblies, or
        filamentous networks.

        These masks are essential for contextual structural biology because
        they allow researchers to focus computational analyses on specific
        cellular regions while excluding irrelevant density or reconstruction
        artifacts. In many modern workflows, segmentation masks also guide
        subtomogram extraction, region-based averaging, and machine learning
        analysis.

        By importing segmentation masks directly into Scipion, the protocol
        enables integrated structural interpretation in which reconstructed
        densities and annotated biological features coexist within the same
        processing environment.

        File Matching and Identification

        The protocol supports flexible matching mechanisms to associate masks
        with tomograms. When regular expressions are enabled, filenames are
        parsed dynamically to identify corresponding tilt-series identifiers.
        Alternatively, direct filename matching can be used when mask names
        already follow a consistent naming convention.

        Additional normalization rules are implemented to support segmentation
        outputs generated by external software packages. For example, suffixes
        such as "_materials" or "_segmented", commonly introduced by membrane
        annotation or neural-network segmentation tools, are automatically
        removed before attempting the association with the tomogram.

        This flexibility is particularly valuable in collaborative projects
        where segmentations may originate from different annotation pipelines
        or external software ecosystems.

        Dimensional Consistency and Validation

        One of the most biologically important aspects of the protocol is the
        validation of dimensional consistency between tomograms and their
        associated masks. For each tomomask, the protocol verifies that the
        X, Y, and Z dimensions exactly match those of the corresponding
        tomogram.

        This validation is critical because segmentation masks must occupy the
        same spatial coordinate system as the tomographic reconstruction. Any
        dimensional mismatch would invalidate the biological interpretation of
        the segmentation and could lead to incorrect localization of structures
        or downstream computational errors.

        Masks with incompatible dimensions are automatically excluded from the
        output dataset, and detailed warning messages are generated to inform
        the user about the inconsistency.

        Sampling Rate and Spatial Consistency

        The protocol propagates the sampling rate from the input tomograms to
        the imported tomomasks. Maintaining voxel-size consistency is essential
        for correct structural interpretation because segmentation masks must
        remain spatially aligned with the tomographic data.

        In biological analyses, incorrect voxel calibration could produce
        misleading measurements of membrane thickness, particle dimensions,
        organelle size, or intermolecular distances. By preserving the original
        tomogram sampling rate, the protocol guarantees geometric consistency
        throughout downstream analyses.

        Data Organization and File Management

        During import, the protocol creates internal symbolic links to the
        segmentation masks inside the Scipion project structure. This approach
        avoids unnecessary duplication of large volumetric datasets while still
        preserving reproducibility and project portability.

        Each imported tomomask is stored with its associated tomogram reference,
        allowing direct integration with visualization and segmentation-aware
        workflows inside Scipion.

        Outputs and Their Interpretation

        After successful execution, the protocol produces a SetOfTomoMasks
        object containing all validated segmentation masks together with their
        corresponding tomogram associations.

        Each tomomask preserves information about the original segmentation file,
        voxel size, and linked tomogram volume. The resulting dataset can then
        be used for structural interpretation, segmentation visualization,
        contextual subtomogram extraction, or downstream quantitative analyses.

        From a biological perspective, the imported masks provide an annotated
        structural map of the tomographic environment, helping researchers
        interpret molecular organization directly inside native cellular
        contexts.

        Error Handling and Warnings

        The protocol includes several validation and warning mechanisms designed
        to prevent biologically inconsistent imports. If no files match the
        provided pattern, execution stops with an explicit error. Similarly,
        masks with incompatible dimensions are excluded automatically.

        Informative warning messages are stored and displayed so users can
        identify which masks failed validation and why. This behavior is
        especially useful in large segmentation projects where mismatches may
        originate from preprocessing differences, cropping operations, or
        inconsistent reconstruction parameters.

        Practical Recommendations

        In practical cryo-ET workflows, users should ensure that segmentation
        masks are generated directly from the same tomograms that will be used
        during import. Even small dimensional differences introduced by binning,
        cropping, or resampling can invalidate the spatial correspondence
        between masks and tomograms.

        Consistent naming conventions are also strongly recommended, especially
        in large datasets containing many tomograms and segmentation outputs.
        Maintaining coherent identifiers greatly simplifies automatic matching
        and reduces the risk of incorrect associations.

        When using segmentation masks generated by machine learning pipelines,
        users should visually inspect several imported masks after execution to
        verify that spatial alignment and biological interpretation remain
        correct.

        Final Perspective

        For cryo-electron tomography users, segmentation masks are not merely
        auxiliary files but biologically meaningful annotations that define the
        structural context of tomographic data. Proper integration of these
        masks is essential for accurate interpretation of cellular organization,
        molecular localization, and structural heterogeneity.

        By validating spatial consistency, preserving tomogram associations,
        and integrating segmentation data into Scipion workflows, this protocol
        provides a reliable foundation for advanced contextual structural
        biology analyses within cryo-electron tomography projects.
    """
