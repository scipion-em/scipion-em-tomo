# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import NEW
from pwem.protocols import EMProtocol

from tomo.constants import BOTTOM_LEFT_CORNER, TR_DYNAMO
from tomo.utils import existsPlugin

EXPORT_TO_TXT = 'txt'
EXPORT_TO_STAR = 'relion'
EXPORT_TO_EMAN = 'eman'
EXPORT_TO_DYNAMO = 'dynamo'
EXPORT_TO_CBOX = 'cbox'


class ProtExportCoordinates3D(EMProtocol):
    """ Export 3D subtomogram coordinates to be used outside Scipion. """


    """
    ProtExportCoordinates3D — Export 3D Subtomogram Coordinates Protocol

    This protocol exports 3D subtomogram coordinates generated inside Scipion into
    external formats commonly used by different cryo-electron tomography software
    packages. Its main purpose is to facilitate interoperability between Scipion
    and external processing environments, allowing users to continue analysis,
    visualization, subtomogram averaging, or particle extraction workflows in
    other specialized software ecosystems.

    AI Generated:

    Export 3D Coordinates (ProtExportCoordinates3D) — User Manual

        Overview

        The Export 3D Coordinates protocol converts a set of subtomogram coordinates
        stored in Scipion into external file formats compatible with several tomography
        processing packages. In practical cryo-ET workflows, this protocol is commonly
        used when coordinates generated during particle picking or annotation inside
        Scipion need to be transferred to external software for downstream analysis.

        Depending on the installed plugins, the protocol supports exporting coordinates
        to plain TXT files, RELION STAR files, EMAN JSON files, Dynamo tables, and
        SPHIRE CBOX files. This flexibility allows users to integrate Scipion into
        heterogeneous cryo-electron tomography pipelines without manually converting
        coordinates between formats.

        Inputs and Workflow

        The protocol requires a SetOfCoordinates3D object as input. These coordinates
        are typically associated with tomograms and contain the spatial positions of
        particles or regions of interest identified during subtomogram analysis.

        During execution, the protocol automatically creates an export directory and
        organizes the output files according to the selected export format. Coordinates
        are grouped by tomogram identifier, ensuring compatibility with downstream
        software packages that expect tomogram-specific coordinate files.

        The export workflow is intentionally lightweight and does not modify the
        original coordinates. Instead, it focuses on transforming the internal Scipion
        representation into the syntax and structure required by external applications.

        Supported Export Formats

        The protocol dynamically detects available export options depending on the
        installed Scipion plugins.

        TXT export produces simple text files containing X, Y, and Z coordinates.
        This format is useful for custom scripts, visualization tools, or lightweight
        processing pipelines where only particle positions are required.

        RELION STAR export generates coordinate STAR files compatible with RelionTomo.
        The protocol preserves tomogram identifiers and sampling information so the
        coordinates can be directly integrated into RELION subtomogram workflows.

        EMAN export creates JSON metadata files compatible with EMAN tomography tools.
        This option facilitates interoperability with EMAN-based subtomogram processing
        environments.

        Dynamo export generates Dynamo table files containing particle coordinates and
        placeholder alignment parameters. Although rotational and translational
        alignment values are initialized to zero in this implementation, the exported
        tables remain compatible with standard Dynamo workflows.

        SPHIRE export produces CBOX-compatible coordinate files for SPHIRE tomography
        pipelines and visualization tools.

        Coordinate Handling and Organization

        The protocol iterates over all coordinates ordered by tomogram identifier.
        For each tomogram, a separate output file is generated, simplifying data
        organization and downstream processing.

        Coordinates are exported using the bottom-left corner reference convention,
        ensuring consistency with tomography coordinate systems commonly used in
        external software packages.

        The export process preserves the original coordinate precision and tomogram
        association while adapting the formatting rules required by each target
        software package.

        Outputs and Interpretation

        After execution, the protocol generates an export directory containing the
        converted coordinate files. The exact content depends on the selected format
        and the number of tomograms present in the input dataset.

        The exported files can be directly imported into external tomography software
        for visualization, particle extraction, subtomogram averaging, or further
        refinement procedures.

        Since the protocol does not alter coordinate geometry or apply transformations,
        the exported coordinates should remain biologically and spatially consistent
        with the original Scipion project.

        Practical Recommendations

        In routine cryo-ET workflows, users should select the export format according
        to the downstream software environment. TXT export is convenient for debugging,
        scripting, or custom analyses, while STAR, Dynamo, EMAN, or CBOX exports are
        better suited for fully integrated processing pipelines.

        Before exporting, it is advisable to verify that the coordinate set is properly
        curated and associated with the correct tomograms, since downstream software
        may assume strict coordinate consistency.

        When transferring coordinates between software packages, users should also
        verify voxel size conventions, coordinate origins, and tomogram orientations
        to avoid downstream alignment inconsistencies.

        Final Perspective

        In cryo-electron tomography workflows, coordinate interoperability is essential
        for combining the strengths of different software ecosystems. This protocol
        provides a simple and reliable bridge between Scipion and external tomography
        platforms, enabling flexible and reproducible subtomogram analysis pipelines
        across multiple computational environments.

    """