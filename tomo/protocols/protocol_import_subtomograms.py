# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *              Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es) [2]
# *
# * [1] EyeSeeTea Ltd, London, UK
# * [2] BCU, Centro Nacional de Biotecnologia, CSIC
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

from os.path import basename

from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.utils.path import createAbsLink

from .protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from ..objects import SubTomogram
from ..utils import _getUniqueFileName


class ProtImportSubTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """Protocol to import a set of tomograms to the project"""

    class ProtImportSubTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
        """
        ProtImportSubTomograms — Import Subtomograms Protocol

        Overview
        --------
        Imports a set of subtomograms from files into the project, associating them
        with acquisition metadata and optional coordinates. This ensures a structured
        and consistent dataset for downstream subtomogram averaging or analysis.

        Inputs and Workflow
        -------------------
        - File Pattern: Defines the filenames of subtomogram volumes to import.
          Practical tips:
            * Ensure the pattern matches the intended set of subtomogram files.
            * Supports .mrc files and single/multi-slice volumes.
        - Sampling Rate: Sets the spatial calibration for the subtomograms.
        - Acquisition Metadata: Optional input describing imaging conditions and parameters.
        - Optional 3D Coordinates: Can be linked to the imported subtomograms (commented out in current version).

        Workflow Steps
        --------------
        1. Iterate through files matching the given pattern.
        2. Determine the dimensions and origin for each subtomogram.
        3. Create a SubTomogram object and assign sampling rate, acquisition data, and file references.
        4. Handle single-slice and multi-slice volumes correctly.
        5. Generate unique filenames and establish links in project structure.
        6. Collect all subtomograms into a SetOfSubTomograms for output.

        Outputs
        -------
        - SetOfSubTomograms: All imported subtomograms with associated sampling rate and acquisition parameters.
        - Maintains provenance with source files and optional acquisition metadata.

        Practical Recommendations
        -------------------------
        - Ensure all input files exist and match the filename pattern.
        - Verify sampling rate is consistent with imaging parameters.
        - Include acquisition metadata for accurate downstream analysis.
        - Optional linking to 3D coordinates should match the number of subtomograms.

        Web-Oriented Variant
        -------------------
        - Not explicitly provided; can be adapted by restricting input pattern and simplifying metadata handling.

        Biological Perspective
        ---------------------
        - Enables structured import of subtomogram datasets for reproducible analysis.
        - Preserves spatial calibration and acquisition context for accurate averaging.
        - Supports efficient downstream structural biology workflows, including subtomogram averaging and classification.
        """