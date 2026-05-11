# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
import os
import re
from glob import glob
from datetime import datetime
from collections import OrderedDict
from os.path import join, basename, exists
from statistics import mean
import numpy as np
from sqlite3 import OperationalError
import pyworkflow as pw
import pyworkflow.protocol
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
import tomo.objects
from pwem.objects import Transform
from pyworkflow.object import Integer, String
from pyworkflow.utils import removeBaseExt
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtImport
from tomo.convert import getAnglesFromHeader, getAnglesFromMdoc, getAnglesAndDosesFromTlt
from tomo.convert.mdoc import normalizeTSId, MDoc
from tomo.objects import TomoAcquisition, SetOfTiltSeries, SetOfTiltSeriesM
from .protocol_base import ProtTomoBase, ProtTomoImportFiles

logger = logging.getLogger(__name__)


class ProtImportTsBase(ProtTomoImportFiles):
    """ Base class for Tilt-Series and Tilt-SeriesMovies import protocols.
    """
    """
    Base protocol for importing tilt-series and tilt-series movies into a tomography workflow.
    The protocol is designed to organize raw tomography acquisitions, extract acquisition
    metadata, generate structured tilt-series objects, and prepare datasets for downstream
    cryo-electron tomography processing.

    AI Generated:

    Import Tilt-Series (ProtImportTsBase, ProtImportTs, ProtImportTsMovies) — User Manual
        Overview

        The Import Tilt-Series protocols provide the entry point for tomography datasets
        within Scipion-based cryo-EM workflows. Their main purpose is to convert raw tilt
        acquisition files into structured tilt-series objects that can be consistently used
        by alignment, reconstruction, CTF estimation, denoising, subtomogram averaging,
        and downstream biological analysis protocols.

        From a biological and experimental perspective, the import step is much more than
        simple file loading. It defines the acquisition geometry, angular organization,
        accumulated electron dose, tilt ordering, microscope metadata, and image sampling
        information that will propagate throughout the entire tomography workflow. Incorrect
        import settings can therefore compromise all subsequent processing stages.

        The protocol supports both tilt-series stacks and tilt-series movies, allowing users
        to import either already-assembled tilt images or raw movie frames acquired during
        tomography experiments. This flexibility makes the protocol suitable for both
        conventional tomography pipelines and modern motion-corrected cryo-ET workflows.

        General Import Workflow

        The protocol begins by identifying input files from a user-defined directory and
        filename pattern. Depending on the dataset organization, metadata can be extracted
        directly from filenames, image headers, mdoc files, tlt files, or manually generated
        angular ranges.

        In practical cryo-ET experiments, different acquisition systems produce datasets in
        different formats. The protocol is therefore designed to accommodate heterogeneous
        acquisition conventions while still generating standardized tilt-series objects.

        Once matching files are detected, the protocol groups them into individual tilt-series,
        assigns acquisition order and tilt angles, calculates acquisition metadata, and creates
        internal Scipion objects representing each tilt-series and its associated tilt images.

        The imported data are then copied or linked into the project structure depending on
        the selected import strategy. This ensures reproducibility while also allowing users
        to optimize disk usage and data management policies.

        Metadata Sources and Experimental Context

        One of the most important aspects of tomography import is defining how angular and
        acquisition metadata are obtained. The protocol supports several complementary strategies.

        When importing using mdoc files, the protocol reads acquisition metadata generated by
        SerialEM or related tomography acquisition software. In this mode, tilt angles,
        acquisition order, accumulated dose, microscope parameters, and movie associations
        are extracted automatically. This is generally the preferred option because it preserves
        the original acquisition information recorded during data collection.

        In many practical situations, however, mdoc information may be incomplete, unavailable,
        or partially unreliable. The protocol therefore allows users to override critical
        acquisition parameters such as voltage, magnification, sampling rate, dose per frame,
        or tilt-axis angle directly from the import form. This is particularly important in
        facility environments where metadata calibration may differ between acquisition and
        processing systems.

        Alternatively, the protocol can infer metadata directly from filename patterns.
        Special tags such as {TS}, {TO}, and {TA} allow the software to reconstruct tilt-series
        organization from file naming conventions. This strategy is especially useful for
        datasets exported from external software packages, manually reorganized experiments,
        or older tomography acquisitions lacking complete metadata support.

        Tilt Angles and Acquisition Geometry

        Correct tilt-angle interpretation is biologically critical because the angular geometry
        directly determines tomographic reconstruction quality and spatial consistency.

        The protocol supports several angle sources. Angles may be read from filenames, image
        headers, mdoc files, tlt/rawtlt files, or generated from a manually specified angular
        range. This flexibility allows the protocol to adapt to many acquisition and processing
        environments commonly encountered in cryo-ET laboratories.

        When importing from manually defined ranges, the protocol generates the expected angular
        sequence from minimum angle, maximum angle, and angular step values. This mode is useful
        for highly standardized acquisitions but requires caution because inconsistencies between
        the number of images and the number of expected angles may lead to unstable downstream
        behavior.

        The protocol also validates angular consistency across imported datasets. This validation
        step is particularly important in heterogeneous collections where different tilt-series
        may contain different numbers of images or incomplete angular coverage.

        Dose Handling and Radiation Damage Considerations

        Electron dose management is a central aspect of cryo-electron tomography because radiation
        damage progressively degrades structural information during acquisition.

        The protocol tracks both incoming dose per tilt image and accumulated dose across the
        entire tilt-series. When mdoc information is available, dose values are extracted directly
        from acquisition metadata. Otherwise, dose accumulation can be estimated from acquisition
        order and user-provided dose-per-frame values.

        From a biological perspective, accurate dose tracking is essential for subsequent
        dose-weighting, motion correction, CTF estimation, and reconstruction quality assessment.
        Incorrect dose information can significantly affect the interpretation of flexible or
        radiation-sensitive biological assemblies.

        Tilt-Series Movies and Motion Correction Workflows

        The ProtImportTsMovies subclass extends the base protocol to support raw tilt-series
        movies instead of already-integrated tilt images. In modern cryo-ET workflows, this
        approach is increasingly important because motion correction is commonly performed
        before tomographic reconstruction.

        In this mode, each tilt image may contain multiple movie frames, and the protocol
        preserves frame-range information required for downstream motion correction protocols.
        Additional acquisition parameters such as gain references and dark references can also
        be associated with the imported movies.

        Biologically, preserving movie-level information improves high-resolution reconstruction
        quality because beam-induced motion correction can substantially enhance signal recovery,
        especially in low-dose cryo-ET experiments.

        File Management and Data Organization

        The protocol provides several strategies for handling imported files. Data may be copied
        directly into the project, linked using absolute symbolic links, or linked using relative
        symbolic links.

        Copying files improves project portability and reproducibility but duplicates raw data
        storage. Symbolic linking minimizes disk usage and is often preferred in large facility
        environments where tomography datasets can occupy many terabytes.

        Relative symbolic links are especially useful when projects may later be moved between
        storage systems or computing infrastructures.

        Validation and Error Handling

        A major strength of the protocol is its extensive validation framework. Before completing
        the import process, the protocol verifies file existence, metadata consistency, angular
        integrity, and acquisition compatibility.

        The protocol detects mismatches between tilt angles and image counts, invalid filename
        patterns, missing metadata files, incorrect gain references, and malformed tilt-series
        identifiers. Invalid datasets are skipped while preserving detailed logging information
        for troubleshooting.

        This validation behavior is particularly important in large cryo-ET facilities where
        datasets often originate from multiple microscopes, operators, or acquisition sessions.

        Acquisition Metadata and Microscope Parameters

        The protocol stores key microscope acquisition parameters including accelerating voltage,
        spherical aberration, amplitude contrast, magnification, sampling rate, and tilt-axis
        orientation.

        These parameters are fundamental for downstream processing accuracy. For example, sampling
        rate consistency is critical for tomographic reconstruction scaling, while tilt-axis
        orientation directly affects alignment and reconstruction geometry.

        The protocol also handles Tomography 5 mdoc conventions by converting tilt-axis angles
        into SerialEM-compatible conventions when required.

        Outputs and Workflow Integration

        After execution, the protocol generates a SetOfTiltSeries or SetOfTiltSeriesM object
        containing all successfully imported datasets. Each tilt-series contains its associated
        tilt images or movies together with complete acquisition metadata.

        The imported datasets become immediately compatible with downstream tomography workflows
        including motion correction, tilt-series alignment, CTF estimation, tomogram reconstruction,
        denoising, segmentation, and subtomogram averaging.

        Importantly, the protocol preserves acquisition relationships and metadata consistency
        throughout the workflow, ensuring reproducibility and biological interpretability.

        Practical Recommendations

        In routine cryo-ET workflows, importing from mdoc metadata is generally recommended
        whenever reliable acquisition files are available. This minimizes manual intervention
        and preserves original experimental information.

        When using filename-based imports, users should carefully verify filename conventions
        and ensure that tilt angles and acquisition order are correctly encoded. Small mistakes
        at this stage can propagate into major reconstruction artifacts later in the workflow.

        For movie-based datasets, providing gain references during import is strongly recommended
        because downstream motion correction quality depends heavily on accurate detector calibration.

        Before large-scale processing, it is advisable to visually inspect imported tilt-series,
        verify angular coverage, confirm dose consistency, and ensure that acquisition geometry
        matches the original experimental setup.

        Final Perspective

        In cryo-electron tomography, data import is not merely a technical preprocessing step
        but the foundation upon which the entire reconstruction workflow is built. Accurate
        interpretation of acquisition geometry, dose accumulation, and microscope metadata is
        essential for obtaining biologically meaningful tomograms.

        The Import Tilt-Series protocols provide a flexible and robust framework capable of
        handling the diversity of modern tomography datasets while maintaining metadata integrity,
        workflow reproducibility, and compatibility with advanced cryo-ET processing pipelines.
    """