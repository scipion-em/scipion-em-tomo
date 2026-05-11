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
    """
    Composes tilt-series in streaming mode from motion-corrected micrographs
    and corresponding mdoc metadata files. The protocol continuously monitors
    incoming data, matches tilt images with metadata information, and generates
    complete tilt-series datasets ready for downstream cryo-electron tomography
    processing workflows.

    AI Generated:

    Compose Tilt Series (ProtComposeTS) — User Manual
        Overview

        The Compose Tilt Series protocol is designed to generate tomography
        tilt-series automatically and incrementally during data acquisition.
        It combines motion-corrected micrographs with metadata extracted from
        mdoc files, creating organized tilt-series objects that can be directly
        used in downstream cryo-ET workflows such as alignment, reconstruction,
        subtomogram averaging, or CTF estimation.

        In practical cryo-electron tomography workflows, tilt images are often
        acquired sequentially while preprocessing steps such as motion correction
        are executed simultaneously. This protocol addresses that streaming
        scenario by continuously monitoring incoming micrographs and metadata,
        composing the tilt-series as soon as sufficient information becomes
        available.

        For biological users working in automated acquisition environments,
        this protocol allows near real-time organization of tomography datasets,
        reducing manual intervention and accelerating processing pipelines.

        Inputs and Streaming Workflow

        The protocol requires two principal inputs: a SetOfMicrographs
        containing the motion-corrected tilt images, and a directory containing
        the corresponding mdoc files generated during acquisition.

        Each mdoc file describes the acquisition metadata of a tilt-series,
        including tilt angles, acquisition order, dose information, and tilt
        axis orientation. The protocol continuously scans the specified folder,
        detects new mdoc files, validates their contents, and attempts to match
        each tilt image described in the metadata with the corresponding
        micrograph available in the streaming dataset.

        The streaming behaviour is controlled through the “Time for the next
        tilt” parameter. This delay determines how long the protocol waits
        after the last modification of an mdoc file before considering that
        tilt-series complete. This mechanism is particularly important in live
        acquisition workflows where tilt images arrive progressively over time.

        In high-throughput tomography facilities or automated acquisition
        sessions such as PACEtomo experiments, carefully adjusting this timeout
        helps avoid premature processing of incomplete tilt-series.

        Mdoc File Detection and Validation

        The protocol searches for mdoc files using user-defined path and pattern
        parameters. Standard wildcard expressions are supported, allowing
        flexible integration with different acquisition folder organizations.

        Before composing a tilt-series, each mdoc file undergoes several
        validation steps. Files containing exclusion keywords may be ignored,
        which is useful when acquisition directories contain temporary,
        corrupted, or unwanted metadata files.

        The protocol also validates the integrity of the mdoc contents. If the
        metadata format is invalid or required information is missing, the
        corresponding tilt-series is skipped to prevent propagation of errors
        into downstream processing.

        A minimum number of tilts can also be enforced. This is biologically
        important because severely incomplete tilt-series generally produce
        unreliable tomographic reconstructions and may compromise later
        subtomogram analysis.

        Matching Micrographs with Metadata

        Once the mdoc information is validated, the protocol attempts to match
        each tilt image described in the metadata with the corresponding
        motion-corrected micrograph.

        The matching process is based on filename consistency between the mdoc
        entries and the imported micrographs. Only successfully matched images
        are incorporated into the final tilt-series.

        During streaming acquisition, it is common for some tilt images to
        still be missing while the protocol is running. In these situations,
        the protocol temporarily delays composition until additional
        micrographs become available.

        If the acquisition stream is already closed, the protocol evaluates
        whether the percentage of available tilts satisfies the minimum
        threshold defined by the user. This behaviour provides flexibility in
        cases where certain movies failed during acquisition or preprocessing,
        while still allowing partially complete tilt-series to be recovered.

        Tilt Ordering and Geometrical Consistency

        Before generating the final tilt-series, the protocol sorts all tilt
        images according to their tilt angle. This guarantees that the output
        tilt-series follows the correct geometrical order required by most
        tomography reconstruction algorithms.

        Correct angular ordering is biologically critical because downstream
        reconstruction software assumes that projections correspond to a
        physically meaningful tilt sequence. Incorrect ordering may lead to
        reconstruction artifacts or completely invalid tomograms.

        Tilt Axis Angle Handling

        The protocol provides several mechanisms to manage the tilt axis angle.
        By default, the value is extracted directly from the mdoc metadata.
        However, users may manually override this value when prior knowledge of
        the acquisition geometry is available.

        An additional correction mode is included for Tomography 5 and related
        acquisition systems that use alternative tilt-axis conventions. In
        these cases, the protocol automatically converts the angle definition
        into the standard geometry expected by Scipion tomography workflows.

        Accurate tilt-axis definition is particularly important because errors
        in this parameter directly affect alignment quality and tomographic
        reconstruction accuracy.

        Generation of Tilt-Series Stacks

        After successful validation and matching, the protocol composes the
        final tilt-series stack by sequentially combining the individual tilt
        images into a single MRC stack file.

        The resulting tilt-series preserves all relevant acquisition metadata,
        including tilt angles, acquisition order, accumulated dose, voltage,
        magnification, and sampling rate.

        This metadata propagation is essential for downstream tomography
        protocols, especially alignment and reconstruction methods that rely on
        accurate acquisition geometry and dose information.

        The protocol also computes global acquisition properties such as
        minimum tilt angle, maximum tilt angle, angular step size, and total
        accumulated dose.

        Odd/Even Tilt-Series Generation

        An optional feature allows the generation of odd and even tilt-series
        stacks when the input micrographs contain odd/even motion-correction
        information.

        This functionality is particularly useful for advanced cryo-ET
        workflows involving noise estimation, validation procedures, or
        resolution assessment strategies based on independent half datasets.

        If odd/even data are requested but unavailable in the metadata, the
        protocol raises a validation error to prevent generation of incomplete
        outputs.

        Acquisition Parameter Estimation

        The protocol automatically generates tomography acquisition objects by
        combining metadata extracted from both the micrographs and the mdoc
        files.

        Parameters such as voltage, magnification, spherical aberration,
        amplitude contrast, dose per frame, angular range, and tilt step are
        propagated into the resulting tilt-series.

        When certain acquisition parameters are already present in the input
        micrographs, these values are prioritized to maintain consistency
        across the processing workflow.

        Streaming Behaviour and Data Persistence

        One of the central aspects of this protocol is its continuous streaming
        execution model. The protocol repeatedly scans the acquisition folders,
        updates the input micrograph set, and dynamically generates processing
        steps for newly detected tilt-series.

        As each tilt-series is completed, it is immediately written into the
        output SetOfTiltSeries object and made available for downstream
        protocols without waiting for the entire acquisition session to finish.

        This behaviour is particularly advantageous in automated cryo-ET
        facilities where reconstruction and quality-control procedures may run
        in parallel with data collection.

        Outputs and Their Interpretation

        The primary output is a SetOfTiltSeries object containing one or more
        composed tilt-series. Each tilt-series includes the ordered tilt-image
        stack together with complete acquisition metadata.

        Individual tilt images retain their associated tilt angle, acquisition
        order, accumulated dose, and optional odd/even information. The output
        is therefore immediately compatible with standard tomography alignment
        and reconstruction protocols.

        From a biological perspective, the quality and completeness of the
        generated tilt-series strongly influence all downstream analyses.
        Missing tilts, incorrect geometry, or inconsistent metadata may reduce
        tomogram quality and limit interpretability of macromolecular
        structures.

        Practical Recommendations

        In routine cryo-ET acquisition workflows, it is recommended to ensure
        consistent naming conventions between mdoc entries and motion-corrected
        micrographs before starting the protocol.

        For streaming acquisitions, the timeout parameter should be adjusted
        according to the acquisition speed and preprocessing latency. Fast
        acquisitions may work well with short waiting times, while slower
        acquisition schemes such as PACEtomo generally require longer delays.

        The minimum percentage of required tilts should also be selected
        carefully. Relaxed thresholds may recover partially incomplete datasets,
        but excessive tilt loss can compromise tomographic reconstruction
        quality.

        When using odd/even datasets, users should verify beforehand that the
        motion-correction protocol generated the required metadata.

        Final Perspective

        For cryo-electron tomography workflows, tilt-series composition is more
        than a simple file organization step. Correct association between tilt
        images, acquisition metadata, angular ordering, and dose information is
        essential for obtaining biologically meaningful tomographic
        reconstructions.

        The Compose Tilt Series protocol provides an automated and streaming-
        oriented solution that integrates naturally into modern cryo-ET
        acquisition pipelines, enabling efficient, scalable, and reliable
        tomography data management.
    """
