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

import os

import pyworkflow as pw
import pyworkflow.protocol.params as params
from pwem.convert.headers import getFileFormat, MRC
from pyworkflow.object import Set, Integer
from pyworkflow.protocol import STATUS_NEW
from pyworkflow.utils.properties import Message
from pwem import emlib

from .protocol_ts_base import ProtTsProcess
from ..objects import SetOfCTFTomoSeries, CTFTomoSeries


class ProtTsEstimateCTF(ProtTsProcess):
    """
    Base class for estimating the CTF on TiltSeries
    """

    """
    Estimates the Contrast Transfer Function (CTF) for each tilt-image
    within a tilt-series dataset. The protocol provides a general
    framework for tomography CTF estimation workflows and is intended
    to be specialized by subclasses implementing specific estimation
    algorithms.

    AI Generated:

    Estimate CTF for Tilt-Series (ProtTsEstimateCTF) — User Manual
        Overview

        The ProtTsEstimateCTF protocol is designed to estimate the
        Contrast Transfer Function (CTF) of tilt-images acquired during
        cryo-electron tomography experiments. In cryo-ET workflows,
        accurate CTF estimation is essential because each tilt-image
        contributes differently to the final tomographic reconstruction,
        and optical distortions introduced by the microscope strongly
        influence the interpretability of structural information.

        The protocol provides a generic and extensible framework for
        tilt-series CTF estimation. Rather than implementing a single
        estimation strategy directly, it defines the common workflow,
        data organization, streaming logic, and metadata management
        required by tomography-oriented CTF estimation procedures.
        Specific estimation algorithms are expected to be implemented
        in subclasses through the `_estimateCtf` method.

        From a biological and structural perspective, CTF estimation
        becomes particularly important when aiming for high-resolution
        subtomogram averaging, accurate tomogram reconstruction, or
        quantitative structural interpretation. Poor CTF estimation
        propagates errors into downstream alignment and reconstruction
        stages, reducing contrast and limiting achievable resolution.

        Inputs and General Workflow

        The protocol accepts as input either a set of tilt-series or
        a previously generated set of tomography CTF series. This
        flexibility allows the protocol to operate both as an initial
        estimation stage and as part of iterative refinement workflows.

        During execution, each tilt-image is processed independently.
        A temporary working directory is created for every image,
        ensuring isolation between estimation tasks and enabling
        efficient streaming execution. Input images are converted into
        a standardized MRC-compatible representation before estimation,
        allowing downstream tools to operate consistently regardless
        of the original image format.

        Once converted, the protocol invokes the estimation method
        implemented in subclasses. This modular design allows the same
        infrastructure to support different CTF estimation engines,
        including Xmipp-based approaches, CTFFIND integrations, or
        custom tomography-specific methods.

        After estimation, the resulting CTF model is attached to the
        corresponding tilt-image and later assembled into a
        CTFTomoSeries object representing the entire tilt-series.

        Image Preparation and Downsampling

        Before estimation, tilt-images may optionally be downsampled
        using Fourier scaling according to the selected downsampling
        factor. Downsampling is frequently used in cryo-ET workflows
        because high-resolution tilt-images can be computationally
        expensive to process, especially in large datasets.

        From a practical perspective, moderate downsampling often
        accelerates estimation substantially while preserving the
        frequency information necessary for reliable defocus
        determination. However, excessive downsampling may reduce
        sensitivity to high-resolution Thon rings and negatively
        impact estimation accuracy.

        The protocol automatically updates the effective sampling rate
        after downsampling so that all downstream calculations remain
        physically consistent.

        CTF Parameters and Microscope Metadata

        The protocol automatically extracts microscope acquisition
        parameters from the input tilt-series. These include voltage,
        spherical aberration, magnification, amplitude contrast,
        scanned pixel size, and sampling rate.

        Additional estimation parameters such as window size,
        low-resolution cutoff, high-resolution cutoff, minimum
        defocus, and maximum defocus are grouped into a global
        parameter dictionary shared across the estimation workflow.

        Biologically and physically, these parameters strongly
        influence the robustness of CTF estimation. Restrictive
        defocus ranges may stabilize estimation for well-behaved
        datasets, whereas broader ranges may be required when
        acquisition conditions vary significantly across tilts.

        Similarly, the selected frequency range determines which
        portions of the power spectrum contribute to fitting.
        Excluding extremely low frequencies helps reduce background
        dominance, while excluding noisy high frequencies improves
        robustness in low-dose tilt-images.

        Streaming and Parallel Processing

        The protocol is designed to support streaming tomography
        workflows. Each tilt-image can be processed independently,
        enabling incremental execution while acquisition or previous
        preprocessing steps are still ongoing.

        Once all images belonging to a tilt-series have been processed,
        the protocol marks the corresponding series as finished and
        updates the output objects dynamically. This design is
        particularly useful in facility-scale environments or
        automated tomography pipelines where large numbers of tilt
        series must be processed continuously.

        Temporary working directories are automatically cleaned after
        processing unless debugging mode is enabled. This minimizes
        storage usage while preserving the ability to inspect
        intermediate files during development or troubleshooting.

        Outputs and Their Interpretation

        The main output of the protocol is a
        `SetOfCTFTomoSeries`, where each entry corresponds to a
        complete tilt-series and contains the estimated CTF models
        for all associated tilt-images.

        Each CTF object preserves the relationship with its original
        tilt-image, including acquisition metadata and ordering within
        the tilt-series. The protocol also maintains streaming-aware
        output states so datasets can be consumed progressively by
        downstream reconstruction or correction protocols.

        From a biological perspective, the estimated CTF parameters
        define the optical transfer characteristics of the microscope
        for every tilt-image. Accurate defocus estimation directly
        impacts the quality of phase correction, tomogram
        reconstruction, and subtomogram averaging.

        Since tilt-images are acquired at different geometric
        orientations and often under varying effective thicknesses,
        defocus values may fluctuate substantially throughout the
        tilt-series. Interpreting these variations carefully is
        important when evaluating data quality.

        Extensibility and Subclass Design

        ProtTsEstimateCTF acts primarily as an abstract framework.
        Several critical methods are intentionally left undefined,
        including `_estimateCtf`, `getCtf`, and
        `_defineProcessParams`. These methods must be implemented by
        subclasses to provide concrete estimation behavior.

        This design allows developers to integrate multiple CTF
        estimation engines while reusing the same tomography-aware
        infrastructure for data handling, streaming management,
        temporary file generation, and output assembly.

        The separation between generic workflow logic and algorithmic
        implementation also improves maintainability and simplifies
        integration of future tomography-specific estimation methods.

        Practical Recommendations

        In routine cryo-ET workflows, it is generally advisable to
        begin with conservative downsampling and moderate frequency
        ranges to ensure robust estimation across all tilt-images.
        Very high tilt angles often exhibit lower signal-to-noise
        ratios, making aggressive high-resolution fitting unstable.

        Careful inspection of estimated defocus trends across the
        tilt-series is recommended. Abrupt inconsistencies may indicate
        poor image quality, inaccurate microscope metadata, or
        unsuitable fitting parameters.

        When processing large datasets, streaming execution combined
        with automatic cleanup provides efficient resource management.
        During development or optimization of new estimation methods,
        enabling debugging mode may help preserve intermediate
        diagnostic files for validation.

        Final Perspective

        In cryo-electron tomography, CTF estimation is not merely a
        preprocessing operation but a fundamental step that determines
        the reliability of downstream structural interpretation.
        Accurate estimation improves contrast restoration, enhances
        tomogram quality, and increases the achievable resolution in
        subtomogram averaging workflows.

        By providing a modular and streaming-oriented framework,
        ProtTsEstimateCTF enables robust integration of tomography CTF
        estimation methods within scalable Scipion processing
        pipelines while remaining flexible enough to accommodate
        future methodological developments.
    """