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
from os.path import abspath, relpath

import numpy as np

import pyworkflow as pw
import pyworkflow.protocol.params as params
from pyworkflow.object import String, CsvList, Set
from pyworkflow.utils import removeBaseExt
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler

from ..objects import TiltSeries, TiltImage, SetOfTiltSeries
from .protocol_ts_base import ProtTsProcess

OUTPUT_TILT_SERIES_ODD = 'TiltSeriesOdd'
OUTPUT_TILT_SERIES_EVEN = 'TiltSeriesEven'
EVEN = 'even'
ODD = 'odd'
OUTPUT_TILT_SERIES_DW = 'TiltSeriesDW'


class ProtTsCorrectMotion(ProtTsProcess):
    """
    Base class for movie alignment protocols such as:
    motioncorr, crosscorrelation and optical flow

    Alignment parameters are defined in common. For example,
    the frames range used for alignment and final sum, the binning factor
    or the cropping options (region of interest)
    """

    """
    Performs motion correction on tilt-series movies generated during cryo-electron tomography acquisition.
    The protocol provides a common framework for different motion-correction strategies such as MotionCor,
    cross-correlation approaches, or optical-flow–based methods. Its purpose is to compensate for beam-induced
    motion across movie frames before generating the final aligned tilt-series.

    AI Generated:

    Tilt-Series Motion Correction (ProtTsCorrectMotion) — User Manual
        Overview

        The Tilt-Series Motion Correction protocol is designed to correct beam-induced motion affecting
        tilt-series movies acquired during cryo-electron tomography experiments. During data acquisition,
        each tilt image is commonly recorded as a movie composed of multiple frames rather than as a
        single exposure. This strategy allows correction of sample drift and beam-induced specimen motion,
        both of which can significantly degrade tomographic reconstruction quality if left uncorrected.

        The protocol acts as a general framework shared by several motion-correction implementations,
        including algorithms based on cross-correlation, MotionCor-style approaches, or optical flow
        strategies. Regardless of the specific backend implementation, the biological objective remains
        the same: producing stable and high-quality aligned tilt images suitable for downstream tomographic
        reconstruction, fiducial alignment, subtomogram analysis, and structural interpretation.

        In practical cryo-ET workflows, motion correction is one of the earliest and most critical
        preprocessing steps. Errors introduced at this stage propagate through the entire reconstruction
        pipeline and may negatively affect alignment accuracy, contrast preservation, and achievable
        resolution.

        Inputs and General Workflow

        The protocol requires as input a set of tilt-series movies. Each tilt image is expected to contain
        multiple movie frames acquired sequentially during exposure. The protocol processes each tilt-image
        movie independently, applies motion correction across frames, and finally reconstructs a corrected
        tilt-series stack preserving the original tomographic acquisition geometry.

        During execution, the protocol iterates through all tilt images belonging to each tilt-series.
        Corrected images are generated individually and later merged into final MRC stacks representing
        the aligned tilt-series. The workflow preserves metadata such as tilt angle, acquisition order,
        sampling rate, and acquisition parameters, ensuring compatibility with downstream tomography tools.

        The protocol also supports optional generation of dose-weighted tilt-series as well as odd/even
        frame splitting workflows commonly used in denoising and validation strategies.

        Frame Alignment and Summation

        One of the most important concepts in this protocol is the distinction between frame alignment
        and frame summation. Alignment determines which movie frames are used to estimate specimen motion,
        while summation determines which aligned frames contribute to the final averaged tilt image.

        In many biological datasets, excluding early frames from the summation can improve image quality
        because initial frames often accumulate the highest radiation damage. At the same time, those
        frames may still contain useful signal for motion estimation. For this reason, the protocol allows
        independent control of alignment and summation frame ranges.

        Users may alternatively choose to use the same frame range for both operations, which is often
        sufficient for routine preprocessing pipelines and simplifies parameter selection.

        From a biological perspective, selecting appropriate frame ranges is particularly important for
        radiation-sensitive specimens, flexible complexes, or low-dose tomographic datasets where signal
        preservation is critical.

        Binning and Sampling Considerations

        The protocol allows application of a binning factor before motion correction. Binning reduces
        image size and computational cost while increasing signal-to-noise ratio. This is especially
        useful for large tomographic datasets or rapid exploratory processing.

        However, excessive binning reduces high-resolution information and may limit downstream structural
        interpretation. In biological workflows aimed at high-resolution subtomogram averaging, users
        typically employ minimal binning or process unbinned data whenever computationally feasible.

        The output sampling rate is automatically adjusted according to the selected binning factor so that
        downstream protocols correctly interpret voxel dimensions.

        Cropping and Region Selection

        Advanced users may define crop offsets and crop dimensions to restrict processing to a specific
        region of the movie frames. This functionality can reduce memory usage and accelerate processing,
        particularly for very large detectors.

        Biologically, cropping is most useful when only a subset of the detector contains relevant signal
        or when edge artifacts, damaged regions, or detector imperfections should be excluded from motion
        estimation.

        Care should be taken to avoid cropping out fiducials, regions of interest, or important structural
        features required for later tomographic alignment.

        Odd and Even Frame Splitting

        The protocol optionally supports generation of odd-frame and even-frame tilt-series. In this mode,
        aligned movie frames are separated into two independent image series after motion correction while
        preserving the same estimated motion trajectory.

        This functionality is particularly valuable for modern denoising approaches, neural-network
        training pipelines, and independent half-set validation strategies. Since odd and even images
        contain statistically independent noise while sharing the same underlying structure, they are
        frequently used for self-supervised learning and noise2noise methodologies in cryo-EM and cryo-ET.

        The protocol automatically creates separate odd and even tilt-series stacks while maintaining
        metadata consistency with the original acquisition.

        Dose Weighting

        The framework also supports generation of dose-weighted tilt-series. Dose weighting compensates
        for radiation damage accumulated throughout exposure by attenuating high-frequency information
        from later frames that are biologically more damaged.

        In cryo-electron tomography, dose weighting is especially important because tilt-series acquisition
        distributes electron dose across many projections, often resulting in low signal conditions.
        Proper dose weighting improves contrast preservation and can significantly benefit subtomogram
        averaging workflows.

        Depending on the subclass implementation, dose-weighted outputs may be generated automatically
        alongside standard corrected tilt-series.

        Internal Processing Architecture

        The protocol is implemented as an extensible base class rather than a standalone correction
        algorithm. The actual motion-correction strategy is delegated to subclasses through the
        _processTiltImageM method, which must be implemented by each specific algorithm.

        This architecture allows different correction engines to share the same data-management workflow,
        parameter definitions, metadata handling, and output generation logic while implementing different
        mathematical approaches for motion estimation.

        The framework also manages temporary working folders, gain and dark reference conversion,
        output stack generation, synchronization of streaming outputs, and tilt-series assembly.

        Outputs and Their Interpretation

        After execution, the protocol generates one or more corrected tilt-series stacks composed of
        aligned tilt images. Each output tilt image preserves its acquisition metadata and geometric
        information, ensuring compatibility with reconstruction and alignment software.

        When enabled, additional outputs may include:
        a dose-weighted tilt-series,
        an odd-frame tilt-series,
        and an even-frame tilt-series.

        The corrected tilt-series generally exhibit improved contrast stability, reduced blurring,
        and better fiducial consistency compared to raw movie averages. These improvements directly
        influence tomographic reconstruction quality and downstream subtomogram analysis.

        Practical Recommendations

        For most routine cryo-ET workflows, using default alignment ranges and moderate binning provides
        a good balance between computational efficiency and image quality. High-resolution studies,
        however, may require reduced binning and careful optimization of frame-selection parameters.

        Odd/even splitting should primarily be enabled when preparing datasets for denoising workflows
        or validation strategies, since it increases storage and processing requirements.

        Dose weighting is generally recommended for biological datasets where radiation damage affects
        high-resolution information, particularly in subtomogram averaging applications.

        Users should visually inspect corrected tilt images and reconstructed stacks whenever possible,
        especially for challenging datasets with low contrast, severe motion, or acquisition artifacts.

        Final Perspective

        Motion correction in cryo-electron tomography is not simply a technical preprocessing operation
        but a biologically critical step that strongly affects the interpretability and quality of the
        final tomographic reconstruction. Proper handling of frame alignment, dose accumulation,
        sampling preservation, and output consistency is essential for obtaining reliable structural
        information from tomographic datasets.

        By providing a unified and extensible framework for multiple correction strategies, the protocol
        enables robust integration of motion correction into both routine and advanced cryo-ET pipelines.
    """