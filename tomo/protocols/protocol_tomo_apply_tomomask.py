# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * National Center of Biotechnology, CSIC, Spain
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
from enum import Enum
from os import remove
from typing import Union
import mrcfile
import numpy as np
from scipy.ndimage import gaussian_filter, binary_dilation
from pwem.emlib.image.image_readers import MRCImageReader, ImageReadersRegistry, ImageStack
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String, Set
from pyworkflow.protocol import PointerParam, BooleanParam, IntParam, STEPS_PARALLEL, GE
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTomograms, SetOfTomoMasks, Tomogram

logger = logging.getLogger(__name__)


class ApplyTomoMaskOutputs(Enum):
    maskedTomograms = SetOfTomograms


class ApplyTomoMaskFormParams(Enum):
    IN_TOMO_SET = 'inTomoSet'
    IN_MASK_SET = 'inMaskSet'
    INVERT_MASK = 'invertMask'
    DILATION_PX = 'dilationPixels'
    SIGMA_GAUSSIAN = 'sigmaGaussian'


class ProtTomoApplyTomoMask(EMProtocol):
    """This protocol applies a set of masks to a given set of tomograms. The protocol
    will try to match the tomograms and the masks by tsId. Once the mask/s are applied.
    Some operations can be applied to the mask: invert, dilate and apply a gaussian filter."""


    """
    Applies one or multiple tomographic masks to a set of tomograms in order
    to isolate regions of interest, suppress unwanted signal, and prepare the
    data for downstream cryo-electron tomography analysis.

    AI Generated:

    Apply Tomomasks to Tomograms (ProtTomoApplyTomoMask) — User Manual

        Overview

        The Apply Tomomasks to Tomograms protocol applies binary or continuous
        masks over tomographic volumes to selectively preserve specific
        structural regions while suppressing background or irrelevant density.
        In cryo-electron tomography workflows, masking is a fundamental step
        for improving visualization, reducing noise, focusing subsequent
        processing steps, and restricting analysis to biologically meaningful
        regions.

        From a biological perspective, masking becomes especially important
        when working with crowded cellular environments, membrane-associated
        complexes, organelles, or heterogeneous intracellular regions where
        unwanted surrounding density may interfere with interpretation or
        downstream computational procedures. By applying masks, the protocol
        helps users emphasize the structures of interest while minimizing the
        contribution of unrelated signal.

        Inputs and General Workflow

        The protocol requires two principal inputs: a set of tomograms and a
        corresponding set of tomographic masks. The matching between tomograms
        and masks is automatically performed through their tsId identifiers.
        Only tomograms and masks sharing the same tsId are processed together.

        During execution, the protocol first identifies the matching pairs
        between both datasets. Tomograms or masks without a corresponding
        partner are excluded from processing, and the protocol reports these
        mismatches to the user. This behavior is especially useful in large
        cryo-ET projects where tomograms and masks may originate from
        different preprocessing pipelines or annotation stages.

        Once matching pairs are identified, the protocol optionally processes
        the masks before applying them to the tomograms. The processed mask is
        then multiplied voxel-by-voxel with the corresponding tomogram,
        generating a new masked tomographic volume.

        Mask Processing and Biological Relevance

        One of the central features of this protocol is the ability to modify
        the masks before application. Several optional operations are
        available, including mask inversion, binary dilation, and Gaussian
        smoothing.

        Mask inversion changes the interpretation of the mask by exchanging
        foreground and background regions. Biologically, this may be useful
        when the original segmentation defines regions to exclude rather than
        regions to preserve. Instead of retaining the segmented structure, the
        inverted mask preserves the surrounding environment.

        Binary dilation expands the masked region by a user-defined number of
        pixels in all spatial directions. This operation is particularly
        important in biological datasets where segmentations may be slightly
        conservative or where neighboring density surrounding the segmented
        structure should also be preserved. For example, membrane proteins,
        ribosome-associated regions, or flexible peripheral domains may
        benefit from moderate dilation to avoid cutting biologically relevant
        signal.

        Gaussian smoothing introduces soft transitions between masked and
        unmasked regions. Instead of abrupt binary boundaries, the protocol
        generates progressively attenuated borders that reduce edge artifacts
        during later processing stages. In cryo-ET workflows, this is
        especially valuable before Fourier-based analyses, subtomogram
        extraction, denoising, or averaging procedures where sharp boundaries
        may introduce undesired mathematical artifacts.

        The sigma parameter controls the degree of smoothing. Small sigma
        values preserve sharper transitions and structural detail, while
        larger values produce softer and blurrier boundaries. Biologically,
        choosing the appropriate smoothing level depends on the balance
        between preserving fine structural features and minimizing boundary
        artifacts.

        Validation and Dataset Consistency

        The protocol performs several validation steps to ensure that masks
        and tomograms are compatible before processing.

        First, it verifies that the sampling rates of the tomograms and masks
        are consistent within a predefined tolerance. Matching voxel size is
        biologically essential because any discrepancy would cause the mask to
        correspond to an incorrect physical region of the tomogram.

        The protocol also validates the dimensions of each tomogram-mask pair.
        If dimensions differ, the corresponding dataset is excluded from
        processing. This protects against applying masks to incompatible
        volumes, which could otherwise generate corrupted or biologically
        meaningless results.

        Importantly, these validations are performed individually for each
        tomogram-mask pair rather than globally at the dataset level. This
        design allows the protocol to operate safely on heterogeneous datasets
        where some entries may be valid while others are not.

        Parallel Execution and Processing Strategy

        The protocol is designed to operate in parallel mode, processing each
        tomogram independently. For every tsId, the workflow is divided into
        three main stages: mask preprocessing, mask application, and output
        registration.

        This strategy improves scalability and efficiency when working with
        large tomographic datasets, which are common in cellular cryo-ET
        projects. Since tomograms can be extremely large in size, independent
        processing reduces bottlenecks and allows more efficient resource
        utilization.

        Temporary smoothed masks are automatically removed after use in order
        to minimize storage consumption during execution. This behavior is
        particularly important for high-throughput cryo-ET pipelines where
        intermediate volumes may occupy substantial disk space.

        Outputs and Interpretation

        The protocol generates a new set of masked tomograms preserving the
        metadata and acquisition information of the original input set. Each
        output tomogram corresponds to the original volume after application
        of the associated processed mask.

        From a biological standpoint, the resulting tomograms contain only
        the regions selected by the mask and are therefore more focused for
        downstream interpretation, segmentation refinement, particle picking,
        subtomogram averaging, or visualization.

        The protocol also reports problematic datasets, including non-matching
        tsIds, incompatible sampling rates, or dimension mismatches. These
        warnings provide transparency and help users diagnose inconsistencies
        in complex cryo-ET projects.

        Practical Recommendations

        In most biological workflows, Gaussian smoothing is strongly
        recommended because it reduces sharp boundary artifacts and produces
        more natural transitions between preserved and suppressed density.
        Moderate sigma values are usually sufficient for routine tomographic
        analyses.

        Dilation should be used carefully. Small expansions are often useful
        for preserving neighboring structural context, but excessive dilation
        may reintroduce unwanted background signal or obscure the intended
        masking effect.

        Users should also verify that masks were generated from tomograms with
        matching voxel size and dimensions before running the protocol. Even
        small inconsistencies in sampling rate can lead to biologically
        inaccurate masking.

        For exploratory visualization tasks, relatively soft masks often
        provide visually pleasing results. In contrast, quantitative analyses
        or segmentation-driven workflows may require tighter and more precise
        masking strategies.

        Final Perspective

        In cryo-electron tomography, masking is not simply a cosmetic
        operation but a biologically meaningful preprocessing step that can
        significantly influence downstream interpretation and computational
        analysis. Properly designed masks help isolate relevant structural
        information, reduce background complexity, and improve the robustness
        of subsequent processing stages.

        The effectiveness of this protocol therefore depends not only on the
        technical correctness of the masks, but also on the biological
        understanding of which regions should be preserved, excluded, or
        softly attenuated in order to best represent the underlying specimen.
    """