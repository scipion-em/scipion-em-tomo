# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
import numpy as np
import math
import csv
import pwem.objects as data
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase

MISALIGNED_TS_NAME = "MisalignedTiltSeries"
INTERPOLATED_TS_NAME = "InterpolatedTiltSeries"


class ProtTomoMisalignTiltSeries(EMProtocol, ProtTomoBase):
    """
    Introduce misalignment in the transformation matrix of a tilt-series.
    NOTE: The Interpolated tilt series in this case resembles a not aligned tilt series
    or an aligned one in case you want to apply the inverse of the misalignment
    transformation matrix.
    """
    """
    Introduces controlled misalignment into tilt-series transformation matrices
    by modifying shifts and angular parameters. The protocol can also generate
    interpolated tilt-series by applying the simulated misalignment matrices
    directly to the image data.

    AI Generated:

    Misalign Tilt-Series (ProtTomoMisalignTiltSeries) — User Manual

        Overview

        The Misalign Tilt-Series protocol is designed to simulate realistic
        alignment errors in cryo-electron tomography tilt-series datasets.
        Instead of reconstructing or correcting tilt-series alignment, the
        protocol intentionally perturbs the transformation matrices associated
        with each tilt image. This allows users to reproduce acquisition or
        alignment inaccuracies commonly observed during tomographic processing.

        In practical cryo-ET workflows, this type of simulation is especially
        useful for benchmarking alignment algorithms, testing reconstruction
        robustness, validating correction methods, or generating synthetic
        datasets for methodological development. By introducing controlled and
        reproducible geometric distortions, users can evaluate how sensitive
        downstream reconstruction pipelines are to different classes of
        alignment errors.

        From a biological perspective, tilt-series misalignment can strongly
        affect tomogram quality, structural interpretability, and subtomogram
        averaging performance. Simulating these imperfections provides a useful
        framework for understanding the limitations and stability of processing
        pipelines under non-ideal experimental conditions.

        Inputs and General Workflow

        The protocol requires a set of input tilt-series containing either
        existing alignment transformations or raw geometric metadata. For each
        tilt image, the protocol modifies the associated transformation matrix
        according to a user-defined mathematical model.

        The workflow operates independently on every tilt-series. During
        execution, the protocol iterates through all tilt images and applies
        controlled perturbations to translational and rotational alignment
        parameters. These modifications are accumulated directly into the
        transformation matrix associated with each image.

        The protocol produces a new set of misaligned tilt-series containing
        the modified transformations. Optionally, it can also generate fully
        interpolated image stacks where the geometric transformations are
        physically applied to the image data itself.

        Shift Misalignment in X and Y

        One of the central features of the protocol is the ability to simulate
        translational misalignment independently along the X and Y axes. These
        perturbations mimic common alignment instabilities observed during
        experimental acquisition, including stage drift, beam-induced motion,
        cumulative alignment inaccuracies, or local registration failures.

        The shift perturbation model combines several independent components.
        A constant offset introduces a systematic displacement affecting all
        images equally. An incremental component propagates a progressive drift
        across the tilt-series, reproducing situations where alignment errors
        accumulate gradually during acquisition.

        In addition, the protocol includes sinusoidal perturbations that model
        cyclic or oscillatory alignment instabilities. Two complementary modes
        are implemented: a half-sine lobe and a full sine cycle. These are
        particularly useful for reproducing mechanical instabilities or
        periodic stage deformations that vary continuously across the angular
        range of the tilt-series.

        Finally, a stochastic Gaussian component can be added to simulate
        random alignment noise. This random contribution reproduces the type
        of local uncertainty commonly observed in low signal-to-noise
        experimental datasets.

        By combining these components, the protocol allows users to reproduce
        highly realistic and biologically plausible alignment distortions.

        Angular Misalignment

        In addition to translational perturbations, the protocol can also
        introduce angular errors directly into the rotational component of the
        transformation matrix. These perturbations simulate inaccuracies in
        tilt-angle estimation or rotational alignment refinement.

        The angular perturbation follows the same mathematical philosophy used
        for translational shifts. Constant rotational offsets simulate global
        calibration errors, while incremental angular drifts reproduce
        progressive orientation inaccuracies accumulated during acquisition.

        Sinusoidal angular perturbations are especially valuable when studying
        systematic rotational oscillations caused by stage instability or
        imperfect microscope mechanics. Random angular noise further allows
        simulation of local orientation uncertainty typical of experimental
        cryo-ET data.

        Biologically, angular inaccuracies are often more damaging than small
        translational shifts because they directly affect projection geometry.
        Even moderate rotational errors can significantly degrade tomogram
        resolution and compromise subtomogram averaging quality.

        Mathematical Model of Misalignment

        The protocol defines misalignment using analytical functions that vary
        across the tilt-series index. Each perturbation is computed as a
        combination of deterministic and stochastic components.

        The translational and angular increments are modeled through offsets,
        linear drift terms, sinusoidal functions, and Gaussian noise
        contributions. This design allows users to generate highly controlled
        synthetic alignment defects ranging from simple systematic offsets to
        complex non-linear perturbation patterns.

        The resulting perturbations are accumulated directly into the affine
        transformation matrix associated with each tilt image. Rotational
        modifications are applied by recomputing the rotation matrix elements,
        while translational perturbations modify the shift coordinates.

        Interpolated Tilt-Series Generation

        The protocol optionally allows the generation of interpolated
        tilt-series by physically applying the modified transformation matrices
        to the image data. This operation produces image stacks that visually
        resemble experimentally misaligned datasets.

        In practical terms, the interpolated output simulates the appearance
        of a tilt-series after geometric distortion has been introduced. This
        is especially useful for testing alignment correction algorithms,
        evaluating reconstruction robustness, or training machine learning
        approaches under realistic acquisition imperfections.

        An optional inverse-matrix mode is also available. In this case, the
        inverse of the generated misalignment matrix is stored in the output
        metadata. This feature is particularly useful for workflows focused on
        alignment recovery or correction benchmarking.

        Transformation Matrix Export

        During execution, the protocol stores the generated transformation
        parameters into external XF matrix files. Two different outputs are
        generated.

        The first file contains only the introduced perturbation increments,
        allowing users to inspect the synthetic misalignment independently of
        the original transformations. The second file stores the final
        transformation matrices after all modifications have been applied.

        These exported matrices are especially useful for debugging,
        benchmarking external reconstruction software, or quantitatively
        comparing alignment correction strategies.

        Outputs and Their Interpretation

        The protocol generates a new set of misaligned tilt-series preserving
        the original image metadata while replacing the transformation
        matrices with the perturbed versions.

        When interpolation is enabled, an additional set of interpolated
        tilt-series is produced. In this output, the geometric transformations
        have already been applied to the image data, generating physically
        distorted image stacks.

        From a methodological perspective, the misaligned output is useful for
        evaluating metadata-driven alignment procedures, while the
        interpolated output is better suited for testing complete image-based
        correction pipelines.

        Practical Recommendations

        In most simulation workflows, it is advisable to begin with small
        translational perturbations and limited angular noise in order to
        evaluate baseline reconstruction sensitivity. Excessively large
        perturbations may generate unrealistic datasets that no longer
        resemble experimentally achievable conditions.

        Incremental drift components are particularly useful for reproducing
        acquisition instabilities observed during long tilt-series collection.
        Sinusoidal perturbations are better suited for studying systematic
        stage oscillations or cyclic alignment artifacts.

        Random noise terms should be introduced carefully because high
        stochastic perturbations can rapidly destabilize reconstruction
        quality. In practice, biologically meaningful simulations usually
        combine moderate systematic drift with limited random noise.

        When generating interpolated datasets, users should visually inspect
        the resulting tilt-series to ensure that the simulated distortions
        remain physically plausible and compatible with the intended
        benchmarking scenario.

        Final Perspective

        The Misalign Tilt-Series protocol provides a flexible framework for
        simulating realistic alignment imperfections in cryo-electron
        tomography datasets. Rather than serving as a correction tool, the
        protocol focuses on controlled degradation of alignment quality in
        order to support methodological validation, robustness testing, and
        synthetic data generation.

        For cryo-ET developers and advanced users, the ability to reproduce
        systematic and stochastic alignment defects represents an important
        resource for understanding how reconstruction pipelines behave under
        imperfect experimental conditions. Careful tuning of translational and
        angular perturbations enables the creation of realistic synthetic
        datasets that closely resemble the variability encountered in real
        tomographic acquisition workflows.
    """