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
from typing import Union

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
from pyworkflow.protocol import PointerParam, STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTiltSeries, TiltSeries, SetOfCTFTomoSeries, CTFTomoSeries

logger = logging.getLogger(__name__)
IN_TS_SET = 'inTsSet'
IN_CTF_SET = 'inCtfSet'


class InvertTiltsOutputs(Enum):
    tiltSeries = SetOfTiltSeries
    ctfs = SetOfCTFTomoSeries


class ProtInvertTiltAngles(EMProtocol):
    """This protocol inverts the physical handedness of the introduced tilt-series by inverting the tilt angles
    in the metadata associated to each tilt-series. Introducing the CTFs will update the pointer from introduced
    CTFs to the tilt-series with the inverted tilt-angles in order to keep the coherence in the
    relationship between both objects after the angle inversion operation."""

    """
    Inverts the tilt angles associated with a set of tilt-series in order
    to change the physical handedness of the tomographic acquisition.
    The protocol generates a new set of tilt-series where each tilt angle
    is multiplied by -1 while preserving the original acquisition structure
    and metadata relationships.

    AI Generated:

    Invert Tilt Angles (ProtInvertTiltAngles) — User Manual
        Overview

        The Invert Tilt Angles protocol modifies the angular metadata of
        a tilt-series dataset by reversing the sign of every tilt angle.
        This operation effectively changes the handedness convention of
        the tomographic acquisition while preserving the original image
        data and acquisition order. The protocol is particularly useful
        in cryo-electron tomography workflows where tilt geometry needs
        to be corrected to maintain consistency between reconstruction,
        visualization, and downstream analysis pipelines.

        In practical biological workflows, handedness inconsistencies may
        appear when datasets are processed with different acquisition
        conventions or software packages. Although the voxel intensities
        remain unchanged, incorrect tilt-angle orientation can lead to
        mirrored reconstructions, incorrect structural interpretation,
        or incompatibilities with subtomogram averaging and segmentation
        procedures. This protocol addresses those issues by generating a
        coherent tilt-series representation with inverted angular metadata.

        Inputs and General Workflow

        The protocol requires as input a set of tilt-series containing
        the angular metadata associated with each tilt image. Optionally,
        a corresponding set of CTF estimations can also be provided.
        When CTF information is included, the protocol automatically
        updates the relationship between the new tilt-series and their
        associated CTF objects to preserve internal consistency across
        the dataset.

        During execution, the protocol iterates through each tilt-series
        independently. A new output tilt-series is created by cloning the
        original metadata and duplicating every tilt image while replacing
        each tilt angle with its negative counterpart. The image ordering,
        acquisition information, and structural organization remain
        unchanged throughout the process.

        Handling of CTF Associations

        One of the most important aspects of this protocol is the
        preservation of coherence between tilt-series and CTF estimation
        objects. In cryo-ET processing pipelines, CTF models are tightly
        associated with the angular geometry of the acquisition. If tilt
        angles are inverted without updating these relationships, later
        reconstruction or refinement stages may become inconsistent.

        When a set of CTF tomo-series is introduced, the protocol first
        validates the compatibility between both datasets using their
        tilt-series identifiers. Only matching tilt-series are processed.
        The protocol then generates a new output CTF set linked to the
        newly generated tilt-series with inverted angles, ensuring that
        all metadata dependencies remain synchronized.

        Validation and Dataset Consistency

        Before processing begins, the protocol checks whether the
        introduced tilt-series and CTF sets share compatible tilt-series
        identifiers. If no common identifiers are found, execution stops
        with an error because the relationship between both datasets
        cannot be established reliably.

        In cases where only part of the datasets match, the protocol
        continues processing the compatible entries while storing a
        warning message describing the non-matching identifiers. This
        behavior allows partially compatible datasets to be reused
        without forcing complete manual curation beforehand.

        Parallel Processing Strategy

        The protocol executes using a parallel step-based strategy in
        which each tilt-series is processed independently. This design
        improves scalability for large cryo-electron tomography projects
        containing many tilt-series acquisitions. Since the operation only
        modifies metadata and does not alter image intensities, execution
        is typically lightweight and computationally efficient.

        Outputs and Their Interpretation

        After execution, the protocol produces a new set of tilt-series
        with inverted tilt angles. The original datasets remain untouched,
        allowing users to preserve both conventions within the same
        project if needed. Each generated tilt-series maintains the same
        image sequence and acquisition structure as the original input,
        differing only in the sign of the angular metadata.

        If CTF information was provided, an additional output set of
        CTF tomo-series is generated and linked directly to the updated
        tilt-series. This ensures compatibility with subsequent
        reconstruction, alignment, or subtomogram analysis workflows.

        Practical Recommendations

        In biological practice, this protocol is most commonly used when
        importing datasets generated under different handedness
        conventions or when correcting inconsistencies discovered during
        tomographic reconstruction. Before applying the protocol, it is
        advisable to confirm that the observed handedness discrepancy is
        truly caused by tilt-angle orientation and not by visualization
        settings or reconstruction artifacts.

        When working with CTF information, users should always introduce
        the associated CTF set together with the tilt-series to preserve
        dataset coherence automatically. Reviewing the generated warning
        messages is also recommended in order to identify missing or
        incompatible tilt-series identifiers before continuing with
        downstream processing.

        Final Perspective

        Although mathematically simple, tilt-angle inversion is a
        biologically significant metadata operation because it directly
        affects the geometric interpretation of tomographic data.
        Maintaining consistency between tilt geometry, CTF estimation,
        and reconstruction conventions is essential for obtaining
        reliable structural interpretations in cryo-electron tomography
        workflows.
    """
