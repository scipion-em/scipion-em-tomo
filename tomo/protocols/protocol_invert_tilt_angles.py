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
    """
    Invert Tilt Angles (ProtInvertTiltAngles) — User Manual

    Overview

    The Invert Tilt Angles protocol reverses the physical handedness of a tilt-series
    dataset by multiplying all tilt angles by -1.

    Its main purpose is to generate a new tilt-series set in which the angular metadata
    are inverted while preserving the original image data and tilt ordering.

    For a biological user, this protocol is useful when a dataset has been imported
    with an incorrect tilt-angle convention or when handedness needs to be corrected
    before tomographic reconstruction or downstream subtomogram analysis.

    In cryo-electron tomography, the sign convention of tilt angles directly affects
    the geometric interpretation of the acquisition. An incorrect handedness may lead
    to mirrored reconstructions or biologically misleading spatial interpretations.

    Inputs and General Workflow

    The protocol requires one principal input:

    - A `SetOfTiltSeries`

    Optionally, it can also receive:

    - A `SetOfCTFTomoSeries`

    During execution, the protocol processes each tilt series independently.

    For every tilt series:

    - a new output tilt series is created,
    - each tilt image is copied,
    - the tilt angle stored in the metadata is multiplied by -1.

    The image files themselves are not modified.

    This means the protocol changes only metadata interpretation and not the actual
    experimental image content.

    Handedness Correction

    The biological relevance of this protocol lies in handedness correction.

    By inverting the tilt angles, the protocol changes the geometric interpretation
    of the tilt acquisition.

    This operation may be required when:

    - data have been imported using the wrong tilt-angle sign convention,
    - external acquisition software uses a different angular convention,
    - reconstruction results suggest a mirrored handedness.

    In tomography workflows, correcting handedness at the metadata level is often
    preferable to manipulating reconstructed volumes afterward.

    Optional CTF Handling

    If a set of tomography CTF estimations is provided, the protocol also generates
    a new CTF output set.

    This is biologically important because CTF metadata remain meaningful only if they
    stay linked to the correct tilt-series geometry.

    During execution:

    - only tilt-series identifiers present in both input sets are considered;
    - matching CTF series are copied to the output;
    - the CTF output is linked to the newly generated tilt-series set.

    This preserves coherence between angular metadata and optical metadata after
    angle inversion.

    Matching Between Tilt Series and CTFs

    When CTFs are provided, the protocol compares the tilt-series identifiers from
    both inputs.

    Three cases are possible.

    Complete match:
    - all tilt series have corresponding CTFs.

    Partial match:
    - only common tilt-series identifiers are processed.

    No match:
    - execution stops with an error.

    If some identifiers are present only in one input set, the protocol reports
    these non-matching tilt-series identifiers.

    This is particularly useful in large tomography datasets where tilt-series and
    CTF estimations may have been generated at different stages of preprocessing.

    Tilt-Series Output Generation

    For every processed tilt series, the protocol creates a new `TiltSeries` object.

    The output tilt series preserves:

    - the original identity,
    - the original metadata,
    - the original image references.

    The only modified parameter is the tilt angle.

    Each output tilt image therefore remains physically linked to the same image data
    but interpreted with inverted angular geometry.

    CTF Output Generation

    When CTFs are provided, the protocol also creates a new
    `SetOfCTFTomoSeries`.

    Each imported CTF series is copied into the output set and linked to the new
    tilt-series output.

    Importantly, the CTF estimations themselves are not recalculated.

    Only the relational consistency between CTFs and tilt-series metadata is updated.

    From a biological perspective, this ensures that downstream CTF-aware processing
    remains coherent after handedness inversion.

    Parallel Execution

    The protocol runs with parallel step execution.

    This means that each tilt series is processed independently.

    In practical terms, this makes the protocol efficient for experiments involving
    many tomograms or large tilt-series collections.

    Since no image recalculation is performed, execution is usually lightweight and
    primarily limited by metadata writing.

    Output Interpretation

    After completion, the protocol produces:

    - a new `SetOfTiltSeries` with inverted tilt angles;
    - optionally, a new `SetOfCTFTomoSeries` linked to the inverted tilt series.

    If partial mismatches between input tilt series and CTFs were detected, the
    summary reports the non-matching tilt-series identifiers.

    This helps users quickly identify incomplete metadata associations.

    Practical Recommendations

    In routine cryo-electron tomography workflows, this protocol is especially useful
    when handedness inconsistencies are detected before reconstruction or subtomogram
    extraction.

    A few practical considerations are important:

    - Verify that handedness inversion is truly required before applying the protocol.
    - If CTF metadata exist, provide them during execution so relational consistency
      is preserved automatically.
    - Carefully inspect reported non-matching tilt-series identifiers when working
      with partially processed datasets.

    In most biological applications, correcting handedness early in the workflow is
    preferable to correcting downstream reconstructed maps.

    Final Perspective

    The Invert Tilt Angles protocol is a metadata-level geometric correction tool.

    Rather than modifying image data, it changes the angular interpretation of the
    tilt acquisition while preserving internal consistency across tilt-series and
    optional CTF metadata.

    In modern cryo-ET workflows, this protocol provides a simple but biologically
    important mechanism for ensuring that tomographic geometry reflects the correct
    physical handedness before downstream analysis.
    """
    _label = 'invert tilt angles'
    _devStatus = BETA
    _possibleOutputs = InvertTiltsOutputs
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.ctfDict = None
        self.nonMatchingTsIdsMsg = String()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS_SET, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-Series')
        form.addParam(IN_CTF_SET, PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label='CTF (opt)',
                      allowsNull=True,
                      help='Introducing the CTFs will update the pointer from introduced CTFs to '
                           'the tilt-series with the inverted tilt-angles in order to keep the coherence in the '
                           'relationship between both objects after the angle inversion operation.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self._invertAnglesStep, tsId,
                                     needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        inTsSet = self._getInTsSet()
        inCtfs = self._getInCtfSet()
        if inCtfs:
            tsIds = set(inTsSet.getTSIds())
            ctfTsIds = set(inCtfs.getTSIds())
            # Check the common elements
            matchingTsIds = tsIds & ctfTsIds
            nonMatchingTsIds = tsIds ^ ctfTsIds
            if not matchingTsIds:
                raise Exception('No matching tsIds were found among the given sets of tilt-series and CTFs.')
            if nonMatchingTsIds:
                msg = f'Some non-matching tsIds were found: {nonMatchingTsIds}'
                self.nonMatchingTsIdsMsg.set(msg)
                logger.info(cyanStr(msg))
                self._store(self.nonMatchingTsIdsMsg)
            self.tsDict = {ts.getTsId(): ts.clone() for ts in inTsSet if ts.getTsId() in matchingTsIds}
            self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inCtfs if ctf.getTsId() in matchingTsIds}
        else:
            self.tsDict = {ts.getTsId(): ts.clone() for ts in inTsSet}

    def _invertAnglesStep(self, tsId: str):
        inTs = self.tsDict[tsId]
        outTsSet = self._getOutputTsSet()
        newTs = TiltSeries()
        newTs.copyInfo(inTs)
        outTsSet.append(newTs)

        for inTi in inTs.iterItems():
            newTi = inTi.clone()
            newTi.setTiltAngle(-1 * inTi.getTiltAngle())
            newTs.append(newTi)

        newTs.write()
        outTsSet.update(newTs)
        outTsSet.write()

        # Generate the output CTFs
        if self.ctfDict:
            inCtf = self.ctfDict[tsId]
            outCtfSet = self._getOutputCtfSet()
            newCtf = inCtf.clone()
            newCtf.setTiltSeries(self.tsDict[tsId])
            outCtfSet.append(newCtf)

            for ctfTomo in inCtf.iterItems():
                newCtfTomo = ctfTomo.clone()
                newCtf.append(newCtfTomo)

            newCtf.write()
            outCtfSet.update(newCtf)
            outCtfSet.write()

    # --------------------------- UTILS functions -----------------------------
    def _getInTsSet(self, returnPointer: bool = False) -> Union[SetOfTiltSeries, Pointer]:
        inTsPointer = getattr(self, IN_TS_SET)
        return inTsPointer if returnPointer else inTsPointer.get()

    def _getInCtfSet(self, returnPointer: bool = False) -> Union[SetOfCTFTomoSeries, Pointer]:
        inCtfsPointer = getattr(self, IN_CTF_SET)
        return inCtfsPointer if returnPointer else inCtfsPointer.get()

    def _getOutputTsSet(self) -> SetOfTiltSeries:
        outSetSetAttrib = self._possibleOutputs.tiltSeries.name
        outTsSet = getattr(self, outSetSetAttrib, None)
        if not outTsSet:
            outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
            outTsSet.copyInfo(self._getInTsSet())
            self._defineOutputs(**{outSetSetAttrib: outTsSet})
            self._defineSourceRelation(self._getInTsSet(returnPointer=True), outTsSet)
        return outTsSet

    def _getOutputCtfSet(self) -> SetOfCTFTomoSeries:
        outSetSetAttrib = self._possibleOutputs.ctfs.name
        outCtfSet = getattr(self, outSetSetAttrib, None)
        if not outCtfSet:
            outCtfSet = SetOfCTFTomoSeries.create(self._getPath(), template='ctfs')
            outCtfSet.copyInfo(self._getInCtfSet())
            outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name)
            outCtfSet.setSetOfTiltSeries(outTsSet)
            self._defineOutputs(**{outSetSetAttrib: outCtfSet})
            self._defineSourceRelation(self._getInCtfSet(returnPointer=True), outCtfSet)
        return outCtfSet

    # --------------------------- INFO functions ------------------------------
    def _summary(self) -> list:
        msgList = []
        nonMatchingTsIdsMsg = self.nonMatchingTsIdsMsg.get()
        if nonMatchingTsIdsMsg:
            msgList.append(f'*{nonMatchingTsIdsMsg}*')
        return msgList
