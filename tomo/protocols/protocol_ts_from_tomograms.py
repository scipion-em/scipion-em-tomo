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
from pyworkflow.object import Pointer
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTiltSeries, SetOfTomograms, TiltSeries, TiltImage

logger = logging.getLogger(__name__)
IN_TOMO_SET = 'inTomoSet'
IN_TS_SET = 'inTsSet'


class OutputsTsFromTomos(Enum):
    tiltSeries = SetOfTiltSeries


class ProtTsFromTomos(EMProtocol):
    """
    Retrieves the tilt-series corresponding to a selected set of tomograms.

    AI Generated:

    Tilt-Series from Tomograms (ProtTsFromTomos) — User Manual
        Overview

        The Tilt-Series from Tomograms protocol extracts from an input set of tilt-series
        only those whose tilt-series identifiers (tsId) match the tomograms provided
        as input.

        This protocol is especially useful in cryo-electron tomography workflows where
        quality inspection is easier to perform at the tomogram level than directly on
        the tilt-series. This situation is particularly common in fiducial-less datasets,
        where reconstruction artifacts or alignment problems become more evident once
        tomograms have already been generated.

        In practical workflows, users often discard low-quality tomograms after visual
        inspection. However, later processing steps may still need to continue from the
        corresponding tilt-series rather than from the tomograms themselves. This protocol
        solves that problem by recovering only the tilt-series associated with the
        selected tomograms.

        Inputs and General Workflow

        The protocol requires two input datasets:

        - A set of tomograms
        - A set of tilt-series

        Both datasets must contain tilt-series identifiers (tsId), which define the
        correspondence between tomograms and tilt-series.

        The protocol compares the identifiers present in both inputs and computes
        their intersection. Only tilt-series whose tsId is present in both datasets
        are preserved in the output.

        Since both inputs may already represent subsets of larger datasets, the protocol
        only considers the common identifiers. Any non-matching identifiers are reported
        in the execution log for user awareness.

        Matching Strategy

        Matching is performed entirely through the tilt-series identifier (tsId).

        This means that:

        - If a tomogram and a tilt-series share the same tsId, they are considered related.
        - If an identifier appears only in one dataset, it is ignored.

        From a practical perspective, this makes the protocol very robust for curated
        workflows in which intermediate subsets are frequently created.

        However, it also means that correctness depends on proper metadata consistency.
        If tsIds are missing, duplicated, or inconsistent across datasets, expected
        matches may not occur.

        Validation and Error Handling

        The protocol performs an important validation step before generating output.

        If no common tsIds exist between the tomogram set and the tilt-series set,
        execution stops with an exception.

        This prevents producing empty outputs caused by accidental mismatch between
        datasets.

        When only partial overlap exists, execution continues normally. The protocol
        simply reports which identifiers were not common between the two inputs.

        Output Generation

        For every matching tsId, the protocol creates a new output tilt-series.

        The output preserves:

        - General metadata from the original tilt-series
        - The original tilt-images
        - The original tilt-image ordering

        Each selected tilt-series is cloned into the new output set and written to disk.

        From a workflow perspective, this means the resulting output is not merely a
        list of references, but a proper Scipion output object that can be used
        directly by downstream protocols.

        Biological Interpretation

        This protocol does not modify alignment, geometry, or image content.

        Its role is purely organizational, but biologically important.

        In cryo-ET workflows, selecting the correct subset of data is often critical.
        Poor-quality tomograms may reflect bad alignment, contamination, missing wedges,
        or reconstruction artifacts. Removing these problematic tomograms and then
        recovering the corresponding tilt-series allows users to continue processing
        only biologically meaningful data.

        Typical use cases include:

        - Discarding low-quality tomograms after visual inspection
        - Recovering the corresponding tilt-series for refinement
        - Restricting downstream analysis to curated subsets
        - Preparing clean datasets for subtomogram averaging pipelines

        Practical Recommendations

        Before running the protocol, verify that both datasets originate from compatible
        processing branches and preserve the same tsId convention.

        This protocol is especially useful after tomogram cleaning or manual curation,
        when users want to return to tilt-series space without manually tracing
        dataset correspondences.

        If the output is unexpectedly empty, the first thing to inspect should be
        metadata consistency rather than image content.

        Final Perspective

        The Tilt-Series from Tomograms protocol is a lightweight but very practical
        utility for cryo-electron tomography workflows.

        Although it performs no image processing itself, it solves a common workflow
        problem: recovering the exact tilt-series corresponding to a curated tomogram
        subset.

        In real biological projects, this can simplify data management considerably
        and helps maintain consistency between tomogram-level curation and subsequent
        tilt-series-based processing.
    """

    _label = 'tilt-series from tomograms'
    _devStatus = BETA
    _possibleOutputs = OutputsTsFromTomos

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET, PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')
        form.addParam(IN_TS_SET, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-Series')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._getTsFromTomosStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _getTsFromTomosStep(self):
        inTomosSet = self._getInTomoSet()
        inTsSet = self._getInTsSet()
        # Compute the matching tsIds among the tilt-series and the tomograms, as they both could be a subset
        tomosTsIds = set(inTsSet.getTSIds())
        tsIds = set(inTomosSet.getTSIds())
        presentTsIds = tomosTsIds & tsIds
        nonMatchingTsIds = (tomosTsIds ^ tsIds) - presentTsIds
        # Validate the intersection
        if len(presentTsIds) <= 0:
            raise Exception("There isn't any common tsIds among the tomograms and the "
                            "tilt-series introduced.")
        if len(nonMatchingTsIds) > 0:
            logger.info(cyanStr(f"TsIds not common in the introduced tomograms and "
                                f"tilt-series are: {nonMatchingTsIds}"))
        tsDict = {ts.getTsId(): ts.clone() for ts in inTsSet.iterItems() if ts.getTsId() in presentTsIds}
        # Create the output set
        outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
        outTsSet.copyInfo(inTsSet)
        self._defineOutputs(**{self._possibleOutputs.tiltSeries.name: outTsSet})
        self._defineSourceRelation(self._getInTsSet(returnPointer=True), outTsSet)
        for tsId in presentTsIds:
            inTs = tsDict[tsId]
            outTs = TiltSeries()
            outTs.copyInfo(inTs)
            outTsSet.append(outTs)
            for ti in inTs.iterItems(orderBy=TiltImage.INDEX_FIELD):
                outTi = TiltImage()
                outTi.copyInfo(ti)
                outTs.append(outTi)
            outTsSet.update(outTs)
            # Data persistence
            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()

        if len(outTsSet) == 0:
            raise Exception(f'No output/s {self._possibleOutputs.tiltSeries.name} were generated. '
                            f'Please check the Output Log > run.stdout and run.stderr')

    # --------------------------- UTILS functions -----------------------------
    def _getInTsSet(self, returnPointer: bool = False) -> Union[SetOfTiltSeries, Pointer]:
        inTsPointer = getattr(self, IN_TS_SET)
        return inTsPointer if returnPointer else inTsPointer.get()

    def _getInTomoSet(self, returnPointer: bool = False) -> Union[SetOfTomograms, Pointer]:
        inTomosPointer = getattr(self, IN_TOMO_SET)
        return inTomosPointer if returnPointer else inTomosPointer.get()
