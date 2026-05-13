# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
import logging
from collections import Counter
from enum import Enum
from typing import List, Union

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
from pyworkflow.protocol import PointerParam, STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr

from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtAssignExcludedViews(EMProtocol):
    """
    Transfer excluded-view annotations (_enabled state) from a source set of
    tilt-series to a target set of tilt-series.

    TiltSeries are matched by tsId. Images within each matched pair are matched
    by acquisition order, so the protocol is robust to re-stacking. Target
    tilt-series without a matching source are preserved unchanged in the output.
    """

    _label = 'assign excluded views'
    _devStatus = BETA
    _possibleOutputs = outputObjects
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matchedMsg = String()
        self.unmatchedMsg = String()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSourceTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Source tilt-series (with exclusions)',
                      help='Set of tilt-series from which the excluded-view '
                           'annotations (_enabled state) will be read.')

        form.addParam('inputTargetTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Target tilt-series (to assign exclusions)',
                      help='Set of tilt-series to which the excluded-view '
                           'annotations will be assigned. The output is a '
                           'copy of this set with updated _enabled metadata '
                           'for matched tilt-series.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.targetTsDict.keys():
            self._insertFunctionStep(self._transferExclusionsStep, tsId,
                                     needsGPU=False)

    # --------------------------- STEPS functions ----------------------------
    def _initialize(self):
        sourceTsSet = self._getSourceTsSet()
        targetTsSet = self._getTargetTsSet()
        sourceTsIds = set(sourceTsSet.getTSIds())
        targetTsIds = set(targetTsSet.getTSIds())

        matchedTsIds = sorted(sourceTsIds & targetTsIds)
        unmatchedTsIds = sorted(targetTsIds - sourceTsIds)

        if matchedTsIds:
            self.matchedMsg.set(", ".join(matchedTsIds))
            self._store(self.matchedMsg)
        if unmatchedTsIds:
            self.unmatchedMsg.set(", ".join(unmatchedTsIds))
            self._store(self.unmatchedMsg)

        self.sourceTsDict = {ts.getTsId(): ts.clone() for ts in sourceTsSet
                             if ts.getTsId() in matchedTsIds}
        self.targetTsDict = {ts.getTsId(): ts.clone() for ts in targetTsSet}

    def _transferExclusionsStep(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - transferring excluded views..."))
        targetTs = self.targetTsDict[tsId]
        sourceTs = self.sourceTsDict.get(tsId, None)
        outTsSet = self._getOutputTsSet()

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(targetTs)
        outTsSet.append(newTs)

        if sourceTs:
            # Build acquisition order -> _enabled mapping from source
            sourceAcqMap = {ti.getAcquisitionOrder(): ti.isEnabled()
                           for ti in sourceTs}

            for tiTarget in targetTs.iterItems():
                newTi = tiTarget.clone()
                acqOrder = tiTarget.getAcquisitionOrder()
                if acqOrder in sourceAcqMap:
                    newTi.setEnabled(sourceAcqMap[acqOrder])
                # If no matching acqOrder in source, preserve target's _enabled
                newTs.append(newTi)
        else:
            # Unmatched target TiltSeries: copy images unchanged
            for tiTarget in targetTs.iterItems():
                newTi = tiTarget.clone()
                newTs.append(newTi)

        newTs.setDim(targetTs.getDim())
        newTs.write()
        outTsSet.update(newTs)
        outTsSet.write()

    # --------------------------- UTILS functions ----------------------------
    def _getSourceTsSet(self, returnPointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        return self.inputSourceTiltSeries if returnPointer else self.inputSourceTiltSeries.get()

    def _getTargetTsSet(self, returnPointer: bool = False) -> Union[Pointer, SetOfTiltSeries]:
        return self.inputTargetTiltSeries if returnPointer else self.inputTargetTiltSeries.get()

    def _getOutputTsSet(self) -> SetOfTiltSeries:
        outAttrib = self._possibleOutputs.tiltSeries.name
        outTsSet = getattr(self, outAttrib, None)
        if not outTsSet:
            outTsSet = SetOfTiltSeries.create(self._getPath(),
                                              template='tiltseries')
            outTsSet.copyInfo(self._getTargetTsSet())
            self._defineOutputs(**{outAttrib: outTsSet})
            self._defineSourceRelation(self._getSourceTsSet(returnPointer=True),
                                       outTsSet)
            self._defineSourceRelation(self._getTargetTsSet(returnPointer=True),
                                       outTsSet)
        return outTsSet

    # --------------------------- INFO functions ----------------------------
    def _validate(self) -> List[str]:
        errors = []
        sourceTsSet = self._getSourceTsSet()
        targetTsSet = self._getTargetTsSet()

        sourceTsIds = set(sourceTsSet.getTSIds())
        targetTsIds = set(targetTsSet.getTSIds())
        matchedTsIds = sourceTsIds & targetTsIds

        if not matchedTsIds:
            errors.append("No matching tsId values found between source "
                          "(%s) and target (%s) sets." %
                          (", ".join(sorted(sourceTsIds)),
                           ", ".join(sorted(targetTsIds))))
            return errors

        # Validate acquisition orders in matched TiltSeries
        for tsId in sorted(matchedTsIds):
            sourceTs = sourceTsSet.getTiltSeriesFromTsId(tsId)
            targetTs = targetTsSet.getTiltSeriesFromTsId(tsId)

            for label, ts in [("source", sourceTs), ("target", targetTs)]:
                acqOrders = [ti.getAcquisitionOrder() for ti in ts]
                if None in acqOrders:
                    errors.append("Missing acquisition order in %s "
                                  "tilt-series %s." % (label, tsId))
                validOrders = [ao for ao in acqOrders if ao is not None]
                counts = Counter(validOrders)
                duplicates = sorted([ao for ao, cnt in counts.items()
                                     if cnt > 1])
                if duplicates:
                    errors.append("Duplicate acquisition orders in %s "
                                  "tilt-series %s: %s" %
                                  (label, tsId, duplicates))

        return errors

    def _summary(self) -> List[str]:
        summary = []
        matchedMsg = self.matchedMsg.get()
        unmatchedMsg = self.unmatchedMsg.get()

        if hasattr(self, self._possibleOutputs.tiltSeries.name):
            if matchedMsg:
                summary.append("Excluded views transferred for tilt-series: "
                               "%s" % matchedMsg)
            else:
                summary.append("No tilt-series were matched between source "
                               "and target sets.")
            if unmatchedMsg:
                summary.append("No matching source tilt-series found for: "
                               "%s" % unmatchedMsg)
            else:
                summary.append("All target tilt-series were matched.")
        else:
            summary.append("Output not ready yet.")

        return summary
