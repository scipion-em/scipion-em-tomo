# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# *****************************************************************************
import logging
import traceback
from collections import Counter
from typing import List, Union, Set
from pwem import (touchHeartbeat, closeStreamJournal, STREAM_HEARTBEAT_TIMEOUT,
                  genExecStatusDir)
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.utils import cyanStr, redStr, yellowStr
from tomo.objects import SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries
from tomo.utils import sleepRandomly

logger = logging.getLogger(__name__)

class ProtocolBaseStreamingTomo(ProtStreamingBase):

    # ------------------------------------------------------------------ #
    # Generic streaming step generator (shared by all tomo streaming      #
    # protocols). The only per-protocol variability is captured by the    #
    # small hook methods below; override those instead of re-implementing #
    # this loop.                                                          #
    # ------------------------------------------------------------------ #
    def stepsGeneratorStep(self) -> None:
        """Discover ready tilt-series from the input producer's append-only
        journal (never its live SQLite set) and inject the per-tilt-series
        processing steps for each new one, until the input stream closes (or
        the producer's heartbeat goes stale). Resumable: previously finished
        tilt-series are recovered by ``_streamingReadingOutput``.

        Protocol-specific behaviour is provided by the hooks:
        ``_streamingInitialize`` (optional setup), ``_getStreamingInputTs``
        (input set), ``_getProcessedTsIds`` (tracking list),
        ``_getStreamingOutputNames`` (output name/s to close at the end) and
        ``_insertCommonSteps`` (the per-tilt-series steps).
        """
        self._streamingInitialize()
        closeSetStepDeps = []
        inTsSet = self._getStreamingInputTs()
        genExecStatusDir(self)
        self._streamingReadingOutput()
        processedTsIds = self._getProcessedTsIds()

        while True:
            try:
                # Discover ready tsIds from the producer's append-only journal
                # (filesystem), not from its live SQLite set.
                inTsIds = set(inTsSet.getTSIds())
                if self._stopGeneratingSteps(inTsSet,
                                             inTsIds=inTsIds,
                                             tsIdReadList=processedTsIds,
                                             outputNames=self._getStreamingOutputNames(),
                                             closeSetStepDeps=closeSetStepDeps):
                    break

                nonProcessedTsIds = inTsIds - set(processedTsIds)
                if nonProcessedTsIds:
                    # Rebuild each new tilt-series in memory from the producer's
                    # JSON sidecar (no producer-DB read).
                    tsToProcessDict = inTsSet.fetchNewTs(nonProcessedTsIds)
                    for tsId, ts in tsToProcessDict.items():
                        self._insertCommonSteps(ts, closeSetStepDeps)
                        logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                        processedTsIds.append(tsId)

                sleepRandomly()

            except Exception as e:
                logger.warning(yellowStr(f'stepsGeneratorStep failed with exception: {e}.'))
                logger.error(traceback.format_exc())
                sleepRandomly()
                continue

    def _streamingReadingOutput(self) -> None:
        """Repopulate the processed-tsId tracking list from the already-created
        output on a continued/resumed run, so previously finished tilt-series
        are not re-inserted. Reads the primary output set (see
        ``_getReadingOutputName``); a tilt-series counts as processed once that
        output exists for it.

        Named ``_streamingReadingOutput`` (not ``readingOutput``) on purpose:
        some bases in the MRO -- e.g. ``ProtImodBase`` -- already define an
        argument-based ``readingOutput`` that would shadow a base method of that
        name (same reasoning as ``_closeStreamOutputsStep``).
        """
        outSet = getattr(self, self._getReadingOutputName(), None)
        processedTsIds = self._getProcessedTsIds()
        if outSet:
            for item in outSet:
                processedTsIds.append(item.getTsId())
            logger.info(cyanStr(f'TsIds processed: {processedTsIds}'))
        else:
            logger.info(cyanStr('No tilt-series have been processed yet'))

    # ----------------------------- hooks ------------------------------- #
    def _streamingInitialize(self) -> None:
        """Optional per-protocol setup run once before the generator loop.
        Default: no-op."""
        pass

    def _getStreamingInputTs(self):
        """Return the input SetOfTiltSeries being consumed (streaming producer)."""
        raise NotImplementedError('_getStreamingInputTs must be implemented by %s'
                                  % self.getClassName())

    def _getProcessedTsIds(self) -> List[str]:
        """Return the mutable list tracking the tsIds already processed. Must
        return the SAME list object on every call (it is appended in place)."""
        raise NotImplementedError('_getProcessedTsIds must be implemented by %s'
                                  % self.getClassName())

    def _getStreamingOutputNames(self) -> Union[List[str], str]:
        """Return the output attribute name (str) or names (list) to
        validate + close when the input stream ends."""
        raise NotImplementedError('_getStreamingOutputNames must be implemented by %s'
                                  % self.getClassName())

    def _getReadingOutputName(self) -> str:
        """Return the single output name scanned by _streamingReadingOutput to
        recover already-processed tsIds. Default: the (first) output name -- for
        multi-output protocols the list is expected to be ordered so its first
        element is the representative 'one output per processed tsId'."""
        names = self._getStreamingOutputNames()
        return names if isinstance(names, str) else names[0]

    def _insertCommonSteps(self, ts, closeSetStepDeps: List[int]) -> None:
        """Insert the per-tilt-series processing steps and append the id of the
        final (output) step to ``closeSetStepDeps``. Implemented per protocol."""
        raise NotImplementedError('_insertCommonSteps must be implemented by %s'
                                  % self.getClassName())

    def _stopGeneratingSteps(self,
                             inSet: Union[SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries],
                             inTsIds: Set[str],
                             tsIdReadList: List[str],
                             outputNames: Union[List[str], str],
                             closeSetStepDeps: List[int]) -> bool:

        closeInputSets = False
        # Refresh this protocol's heartbeat so its own consumers can tell it
        # is alive even during long gaps with no new tilt-series.
        touchHeartbeat(self)

        if inSet.isStreamClosed() and Counter(tsIdReadList) == Counter(inTsIds):
            logger.info(cyanStr('Input set closed.\n'))
            self._insertFunctionStep(self._closeStreamOutputsStep,
                                     outputNames,
                                     prerequisites=closeSetStepDeps,
                                     needsGPU=False)
            closeInputSets = True

        # Producer-liveness: if the stream was never closed but the producer's
        # heartbeat is stale, it likely died. Close gracefully with whatever
        # was processed instead of looping forever.
        if not inSet.isStreamClosed():
            hbAge = inSet.getProducerHeartbeatAge()
            if hbAge is not None and hbAge > STREAM_HEARTBEAT_TIMEOUT:
                logger.error(redStr(
                    f'Producer heartbeat stale ({hbAge:.0f}s) and stream not '
                    f'closed; closing with partial outputs.'))
                self._insertFunctionStep(self._closeStreamOutputsStep,
                                         outputNames,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                closeInputSets= True

        return closeInputSets

    def _closeStreamOutputsStep(self, outputNames: Union[List[str], str]) -> None:
        """ Centralized stream shutdown for tomo streaming protocols.

        Validates + closes the named output sets (Protocol._closeOutputSet raises
        if any named output is empty) AND publishes the terminal record to this
        protocol's own stream journal, so its downstream streaming consumers can
        detect completion via isStreamClosed(). A uniquely named step (not
        'closeOutputSetsStep') so it never collides via MRO with a non-streamified
        base such as ProtImodBase that defines its own closeOutputSetsStep.
        """
        self._closeOutputSet(outputNames)
        closeStreamJournal(self)
