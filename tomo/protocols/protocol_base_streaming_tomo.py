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
        inputSets = self._getStreamingInputSets()
        genExecStatusDir(self)
        self._streamingReadingOutput()
        processedTsIds = self._getProcessedTsIds()

        while True:
            try:
                # Ready ids come from the producers' append-only journals (a file
                # read, never the live SQLite set). For multi-input protocols this
                # is the intersection across inputs, so an id is "ready" only once
                # ALL of its inputs have published it.
                readyTsIds = self._getReadyTsIds(inputSets)
                if self._stopGeneratingSteps(inputSets,
                                             inTsIds=readyTsIds,
                                             tsIdReadList=processedTsIds,
                                             outputNames=self._getStreamingOutputNames(),
                                             closeSetStepDeps=closeSetStepDeps):
                    break

                nonProcessedTsIds = readyTsIds - set(processedTsIds)
                if nonProcessedTsIds:
                    # Rebuild each new item in memory from the producer(s') JSON
                    # sidecar(s) (no producer-DB read). _discoverReadyWork returns
                    # {tsId: payload}; payload is the item itself for single-input
                    # protocols and a per-protocol tuple/dict for multi-input ones.
                    newWork = self._discoverReadyWork(nonProcessedTsIds, inputSets)
                    for tsId, payload in newWork.items():
                        # payload is the single item for single-input protocols, or
                        # a tuple of items for multi-input ones (from an overridden
                        # _discoverReadyWork). Normalize to a positional-args tuple
                        # and pass closeSetStepDeps BY KEYWORD, matching the
                        # `_insertCommonSteps(self, *stepsInputs, closeSetStepDeps)`
                        # contract. Passing closeSetStepDeps positionally made it land
                        # in *stepsInputs and left the keyword-only parameter unfilled
                        # -> "TypeError: _insertCommonSteps() missing 1 required
                        # keyword-only argument: 'closeSetStepDeps'".
                        stepInputs = payload if isinstance(payload, tuple) else (payload,)
                        self._insertCommonSteps(*stepInputs, closeSetStepDeps=closeSetStepDeps)
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

    def _getStreamingInputSets(self) -> List:
        """Return the list of input sets to synchronise. Used both for the
        ready-id intersection and for termination/heartbeat across all producers."""
        raise NotImplementedError('_getStreamingInputSets must be implemented by %s'
                                  % self.getClassName())

    def _getReadyTsIds(self, inputSets: List) -> Set[str]:
        """Ids ready to process this cycle: the INTERSECTION of the tsIds each
        input producer has published in its journal. For a single input this is
        just that set's ids; for multiple inputs an id becomes ready only once
        EVERY input has it, so paired work is never scheduled half-ready. Override
        only if a protocol needs a different join (e.g. union)."""
        return set.intersection(*(set(s.getTSIds()) for s in inputSets))

    def _discoverReadyWork(self, tsIds: Set[str], inputSets: List) -> dict:
        """Rebuild the ready items in memory (from sidecars, no producer-DB read)
        and return ``{tsId: payload}``. Default (single input): the item itself,
        via ``fetchNewItems``. Multi-input protocols override to fetch from each
        input and join by tsId into a payload tuple/dict, skipping a tsId whose
        partner is not materialisable yet (it is retried next cycle)."""
        return inputSets[0].fetchNewItems(tsIds)

    def _getProcessedTsIds(self) -> List[str]:
        """Mutable list tracking the tsIds already processed. It is seeded once by
        ``_streamingReadingOutput`` (resume) and appended per scheduled tsId by the
        loop, so it MUST return the SAME list object on every call. Base-owned by
        default (a lazily-created per-run list); protocols need not override it."""
        if not hasattr(self, '_streamProcessedTsIds'):
            self._streamProcessedTsIds = []
        return self._streamProcessedTsIds

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

    def _releaseOutputWriteLock(self, output, tsId: str = None) -> None:
        """Release the write lock on a streaming output set after a SQLite lock
        error, so a ``@retry_on_sqlite_lock`` retry is a clean, non-hogging redo.

        Call this from the ``except sqlite3.OperationalError`` block of a
        protocol's output-registration method, then re-raise the original error
        so the decorator retries. Rolling back here means the producer does NOT
        keep the write transaction open across the retry/backoff window --
        otherwise it starves the very readers it is waiting for (a downstream
        streaming consumer reading the shared sqlite under journal_mode=DELETE).

        The rollback + duplicate-guard-cache reset lives on the set itself:
        SetOfTiltSeriesBase / SetOfCTFTomoSeries / SetOfTomograms /
        SetOfLandmarkModels all provide ``rollbackFailedAppend`` via the shared
        ``_AppendRollbackMixin`` (each knows its own cache, if any). This wrapper
        only guarantees it never raises -- it runs inside an except handler, so a
        failure here must not mask the original lock error the caller re-raises
        to drive the retry.

        :param output: the output set being written.
        :param tsId: tsId of the item whose append failed; used to clear any
            per-tsId duplicate-guard cache. Optional for cache-less leaf sets.
        """
        try:
            output.rollbackFailedAppend(tsId)
        except Exception as e:
            logger.error(yellowStr(f'_releaseOutputWriteLock failed for '
                                   f'tsId={tsId}: {e}'))

    def _insertCommonSteps(self, *stepsInputs, closeSetStepDeps: List[int]) -> None:
        """Insert the per-tilt-series processing steps and append the id of the
        final (output) step to ``closeSetStepDeps``. Implemented per protocol.

        ``payload`` is whatever ``_discoverReadyWork`` produced for that tsId: the
        single item for single-input protocols, or the joined tuple/dict for
        multi-input ones (which unpack it, e.g. ``ts, ctf, orders = payload``)."""
        raise NotImplementedError('_insertCommonSteps must be implemented by %s'
                                  % self.getClassName())

    def _stopGeneratingSteps(self,
                             inSets: Union[List, SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries],
                             inTsIds: Set[str],
                             tsIdReadList: List[str],
                             outputNames: Union[List[str], str],
                             closeSetStepDeps: List[int]) -> bool:
        # Accept either a single input set (legacy callers) or a list of input
        # sets (multi-input protocols). Termination requires ALL inputs closed.
        if not isinstance(inSets, (list, tuple)):
            inSets = [inSets]

        closeInputSets = False
        # Refresh this protocol's heartbeat so its own consumers can tell it
        # is alive even during long gaps with no new tilt-series.
        touchHeartbeat(self)

        allClosed = all(s.isStreamClosed() for s in inSets)
        if allClosed and Counter(tsIdReadList) == Counter(inTsIds):
            logger.info(cyanStr('Input set(s) closed.\n'))
            self._insertFunctionStep(self._closeStreamOutputsStep,
                                     outputNames,
                                     prerequisites=closeSetStepDeps,
                                     needsGPU=False)
            closeInputSets = True

        # Producer-liveness: if some input stream was never closed but its
        # producer's heartbeat is stale, it likely died. Close gracefully with
        # whatever was processed instead of looping forever.
        if not allClosed and not closeInputSets:
            for s in inSets:
                if s.isStreamClosed():
                    continue
                hbAge = s.getProducerHeartbeatAge()
                if hbAge is not None and hbAge > STREAM_HEARTBEAT_TIMEOUT:
                    logger.error(redStr(
                        f'Producer heartbeat stale ({hbAge:.0f}s) and stream not '
                        f'closed; closing with partial outputs.'))
                    self._insertFunctionStep(self._closeStreamOutputsStep,
                                             outputNames,
                                             prerequisites=closeSetStepDeps,
                                             needsGPU=False)
                    closeInputSets = True
                    break

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
