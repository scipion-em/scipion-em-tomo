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
from collections import Counter
from typing import List, Union, Set
from pwem import touchHeartbeat, closeStreamJournal, STREAM_HEARTBEAT_TIMEOUT
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.utils import cyanStr, redStr
from tomo.objects import SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries

logger = logging.getLogger(__name__)

class ProtocolBaseStreamingTomo(ProtStreamingBase):

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
