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
    """This protocol gets the tilt-series that corresponds to a given set of tomograms. It is very useful for
    fiducial-less samples, when the quality of alignment results is difficult to be observed from
    the tilt-series, but the tomograms. In that case, the undesired data objects would be discarded
    at tomogram level, but further processing may be desired to be carried out with the corresponding
    tilt-series before getting the final tomograms."""

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
