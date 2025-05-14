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
from enum import Enum
from typing import Union

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries, TiltSeries

IN_TS_SET = 'inTsSet'


class InvertTiltsOutputs(Enum):
    tiltSeries = SetOfTiltSeries


class ProtInverTiltAngles(EMProtocol):
    """This protocol inverts the physical handedness of the introduced tilt-series by inverting the tilt angles
    in the metadata associated to each tilt.-series."""
    _label = 'invert tilt angles'
    _devStatus = BETA
    _possibleOutputs = InvertTiltsOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS_SET, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-Series')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self._invertAnglesStep, tsId, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone() for ts in self._getInTsSet()}

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

    # --------------------------- UTILS functions -----------------------------
    def _getInTsSet(self, returnPointer: bool = False) -> Union[SetOfTiltSeries, Pointer]:
        inTsPointer = getattr(self, IN_TS_SET)
        return inTsPointer if returnPointer else inTsPointer.get()

    def _getOutputTsSet(self) -> SetOfTiltSeries:
        outSetSetAttrib = self._possibleOutputs.tiltSeries.name
        outTsSet = getattr(self, outSetSetAttrib, None)
        if not outTsSet:
            outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
            outTsSet.copyInfo(self._getInTsSet())
            self._defineOutputs(**{outSetSetAttrib: outTsSet})
            self._defineSourceRelation(self._getInTsSet(returnPointer=True), outTsSet)
        return outTsSet
