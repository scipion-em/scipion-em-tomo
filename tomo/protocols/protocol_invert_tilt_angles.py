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


class ProtInverTiltAngles(EMProtocol):
    """This protocol inverts the physical handedness of the introduced tilt-series by inverting the tilt angles
    in the metadata associated to each tilt.-series. Introducing the CTFs will update the pointer from introduced
    CTFs to the tilt-series with the inverted tilt-angles in order to keep the coherence in the
    relationship between both objects after the angle inversion operation."""
    _label = 'invert tilt angles'
    _devStatus = BETA
    _possibleOutputs = InvertTiltsOutputs
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.ctfDict = None
        self.nonMatchingTsIdsMsg = String()

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
