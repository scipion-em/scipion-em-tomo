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

from pwem.objects import VolumeMask
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
from pyworkflow.protocol import PointerParam, BooleanParam
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTomograms, SetOfTomoMasks, TomoMask

logger = logging.getLogger(__name__)
IN_TOMO_SET = 'inTomoSet'
IN_MASK_SET = 'inMaskSet'


class TomoApplyMaskOutputs(Enum):
    maskedTomograms = SetOfTomograms


class ProtTomoApplyMask(EMProtocol):
    """This protocol applies a mask or a set of masks to a given set of tomograms. If a single mask is given, it will
    be applied to all the tomograms that compose the given set of tomograms. If a set of masks is provided, the
    protocol will try to match the tomograms and the masks by tsId. Once the mask/s are applied, the voxels whose
    value is zero may be optionally filled with noise."""

    _label = 'apply mask/s to tomograms'
    _devStatus = BETA
    _possibleOutputs = TomoApplyMaskOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomosDict = None
        self.tomoMaskDict = None
        self.mask = None
        self.nonMatchingTsIdsMsg = String()

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET, PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')
        form.addParam(IN_MASK_SET, PointerParam,
                      pointerClass="VolumeMask, SetOfTomoMasks",
                      important=True,
                      label='Masks',
                      help='If a set of masks is provided, the protocol will try to match the '
                           'tomograms and the masks by tsId')
        form.addParam('fillWithNoise', BooleanParam,
                      default=False,
                      label='Fill blank spaces with noise?')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tomosDict.keys():
            self._insertFunctionStep(self.applyMaskStep, tsId, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        inTomos = self._getInTomoSet()
        inTomoMasks = self._getInTomoMasks()
        if type(inTomoMasks) == SetOfTomoMasks:
            logger.info(cyanStr('A set of tomo masks was provided'))
            tomosTsIds = set(inTomos.getTSIds())
            tomoMaskTsIds = set(inTomoMasks.getTSIds())
            # Check the common elements
            matchingTsIds = tomosTsIds & tomoMaskTsIds
            nonMatchingTsIds = tomosTsIds ^ tomoMaskTsIds
            if not nonMatchingTsIds:
                raise Exception('No matching tsIds were found among the given sets of tomogramas and tomo masks.')
            if nonMatchingTsIds:
                msg = f'Some non-matching tsIds were found: {nonMatchingTsIds}'
                self.nonMatchingTsIdsMsg.set(msg)
                logger.info(cyanStr(msg))
                self._store(self.nonMatchingTsIdsMsg)
            self.tomosDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos
                              if tomo.getTsId() in matchingTsIds}
            self.tomoMaskDict = {tomoMask.getTsId(): tomoMask.clone() for tomoMask in inTomoMasks
                                 if tomoMask.getTsId() in matchingTsIds}
        else:
            logger.info(cyanStr('A single mask was provided'))
            self.tomosDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos}
            self.mask = inTomoMasks

    def applyMaskStep(self, tsId: str):
        mask = self._getCurrentMask(tsId)
        # TODO: finish this

    # --------------------------- UTILS functions -----------------------------
    def _getInTomoSet(self, returnPointer: bool = False) -> Union[SetOfTomograms, Pointer]:
        inTomoSetPointer = getattr(self, IN_TOMO_SET)
        return inTomoSetPointer if returnPointer else inTomoSetPointer.get()

    def _getInTomoMasks(self) -> Union[SetOfTomoMasks, VolumeMask]:
        return getattr(self, IN_MASK_SET).get()

    def _getCurrentMask(self, tsId: str) -> Union[TomoMask, VolumeMask]:
        if self.tomoMaskDict:
            return self.tomoMaskDict[tsId]
        else:
            return self.mask

    # --------------------------- INFO functions ------------------------------
    def _summary(self) -> list:
        msgList = []
        nonMatchingTsIdsMsg = self.nonMatchingTsIdsMsg.get()
        if nonMatchingTsIdsMsg:
            msgList.append(f'*{nonMatchingTsIdsMsg}*')
        return msgList

    def _validate(self) -> list:
        errorList = []
        inTomos = self._getInTomoSet()
        inMasks = self._getInTomoMasks()
        # Check the sampling rate
        sRateTol = 1e-3
        inTomosSRate = inTomos.getSamplingRate()
        inMasksSRate = inMasks.getSamplingRate()
        if abs(inTomosSRate - inMasksSRate) > sRateTol:
            errorList.append(f'The sampling rate of the given tomograms and mask/s are different within tolerance: '
                             f'abs({inTomosSRate:.3f} - {inMasksSRate:.3f} > {sRateTol:3f}')
        # TODO: check tomo by tomo (to cover heterogeneous sets) the dimensions respecting the masks

        return errorList




