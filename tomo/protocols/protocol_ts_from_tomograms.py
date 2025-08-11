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
from pyworkflow.utils import Message
from tomo.objects import SetOfTiltSeries, SetOfTomograms

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
                      label='Tilt-Series (opt.)',
                      help='If not provided, the protocol will try to find the corresponding '
                           'tilt-series via data relations.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self._getTsFromTomosStep, tsId,
                                     needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        inTsSet = self._getInTsSet()

    def _getTsFromTomosStep(self, tsId: str):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _getInTsSet(self, returnPointer: bool = False) -> Union[SetOfTiltSeries, Pointer]:
        inTsPointer = getattr(self, IN_TS_SET)
        return inTsPointer if returnPointer else inTsPointer.get()

    def _getInTomoSet(self, returnPointer: bool = False) -> Union[SetOfTomograms, Pointer]:
        inTomosPointer = getattr(self, IN_TOMO_SET)
        return inTomosPointer if returnPointer else inTomosPointer.get()

    # def _getTsFromRelations(self) -> Union[SetOfTiltSeries, None]:
    #     inTomos = self._getFormAttrib(IN_CTF_SET)
    #     return getObjFromRelation(inTomos, self, SetOfTiltSeries)
