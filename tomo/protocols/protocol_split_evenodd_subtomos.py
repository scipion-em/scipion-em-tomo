# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam
from tomo.protocols import ProtTomoBase


class ProtSplitEvenOddSubtomoSet(EMProtocol, ProtTomoBase):
    """ Protocol to split even/odd subtomograms.
    """
    _label = 'split subtomos even/odd'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Set of subtomograms",
                      help='Select the set of subtomograms that you '
                           'want to split in even/odd sets.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        evenSet = self._createSetOfSubTomograms(suffix='_even')
        evenSet.copyInfo(inputSet)
        oddSet = self._createSetOfSubTomograms(suffix='_odd')
        oddSet.copyInfo(inputSet)

        for subtomo in inputSet:
            if subtomo.getObjId() % 2 == 0:
                evenSet.append(subtomo)
            else:
                oddSet.append(subtomo)

        self._defineOutputs(outputset_even=evenSet)
        self._defineSourceRelation(inputSet, evenSet)
        self._defineOutputs(outputset_odd=oddSet)
        self._defineSourceRelation(inputSet, oddSet)

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        if not self.isFinished():
            return["Output sets not ready yet."]
        else:
            return["We have split the input subtomogram set in even and odd sets."]
