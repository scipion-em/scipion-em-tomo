# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import os
import time

import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow import BETA


class ProtCTFTomoSeriesValidate(EMProtocol):
    """
     Validate a set of CTF tomo series and separate into two sets (good and
     bad tomo series )
    """
    _label = 'ctf validate'
    _devStatus = BETA
    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('inputCtfTomoSeries', params.PointerParam, important=True,
                      pointerClass='SetOfCTFTomoSeries',
                      label='Input ctf tomo series')
        form.addParam('criteria',  params.EnumParam,
                      choices=["defocus tolerance"], label="Criteria",
                      default=0, help="Criteria to validate")

    def _insertAllSteps(self):
        self._insertFunctionStep(self.ctfValidateStep)
        self._insertFunctionStep(self.createOutputStep)

    def ctfValidateStep(self):
        """
        Validate all ctf tomo series and separate into two sets(good and bad
        following the selected criteria)
        """
        ctfSeries = self.inputCtfTomoSeries.get()
        self.goodCTFName = 'goodSetOfCTFTomoSeries'
        self.badCTFName = 'badSetOfCTFTomoSeries'

        self.outputSetOfgoodCTFTomoSeries = ctfSeries.createCopy(self._getExtraPath(),
                                                                 prefix=self.goodCTFName,
                                                                 copyInfo=True)
        self.outputSetOfbadCTFTomoSeries = ctfSeries.createCopy(self._getExtraPath(),
                                                                prefix=self.badCTFName,
                                                                copyInfo=True)
        if self.criteria.get() == 0:  # Defocus angles case
            self._validateCtfDefocusDeviation(ctfSeries)

    def _validateCtfDefocusDeviation(self, ctfSeries):
        """
        Validate the set of ctf tomo series taking into account de defocus angle
        deviation
        """

        for ctfSerie in ctfSeries:
            ctfSerieClon = ctfSerie.clone()
            # Store good and bad ctf series
            if (ctfSerie.getIsDefocusVDeviationInRange() and
                    ctfSerie.getIsDefocusUDeviationInRange()):
                self.outputSetOfgoodCTFTomoSeries.append(ctfSerieClon)
            else:
                self.outputSetOfbadCTFTomoSeries.append(ctfSerieClon)
            for item in ctfSerie.iterItems():
                ctfEstItem = item.clone()
                ctfSerieClon.append(ctfEstItem)

    def createOutputStep(self):
        if len(self.outputSetOfgoodCTFTomoSeries) > 0:
            self._defineOutputs(**{self.goodCTFName: self.outputSetOfgoodCTFTomoSeries})
        else:
            os.remove(self._getExtraPath(self.goodCTFName+'..sqlite'))

        if len(self.outputSetOfbadCTFTomoSeries) > 0:
            self._defineOutputs(**{self.badCTFName: self.outputSetOfbadCTFTomoSeries})
        else:
            os.remove(self._getExtraPath(self.badCTFName + '..sqlite'))

    def allowsDelete(self, obj):
        return True

    def _summary(self):
        summary = []
        if self.criteria.get() == 0:

            if hasattr(self, 'goodSetOfCTFTomoSeries'):
                summary.append("Number of good ctf series: %d." % (self.goodSetOfCTFTomoSeries.getSize()))
            if hasattr(self, 'badSetOfCTFTomoSeries'):
                summary.append("Number of bad ctf series: %d." % (self.badSetOfCTFTomoSeries.getSize()))

        return summary

