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

import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from tomo.objects import SetOfCTFTomoSeries


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
                      default=0, help="Criteria to validate: \n\n"
                                      "defocus tolerance: Calculate the defocus deviation taking into account a threshold(percent tolerance) respect to the mean.")
        form.addParam('tolerance', params.FloatParam, label='Tolerance percent',
                      condition='criteria==0', default=20,
                      help="Percent Tolerance")

    def _insertAllSteps(self):
        self._insertFunctionStep(self.ctfValidateStep)
        self._insertFunctionStep(self.createOutputStep)

    def ctfValidateStep(self):
        """
        Validate all ctf tomo series and separate into two sets(good and bad
        following the selected criteria)
        """
        ctfSeries = self.inputCtfTomoSeries.get()
        self.goodCtfName = 'goodSetOfCTFTomoSeries'
        self.badCtfName = 'badSetOfCTFTomoSeries'
        self.outputCtfName = 'outputSetOfCTFTomoSeries'

        self.outputSetOfCtfTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                  template='CTFmodels%s.sqlite')
        self.outputSetOfCtfTomoSeries.setSetOfTiltSeries(ctfSeries.getSetOfTiltSeries(pointer=True).get())

        if self.criteria.get() == 0:  # Defocus angles case
            self._validateCtfDefocusDeviation(ctfSeries, self.tolerance.get())

    def _validateCtfDefocusDeviation(self, ctfSeries, tolerance):
        """
        Validate the set of ctf tomo series taking into account de defocus angle
        deviation
        """
        # Creating a new copy of ctf tomo series and recalculate the defocus
        # deviation
        self.outputSetOfgoodCtfTomoSeries = None
        self.outputSetOfbadCtfTomoSeries = None

        for ctfSerie in ctfSeries:
            newCTFTomoSeries = ctfSerie.clone()
            newCTFTomoSeries._isDefocusUDeviationInRange.set(True)
            newCTFTomoSeries._isDefocusVDeviationInRange.set(True)
            self.outputSetOfCtfTomoSeries.append(newCTFTomoSeries)
            for item in ctfSerie.iterItems():
                ctfEstItem = item.clone()
                newCTFTomoSeries.append(ctfEstItem)
            newCTFTomoSeries.calculateDefocusUDeviation(defocusUTolerance=tolerance)
            newCTFTomoSeries.calculateDefocusVDeviation(defocusVTolerance=tolerance)
            newCTFTomoSeries.write()
            self.outputSetOfCtfTomoSeries.update(newCTFTomoSeries)

        for ctfSerie in self.outputSetOfCtfTomoSeries:
            ctfSerieClon = ctfSerie.clone()
            # Store good and bad ctf series
            if ctfSerie.getIsDefocusUDeviationInRange():
                if self.outputSetOfgoodCtfTomoSeries is None:
                    self.outputSetOfgoodCtfTomoSeries = self.outputSetOfCtfTomoSeries.createCopy(
                        self._getExtraPath(),
                        prefix=self.goodCtfName,
                        copyInfo=True)
                else:
                    self.outputSetOfgoodCtfTomoSeries.enableAppend()
                self.outputSetOfgoodCtfTomoSeries.append(ctfSerieClon)
            else:
                if self.outputSetOfbadCtfTomoSeries is None:
                    self.outputSetOfbadCtfTomoSeries = self.outputSetOfCtfTomoSeries.createCopy(
                        self._getExtraPath(),
                        prefix=self.badCtfName,
                        copyInfo=True)
                else:
                    self.outputSetOfbadCtfTomoSeries.enableAppend()
                self.outputSetOfbadCtfTomoSeries.append(ctfSerieClon)

            for item in ctfSerie.iterItems():
                ctfEstItem = item.clone()
                ctfSerieClon.append(ctfEstItem)

    def createOutputStep(self):
        if self.outputSetOfgoodCtfTomoSeries is not None:
            self._defineOutputs(**{self.goodCtfName: self.outputSetOfgoodCtfTomoSeries})
        if self.outputSetOfbadCtfTomoSeries is not None:
            self._defineOutputs(**{self.badCtfName: self.outputSetOfbadCtfTomoSeries})

    def allowsDelete(self, obj):
        return True

    def _summary(self):
        summary = []
        if hasattr(self, 'goodSetOfCTFTomoSeries'):
            summary.append("Number of good ctf series: %d." % (self.goodSetOfCTFTomoSeries.getSize()))
        if hasattr(self, 'badSetOfCTFTomoSeries'):
            summary.append("Number of bad ctf series: %d." % (self.badSetOfCTFTomoSeries.getSize()))

        return summary

