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
import enum

from pwem.objects import SetOfClasses2D, SetOfParticles
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtTomoBase


class importSubTomograms(enum.Enum):
    outputSetOfSubtomograms = SetOfSubTomograms()


class Prot2DParticlesToSubtomograms(EMProtocol, ProtTomoBase):
    """ Protocol to create a set of subtomograms from a selected 2D particles.
    """
    _label = '2D particles to subtomograms'
    _devStatus = BETA
    _possibleOutputs = importSubTomograms

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSubtomogramSet', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Set of subtomograms",
                      help='Select the set of subtomograms ')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfClasses2D, SetOfParticles',
                      label="Input set",
                      help='Select the 2D classes or a set of particles')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------

    def createOutputStep(self):
        subtomogramsSet = self.inputSubtomogramSet.get()
        inputParticles = self.inputSet.get()
        self.outputSubtomograms = 'outputSubtomograms'

        self.outputSetOfSubtomograms = subtomogramsSet.createCopy(self._getExtraPath(),
                                                                  prefix=self.outputSubtomograms,
                                                                  copyInfo=True)
        if isinstance(inputParticles, SetOfClasses2D):
            for clazz in inputParticles.iterItems():
                self._appendSubtomograms(clazz)
        elif isinstance(inputParticles, SetOfParticles):
            self._appendSubtomograms(inputParticles)

        self._defineOutputs(**{self.outputSubtomograms: self.outputSetOfSubtomograms})

    def _appendSubtomograms(self, inputSet):
        subTomogramsIds = inputSet.aggregate(["COUNT"], "_subtomogramID", ["_subtomogramID"])
        subTomogramsIds = [d['_subtomogramID'] for d in subTomogramsIds]

        for item in subTomogramsIds:
            subtomogramId = int(item)
            subtomogram = self.inputSubtomogramSet.get()[subtomogramId].clone()
            self.outputSetOfSubtomograms.append(subtomogram)

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSubtomograms'):
            newSize = self.outputSubtomograms.getSize()
            excluded = self.inputSubtomogramSet.get().getSize() - newSize
            summary.append("Number of subtomogram: %d" % newSize)
            summary.append("Number of excluded subtomogram: %d" % excluded)

        return summary


