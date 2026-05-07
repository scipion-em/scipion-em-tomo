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
    """
    Creates a subset of subtomograms based on selected 2D particles or classes,
    preserving the correspondence between 2D projections and their originating
    3D subtomograms.

    AI Generated:

    2D Particles to Subtomograms (Prot2DParticlesToSubtomograms) — User Manual
        Overview

        The 2D Particles to Subtomograms protocol generates a subset of
        subtomograms based on a selection performed in 2D space. Its main purpose
        is to bridge 2D classification or particle selection with the original
        3D subtomogram data, enabling users to recover the corresponding 3D
        volumes associated with specific 2D particles or classes.

        In cryo-electron tomography workflows, this protocol is especially useful
        after performing 2D classification or particle cleaning. It allows users
        to translate decisions made in 2D—such as selecting high-quality classes
        or removing junk particles—into a refined set of 3D subtomograms.

        For biological users, this step is important when filtering data prior to
        subtomogram averaging, structural classification, or downstream analysis,
        ensuring that only biologically meaningful particles are retained.

        Inputs and General Workflow

        The protocol requires two main inputs: a reference set of subtomograms
        and a set of 2D particles or 2D classes derived from those subtomograms.

        Each 2D particle contains a reference to its originating subtomogram
        through an internal identifier. The protocol uses this relationship to
        determine which subtomograms should be included in the output.

        The workflow consists of identifying all unique subtomogram identifiers
        present in the selected 2D dataset and retrieving the corresponding
        subtomograms from the original set.

        If the input is a set of 2D classes, the protocol iterates over all
        classes and gathers subtomograms from each class. If the input is a
        direct set of particles, the extraction is performed directly from the
        particle set.

        Mapping Between 2D and 3D Data

        The key concept behind this protocol is the mapping between 2D particles
        and their parent 3D subtomograms.

        Each particle carries metadata linking it to a specific subtomogram.
        By aggregating these identifiers, the protocol reconstructs the subset
        of 3D volumes that contributed to the selected 2D dataset.

        This ensures that the resulting subtomogram set is fully consistent
        with the selection performed in 2D, preserving biological relevance.

        Handling of Duplicates and Selection Logic

        The protocol extracts unique subtomogram identifiers from the input
        dataset. This means that even if multiple particles correspond to the
        same subtomogram, that subtomogram will only be included once in the
        output set.

        This behavior is important for avoiding redundancy and ensuring that
        the resulting dataset reflects unique 3D volumes rather than repeated
        projections.

        Output Generation

        The output is a new SetOfSubTomograms object, which is a filtered copy
        of the original subtomogram set.

        Each selected subtomogram is cloned from the input set and appended to
        the new output set, preserving metadata and acquisition information.

        The output maintains consistency with the original dataset while
        containing only the subtomograms associated with the selected 2D
        particles or classes.

        Interpretation of Results

        The resulting subset represents the 3D volumes corresponding to the
        selected 2D dataset.

        Biologically, this allows users to refine their dataset based on visual
        or statistical criteria applied in 2D, such as selecting well-defined
        classes, removing noisy particles, or isolating specific structural
        features.

        This step is often used before subtomogram averaging, where dataset
        quality has a direct impact on resolution and interpretability.

        Summary Information

        The protocol reports the number of subtomograms included in the output
        as well as the number of excluded subtomograms relative to the original
        dataset.

        This provides a quick overview of how restrictive the 2D-based selection
        was and helps users assess the impact of their filtering strategy.

        Practical Recommendations

        When using 2D classification results, it is advisable to carefully
        select classes that represent meaningful biological features while
        excluding poorly defined or noisy classes.

        If working with individual particles, prior cleaning steps—such as
        removing outliers or low-quality picks—will improve the quality of
        the resulting subtomogram set.

        Since the protocol removes duplicate references to the same
        subtomogram, users should be aware that the output reflects unique
        volumes rather than particle counts.

        Final Perspective

        This protocol provides a direct and practical link between 2D analysis
        and 3D tomographic data.

        By enabling users to propagate 2D-based decisions back to the original
        subtomograms, it plays a key role in refining datasets and improving
        the biological interpretability of downstream cryo-EM and cryo-ET
        analyses.
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


