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

    """
    Creates a new set of subtomograms from selected 2D particles or 2D classes.
    The protocol establishes the relationship between 2D particle selections and
    their originating subtomograms, allowing users to recover only the subtomograms
    associated with specific particles or classes of interest.

    AI Generated:

    2D Particles to Subtomograms (Prot2DParticlesToSubtomograms) — User Manual
        Overview

        The 2D Particles to Subtomograms protocol generates a filtered subset of
        subtomograms based on selections performed at the 2D particle level. Its
        main purpose is to reconnect particle-based classification or selection
        results with the original subtomograms from which those particles were
        extracted. In cryo-electron tomography workflows, this protocol becomes
        especially useful after 2D classification, particle cleaning, or quality
        assessment steps, where users wish to retain only the subtomograms linked
        to biologically meaningful or high-quality particle projections.

        From a biological perspective, the protocol provides a bridge between
        2D analysis and 3D subtomogram processing. Researchers commonly perform
        extensive cleaning or classification in 2D because it is computationally
        efficient and visually intuitive. Once relevant particles or classes have
        been identified, this protocol allows the recovery of the corresponding
        subtomograms for downstream tomographic refinement, averaging, or structural
        interpretation.

        Inputs and General Workflow

        The protocol requires two primary inputs. The first input is a set of
        subtomograms, representing the original 3D particles extracted from
        tomographic data. The second input is either a set of 2D classes or a
        direct set of particles. These particles are expected to contain references
        to the subtomograms from which they originated.

        During execution, the protocol analyzes the selected particles or classes
        and retrieves the identifiers associated with their parent subtomograms.
        It then creates a new output set containing only the corresponding
        subtomograms. The resulting subset preserves the original metadata and
        structural information while excluding subtomograms unrelated to the
        selected particles.

        Biological Relevance of Particle Selection

        In many cryo-ET workflows, 2D classification acts as an important quality
        control step. Poorly aligned particles, contaminants, damaged particles,
        or projections dominated by noise are frequently removed at this stage.
        By transferring the selection back to the subtomogram level, this protocol
        ensures that only biologically relevant or structurally consistent
        subtomograms continue into downstream analyses.

        This strategy is particularly important when working with heterogeneous
        datasets, low signal-to-noise ratios, or large in situ datasets where
        manual inspection of subtomograms would be impractical. Selecting
        subtomograms through validated 2D particle populations often improves the
        quality of subsequent subtomogram averaging and classification steps.

        Working with 2D Classes

        When a set of 2D classes is provided, the protocol iterates through each
        class and retrieves all associated subtomogram identifiers. This workflow
        is commonly used after 2D classification procedures where users manually
        select classes displaying clear structural features while excluding classes
        dominated by noise or artifacts.

        From a practical perspective, this enables efficient dataset cleaning.
        Instead of manually tracking subtomogram identifiers, the protocol
        automatically reconstructs the corresponding subset directly from the
        selected 2D classes.

        Working with Particle Sets

        The protocol can also operate directly on a set of particles without the
        need for classification. In this scenario, the selected particles define
        the subtomograms that will be retained in the output.

        This approach is useful in workflows involving automated scoring,
        particle filtering, or selection methods based on confidence metrics,
        template matching scores, or deep-learning predictions. Any subset of
        particles carrying valid subtomogram references can be converted into
        a corresponding subtomogram subset.

        Outputs and Their Interpretation

        After execution, the protocol produces a new set of subtomograms
        containing only the entries associated with the selected particles or
        classes. The output preserves the original subtomogram metadata and can
        be used directly in subsequent tomography workflows.

        Biologically, the resulting subset represents a cleaner and more focused
        population of particles. Depending on the upstream selection criteria,
        this may correspond to a specific conformational state, a structurally
        homogeneous population, or simply a higher-quality subset suitable for
        refinement and averaging.

        The protocol also reports the number of retained subtomograms and the
        number excluded from the original dataset. This information provides a
        simple but useful indication of the stringency of the selection process.

        Practical Recommendations

        In routine cryo-ET workflows, this protocol is most effective after
        careful 2D classification and visual inspection of classes. Retaining
        only classes with recognizable structural features generally leads to
        improved subtomogram averages and more stable downstream refinements.

        Users should also ensure that particle metadata correctly preserves the
        subtomogram identifiers. If the linkage between particles and their
        originating subtomograms has been lost during previous processing steps,
        the protocol will not be able to reconstruct the correct subset.

        For heterogeneous biological systems, iterative cycles of 2D cleaning
        followed by subtomogram recovery can substantially improve the quality
        of the final dataset while reducing computational cost in later stages
        of analysis.

        Final Perspective

        For cryo-electron tomography users, this protocol serves as a practical
        connection between particle-based 2D analysis and subtomogram-level 3D
        processing. Although computationally simple, it plays an important role
        in maintaining dataset consistency and transferring biologically relevant
        selections across different stages of the workflow. Careful particle or
        class selection at the 2D level can significantly improve the quality,
        interpretability, and biological relevance of subsequent subtomogram
        analyses.
    """
