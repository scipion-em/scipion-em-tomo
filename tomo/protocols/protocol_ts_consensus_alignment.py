# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
logger = logging.getLogger(__name__)

import math
import numpy as np

from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.protocol.params import MultiPointerParam, FloatParam, EnumParam
from pyworkflow.object import Set, Float

from pwem.protocols import EMProtocol

from tomo.objects import TiltSeries, TiltImage
from tomo.protocols import ProtTomoBase


class ProtConsensusAlignmentTS(EMProtocol, ProtTomoBase):
    """
    Perform a consensus of a set of alignments for the same tilt series. Returns the average alignment matrix of the
    consensus alignments and its standard deviation of shift and angle.
    """

    _label = 'Tilt-series consensus alignment'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputMultiSoTS',
                      MultiPointerParam,
                      important=True,
                      label="Input tilt series",
                      pointerClass='SetOfTiltSeries',
                      help='Select several sets of tilt-series where to evaluate the consensus in their alignment. '
                           'Output set will bring the information from the first selected set.')

        form.addParam('shiftTolerance',
                      FloatParam,
                      label="Shift tolerance (A)",
                      default=20,
                      help='Maximum shift difference between alignments to consider them as equal. it is measured in '
                           'Angstroms.')

        form.addParam('angleTolerance',
                      FloatParam,
                      label="Angle tolerance (degrees)",
                      default=3,
                      help='Maximum angle difference between alignments to consider them as equal. It is measured in '
                           'degrees')

        form.addParam('globalAlignment',
                      EnumParam,
                      choices=["Local", "Global"],
                      label="Consensus",
                      default=1,
                      help="Criteria for consensus alignment. If local, consensus will be performed to the tilt-images "
                           "individually (different consensus applied individually) and if there is no consensus in "
                           "some image the series will be discarded. If global, consensus will be performed at the "
                           "level o tilt-series (the whole pair of alignments must agree).")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.consensusAlignment)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def consensusAlignment(self):
        tsIdList = self.generateTsIdList()

        for tsId in tsIdList:
            Mset = []
            SRset = []

            for sots in self.inputMultiSoTS:
                ts = sots.get().getTiltSeriesFromTsId(tsId)

                M = []

                for ti in ts:
                    M.append(ti.getTransform().getMatrix())

                Mset.append(M)
                SRset.append(ts.getSamplingRate())

            self.info("\nAnalyzing tilt series " + tsId + "...")

            shiftTolPx = round(self.shiftTolerance.get() / SRset[0])

            if self.globalAlignment.get() == 1:
                averageAlignmentV, angleSDV, shiftSDV = \
                    self.compareTransformationMatricesGlobal(Mset,
                                                             shiftTolPx,
                                                             self.angleTolerance.get(),
                                                             SRset)
            else:
                averageAlignmentV, angleSDV, shiftSDV = \
                    self.compareTransformationMatricesLocal(Mset,
                                                            shiftTolPx,
                                                            self.angleTolerance.get(),
                                                            SRset)

            # Consensus achieved
            if averageAlignmentV is not None:
                self.getOutputAliConsensusSoTS()

                ts = self.inputMultiSoTS[0].get().getTiltSeriesFromTsId(tsId)
                newTs = TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                self.outputAliConsensusSoTS.append(newTs)

                for i, ti in enumerate(ts):
                    newTi = TiltImage()
                    newTi.copyInfo(ti, copyId=True)
                    newTi.setLocation(ti.getLocation())

                    transform = Transform()
                    if averageAlignmentV[i] is None:
                        newTi.setEnabled(False)
                    else:
                        transform.setMatrix(averageAlignmentV[i])
                    newTi.setTransform(transform)

                    newTi._angleStd = Float(angleSDV[i])
                    newTi._shiftStd = Float(shiftSDV[i])

                    newTs.append(newTi)

                newTs.setDim(ts.getDim())
                newTs.write()

                self.outputAliConsensusSoTS.update(newTs)
                self.outputAliConsensusSoTS.updateDim()
                self.outputAliConsensusSoTS.write()

            # No consensus achieved (return tilt-series with no alignment information)
            else:
                self.getOutputNoAliConsensusSoTS()

                ts = self.inputMultiSoTS[0].get().getTiltSeriesFromTsId(tsId)
                newTs = TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                self.outputNoAliConsensusSoTS.append(newTs)

                for i, ti in enumerate(ts):
                    newTi = TiltImage()
                    newTi.copyInfo(ti, copyId=True)
                    newTi.setLocation(ti.getLocation())

                    newTs.append(newTi)

                newTs.setDim(ts.getDim())
                newTs.write()

                self.outputNoAliConsensusSoTS.update(newTs)
                self.outputNoAliConsensusSoTS.updateDim()
                self.outputNoAliConsensusSoTS.write()

        self._store()

    def closeOutputSetsStep(self):
        if hasattr(self, "outputAliConsensusSoTS"):
            self.outputAliConsensusSoTS.setStreamState(Set.STREAM_CLOSED)
            self.outputAliConsensusSoTS.write()

        if hasattr(self, "outputNoAliConsensusSoTS"):
            self.outputNoAliConsensusSoTS.setStreamState(Set.STREAM_CLOSED)
            self.outputNoAliConsensusSoTS.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def compareTransformationMatricesGlobal(Mset, shiftTol, angleTol, SRset):

        logger.info("Running global consensus alignment...")

        Nts = len(Mset)

        # If there is only one matrix in Mset then there has been no consensus in the recursion.
        if Nts < 2:
            logger.info("No consensus achieved for this tilt-series.")
            return None, None, None

        Nti = len(Mset[0])

        logger.info("Number of tilt-series analyzed: " + str(Nts))
        logger.info("Number of tilt-images per tilt-series analyzed: " + str(Nti))

        # Matrix reporting consensus between each pair of tilt series
        consensusAlignmentMatrix = np.zeros((Nts, Nts))

        # Compare each pair of matrices
        for j in range(Nts):
            for k in range(j + 1, Nts):
                # Calculate the sampling factor between the two matrices (for comparing shifts)
                samplingFactor = SRset[k] / SRset[j]

                # Calculate p matrix for the pair of tilt-series
                p = np.zeros((3, 3))

                for i in range(Nti):  # Iterate each tilt-image
                    sampledMatrix = Mset[k][i].copy()
                    sampledMatrix[0, 2] *= samplingFactor
                    sampledMatrix[1, 2] *= samplingFactor

                    p += np.matmul(Mset[j][i], np.linalg.inv(sampledMatrix))

                # Calculate error matrix given a calculated p matrix (for a pair of matrices)
                p /= Nti  # Normalized by the number of comparisons performed

                detectedMisali = False

                # Calculate error matrix for the pair of tilt-series.
                for i in range(Nti):  # Iterate each tilt-image
                    # Only use p matrix to correct for shiftY in case there exist an Y offset in the whole series
                    matrixShiftYCorrected = Mset[k][i].copy()
                    matrixShiftYCorrected[0, 2] *= samplingFactor
                    matrixShiftYCorrected[1, 2] *= samplingFactor
                    matrixShiftYCorrected[1, 2] += p[1][2]

                    pError = Mset[j][i] - matrixShiftYCorrected

                    # Angle from pTotalError matrix
                    cosRotationAngle = pError[0][0]
                    sinRotationAngle = pError[1][0]

                    if math.isnan(sinRotationAngle / cosRotationAngle):
                        angleError = 0
                    else:
                        angleError = math.degrees(math.asin(sinRotationAngle))

                    # Shifts from pError matrix
                    shiftX = pError[0][2]
                    shiftY = pError[1][2]

                    if abs(shiftX) > shiftTol or abs(shiftY) > shiftTol or abs(angleError) > angleTol:
                        detectedMisali = True

                if detectedMisali:
                    logger.info(
                        "No consensus achieved between tilt-series " + str(j) + " and tilt-series " + str(k))
                    consensusAlignmentMatrix[j][k] = 1

        colOnes = []
        rowOnes = []  # This should be a star wars Easter Egg

        indexFailed = np.sum(consensusAlignmentMatrix, axis=0)

        for i in range(Nts):
            for j in range(Nts):
                if consensusAlignmentMatrix[i][j] == 1:
                    rowOnes.append(i)
                    colOnes.append(j)

        if indexFailed[np.argmax(indexFailed)] != 0:
            discardedIndexes = [np.argmax(indexFailed)]
        else:
            discardedIndexes = []

        discardedIndexes.sort(reverse=True)

        if len(discardedIndexes) != 0:
            for n in discardedIndexes:
                del Mset[n]
            return ProtConsensusAlignmentTS.compareTransformationMatricesGlobal(Mset, shiftTol, angleTol, SRset)
        else:
            logger.info("\nConsensus achieved for this tilt-series.")
            averageAlignmentV = []
            angleSDV = []
            shiftSDV = []

            for i in range(Nti):
                averageMatrix = np.zeros((3, 3))
                angleV = []
                shiftV = []

                for j in range(Nts):
                    # Calculate the sampling factor between the two matrices (for comparing shifts)
                    samplingFactor = SRset[j] / SRset[0]

                    sampledMatrix = Mset[j][i].copy()
                    sampledMatrix[0, 2] *= samplingFactor
                    sampledMatrix[1, 2] *= samplingFactor

                    m = sampledMatrix.copy()

                    # Calculate average alignment
                    averageMatrix += m

                    # Angle from average alignment matrix
                    sinRotationAngle = m[1][0]
                    angle = math.degrees(math.asin(sinRotationAngle))

                    # Shifts from average alignment matrix
                    shiftX = m[0][2]
                    shiftY = m[1][2]

                    angleV.append(angle)
                    shiftV.append((shiftX + shiftY) / 2)

                angleSDV.append(np.std(angleV))
                shiftSDV.append(np.std(shiftV))

                averageMatrix /= Nts
                averageAlignmentV.append(averageMatrix)

            return averageAlignmentV, angleSDV, shiftSDV

    @staticmethod
    def compareTransformationMatricesLocal(Mset, shiftTol, angleTol, SRset):

        logger.info("Running local consensus alignment...")

        Nts = len(Mset)
        Nti = len(Mset[0])

        logger.info("Number of tilt-series analyzed: " + str(Nts))
        logger.info("Number of tilt-images per tilt-series analyzed: " + str(Nti))

        # 2D-array saving all p matrices for each pair of comparisons. Keys follow #tsX_#tsY being p the matrix that
        # spams X set onto Y set
        P_dict = {}

        # Array saving the final consensus alignment
        averageAlignmentV = []
        angleSDV = []
        shiftSDV = []

        # Calculate p matrix for each
        for j in range(Nts):
            for k in range(j + 1, Nts):
                # Calculate the sampling factor between the two matrices (for comparing shifts)
                samplingFactor = SRset[k] / SRset[j]

                # Calculate p matrix for the pair of tilt-series
                p = np.zeros((3, 3))

                for i in range(Nti):  # Iterate each tilt-image
                    sampledMatrix = Mset[k][i].copy()
                    sampledMatrix[0, 2] *= samplingFactor
                    sampledMatrix[1, 2] *= samplingFactor

                    p += np.matmul(Mset[j][i], np.linalg.inv(sampledMatrix))

                # Calculate error matrix given a calculated p matrix (for a pair of matrices)
                p /= Nti  # Normalized by the number of comparisons performed

                key = "%d_%d" % (j, k)
                P_dict.update({key: p})

        # Apply consensus over each tilt-image
        for i in range(Nti):  # Iterate each tilt-image
            logger.info("Analyzing image " + str(i))

            tiConsensusAlignment = np.zeros((3, 3))

            angleV = []
            shiftV = []

            # Save all indexes that consensus for this tilt-image
            consensusIndexes = []

            for j in range(Nts):
                for k in range(j + 1, Nts):
                    logger.info("Comparing matrices from series " + str(j) + " and " + str(k))
                    # Matrix p for each pair of sets
                    key = "%d_%d" % (j, k)
                    p = P_dict[key]

                    # Calculate the sampling factor between the two matrices (for comparing shifts)
                    samplingFactor = SRset[k] / SRset[j]

                    # Only use p matrix to correct for shiftY in case there exist an Y offset in the whole series
                    matrixShiftYCorrected = Mset[k][i].copy()
                    matrixShiftYCorrected[0, 2] *= samplingFactor
                    matrixShiftYCorrected[1, 2] *= samplingFactor
                    matrixShiftYCorrected[1, 2] += p[1][2]

                    pError = Mset[j][i] - matrixShiftYCorrected

                    # Angle from pTotalError matrix
                    cosRotationAngle = pError[0][0]
                    sinRotationAngle = pError[1][0]

                    if math.isnan(sinRotationAngle / cosRotationAngle):
                        angleError = 0
                    else:
                        angleError = math.degrees(math.asin(sinRotationAngle))

                    # Shifts from pError matrix
                    shiftX = pError[0][2]
                    shiftY = pError[1][2]

                    if abs(shiftX) < shiftTol and abs(shiftY) < shiftTol and abs(angleError) < angleTol:
                        logger.info("Consensus achieved for " + str(j) + " and " + str(k))
                        consensusIndexes.append(j)
                        consensusIndexes.append(k)

            # Remove repeated indexes
            consensusIndexes = list(dict.fromkeys(consensusIndexes))

            if not consensusIndexes:
                logger.info("No consensus achieved for tilt-image " + str(i))
                averageAlignmentV.append(None)
                angleSDV.append(None)
                shiftSDV.append(None)

            else:
                logger.info("Consensus achieved for tilt-image " + str(i))

                tiConsensusAlignment += Mset[0][i]
                numberOfConsensusAlignment = len(consensusIndexes)

                # Angle from alignment matrix
                if math.isnan(Mset[0][i][1][0] / Mset[0][i][0][0]):
                    angleV.append(0)
                else:
                    angleV.append(math.degrees(math.asin(Mset[0][i][1][0])))

                # Shifts from average alignment matrix
                shiftX = Mset[0][i][0][2]
                shiftY = Mset[0][i][1][2]

                shiftV.append((shiftX + shiftY) / 2)

                for m in range(1, numberOfConsensusAlignment):
                    # p matrix for this pair of matrices
                    key = "%d_%d" % (consensusIndexes[0], consensusIndexes[m])
                    p = P_dict[key]

                    # Calculate the sampling factor between the two matrices (for comparing shifts)
                    samplingFactor = SRset[consensusIndexes[0]] / SRset[consensusIndexes[m]]

                    # Only use p matrix to correct for shiftY in case there exist an Y offset in the whole series
                    matrixShiftYCorrected = Mset[consensusIndexes[m]][i].copy()
                    matrixShiftYCorrected[0, 2] *= samplingFactor
                    matrixShiftYCorrected[1, 2] *= samplingFactor
                    matrixShiftYCorrected[1, 2] += p[1][2]

                    tiConsensusAlignment += matrixShiftYCorrected

                    # Angle from alignment matrix
                    if math.isnan(matrixShiftYCorrected[1][0] / matrixShiftYCorrected[0][0]):
                        angleV.append(0)
                    else:
                        angleV.append(math.degrees(math.asin(matrixShiftYCorrected[1][0])))

                    # Shifts from average alignment matrix
                    shiftX = matrixShiftYCorrected[0][2]
                    shiftY = matrixShiftYCorrected[1][2]

                    shiftV.append((shiftX + shiftY) / 2)

                # Calculate STD angles
                angleSDV.append(np.std(angleV))

                # Calculate STD shifts
                shiftSDV.append(np.std(shiftV))

                tiConsensusAlignment /= numberOfConsensusAlignment
                averageAlignmentV.append(tiConsensusAlignment)

        # Check if all images are misaligned
        for a in averageAlignmentV:
            if a is not None:
                return averageAlignmentV, angleSDV, shiftSDV

        return None, None, None

    def generateTsIdList(self):
        resultSet = set(self.inputMultiSoTS[0].get().getUniqueValues('_tsId'))

        for i in range(1, len(self.inputMultiSoTS)):
            resultSet.intersection_update(set(self.inputMultiSoTS[i].get().getUniqueValues('_tsId')))

        return list(resultSet)

    def getOutputAliConsensusSoTS(self):
        if not hasattr(self, "outputAliConsensusSoTS"):
            outputAliConsensusSoTS = self._createSetOfTiltSeries(suffix='AliConsensus')
            outputAliConsensusSoTS.copyInfo(self.inputMultiSoTS[0].get())
            outputAliConsensusSoTS.setDim(self.inputMultiSoTS[0].get().getDim())

            self._defineOutputs(outputAliConsensusSoTS=outputAliConsensusSoTS)
            self._defineSourceRelation(self.inputMultiSoTS, outputAliConsensusSoTS)
        return self.outputAliConsensusSoTS

    def getOutputNoAliConsensusSoTS(self):
        if not hasattr(self, "outputNoAliConsensusSoTS"):
            outputNoAliConsensusSoTS = self._createSetOfTiltSeries(suffix='NoAliConsensus')
            outputNoAliConsensusSoTS.copyInfo(self.inputMultiSoTS[0].get())
            outputNoAliConsensusSoTS.setDim(self.inputMultiSoTS[0].get().getDim())

            self._defineOutputs(outputNoAliConsensusSoTS=outputNoAliConsensusSoTS)
            self._defineSourceRelation(self.inputMultiSoTS, outputNoAliConsensusSoTS)
        return self.outputNoAliConsensusSoTS

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = [] if len(self.inputMultiSoTS) > 1 else \
            ["More than one input set of tilt-series is needed to compute the consensus."]

        for i, sots in enumerate(self.inputMultiSoTS):
            ts = sots.get().getFirstItem()
            if not ts.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set of tilt-series %d does not have a "
                                    "transformation matrix assigned." % (i + 1))

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputAliConsensusSoTS') or hasattr(self, 'outputNoAliConsensusSoTS'):
            summary.append("Input Tilt-Series: %d."
                           % (len(self.inputMultiSoTS)))

        if hasattr(self, 'outputAliConsensusSoTS'):
            summary.append("Output tilt series with consensus achieved: %d."
                           % (self.outputAliConsensusSoTS.getSize()))

        if hasattr(self, 'outputNoAliConsensusSoTS'):
            summary.append("Output tilt-series with no consensus achieved: %d"
                           % (self.outputNoAliConsensusSoTS.getSize()))

        if not (hasattr(self, 'outputAliConsensusSoTS') or hasattr(self, 'outputNoAliConsensusSoTS')):
            summary.append("Output classes not ready yet.")

        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputAliConsensusSoTS') or hasattr(self, 'outputNoAliConsensusSoTS'):
            methods.append("Input Tilt-Series: %d."
                           % (len(self.inputMultiSoTS)))

        if hasattr(self, 'outputAliConsensusSoTS'):
            methods.append("Output tilt series with consensus achieved: %d."
                           % (self.outputAliConsensusSoTS.getSize()))

        if hasattr(self, 'outputNoAliConsensusSoTS'):
            methods.append("Output tilt-series with no consensus achieved: %d"
                           % (self.outputNoAliConsensusSoTS.getSize()))

        if not (hasattr(self, 'outputAliConsensusSoTS') or hasattr(self, 'outputNoAliConsensusSoTS')):
            methods.append("Output classes not ready yet.")

        return methods
