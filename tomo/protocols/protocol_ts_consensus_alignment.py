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
import math

import numpy as np

from pyworkflow import BETA
from pyworkflow.protocol.params import MultiPointerParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.object import Set

from pwem.protocols import EMProtocol

from tomo.objects import TiltSeries, TiltImage
from tomo.protocols import ProtTomoBase


class ProtConsensusAlignmentTS(EMProtocol, ProtTomoBase):
    """
    Perform a consensus of two alignments for the same tilt series.
    """

    _label = 'Tilt-series consensus alignment'
    _devStatus = BETA

    _tsIdList = []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        # form.addParam('setOfTiltSeries1',
        #               PointerParam,
        #               pointerClass='SetOfTiltSeries',
        #               important=True,
        #               help='First set of tilt-series to be analyzed in the consensus alignment. The unrelated '
        #                    'alignment information will be taken form this set.',
        #               label='First set of tilt-series')
        #
        # form.addParam('setOfTiltSeries2',
        #               PointerParam,
        #               pointerClass='SetOfTiltSeries',
        #               important=True,
        #               help='Second set of tilt-series to be analyzed in the consensus alignment.',
        #               label='Second set of tilt-series')

        form.addParam('inputMultiSoTS',
                      MultiPointerParam,
                      important=True,
                      label="Input tilt series",
                      pointerClass='SetOfTiltSeries',
                      help='Select several sets of tilt-series where to evaluate the consensus in their alignment. '
                           'Output set will bring the information from the first selected set.')

        form.addParam('shiftTolerance',
                      FloatParam,
                      label="Shift tolerance (px)",
                      expertLevel=LEVEL_ADVANCED,
                      help='Maximum shift difference between alignments to consider them as equal. it is measured in '
                           'pixels.')

        form.addParam('angleTolerance',
                      FloatParam,
                      label="Angle tolerance (degrees)",
                      expertLevel=LEVEL_ADVANCED,
                      help='Maximum angle difference between alignments to consider them as equal. It is measured in '
                           'degrees')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.consensusAlignment)

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def consensusAlignment(self):
        tsIdList = self.generateTsIdList()

        self.getOutputAlignmentConsensusSetOfTiltSeries()

        for tsId in tsIdList:
            Mset = []

            for sots in self.inputMultiSoTS:
                ts = sots.get().getTiltSeriesFromTsId(tsId)

                M = []

                for ti in ts:
                    M.append(ti.getTransform().getMatrix())

                Mset.append(M)

            # anglesAvgV, anglesStdV, shiftXAvgV, shiftXStdV, shiftYAvgV, shiftYStdV = self.calculateAngleAndShiftDistribution(Mset)
            # anglesStdV, shiftXStdV, shiftYStdV = self.calculateAngleAndShiftDistribution(Mset)

            # enable, anglesStdV, shiftXStdV, shiftYStdV = self.compareTransformationMatrices(Mset)
            enable = self.compareTransformationMatrices(Mset)

            ts = self.inputMultiSoTS[0].get().getTiltSeriesFromTsId(tsId)
            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            newTs._enable = enable
            self.outputAlignmentConsensusSetOfTiltSeries.append(newTs)

            for i, ti in enumerate(ts):
                newTi = TiltImage()
                newTi.copyInfo(ti, copyId=True)
                newTi.setLocation(ti.getLocation())
                newTi.setTransform(ti.getTransform())

                # # newTi._angleAvg = Float(anglesAvgV[i])
                # newTi._angleStd = Float(anglesStdV[i])
                # # newTi._shiftXAvg = Float(shiftXAvgV[i])
                # newTi._shiftXStd = Float(shiftXStdV[i])
                # # newTi._shiftYAvg = Float(shiftYAvgV[i])
                # newTi._shiftYStd = Float(shiftYStdV[i])

                newTs.append(newTi)

            newTs.setDim(ts.getDim())
            newTs.write()

            self.outputAlignmentConsensusSetOfTiltSeries.update(newTs)
            self.outputAlignmentConsensusSetOfTiltSeries.updateDim()
            self.outputAlignmentConsensusSetOfTiltSeries.write()

        self._store()

    # def consensusAlignment(self):
    #     tsIdList = self.generateTsIdList()
    #
    #     self.getOutputAlignmentConsensusSetOfTiltSeries()
    #
    #     for tsId in tsIdList:
    #         ts1 = self.setOfTiltSeries1.get().getTiltSeriesFromTsId(tsId)
    #         ts2 = self.setOfTiltSeries2.get().getTiltSeriesFromTsId(tsId)
    #
    #         newTs = TiltSeries(tsId=tsId)
    #         newTs.copyInfo(ts1)
    #         self.outputAlignmentConsensusSetOfTiltSeries.append(newTs)
    #
    #         # First method
    #         M1set = []
    #         M2set = []
    #
    #         # Second method
    #         Mset = []
    #
    #         for ti1, ti2 in zip(ts1, ts2):
    #             newTi = TiltImage()
    #             newTi.copyInfo(ti1, copyId=True)
    #             newTi.setLocation(ti1.getLocation())
    #             newTi.setTransform(ti1.getTransform())
    #
    #             ra1, sh1 = utils.getRotationAngleAndShiftFromTM(ti1)
    #             ra2, sh2 = utils.getRotationAngleAndShiftFromTM(ti2)
    #
    #             newTi._angleDiff = Float(abs(ra1 - ra2))
    #             newTi._shiftDiff = Float(
    #                 (abs(sh1[0] - sh2[0]) + abs(sh1[1] - sh2[1])) / 2)  # Shift average in both directions
    #
    #             newTs.append(newTi)
    #
    #             # First method
    #             M1set.append(ti1.getTransform().getMatrix())
    #             M2set.append(ti2.getTransform().getMatrix())
    #
    #         # Second method
    #         Mset.append(M1set)
    #         Mset.append(M2set)
    #
    #         # First method
    #         self.transformationMatrix(M1set, M2set)
    #
    #         # Second method
    #         self.compareTransformationMatrices(Mset)
    #
    #         newTs.setDim(ts1.getDim())
    #         newTs.write()
    #
    #         self.outputAlignmentConsensusSetOfTiltSeries.update(newTs)
    #         self.outputAlignmentConsensusSetOfTiltSeries.updateDim()
    #         self.outputAlignmentConsensusSetOfTiltSeries.write()
    #
    #     self._store()

    def closeOutputSetsStep(self):
        self.outputAlignmentConsensusSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.outputAlignmentConsensusSetOfTiltSeries.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def compareTransformationMatrices(Mset):
        Nts = len(Mset)
        Nti = len(Mset[0])

        # If there is only one matrix in Mset then there has been no consensus in the recursion.
        if Nts < 2:
            print("No consensus achieved for this set of tilt-series.")
            return False

        # Number of possible combinations picking 2 tilt-series out of the set
        numberOfCombinations = math.factorial(Nts) / (math.factorial(2) * math.factorial(Nts - 2))

        print("Number of tilt-series analyzed: " + str(Nts))
        print("Number of tilt-images per tilt-series analyzed: " + str(Nti))

        # Matrix for saving the final accumulated error in the whole set
        pTotalError = np.zeros((3, 3))

        # Matrix saving the module error from comparing matrices
        # moduleErrorMatrix = np.zeros((Nts, Nts))

        # Vectors for angle and shift average errors between tilt-series
        avgAngleErrorsV = np.zeros(Nts)
        avgShiftErrorsV = np.zeros(Nts)

        # Calculate p matrix
        for k in range(1, Nts):  # Compare each matrix with the following

            print("\nComparing matrices " + str(0) + " and " + str(k))

            p = np.zeros((3, 3))

            # Calculate p matrix for the pair of tilt-series
            for i in range(Nti):  # Iterate each tilt-image
                p += np.matmul(Mset[0][i], np.linalg.inv(Mset[k][i]))

            # Calculate error matrix given a calculated p matrix (for a pair of matrices)
            p /= Nti  # Normalized by the number of comparisons performed
            print("p"+str(k))
            print(np.matrix.round(p, 2))

            pError = np.zeros((3, 3))

            # Calculate error matrix for the pair of tilt-series.
            for i in range(Nti):  # Iterate each tilt-image
                # pError += np.absolute(p - np.matmul(Mset[0][i], np.linalg.inv(Mset[k][i])))

                # Only use p matrix to correct for shiftY in case there exist an offset in the whole series
                matrixShiftYCorrected = Mset[k][i]
                matrixShiftYCorrected[1, 2] = matrixShiftYCorrected[1][2] + p[1][2]
                pError += np.absolute(Mset[0][i] - matrixShiftYCorrected)

            pError /= Nti  # Normalized by the number of tilt-images analyzed
            print("pError"+str(k))
            print(np.matrix.round(pError, 2))

            # Angle from pError matrix
            cosRotationAngle = pError[0][0]
            sinRotationAngle = pError[1][0]

            if math.isnan(sinRotationAngle / cosRotationAngle):
                avgAngleErrorsV[k] = 0
            else:
                avgAngleErrorsV[k] = math.degrees(math.asin(sinRotationAngle))

            # Shifts from pError matrix
            shiftX = pError[0][2]
            shiftY = pError[1][2]
            avgShiftErrorsV[k] = (shiftX + shiftY) / 2

            # Add the matrix error for a given pair of matrices to the total error matrix
            pTotalError += pError

        # Normalize the total error matrix
        pTotalError /= (Nts - 1)

        print("pTotalError")
        print(np.matrix.round(pTotalError, 2))

        # Angle from pTotalError matrix
        cosRotationAngle = pTotalError[0][0]
        sinRotationAngle = pTotalError[1][0]

        if math.isnan(sinRotationAngle / cosRotationAngle):
            avgAngleError = 0
        else:
            avgAngleError = math.degrees(math.asin(sinRotationAngle))

        # Shifts from pTotalError matrix
        shiftX = pTotalError[0][2]
        shiftY = pTotalError[1][2]
        avgShiftError = (shiftX + shiftY) / 2

        print("\nMean shift error vector: ")
        print(avgShiftErrorsV)
        print("Mean angle error vector: ")
        print(avgAngleErrorsV)

        print("Mean angle error: " + str(avgAngleError))
        print("Mean shift error: " + str(avgShiftError))

        # Discard series by MAD
        shiftMedian = np.median(avgShiftErrorsV)
        angleMedian = np.median(avgAngleErrorsV)

        shiftMAD = 0
        angleMAD = 0

        for i in range(Nts):
            shiftMAD += abs(avgShiftErrorsV[i] - shiftMedian)
            angleMAD += abs(avgAngleErrorsV[i] - angleMedian)

        shiftMAD /= Nts
        angleMAD /= Nts

        print("Angle MAD: " + str(angleMAD))
        print("Shift MAD: " + str(shiftMAD))

        discardedIndexes = []

        for i in range(Nts):
            if (abs(avgShiftErrorsV[i] - shiftMedian) > 1 * shiftMAD or
                    abs(avgAngleErrorsV[i] - angleMedian) > 1 * angleMAD):
                discardedIndexes.append(i)

        print("Discarded indexes")
        print(discardedIndexes)

        if len(discardedIndexes) != 0:
            for n in discardedIndexes:
                del Mset[n]
            ProtConsensusAlignmentTS.compareTransformationMatrices(Mset)
        else:
            return True

        return True  # enable

    @staticmethod
    def compareTransformationMatricesBis(Mset):
        Nts = len(Mset)
        Nti = len(Mset[0])

        # Number of possible combinations picking 2 tilt-series out of the set
        numberOfCombinations = math.factorial(Nts) / (math.factorial(2) * math.factorial(Nts - 2))

        print("Number of tilt-series analyzed: " + str(Nts))
        print("Number of tilt-images per tilt-series analyzed: " + str(Nti))

        # Vector for posterior calculation of std and avg of both shifts and angles.
        sumAnglesV = np.zeros(Nti)
        sum2AnglesV = np.zeros(Nti)
        sumShiftXV = np.zeros(Nti)
        sum2ShiftXV = np.zeros(Nti)
        sumShiftYV = np.zeros(Nti)
        sum2ShiftYV = np.zeros(Nti)

        # Calculate shifts and vectors values of the first matrix (since it is the only one not compared in the
        # posterior loop). Also, it is the only one whose shiftY does not have to be corrected by the matrix p (all the
        # matrices in set are compared to the first one).
        for i in range(Nti):
            # Angle from TM
            cosRotationAngle = Mset[0][i][0][0]
            sinRotationAngle = Mset[0][i][1][0]
            angle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

            # Shifts from TM
            shiftX = Mset[0][i][0][2]
            shiftY = Mset[0][i][1][2]

            sumAnglesV[i] += angle
            sum2AnglesV[i] += angle * angle
            sumShiftXV[i] += shiftX
            sum2ShiftXV[i] += shiftX * shiftX
            sumShiftYV[i] += shiftY
            sum2ShiftYV[i] += shiftY * shiftY

        # Matrix for saving the final accumulated error in the whole set
        pTotalError = np.zeros((3, 3))

        for j in range(1):  # Iterate each tilt-series

            # Calculate p matrix
            for k in range(j + 1, Nts):  # Compare each matrix with the following

                p = np.zeros((3, 3))

                # Calculate p matrix for the pair of tilt-series
                for i in range(Nti):  # Iterate each tilt-image
                    p += np.matmul(Mset[j][i], np.linalg.inv(Mset[k][i]))

                    print(Mset[j][i])
                    print(Mset[k][i])
                    print("-----------------------------------")

                # Calculate error matrix given a calculated p matrix (for a pair of matrices)
                p /= Nti  # Normalized by the number of comparisons performed
                print("p")
                print(np.matrix.round(p, 2))

                pError = np.zeros((3, 3))

                # Calculate shifts and vectors values
                for i in range(Nti):  # Iterate each tilt-image
                    # Angle from TM
                    cosRotationAngle = Mset[k][i][0][0]
                    sinRotationAngle = Mset[k][i][1][0]
                    angle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

                    # Shifts from TM
                    shiftX = Mset[k][i][0][2]
                    shiftY = Mset[k][i][1][2] + p[1][2]  # Correct the shiftY from calculated p matrix

                    sumAnglesV[i] += angle
                    sum2AnglesV[i] += angle * angle
                    sumShiftXV[i] += shiftX
                    sum2ShiftXV[i] += shiftX * shiftX
                    sumShiftYV[i] += shiftY
                    sum2ShiftYV[i] += shiftY * shiftY

                print("\nComparing matrices " + str(j) + " and " + str(k))

                # Calculate error matrix for the pair of tilt-series.
                for i in range(Nti):  # Iterate each tilt-image
                    # pError += np.absolute(p - np.matmul(Mset[j][i], np.linalg.inv(Mset[k][i])))

                    # Only use p matrix to correct for shiftY in case there exist an offset in the whole series
                    matrixShiftYCorrected = Mset[k][i]
                    matrixShiftYCorrected = matrixShiftYCorrected[0][1] - p[0][1]
                    pError += np.absolute(Mset - matrixShiftYCorrected)

                pError /= Nti  # Normalized by the number of tilt-series analyzed
                print("pError")
                print(np.matrix.round(pError, 2))

                # Add the matrix error for a given pair of matrices to the total error matrix
                pTotalError += pError

        # Calculate the final avg and stf for each tilt-image
        anglesAvgV = np.zeros(Nti)
        anglesStdV = np.zeros(Nti)
        shiftXAvgV = np.zeros(Nti)
        shiftXStdV = np.zeros(Nti)
        shiftYAvgV = np.zeros(Nti)
        shiftYStdV = np.zeros(Nti)

        for i in range(Nti):
            anglesAvgV[i] = sumAnglesV[i] / (
                    numberOfCombinations + 1)  # Plus one because we added the first tilt-series
            shiftXAvgV[i] = sumShiftXV[i] / (numberOfCombinations + 1)
            shiftYAvgV[i] = sumShiftYV[i] / (numberOfCombinations + 1)

            # Abs to avoid 0 as small negative number
            anglesStdV[i] = math.sqrt(abs(sum2AnglesV[i] / (numberOfCombinations + 1) - anglesAvgV[i] * anglesAvgV[i]))
            shiftXStdV[i] = math.sqrt(abs(sum2ShiftXV[i] / (numberOfCombinations + 1) - shiftXAvgV[i] * shiftXAvgV[i]))
            shiftYStdV[i] = math.sqrt(abs(sum2ShiftYV[i] / (numberOfCombinations + 1) - shiftYAvgV[i] * shiftYAvgV[i]))

        print("anglesAvgV")
        print(anglesAvgV)
        print("anglesStdV")
        print(anglesStdV)
        print("shiftXAvgV")
        print(shiftXAvgV)
        print("shiftXStdV")
        print(shiftXStdV)
        print("shiftYAvgV")
        print(shiftYAvgV)
        print("shiftYStdV")
        print(shiftYStdV)

        # Normalize the total error matrix
        pTotalError /= numberOfCombinations

        print("\n\npTotalError")
        print(np.matrix.round(pTotalError, 2))

        # Calculate the determinant of the adjugate matrix (transpose of the cofactor matrix)
        m = np.linalg.det(pTotalError)
        c = [[i for i in range(3)] for j in range(3)]

        for i in range(3):
            for j in range(3):
                c[i][j] = (-1) * (i + j) * m

        print("Determinant of the adjugate matrix")
        print(np.linalg.det(c))
        print("+++++++++++++++++++++++++++++++++++")

        # *** add condition based on pError
        enable = True

        #     # Calculate the module of the pError
        #     pError[2][2] = 1
        #
        #     m = np.linalg.det(p)
        #     c = [[i for i in range(3)] for j in range(3)]
        #
        #     print("m")
        #     print(m)
        #
        #     for ii in range(3):
        #         for jj in range(3):
        #             c[ii][jj] = (-1) * (ii + jj) * m
        #
        #     print("c")
        #     print(np.matrix(c))
        #
        #     modError = np.linalg.det(np.matrix(c))
        #     print("modError")
        #     print(modError)
        #     moduleErrorMatrix[1][k] = modError
        #
        # print(moduleErrorMatrix)

        return enable, anglesStdV, shiftXStdV, shiftYStdV

    # def generateTsIdList(self):
    #     tsIdList = []
    #     tmpTsIdList = []
    #
    #     for ts in self.setOfTiltSeries1.get():
    #         tsIdList.append(ts.getTsId())
    #
    #     for ts in self.setOfTiltSeries2.get():
    #         tmpTsIdList.append(ts.getTsId())
    #
    #     for tsId in tsIdList:
    #         if tsId not in tmpTsIdList:
    #             tsIdList.remove(tsId)
    #
    #     if len(tsIdList) == 0:
    #         raise Exception("None matching tilt-series between two sets.")
    #
    #     return tsIdList

    def generateTsIdList(self):
        tsIdList = []
        tmpTsIdList = []

        for ts in self.inputMultiSoTS[0].get():
            tsIdList.append(ts.getTsId())

        for tsId in tsIdList:
            for i in range(1, len(self.inputMultiSoTS)):
                for ts in self.inputMultiSoTS[i].get():
                    tmpTsIdList.append(ts.getTsId())
                if tsId not in tmpTsIdList:
                    tsIdList.remove(tsId)
                    break
                tmpTsIdList = []

        if len(tsIdList) == 0:
            raise Exception("None matching tilt-series between two sets.")

        return tsIdList

    def getOutputAlignmentConsensusSetOfTiltSeries(self):
        if not hasattr(self, "outputAlignmentConsensusSetOfTiltSeries"):
            outputAlignmentConsensusSetOfTiltSeries = self._createSetOfTiltSeries(suffix='AliConsensus')
            outputAlignmentConsensusSetOfTiltSeries.copyInfo(self.inputMultiSoTS[0].get())
            outputAlignmentConsensusSetOfTiltSeries.setDim(self.inputMultiSoTS[0].get().getDim())

            self._defineOutputs(outputAlignmentConsensusSetOfTiltSeries=outputAlignmentConsensusSetOfTiltSeries)
            self._defineSourceRelation(self.inputMultiSoTS, outputAlignmentConsensusSetOfTiltSeries)
        return self.outputAlignmentConsensusSetOfTiltSeries

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

    # def _validate(self):
    #     validateMsgs = []
    #
    #     for ts1, ts2 in zip(self.setOfTiltSeries1.get(), self.setOfTiltSeries2.get()):
    #         if not ts1.getFirstItem().hasTransform():
    #             validateMsgs.append("Some tilt-series from the input set of tilt-series 1 does not have a "
    #                                 "transformation matrix assigned.")
    #
    #         if not ts2.getFirstItem().hasTransform():
    #             validateMsgs.append("Some tilt-series from the input set of tilt-series 2 does not have a "
    #                                 "transformation matrix assigned.")
    #
    #         if ts1.getSize() != ts2.getSize():
    #             validateMsgs.append("Some tilt-series from the input set of tilt-series 1 and its target in the "
    #                                 "set of tilt-series 2 sizes' do not match. Every input tilt-series "
    #                                 "and its target must have the same number of tilt-images")
    #
    #     return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputAlignmentConsensusSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nOutput tilt series with consensus applied : %d.\n"
                           % (len(self.inputMultiSoTS),
                              self.outputAlignmentConsensusSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    # def _summary(self):
    #     summary = []
    #     if hasattr(self, 'outputAlignmentConsensusSetOfTiltSeries'):
    #         summary.append("Input Tilt-Series: %d.\nOutput tilt series with consensus applied : %d.\n"
    #                        % (self.setOfTiltSeries1.get().getSize(),
    #                           self.setOfTiltSeries2.get().getSize()))
    #     else:
    #         summary.append("Output classes not ready yet.")
    #     return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputAlignmentConsensusSetOfTiltSeries'):
            methods.append("Consensus have been performed over %d transformation matrices from the input sets of "
                           "tilt-series. New consensus set of tilt-series generated with %d elements.\n"
                           % (len(self.inputMultiSoTS),
                              self.outputAlignmentConsensusSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods

    # def _methods(self):
    #     methods = []
    #     if hasattr(self, 'outputAlignmentConsensusSetOfTiltSeries'):
    #         methods.append("Consensus have been performed over %d transformation matrices from the input sets of "
    #                        "tilt-series.\n"
    #                        % (self.outputAlignmentConsensusSetOfTiltSeries.getSize()))
    #     else:
    #         methods.append("Output classes not ready yet.")
    #     return methods

    # @staticmethod
    # def transformationMatrix(M1set, M2set):
    #     p = np.zeros((3, 3))
    #
    #     for m1, m2 in zip(M1set, M2set):
    #         m1 = np.matrix(m1)
    #         m2 = np.matrix(m2)
    #
    #         p += np.matmul(m1, np.linalg.inv(m2))
    #
    #     p /= len(M1set)
    #     print("p")
    #     print(np.matrix.round(p, 2))
    #
    #     mError = np.zeros((3, 3))
    #     for m1, m2 in zip(M1set, M2set):
    #         mError += np.absolute(p - np.matmul(m1, np.linalg.inv(m2)))
    #
    #     print("mError")
    #     print(np.matrix.round(mError, 2))
    #
    #     print("----------------")

    # @staticmethod
    # def calculateAngleAndShiftDistribution(Mset):
    #     Nts = len(Mset)
    #     Nti = len(Mset[0])
    #
    #     anglesAvgV = []
    #     anglesStdV = []
    #     shiftXAvgV = []
    #     shiftXStdV = []
    #     shiftYAvgV = []
    #     shiftYStdV = []
    #
    #     for i in range(Nti):  # Iterate each tilt-image
    #         sumAngle = 0
    #         sum2Angle = 0
    #         sumShiftX = 0
    #         sum2ShiftX = 0
    #         sumShiftY = 0
    #         sum2ShiftY = 0
    #
    #         for j in range(Nts):  # Iterate each tilt-series
    #             # Angle from TM
    #             cosRotationAngle = Mset[j][i][0][0]
    #             sinRotationAngle = Mset[j][i][1][0]
    #             angle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))
    #
    #             # Shifts from TM
    #             shifts = [Mset[j][i][0][2], Mset[j][i][0][2]]
    #
    #             sumAngle += angle
    #             sum2Angle += angle*angle
    #             sumShiftX += shifts[0]
    #             sum2ShiftX += shifts[0] * shifts[0]
    #             sumShiftY += shifts[1]
    #             sum2ShiftY += shifts[1] * shifts[1]
    #
    #         # Append angle information
    #         anglesAvg = sumAngle / Nts
    #         anglesStd = math.sqrt(abs(sum2Angle/Nts - anglesAvg*anglesAvg))  # Abs avoid 0 as small negative number
    #         anglesAvgV.append(anglesAvg)
    #         anglesStdV.append(anglesStd)
    #
    #         # Append shift X information
    #         shiftXAvg = sumShiftX / Nts
    #         shiftXStd = math.sqrt(abs(sum2ShiftX/Nts - shiftXAvg*shiftXAvg))  # Abs avoid 0 as small negative number
    #         shiftXAvgV.append(shiftXAvg)
    #         shiftXStdV.append(shiftXStd)
    #
    #         # Append shift Y information
    #         shiftYAvg = sumShiftY / Nts
    #         shiftYStd = math.sqrt(abs(sum2ShiftY/Nts - shiftYAvg*shiftYAvg))  # Abs avoid 0 as small negative number
    #         shiftYAvgV.append(shiftYAvg)
    #         shiftYStdV.append(shiftYStd)
    #
    #     # return anglesAvgV, anglesStdV, shiftXAvgV, shiftXStdV, shiftYAvgV, shiftYStdV
    #     return anglesStdV, shiftXStdV, shiftYStdV
