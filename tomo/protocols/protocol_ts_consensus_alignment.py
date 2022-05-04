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
from pyworkflow.protocol.params import MultiPointerParam, PointerParam
from pyworkflow.object import Set

from pwem.protocols import EMProtocol
from pwem.objects import Float

from tomo.objects import TiltSeries, TiltImage
from tomo.protocols import ProtTomoBase
from tomo import utils


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

            # enable, anglesAvgV, anglesStdV, shiftXAvgV, shiftXStdV, shiftYAvgV, shiftYStdV = self.compareTransformationMatrices(Mset)
            enable, anglesStdV, shiftXStdV, shiftYStdV = self.compareTransformationMatrices(Mset)

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

                # newTi._angleAvg = Float(anglesAvgV[i])
                newTi._angleStd = Float(anglesStdV[i])
                # newTi._shiftXAvg = Float(shiftXAvgV[i])
                newTi._shiftXStd = Float(shiftXStdV[i])
                # newTi._shiftYAvg = Float(shiftYAvgV[i])
                newTi._shiftYStd = Float(shiftYStdV[i])

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

    @staticmethod
    def compareTransformationMatrices(Mset):
        Nts = len(Mset)
        Nti = len(Mset[0])

        # Number od possible combinations picking 2 tilt-series out of the set
        numberOfCombinations = math.factorial(Nts)/(math.factorial(2)*math.factorial(Nts-2))

        print("Number of tilt-series analyzed: " + str(Nts))
        print("Number of tilt-images per tilt-series analyzed: " + str(Nti))

        p = np.zeros((3, 3))

        anglesAvgV = []
        anglesStdV = []
        shiftXAvgV = []
        shiftXStdV = []
        shiftYAvgV = []
        shiftYStdV = []

        for i in range(Nti):  # Iterate each tilt-image
            sumAngle = 0
            sum2Angle = 0
            sumShiftX = 0
            sum2ShiftX = 0
            sumShiftY = 0
            sum2ShiftY = 0

            for j in range(Nts):  # Iterate each tilt-series
                # Angle from TM
                cosRotationAngle = Mset[j][i][0][0]
                sinRotationAngle = Mset[j][i][1][0]
                angle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

                # Shifts from TM
                shifts = [Mset[j][i][0][2], Mset[j][i][0][2]]

                sumAngle += angle
                sum2Angle += angle*angle
                sumShiftX += shifts[0]
                sum2ShiftX += shifts[0] * shifts[0]
                sumShiftY += shifts[1]
                sum2ShiftY += shifts[1] * shifts[1]

                for k in range(j+1, Nts):  # Compare each matrix with the following
                    p += np.matmul(Mset[j][i], np.linalg.inv(Mset[k][i]))

            # Append angle information
            anglesAvg = sumAngle / Nts
            anglesStd = math.sqrt(sum2Angle/Nts - anglesAvg*anglesAvg)
            anglesAvgV.append(anglesAvg)
            anglesStdV.append(anglesStd)

            # Append shift X information
            shiftXAvg = sumShiftX / Nts
            shiftXStd = math.sqrt(sum2ShiftX/Nts - shiftXAvg*shiftXAvg)
            shiftXAvgV.append(shiftXAvg)
            shiftXStdV.append(shiftXStd)

            # Append shift Y information
            shiftYAvg = sumShiftY / Nts
            shiftYStd = math.sqrt(sum2ShiftY/Nts - shiftYAvg*shiftYAvg)
            shiftYAvgV.append(shiftYAvg)
            shiftYStdV.append(shiftYStd)

        p /= Nti * numberOfCombinations  # Normalized by the number of comparisons performed
        print("p")
        print(np.matrix.round(p, 2))

        pError = np.zeros((3, 3))

        for i in range(Nti):  # Iterate each tilt-image
            for j in range(Nts):  # Iterate each tilt-series
                for k in range(j + 1, Nts):  # Compare each matrix with the following
                    pError += np.absolute(p - np.matmul(Mset[j][i], np.linalg.inv(Mset[k][i])))

        pError /= numberOfCombinations  # Normalized by the number of tilt-series analyzed

        print("pError")
        print(np.matrix.round(pError, 2))

        print("+++++++++++++++++++++++++++++++++++")

        # *** add condition based on pError
        enable = True

        # return enable, anglesAvgV, anglesStdV, shiftXAvgV, shiftXStdV, shiftYAvgV, shiftYStdV
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
                    tmpTsIdList.append(ts.getTsId)
                if tsId not in tmpTsIdList:
                    tsIdList.remove(tsId)
                    break

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
                                    "transformation matrix assigned." % (i+1))

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
