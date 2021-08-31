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


import numpy as np
from scipy.spatial import ConvexHull

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import path
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase


class ProtResidualMisalignmentDetection(EMProtocol, ProtTomoBase):
    """
    Detection of misalignment in a set of tilt-series according to the obtained fiducial residuals.
    """

    _label = 'Residual misalignment detection'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfLandmarkModels',
                      params.PointerParam,
                      pointerClass='SetOfLandmarkModels',
                      important=True,
                      label='Input set of fiducial models containing residual information.')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for lm in self.inputSetOfLandmarkModels.get():
            lmObjId = lm.getObjId()
            self._insertFunctionStep(self.generateResidualList, lmObjId)

        # self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ----------------------------
    def generateResidualList(self, lmObjId):
        """ This method generates a residual list from each landmark model input set"""
        lm = self.inputSetOfLandmarkModels.get()[lmObjId]

        sr = lm.getTiltSeries().getSamplingRate()

        tsId = lm.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        listOfLMChainsMatrix = lm.retrieveLandmarkModelChains()

        listOfLMChainsMatrixNorm = self.normalizeResiduals(listOfLMChainsMatrix, sr)

        for listOfLMChainsNorm in listOfLMChainsMatrixNorm:

            listOfLMResid = []

            for lm in listOfLMChainsNorm:
                listOfLMResid.append([lm[4], lm[5]])

            chArea, chPerimeter = self.getResidualStatistics(listOfLMResid)

            maxDistance = self.getMaximumDistance(listOfLMResid)

            print(chArea)
            print(chPerimeter)
            print(maxDistance)
            print("-------------------------")

    def getResidualStatistics(self, listOfLMResid):
        """ This method calculate the statistics from each chain of residual to filter the properly aligned set of
        tilt series. """

        convexHullArea, convexHullPerimeter = self.getConvexHullAreaAndPerimeter(listOfLMResid)

        return convexHullArea, convexHullPerimeter

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def normalizeResiduals(listOfLMChainsMatrix, sr):
        """ This method multiplies the residual values by the sampling rate of the tilt-series to normalize its values
        to a common scale.
        :param listOfLMChainsMatrix: input set of landmark models separated by chainID.
        :param sr: sampling rate.
        """

        for listOfLMChains in listOfLMChainsMatrix:
            for lm in listOfLMChains:
                lm[4] = float(lm[4]) * sr
                lm[5] = float(lm[5]) * sr

        return listOfLMChainsMatrix

    def getConvexHullAreaAndPerimeter(self, listOfLMResid):
        """ Method to calculate the convex hull area and perimeter containing the residuals """

        hull = ConvexHull(listOfLMResid)

        # For 2-Dimensional convex hulls volume attribute equals to the area.
        area = hull.volume
        perimeter = 0

        hullVertices = []

        for position in hull.vertices:
            hullVertices.append(hull.points[position])

        for i in range(len(hullVertices)):
            shiftedIndex = (i + 1) % len(hullVertices)
            perimeter += self.getDistance2D(np.array(hullVertices[i]), np.array(hullVertices[shiftedIndex]))

        return area, perimeter

    def getMaximumDistance(self, listOfLMResid):
        """ Method to calculate the maximum residual distance belonging to the same landmark chain. """

        maxDistance = 0

        for resid in listOfLMResid:
            distance = self.getDistance2D(np.array([0, 0]), resid)

            if distance > maxDistance:
                maxDistance = distance

        return maxDistance

    @staticmethod
    def getDistance2D(coordinate2Da, coordinate2Db):
        """ Method to calculate the distance between two 2D coordinates. """

        distanceVector = coordinate2Da - coordinate2Db
        distanceVector = [i ** 2 for i in distanceVector]
        distance = np.sqrt(sum(distanceVector))

        return distance

    # @staticmethod
    # def getLandMarkModelFromTs(SoLM, tsId):
    #     """ This method inputs a set of Landmark Models and the TsId and search for a Landmark Model with a coincident
    #     tsId.
    #     :param SoLM: input set of landmark models.
    #     :param tsId: is of the landmark to search.
    #     """
    #
    #     for lm in SoLM:
    #         if lm.getTsId() == tsId:
    #             return lm
    #
    # @staticmethod
    # def getTiltSeriesFromTs(SoTS, tsId):
    #     """ This method inputs a set of Landmark Models and the TsId and search for a Landmark Model with a coincident
    #     tsId.
    #     :param SoTS: input set of landmark models.
    #     :param tsId: is of the landmark to search.
    #     """
    #
    #     for ts in SoTS:
    #         if ts.getTsId() == tsId:
    #             return ts
