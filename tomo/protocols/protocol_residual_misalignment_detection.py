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
from sklearn.decomposition import PCA

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

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for lm in self.inputSetOfLandmarkModels.get():
            lmObjId = lm.getObjId()
            self._insertFunctionStep(self.filterLandmarkModelsByResidualStatistics, lmObjId)

        # self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ----------------------------
    def filterLandmarkModelsByResidualStatistics(self, lmObjId):
        lm = self.inputSetOfLandmarkModels.get()[lmObjId]
        ts = lm.getTiltSeries()

        tsId = lm.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        residualStatistics = self.getResidualStatistics(lm, ts)

    # --------------------------- UTILS functions ----------------------------
    def generateResidualList(self, lm, ts):
        """ This method generates a residual list from each landmark model input set"""
        firstItem = ts.getFirstItem()

        listOfLMChainsMatrix = lm.retrieveLandmarkModelChains()

        xDim = firstItem.getXDim()
        yDim = firstItem.getYDim()

        listOfLMChainsMatrixNorm = self.normalizeResiduals(listOfLMChainsMatrix, xDim, yDim)

        # Vector containing the statistical features calculated for the set of residuals from each landmark chain
        # belonging to the landmark model
        residualStatistics = []

        for listOfLMChainsNorm in listOfLMChainsMatrixNorm:

            listOfLMResid = []

            for lm in listOfLMChainsNorm:
                listOfLMResid.append([lm[4], lm[5]])

            chArea, chPerimeter = self.getResidualStatistics(listOfLMResid)

            maxDistance = self.getMaximumDistance(listOfLMResid)

            totalDistance = self.getTotalDistance(listOfLMResid)

            pca = self.getPCA(listOfLMResid)

            residualStatistics.append([totalDistance, maxDistance, chArea, chPerimeter, pca[0], pca[1]])

        return residualStatistics

    def getResidualStatistics(self, listOfLMResid):
        """ This method calculate the statistics from each chain of residual to filter the properly aligned set of
        tilt series. """

        convexHullArea, convexHullPerimeter = self.getConvexHullAreaAndPerimeter(listOfLMResid)

        return convexHullArea, convexHullPerimeter

    @staticmethod
    def normalizeResiduals(listOfLMChainsMatrix, xDim, yDim):
        """ This method multiplies the residual values by the sampling rate of the tilt-series to normalize its values
        to a common scale.
        :param listOfLMChainsMatrix: input set of landmark models separated by chainID.
        :param xDim: X dimension of the input tilt-image used as a normalization factor.
        :param yDim: Y dimension of the input tilt-image used as a normalization factor.
        """

        for listOfLMChains in listOfLMChainsMatrix:
            for lm in listOfLMChains:
                lm[4] = (float(lm[4]) / xDim) * 100
                lm[5] = (float(lm[5]) / yDim) * 100

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

    def getTotalDistance(self, listOfLMResid):
        """ Method to calculate the total distance of the trajectory """

        totalDistance = 0

        for resid in listOfLMResid:
            totalDistance += self.getDistance2D(np.array([0, 0]), resid)

        return totalDistance

    @staticmethod
    def getPCA(listOfLMResid):
        """ Method to calculate the PCA of the scatter cloud formed by the residual vectors. """

        pca = PCA(n_components=2)
        pca.fit(listOfLMResid)

        # Return only the first component (remove redundant information)
        return pca.components_[0]

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
