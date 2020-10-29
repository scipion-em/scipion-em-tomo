# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             Tomas Majtner (tmajtner@cnb.csic.es)  -- streaming version
# *             David Herreros Calero (dherreros@cnb.csic.es) -- Tomo version
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
from scipy.spatial.distance import pdist

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

from pwem.emlib.image import ImageHandler

from xmipp3.convert import openMd

from tomo.protocols import ProtTomoPicking
from tomo.utils import delaunayTriangulation, normalFromMatrix, extractVesicles


class ProtTomoScoreCoordinates(ProtTomoPicking):
    '''Scoring and (optional) filtering of coordinates based on different scoring
    functions (normals angle, carbon distance, neighbour distance)'''

    _label = 'score/filter coordinates'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input 3D coordinates", important=True,
                      help='Select the set of 3D coordinates to compare')
        form.addParam('outliers', params.BooleanParam, default=True,
                      label="Score outluiers?")
        form.addParam('outliersThreshold', params.FloatParam, default=10,
                      label="Outliers distance threshold", condition='outliers == True',
                      help='Distance expressed in pixels')
        form.addParam('angle', params.BooleanParam, default=True,
                      label="Score normals?")
        form.addParam('angleThreshold', params.FloatParam, default=0,
                      label="Angle threshold", condition='angle == True',
                      help='Angle expressed in ***')
        form.addParam('carbon', params.BooleanParam, default=True,
                      label="Score carbon closeness?")
        form.addParam('carbonThreshold', params.FloatParam, default=10,
                      label="Carbon distance threshold", condition='outliers == True',
                      help='Distance expressed in pixels')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        pwutils.makePath(self._getExtraPath('inputCoords'))
        pwutils.makePath(self._getExtraPath('outputCoords'))
        coordinates = self.inputCoordinates.get()
        self.tomos = coordinates.getPrecedents()
        self.tomoNames = [pwutils.removeBaseExt(tomo.getFileName()) for tomo in self.tomos]
        self._insertFunctionStep('computeParams', coordinates)
        if self.outliers.get():
            self._insertFunctionStep('detectOutliers')
        if self.carbon.get():
            self._insertFunctionStep('detectCarbonCloseness', coordinates)
        self._insertFunctionStep('createOutputStep')

    def computeParams(self, coordinates):
        self.tomo_vesicles = extractVesicles(coordinates)

    def detectOutliers(self):
        self.scoreOutliers = []
        self.maskOutliers = []
        threshold = np.exp(-self.outliersThreshold.get())
        for tomo in self.tomoNames:
            for vesicle in self.tomo_vesicles[tomo]['vesicles']:
                distance = pdist(vesicle)
                self.scoreOutliers.append(np.exp(-distance))
                self.maskOutliers.append(self.scoreOutliers[-1] >= threshold)

    def detectCarbonCloseness(self, coordinates):
        for tomo, tomoName in zip(self.tomos.iterItems(), self.tomoNames):
            projFile, dfactor = self.projectTomo(tomo)
            coordList = self.generateCoordList(tomo, coordinates, dfactor)
            self.writeTomoCoordinates(tomo, coordList, self._getTomoPos(projFile))
            args = '-i %s -c %s -o %s -b %d' \
                   % (projFile, self._getExtraPath('inputCoords'),
                      self._getExtraPath('outputCoords'), tomo.getSamplingRate())
            self.runJob('xmipp_deep_micrograph_cleaner', args)


    def createOutputStep(self):
        pass

    # --------------------------- UTILS functions ----------------------
    def projectTomo(self, tomo):
        dfactor = None
        outFile = pwutils.removeBaseExt(tomo.getFileName()) + '_projected.mrc'
        ih = ImageHandler()
        outProjection = ih.createImage()
        tomoData = np.squeeze(ih.read(tomo.getFileName()).getData())
        projection = np.sum(tomoData, axis=0)
        outProjection.setData(projection)
        ih.write(outProjection, self._getExtraPath(outFile))
        bgDim = max(projection.shape)
        if bgDim > 1000:
            dfactor = bgDim / 1000
            args = '-i %s --step %f --method fourier' % \
                   (self._getExtraPath(outFile), dfactor)
            self.runJob("xmipp_transform_downsample", args)
        return self._getExtraPath(outFile), dfactor

    def square_to_condensed(self, i, j, n):
        assert i != j, "no diagonal elements in condensed matrix"
        if i < j:
            i, j = j, i
        return n * j - j * (j + 1) // 2 + i - 1 - j

    def _getTomoPos(self, fileName):
        """ Return the corresponding .pos file for a given tomogram. """
        baseName = pwutils.removeBaseExt(fileName)
        return self._getExtraPath('inputCoords', baseName + ".pos")

    def writeTomoCoordinates(self, tomo, coordList, outputFn, isManual=True,
                             getPosFunc=None):
        """ Write the pos file as expected by Xmipp with the coordinates
        of a given tomogram.
        Params:
            tomo: input tomogram.
            coordList: list of (x, y) pairs of the mic coordinates.
            outputFn: output filename for the pos file .
            isManual: if the coordinates are 'Manual' or 'Supervised'
            getPosFunc: a function to get the positions from the coordinate,
                it can be useful for scaling the coordinates if needed.
        """
        if getPosFunc is None:
            getPosFunc = lambda coord: coord.getPosition()

        state = 'Manual' if isManual else 'Supervised'
        f = openMd(outputFn, state)

        for coord in coordList:
            x, y, z = getPosFunc(coord)
            f.write(" %06d   1   %d  %d  %d   %06d \n"
                    % (coord.getObjId(), x, y, 1, tomo.getObjId()))

        f.close()

    def generateCoordList(self, tomo, coordinates, dfactor):
        if dfactor is not None:
            coordList = []
            for coord in coordinates.iterCoordinates(volume=tomo):
                cloned_coord = coord.clone()
                x, y, z = cloned_coord.getPosition()
                cloned_coord.setPosition(x / dfactor,
                                         y / dfactor,
                                         z / dfactor)
                coordList.append(cloned_coord)
            return coordList
        else:
            return [coord.clone() for coord in coordinates.iterCoordinates(volume=tomo)]

    # --------------------------- INFO functions ----------------------

