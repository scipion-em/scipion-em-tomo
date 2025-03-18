# **************************************************************************
# *
# * Authors:     Oier Lauzirika Zarrabeitia (oierlauzi@bizkaia.eu) [1]
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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol

from tomo.objects import (SetOfTiltSeries, TiltImage,
                          SetOfCoordinates3D, Coordinate3D,
                          SetOfLandmarkModels, LandmarkModel )
from tomo.protocols import ProtTomoBase
import tomo.constants as constants
from tomo.utils import getObjFromRelation

import enum
import numpy as np

class OutputProjectCoordinates(enum.Enum):
    landmarkModels = SetOfLandmarkModels
    
class ProtProjectCoordinates(EMProtocol, ProtTomoBase):
    """
    Project 3D coordinates into a set of landmarks.
    """

    _label = 'project coordinates'
    _devStatus = BETA
    _possibleOutputs = OutputProjectCoordinates

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputCoordinates',
                      params.PointerParam,
                      label="Coordinates",
                      pointerClass=SetOfCoordinates3D,
                      important=True,
                      help='Coordinates to be projected onto the tilt-series' )

        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      label="Tilt-series",
                      pointerClass=SetOfTiltSeries,
                      allowsNull=True,
                      help='Tilt series on which coordinates are projected. '
                      'If not specified, it will be deduced from the coordinates.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ----------------------------
    def createOutputStep(self):
        inputTiltSeries = self._getInputSetOfTiltSeries()
        inputCoodinates = self._getInputSetOfCoordinates3d()
        landmarkSize = inputCoodinates.getBoxSize() * inputTiltSeries.getSamplingRate()
        offset = np.array(inputTiltSeries.getDim()[:2]) / 2
        scale = inputCoodinates.getSamplingRate() / inputTiltSeries.getSamplingRate()
        
        outputSetOfLandmarkModels: SetOfLandmarkModels = self._createSetOfLandmarkModels()
        outputSetOfLandmarkModels.copyInfo(inputTiltSeries)
        outputSetOfLandmarkModels.setSetOfTiltSeries(inputTiltSeries)
        
        for tiltSeries in inputTiltSeries:
            tsId = tiltSeries.getTsId()
            where = '%s="%s"' % (Coordinate3D.TOMO_ID_ATTR, tsId)
            landmarkModel = LandmarkModel(
                tsId=tsId,
                fileName=self._getExtraPath(tsId + '.sfid'),
                size=landmarkSize,
                applyTSTransformation=False
            )
            landmarkModel.setTiltSeries(tiltSeries)
            
            for coordinate3d in inputCoodinates.iterItems(where=where):
                position3d = np.array(coordinate3d.getPosition(constants.SCIPION) + (1, ))
                position3d[:3] *= scale # Scale to match the binning of the TS
                chainId = coordinate3d.getObjId()
                
                for tiltImage in tiltSeries:
                    position2d = self._projectCoordinate(tiltImage, position3d)
                    position2d = position2d[:2]
                    position2d += offset
                    
                    landmarkModel.addLandmark(
                        xCoor=position2d[0],
                        yCoor=position2d[1],
                        tiltIm=tiltImage.getIndex(),
                        chainId=chainId,
                        xResid=0,
                        yResid=0
                    )

            outputSetOfLandmarkModels.append(landmarkModel)
                    
        self._defineOutputs(**{OutputProjectCoordinates.landmarkModels.name: outputSetOfLandmarkModels})
        self._defineSourceRelation(self.inputCoordinates, outputSetOfLandmarkModels)
        if self.inputTiltSeries.get() is not None:
            self._defineSourceRelation(self.inputTiltSeries, outputSetOfLandmarkModels)
    
    # --------------------------- UTILS functions ----------------------------
    def _getInputSetOfCoordinates3d(self) -> SetOfCoordinates3D:
        return self.inputCoordinates.get()

    def _getInputSetOfTiltSeries(self) -> SetOfTiltSeries:
        result = self.inputTiltSeries.get()
        if result is None:
            coordinates = self._getInputSetOfCoordinates3d()
            result = getObjFromRelation(coordinates, self, SetOfTiltSeries) 
            
        return result

    def _projectCoordinate(self,
                           tiltImage: TiltImage,
                           position3d: np.ndarray ) -> np.ndarray:
        projection = self._getProjectionMatrix(tiltImage)
        projected = projection @ position3d
        
        if tiltImage.hasTransform():
            transform = tiltImage.getTransform().getMatrix()
            projected = np.linalg.inv(transform) @ projected

        return projected    
    
    def _getProjectionMatrix(self, tiltImage: TiltImage) -> np.ndarray:
        tiltAngle = tiltImage.getTiltAngle()
        tiltAngle = np.deg2rad(tiltAngle)
        return np.array([
            [np.cos(tiltAngle), 0, np.sin(tiltAngle), 0],
            [0,                 1, 0,                 0],
            [0,                 0, 0,                 1]
        ])
        
    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        result = []
        
        tiltSeries = self._getInputSetOfTiltSeries()
        if tiltSeries is None:
            result.append(
                'Could not deduce tilt series from the coordinates. '
                'Please, specify the tilt series manually. '
            )
        
        return result
