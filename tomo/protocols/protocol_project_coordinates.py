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
    Projects a set of 3D coordinates onto a tilt-series and converts them
    into landmark models.

    AI Generated:

    Project Coordinates (ProtProjectCoordinates) — User Manual
        Overview

        The Project Coordinates protocol projects 3D coordinates onto the
        images of a tilt-series and generates a set of landmark models.
        Its main purpose is to transform volumetric coordinate information
        into 2D landmark positions that are geometrically consistent with
        each tilt image.

        In cryo-electron tomography workflows, this protocol is especially
        useful when 3D positions—such as particle centers, annotated
        features, or reference points—must be mapped back onto the
        original acquisition images. This enables direct visualization of
        the projected positions and provides input for downstream
        landmark-based processing.

        Inputs and General Workflow

        The protocol requires a set of 3D coordinates as the main input.
        Optionally, the user may also provide a tilt-series. If the
        tilt-series is not explicitly provided, the protocol attempts to
        deduce it automatically from the coordinate metadata.

        During execution, the protocol first retrieves the input
        tilt-series and computes the sampling relationship between the
        coordinate set and the tilt-series. This scaling step ensures
        that the coordinates are expressed in the same spatial sampling
        system as the tilt images.

        For each tilt-series, the protocol creates one landmark model.
        Then, every 3D coordinate belonging to that tilt-series is
        projected onto every tilt image, producing a set of 2D landmark
        positions across the angular series.

        Coordinate Scaling and Centering

        Since coordinates and tilt-series may have different sampling
        rates, the protocol rescales each coordinate before projection.
        This guarantees that the projected positions match the binning of
        the tilt-series.

        In addition, the projected 2D coordinates are shifted by the
        image center offset. This is required because the geometric
        projection is computed around the origin, whereas image
        coordinates are expressed relative to the pixel grid.

        From a practical perspective, this means that the final landmark
        positions are directly usable in the coordinate system of the
        tilt images.

        Projection Geometry

        The protocol applies a simple projection model based on the tilt
        angle of each image.

        For every tilt image, a projection matrix is generated from the
        corresponding tilt angle. This matrix rotates the 3D coordinate
        according to the acquisition geometry and projects it into 2D
        image space.

        Biologically, this reproduces the expected apparent position of a
        3D object when observed under the viewing angle of each tilt
        image.

        Tilt-Image Transformations

        If a tilt image contains an associated geometric transformation,
        the protocol applies the inverse transformation after projection.

        This is important because tilt images may already contain prior
        alignment corrections. By applying the inverse transformation,
        the projected landmarks remain consistent with the actual image
        coordinates.

        In practical tomography workflows, this allows landmarks to be
        correctly placed even when the tilt-series has undergone previous
        alignment or motion correction steps.

        Landmark Construction

        For each projected coordinate, the protocol adds one landmark per
        tilt image.

        Each landmark stores:
            - the projected X coordinate,
            - the projected Y coordinate,
            - the tilt-image index,
            - the chain identifier corresponding to the original 3D
              coordinate.

        As a result, every 3D coordinate becomes a landmark trajectory
        across the tilt-series.

        This structure is particularly useful for visualization,
        refinement, or subsequent landmark-based alignment procedures.

        Outputs and Their Interpretation

        The main output of the protocol is a SetOfLandmarkModels.

        One landmark model is created for each tilt-series. Each model
        contains all projected landmark tracks corresponding to the 3D
        coordinates associated with that series.

        From a biological and experimental point of view, the output
        allows users to verify whether reconstructed 3D positions are
        geometrically consistent with the raw tilt images.

        This can be especially useful for:
            - validating particle localization,
            - checking reconstruction consistency,
            - preparing landmarks for downstream tomography processing.

        Validation and Input Consistency

        Before execution, the protocol verifies that a tilt-series can be
        obtained.

        If no tilt-series is explicitly provided and none can be deduced
        from the coordinates, execution cannot proceed and the user is
        asked to specify the tilt-series manually.

        This validation step is essential because the projection geometry
        depends entirely on the tilt-series acquisition parameters.

        Practical Recommendations

        In routine tomography workflows, this protocol is most reliable
        when the input coordinates originate from the same tilt-series
        that will be used for projection.

        If coordinates come from data that has been resampled, rebinned,
        or transformed independently, special care should be taken to
        ensure that the sampling rates remain consistent.

        It is often good practice to visually inspect the resulting
        projected landmarks after execution. Incorrect scaling,
        mismatched tilt-series, or inconsistent geometry usually become
        immediately visible at this stage.

        Final Perspective

        For most cryo-ET users, this protocol provides a direct link
        between reconstructed 3D information and the original tilt
        images.

        Although mathematically simple, this projection step is
        biologically important because it enables validation of spatial
        interpretations directly against the experimental acquisition
        data.
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
