# **************************************************************************
# *
# * Authors:     Oier Lauzirika Zarrabeitia (olauzirika@cnb.csic.es)
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
import numpy as np
import math

from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam, EnumParam, 
                                        BooleanParam, GT)

from pwem.objects import Set
from pwem.protocols import EMProtocol

import tomo.constants as const
from tomo.objects import (SetOfTomograms, Tomogram, SetOfMeshes, MeshPoint,
                          SetOfCoordinates3D, Coordinate3D)
from tomo.protocols import ProtTomoBase
from tomo.utils import fit_ellipsoid, fit_sphere

COORDINATES = 'Coordinates'

class ProtGeometricPickingFilter(EMProtocol, ProtTomoBase):
    """ Fit a geometric shape to a set of meshes and filter particles
    """
    _label = 'geometric picking filter'
    _devStatus = BETA
    _possibleOutputs = {COORDINATES: SetOfCoordinates3D}

    SHAPE_SPHERE = 0
    SHAPE_ELLIPSOID = 1

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', 
                      PointerParam, pointerClass=SetOfCoordinates3D,
                      label='Input coordinates', important=True,
                      help='Select the coordinates to be filtered')
        form.addParam('inputMeshes', 
                      PointerParam, pointerClass=SetOfMeshes,
                      label='Input meshes', important=True,
                      help='Select the meshes to be used as template')
        
        form.addParam('shape', EnumParam, label='Shape',
                      choices=['Sphere', 'Ellipsoid'],
                      default=0,
                      help='Geometric shape to be fitted to the meshes')
        form.addParam('minimumFitQuality', FloatParam, label='Minimum fit quality',
                      default=16.0, validators=[GT(0.0)],
                      help='Minimum fit quality for a shape to be picked. '
                      'Measured as RMS error in pixels. Use inf to disable')
        
        line = form.addLine('Fit radius range',
                             help='Only shapes in this rage will be considered. '
                                  'Radius are masured in pixels')
        line.addParam('minimumFitRadius', FloatParam, default=0.1, label='Minimum')
        line.addParam('maximumFitRadius', FloatParam, default=1000.0, label='Maximum')
        
        line = form.addLine('Offset range',
                             help='A band of these dimensions will be selected '
                                  'around the fitted shape. In pixels.')
        line.addParam('minimumOffset', FloatParam, default=-10.0, label='Minimum')
        line.addParam('maximumOffset', FloatParam, default=10.0, label='Maximum')
        
        form.addParam('writeGroupId', BooleanParam, label='Write group id',
                      default=True,
                      help='Write the group identifier from the assotiated mesh')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        coodinates = self._getInputCoordintates()
        tomograms: SetOfTomograms = coodinates.getPrecedents()

        self._insertFunctionStep(self.createOutputStep)
        for tomogram in tomograms:
            tsId: str = tomogram.getTsId()
            self._insertFunctionStep(self.processTomogramStep, tsId)
        self._insertFunctionStep(self.closeOutputStep)
        
    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        inputCoordinates = self._getInputCoordintates()
        tomograms: SetOfTomograms = inputCoordinates.getPrecedents()
        outputCoordinates = self._createSetOfCoordinates3D(tomograms)
        outputCoordinates.copyInfo(inputCoordinates)
        outputCoordinates.setStreamState(Set.STREAM_OPEN)
        
        self._defineOutputs(**{COORDINATES: outputCoordinates})
        self._defineSourceRelation(self.inputMeshes, outputCoordinates)
        self._defineSourceRelation(self.inputCoordinates, outputCoordinates)

    def processTomogramStep(self, tsId: str):
        meshes = self._getInputMeshes()
        coodinates = self._getInputCoordintates()
        
        processed = set()
        query = "%s='%s'" % (MeshPoint.TOMO_ID_ATTR, tsId)
        groupIds = meshes.getUniqueValues(MeshPoint.GROUP_ID_ATTR, where=query)
        for groupId in groupIds:
            self._processMesh(coodinates, meshes, tsId=tsId, groupId=groupId, processed=processed)
        
    def closeOutputStep(self):
        self._closeOutputSet()
    
    # --------------------------- UTILS functions ------------------------------
    def _getInputCoordintates(self) -> SetOfCoordinates3D:
        return self.inputCoordinates.get()
    
    def _getInputMeshes(self) -> SetOfMeshes:
        return self.inputMeshes.get()

    def _processMesh(self, 
                     coordinates: SetOfCoordinates3D,
                     meshes: SetOfMeshes,
                     tsId: str,
                     groupId: int,
                     processed: set ):
        points = self._getGroupPoints(meshes=meshes, tsId=tsId, groupId=groupId)
        
        if self.shape == self.SHAPE_SPHERE:
            self._processSphere(points, coordinates, tsId=tsId, groupId=groupId, processed=processed)
        elif self.shape == self.SHAPE_ELLIPSOID:
            self._processEllipsoid(points, coordinates, tsId=tsId, groupId=groupId, processed=processed)
            
    def _getGroupPoints(self, 
                        meshes: SetOfMeshes, 
                        groupId: int,
                        tsId: str ) -> np.ndarray:
        points = []
        query = "%s='%s' AND %s='%s'" % (MeshPoint.GROUP_ID_ATTR, groupId, 
                                         MeshPoint.TOMO_ID_ATTR, tsId)
        for coordinate in meshes.iterItems(where=query):
            points.append(coordinate.getPosition(const.SCIPION))
        return np.array(points)

    def _append(self, coordinate: Coordinate3D, groupId: int, processed: set):
        coordinates: SetOfCoordinates3D = getattr(self, COORDINATES)
        if coordinate.getObjId() not in processed:
            if self.writeGroupId:
                coordinate.setGroupId(groupId)
            coordinates.append(coordinate)
            processed.add(coordinate.getObjId())

    def _filterSphere(self,
                      coordinates: SetOfCoordinates3D,
                      tsId: str,
                      groupId: int,
                      centre: np.ndarray,
                      radius: float,
                      processed: set ):
        minDistance2 = (radius + self.minimumOffset.get())**2
        maxDistance2 = (radius + self.maximumOffset.get())**2
        query = "%s='%s'" % (Coordinate3D.TOMO_ID_ATTR, tsId)
        for coordinate in coordinates.iterItems(where=query):
            position = coordinate.getPosition(const.SCIPION)
            distance2 = np.sum(np.square(centre - position))
            if minDistance2 < distance2 and distance2 < maxDistance2:
                self._append(coordinate, groupId, processed)

    def _processSphere(self, 
                       points: np.ndarray, 
                       coordinates: SetOfCoordinates3D,
                       tsId: str,
                       groupId: int,
                       processed: set ):
        centre, radius, chi2 = fit_sphere(points)
        rms = np.sqrt(chi2 / len(points))
        
        minRadius = self.minimumFitRadius.get()
        maxRadius = self.maximumFitRadius.get()
        if rms <= self.minimumFitQuality and minRadius <= radius and radius <= maxRadius:
            self._filterSphere(coordinates, tsId, groupId, centre, radius, processed)

    def _filterEllipsoid(self,
                         coordinates: SetOfCoordinates3D,
                         tsId: str,
                         groupId: int,
                         centre: np.ndarray,
                         radii: np.ndarray,
                         evecs: np.ndarray,
                         processed: set ):
        minRadii2 = np.square(radii + self.minimumOffset.get())
        maxRadii2 = np.square(radii + self.maximumOffset.get())
        query = "%s='%s'" % (Coordinate3D.TOMO_ID_ATTR, tsId)
        for coordinate in coordinates.iterItems(where=query):
            position = coordinate.getPosition(const.SCIPION)
            centered = position - centre # Center the coordinate system in the ellipsoid
            aligned = evecs.T @ centered # Axis-align the coordinate system
            
            aligned2 = np.square(aligned)
            minEllipsoid = np.sum(aligned2 / minRadii2) - 1
            maxEllipsoid = np.sum(aligned2 / maxRadii2) - 1
            if minEllipsoid >= 0 and maxEllipsoid <= 0:
                self._append(coordinate, groupId, processed)
        
    def _processEllipsoid(self, 
                          points: np.ndarray,
                          coordinates: SetOfCoordinates3D,
                          tsId: str,
                          groupId: int,
                          processed: set ):
        center, radii, _, evecs, chi2 = fit_ellipsoid(points)
        rms = np.sqrt(chi2 / len(points))
        
        minRadius = self.minimumFitRadius.get()
        maxRadius = self.maximumFitRadius.get()
        if rms <= self.minimumFitQuality and minRadius <= min(radii) and max(radii) <= maxRadius:
            self._filterEllipsoid(coordinates, tsId, groupId, center, radii, evecs, processed)
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        result = []
        coordinates = self._getInputCoordintates()
        meshes = self._getInputMeshes()
        
        if self.minimumFitRadius > self.maximumFitRadius:
            result.append('Minimum fit radius must be smaller than the maximum radius')
        
        if self.minimumOffset > self.maximumOffset:
            result.append('Minimum offset must be smaller than the maximum offset')

        if coordinates.getSamplingRate() != meshes.getSamplingRate():
            result.append('Coordinates and meshes must have the same sampling')
        
        return result
    