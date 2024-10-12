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
                                        IntParam, GE, GT)

from pwem.protocols import EMProtocol

import tomo.constants as const
from tomo.objects import (SetOfCoordinates3D, Coordinate3D,
                          SetOfMeshes, MeshPoint, Transform)
from tomo.protocols import ProtTomoBase
from tomo.utils import fit_ellipsoid, fit_sphere

COORDINATES = 'Coordinates'

class ProtGeometricPicking(EMProtocol, ProtTomoBase):
    """ Pick coordinates around a geometric shape
    """
    _label = 'geometric picking'
    _devStatus = BETA
    _possibleOutputs = {COORDINATES: SetOfCoordinates3D}

    SHAPE_SPHERE = 0
    SHAPE_ELLIPSOID = 1

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMeshes', 
                      PointerParam, pointerClass=SetOfMeshes,
                      label='Input meshes', important=True,
                      help='Select the Coordinates3D to be filtered')
        
        form.addParam('shape', EnumParam, label='Shape',
                      choices=['Sphere', 'Ellipsoid'],
                      default=0,
                      help='Geometric shape to be fitted to the meshes')
        form.addParam('minimumFitQuality', FloatParam, label='Minimum fit quality',
                      default=16.0, validators=[GT(0.0)],
                      help='Minimum fit quality for a shape to be picked. '
                      'Measured as RMS error in pixels. Use inf to disable')
        
        line = form.addLine('Radius range',
                             help='Only shapes in this rage will be considered. '
                                  'Radius are masured in pixels')
        line.addParam('minimumRadius', FloatParam, default=0.1, label='Minimum')
        line.addParam('maximumRadius', FloatParam, default=1000.0, label='Maximum')
        
        form.addParam('density', FloatParam, label='Picking density',
                      default = 2.0, validators=[GE(0.0)],
                      help='Higher values involve higher picking densities')
        form.addParam('radialOffset', FloatParam, label='Radial offset (px)',
                      default=0.0, 
                      help='Increment or decrement the radius of the ellipsoid '
                           'when picking')
        form.addParam('boxSize', IntParam, label='Box size',
                      default = 32, validators=[GT(0)],
                      help='Box size of the picked particles')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)
        
    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        meshes = self._getInputMeshes()
        tomograms = meshes.getPrecedents()
        outputCoordinates = self._createSetOfCoordinates3D(tomograms)
        outputCoordinates.copyInfo(meshes)
        outputCoordinates.setBoxSize(self.boxSize.get())
        
        groupIds =  meshes.getUniqueValues(MeshPoint.GROUP_ID_ATTR)
        for groupId in groupIds:
            self._processMesh(outputCoordinates, meshes, groupId)
            
        self._defineOutputs(**{COORDINATES: outputCoordinates})
        self._defineSourceRelation(self.inputMeshes, outputCoordinates)
    
    # --------------------------- UTILS functions ------------------------------
    def _getInputMeshes(self) -> SetOfMeshes:
        return self.inputMeshes.get()
    
    def _processMesh(self, 
                     outputCoordinates: SetOfCoordinates3D, 
                     meshes: SetOfMeshes,
                     groupId: int):
        points = self._getGroupPoints(meshes=meshes, groupId=groupId)
        precedent: MeshPoint = meshes[{MeshPoint.GROUP_ID_ATTR: groupId}]
        
        coordinates = None 
        normals = None
        rms = None       
        if self.shape == self.SHAPE_SPHERE:
            coordinates, normals, rms = self._processSphere(points)
        elif self.shape == self.SHAPE_ELLIPSOID:
            coordinates, normals, rms = self._processEllipsoid(points)

        if coordinates is not None:
            self._pointsToCoordinates(outputCoordinates, precedent, 
                                      coordinates, normals, rms)
            
    def _getGroupPoints(self, meshes: SetOfMeshes, groupId: int) -> np.ndarray:
        points = []
        query = '%s=%d' % (MeshPoint.GROUP_ID_ATTR, groupId)
        for coord in meshes.iterItems(where=query):
            points.append(coord.getPosition(const.SCIPION))
        return np.array(points)
    
    def _computeTransformMatricesFromNormals(self, normals: np.ndarray) -> np.ndarray:
        result = np.zeros((len(normals), 4, 4))
        l = np.linalg.norm(normals[...,0:2], axis=-1)
        r = 1.0 / l
        result[...,0,0] = normals[...,1]*r
        result[...,0,1] = -normals[...,0]*r
        result[...,0,2] = 0.0
        result[...,1,0] = normals[...,0]*normals[...,2]*r
        result[...,1,1] = normals[...,1]*normals[...,2]*r
        result[...,1,2] = -l
        result[...,2,:3] = normals
        result[...,3,3] = 1 # Affine
        return result
    
    def _pointsToCoordinates(self,
                             outputCoordinates: SetOfCoordinates3D,
                             precedent: MeshPoint,
                             points: np.ndarray,
                             normals: np.ndarray,
                             rms: float):
        
        matrices = self._computeTransformMatricesFromNormals(normals)
        for (x, y, z), matrix in zip(points, matrices):
            coordinate = Coordinate3D()
            coordinate.copy(precedent, copyId=False)
            coordinate.setPosition(x, y, z, const.SCIPION)
            coordinate.setScore(rms)
            coordinate.setMatrix(matrix) # TODO convention
            outputCoordinates.append(coordinate)
    
    def _getSampleCoundForArea(self, area) -> int:
        boxSectionArea = self.boxSize.get() ** 2
        return math.ceil(self.density.get()*(area / boxSectionArea))
    
    def _sphereArea(self, r):
        return 4*np.pi*r*r
            
    def _ellipsoidArea(self, a, b, c):
        P = 1.6075
        return 4*np.pi*(((a*b)**P + (a*c)**P + (b*c)**P)/3)**(1/P)
    
    def _sampleUnitSphere(self, nPoints: int):
        theta = 2*np.pi*np.random.rand(nPoints)
        phi = np.arccos(2*np.random.rand(nPoints)-1)
        x = np.cos(theta)*np.sin(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(phi)
        return np.stack((x,y,z), axis=1)
    
    def _processSphere(self, points: np.ndarray) -> np.ndarray:
        center, radius, chi2 = fit_sphere(points)
        rms = np.sqrt(chi2 / len(points))
        
        coordinates = None
        normals = None
        minRadius = self.minimumRadius.get()
        maxRadius = self.maximumRadius.get()
        if rms <= self.minimumFitQuality and minRadius <= radius and radius <= maxRadius:
            radius += self.radialOffset.get()
            
            area = self._sphereArea(radius)
            nPoints = self._getSampleCoundForArea(area)
            
            sphere = self._sampleUnitSphere(nPoints)
            coordinates = radius*sphere
            coordinates += center

        return coordinates, sphere, rms
    
    def _processEllipsoid(self, points: np.ndarray) -> np.ndarray:
        center, radii, _, evecs, chi2 = fit_ellipsoid(points)
        rms = np.sqrt(chi2 / len(points))
        
        coordinates = None
        normals = None
        minRadius = self.minimumRadius.get()
        maxRadius = self.maximumRadius.get()
        if rms <= self.minimumFitQuality and minRadius <= min(radii) and max(radii) <= maxRadius:
            radii += self.radialOffset.get()
            radii = np.maximum(radii, 0.0)
            vertices = evecs * radii
            
            # Estimate the sample count
            area = self._ellipsoidArea(*radii)
            nPoints = self._getSampleCoundForArea(area)
            
            # Compute coordinates
            sphere = self._sampleUnitSphere(nPoints)
            coordinates = sphere @ vertices.T # (vertices @ sphere.T).T
            coordinates += center
            
            # Compute normals
            normVertices = np.empty((3, 3))
            normVertices[:,0] = np.cross(vertices[:,1], vertices[:,2])
            normVertices[:,1] = np.cross(vertices[:,2], vertices[:,0])
            normVertices[:,2] = np.cross(vertices[:,0], vertices[:,1])
            normals = sphere @ normVertices.T # (normVertices @ sphere.T).T
            normals /= np.linalg.norm(normals, axis=-1, keepdims=True)

        return coordinates, normals, rms
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        result = []
        
        if self.minimumRadius > self.maximumRadius:
            result.append('Minimum radius must be smaller than the maximum radius')
        
        return result
    