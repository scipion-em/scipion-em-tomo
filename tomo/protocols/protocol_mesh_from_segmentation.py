# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
# *              Oier Lauzirika Zarrabeita (oierlauzi@bizkaia.eu)
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

from enum import Enum
import math
import numpy as np
import skimage.morphology

from pwem.protocols import EMProtocol
from pwem.objects import Set, Integer
from pwem.emlib.image import ImageHandler

from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam, 
                                        BooleanParam, IntParam)

from tomo.objects import (SetOfMeshes, MeshPoint, SetOfTomoMasks, TomoMask,
                          SetOfTomograms, Tomogram)
from tomo.protocols import ProtTomoBase
import tomo.constants as const

class OutputMeshesFromSegmentation(Enum):
    Meshes = SetOfMeshes

class ProtMeshFromSegmentation(EMProtocol, ProtTomoBase):
    """
    Creates meshes based on segmentations or voxels values (TomoMasks).
    """
    _label = 'meshes from tomo mask'
    _devStatus = BETA
    _possibleOutputs = OutputMeshesFromSegmentation

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMasks', PointerParam, pointerClass=SetOfTomoMasks,
                      label='Tomo Masks', important=True,
                      help='Set of tomo mask from which the meshes will be created')
        form.addParam('inputTomograms', PointerParam, pointerClass=SetOfTomograms,
                      label='Tomograms',
                      help='Tomograms to which the meshes will be asotiated')
        
        form.addSection(label='Parameters')
        form.addParam('smoothMask', BooleanParam, label='Smooth mask',
                      important=True, default=False,
                      help='Wether the input is a segmentation or a smooth mask')
        form.addParam('backgroundLabel', IntParam, label='Background label',
                      condition='not smoothMask', default=0,
                      help='Label in the segmentation to be considered as '
                           'background')

        line = form.addLine('Threshold to keep',
                            condition='smoothMask',
                            help="Only the voxels between these two values "
                            "will be considered to create the meshes.")
        line.addParam('lowLimit', FloatParam, default=0.1, label='Lowest')
        line.addParam('highLimit', FloatParam, default=1, label='Highest')

        group = form.addGroup('Morphological operations', 
                              help='Operations are performed in the same order '
                              'as in the form')
        group.addParam('applyDilation', IntParam, label='Dilation',
                      default=0,
                      help='When set to a positive number, a dilation '
                      'operation is applied with the provided pixel count.')
        group.addParam('applySkeletonization', BooleanParam, label='Skeletonization',
                      default=True,
                      help='When set to yes, a skeletonization operation is applied.')
        
        form.addParam('density',
                      FloatParam,
                      label='Percentage of density ',
                      default=5.0,
                      help='This parameter goes from 0 - 100 and defines the '
                           'percentage of voxel of the tomoMask that '
                           'will be considered as points of the mesh.')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        
        for tsId in self.masks.keys():
            self._insertFunctionStep(self.processTomogramStep, tsId)
        self._insertFunctionStep(self._closeOutputSet)

    # --------------------------- STEPS functions ------------------------------
    def processTomogramStep(self, tsId):
        mask = self._getInputMask(tsId)
        tomogram = self._getInputTomogram(tsId)
        
        if mask is not None and tomogram is not None:
            with self._lock:
                outputMeshes = self._getOutputMeshes()
                maskData = self._loadMask(mask)
        
                if self.smoothMask:
                    self._processSmoothMask(
                        mesh=outputMeshes,
                        tomogram=tomogram,
                        mask=maskData
                    )
                else:
                    self._processSegmentation(
                        mesh=outputMeshes,
                        tomogram=tomogram,
                        segmentation=maskData
                    )

            
                outputMeshes.write()
                self._store()
    
    # --------------------------- UTILS functions ------------------------------
    def _initialize(self):
        self.masks = { mask.getTsId(): mask.clone() for mask in  self.inputMasks.get() }
        self.tomos = { tomo.getTsId(): tomo.clone() for tomo in  self.inputTomograms.get() } 
        self.baseGroupId = Integer(1)
    
    def _getOutputMeshes(self) -> SetOfMeshes:
        output: SetOfMeshes = getattr(self, self._possibleOutputs.Meshes.name, None)
        if output is not None:
            output.enableAppend()
        else:
            output = self._createSetOfMeshes(self.inputTomograms)
            output.setSamplingRate(self.inputMasks.get().getSamplingRate())
            output.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{self._possibleOutputs.Meshes.name: output})
            self._defineSourceRelation(self.inputMasks, output)

        return output
    
    def _getInputMask(self, tsId: str) -> TomoMask:
        return self.masks[tsId]

    def _getInputTomogram(self, tsId: str) -> Tomogram:
        return self.tomos[tsId]
    
    def _loadMask(self, mask: TomoMask) -> np.ndarray:
        ih = ImageHandler()
        image = ih.read(mask)
        return image.getData()

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        
        if self.density.get()>100:
            errors.append('The density must be lesser than 100')
        if self.density.get()<0.0:
            errors.append('The density must be greater than zero')
            
        return errors

    def _processSegmentation(self, 
                             mesh: SetOfMeshes, 
                             tomogram: Tomogram, 
                             segmentation: np.ndarray ) -> None:
        labels = np.unique(segmentation)
        for label in labels:
            if label != self.backgroundLabel.get():
                mask = (segmentation == label)
                self._processBinaryMask(
                    mesh=mesh, 
                    tomogram=tomogram,
                    mask=mask
                )
    
    def _processSmoothMask(self,
                           mesh: SetOfMeshes,
                           tomogram: Tomogram,
                           mask: np.ndarray ) -> None:
        mask = (self.lowLimit.get() <= mask) & (mask <= self.highLimit.get())
        self._processBinaryMask(
            mesh=mesh,
            tomogram=tomogram,
            mask=mask
        )
                
    def _processBinaryMask(self,
                           mesh: SetOfMeshes,
                           tomogram: Tomogram,
                           mask: np.ndarray ):
        if self.applyDilation.get() > 0:
            footprint =skimage.morphology.ball(self.applyDilation.get()) 
            mask = skimage.morphology.binary_dilation(mask, footprint)
            
        if self.applySkeletonization.get():
            mask = skimage.morphology.skeletonize_3d(mask)
        
        
        probability = self.density.get() / 100.0
        coordinates = np.argwhere(mask)
        indices = np.arange(0, len(coordinates))
        nPoints = math.floor(probability*len(indices))
        selection = np.random.choice(indices, size=nPoints, replace=False)
        
        for z, y, x in coordinates[selection,:]:
            point = MeshPoint()
            point.setVolume(tomogram)
            point.setGroupId(self.baseGroupId)
            point.setPosition(x, y, z, const.BOTTOM_LEFT_CORNER)
            mesh.append(point)
        
        self.baseGroupId.increment()
