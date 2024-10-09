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
import scipy.ndimage

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, EnumParam

from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.objects import Set

from tomo.objects import (SetOfCoordinates3D, Coordinate3D, SetOfTomograms,
                          SetOfTomoMasks, TomoMask, Tomogram )
from tomo.protocols import ProtTomoBase
import tomo.constants as const

COORDINATES = 'Coordinates'

MORPHOLOGICAL_OP_NONE = 0
MORPHOLOGICAL_OP_DILATION = 1
MORPHOLOGICAL_OP_EROSION = 2
MORPHOLOGICAL_OP_OPENING = 3
MORPHOLOGICAL_OP_CLOSING = 4

class ProtMaskCoordinates(EMProtocol, ProtTomoBase):
    """ TODO
    """
    _label = 'mask 3d coordinates'
    _devStatus = BETA
    _possibleOutputs = {COORDINATES: SetOfCoordinates3D}

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', 
                      PointerParam, pointerClass=SetOfCoordinates3D,
                      label='Input coordinates', important=True,
                      help='Select the Coordinates3D to be filtered')
        form.addParam('inputSegmentations', 
                      PointerParam, pointerClass=SetOfTomoMasks,
                      label='Input segmentation', important=True,
                      help='Select the tomo mask used for filtering coordinates')
        form.addParam('segmentationLabel', IntParam, label='Segmentation label',
                      default=-1,
                      help='Labels to consider. If negative, it will consider '
                           'all non-zero areas in the input segmentation')
        form.addParam('morphologicalOperation', EnumParam, 
                      label='Morphological operation',
                      choices=['None', 'Dilation', 'Erosion', 'Opening', 'Closing'],
                      help='Apply a morphological operation to the binarized mask')
        form.addParam('morphologicalOperationSphRadius', IntParam,
                      label='Sphere size', default = 10,
                      condition='morphologicalOperation != MORPHOLOGICAL_OP_NONE',
                      help='Radiuis of the sphere used for the morphological '
                           'operation. In pixels')
        
    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        coordinates = self._getInputCoordinates()
        tomograms: SetOfTomograms = coordinates.getPrecedents()
        
        self._insertFunctionStep(self.createOutputStep)
        for tomogram in tomograms.iterItems():
            tsId = tomogram.getTsId()
            self._insertFunctionStep(self.filterTomogramCoordinatesStep, tsId)
        self._insertFunctionStep(self.finalizeStep)
        

    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        inputCoordintates = self._getInputCoordinates()
        tomograms = inputCoordintates.getPrecedents()
        outputCoordinates = self._createSetOfCoordinates3D(tomograms)
        outputCoordinates.copyInfo(inputCoordintates)
        outputCoordinates.setStreamState(Set.STREAM_OPEN)
        
        self._defineOutputs(**{COORDINATES: outputCoordinates})
        self._defineSourceRelation(self.inputCoordinates, outputCoordinates)
        self._defineSourceRelation(self.inputSegmentations, outputCoordinates)
    
    def filterTomogramCoordinatesStep(self, tsId: str):
        mask = self._calculateMask(tsId)
        inputCoordinates3d = self._getInputCoordinates()
        outputCoordinates: SetOfCoordinates3D = self.outputCoordinates
        tomogram: Tomogram = inputCoordinates3d.getPrecedent(tsId)
        
        for item in inputCoordinates3d.iterCoordinates(tomogram):
            if self._checkCoordinate(item, mask):
                outputCoordinates.append(item)

    def closeOuputStep(self):
        self._closeOutputSet()

    # --------------------------- UTILS functions ------------------------------
    def _getInputCoordinates(self) -> SetOfCoordinates3D:
        return self.inputCoordinates.get()
    
    def _getInputSegmentation(self, tsId: str) -> TomoMask:
        segmentations: SetOfTomoMasks = self.inputSegmentations.get()
        for item in segmentations.iterItems():
            if item.getTsId() == tsId:
                return item
            
        raise RuntimeError(f'No segmentation found for {tsId}')
    
    def _loadSegmentation(self, tsId) -> np.ndarray:
        segmentation = self._getInputSegmentation(tsId)
        ih = ImageHandler()
        image = ih.read(segmentation)
        return image.getData()
    
    def _applyMorphologicalOperation(self, mask: np.ndarray) -> np.ndarray:
        opCode = self.morphologicalOperation.get()
        if opCode != MORPHOLOGICAL_OP_NONE:
            radius: int = self.morphologicalOperationSphRadius.get()
            struct = scipy.ndimage.generate_binary_structure(3, 1)
            struct = scipy.ndimage.iterate_structure(struct, iterations=radius)
            
            if opCode == MORPHOLOGICAL_OP_DILATION:
                mask = scipy.ndimage.binary_dilation(mask, struct)
            elif opCode == MORPHOLOGICAL_OP_EROSION:
                mask = scipy.ndimage.binary_erosion(mask, struct)
            elif opCode == MORPHOLOGICAL_OP_OPENING:
                mask = scipy.ndimage.binary_opening(mask, struct)
            elif opCode == MORPHOLOGICAL_OP_CLOSING:
                mask = scipy.ndimage.binary_closing(mask, struct)
            else:
                raise RuntimeError('Unknown morphological operation')
            
        return mask
            
    def _calculateMask(self, tsId: str) -> np.ndarray:
        segmentation = self._loadSegmentation(tsId)
        label: int = self.segmentationLabel.get()
        
        if label < 0:
            mask = (segmentation > 0)
        else:
            mask = (segmentation == label)
            
        return self._applyMorphologicalOperation(mask)
        
    def _checkCoordinate(self, coordinate: Coordinate3D, mask: np.ndarray) -> bool:
        x, y, z = coordinate.getPosition(const.SCIPION)
        x = round(x - mask.shape[-1]/2)
        y = round(y - mask.shape[-2]/2)
        z = round(z - mask.shape[-3]/2)
        return bool(mask[z,y,x] == True)
    
    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        return [] # TODO