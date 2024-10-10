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

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam

from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.objects import Set

from tomo.objects import (SetOfCoordinates3D, Coordinate3D, SetOfTomograms,
                          SetOfTomoMasks, TomoMask, Tomogram)
from tomo.protocols import ProtTomoBase
import tomo.constants as const

COORDINATES = 'Coordinates'

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
        form.addParam('excludeUnsegmented', BooleanParam, 
                      label='Exclude unsegmented coordinates', default=True,
                      help='Determines behaviour when encountering coordinates '
                           'without a segmentation. When true, those coordintates ' 
                           'are not outputed. If false, all coordintates from '
                           'those coordinates are outputed.')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        coordinates = self._getInputCoordinates()
        tomograms: SetOfTomograms = coordinates.getPrecedents()
        
        self._insertFunctionStep(self.createOutputStep)
        for tomogram in tomograms.iterItems():
            tsId = tomogram.getTsId()
            self._insertFunctionStep(self.filterTomogramCoordinatesStep, tsId)
        self._insertFunctionStep(self.closeOuputStep)
        
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
        inputCoordinates3d = self._getInputCoordinates()
        outputCoordinates: SetOfCoordinates3D = getattr(self, COORDINATES)
        tomogram: Tomogram = inputCoordinates3d.getPrecedent(tsId)
        segmentation = self._getInputSegmentation(tsId)
        
        if segmentation is not None:
            mask = self._calculateMask(segmentation)
            
            for item in inputCoordinates3d.iterCoordinates(tomogram):
                if self._checkCoordinate(item, mask):
                    outputCoordinates.append(item)
        
        elif not self.excludeUnsegmented:
            for item in inputCoordinates3d.iterCoordinates(tomogram):
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
        return None
    
    def _loadSegmentation(self, segmentation: TomoMask) -> np.ndarray:
        ih = ImageHandler()
        image = ih.read(segmentation)
        return image.getData()
    
    def _calculateMask(self, segmentation: TomoMask) -> np.ndarray:
        segmentationData = self._loadSegmentation(segmentation)
        label: int = self.segmentationLabel.get()
        
        if label < 0:
            mask = (segmentationData > 0)
        else:
            mask = (segmentationData == label)
            
        return mask
        
    def _checkCoordinate(self, coordinate: Coordinate3D, mask: np.ndarray) -> bool:
        position = coordinate.getPosition(const.BOTTOM_LEFT_CORNER)
        position = tuple(map(round, position))
        x, y, z = position
        return bool(mask[z, y, x] == True)
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        result = []
        
        coordinates = self._getInputCoordinates()
        tomograms: SetOfTomograms = coordinates.getPrecedents()
        segmentations: SetOfTomoMasks = self.inputSegmentations.get()

        if coordinates.getSamplingRate() != segmentations.getSamplingRate():
            result.append('Sampling rate of the segmentation does not '
                          'match the sampling rate of the segmentation')
        
        if tomograms.getDim() != segmentations.getDim():
            result.append('Dimensions of the segmentation does not match '
                          'the dimensions of the tomogram used for picking')
        
        return result
    