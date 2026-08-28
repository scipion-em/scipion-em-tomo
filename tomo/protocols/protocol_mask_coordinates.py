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
    """
    Filters a set of 3D coordinates using segmentation masks.

    AI Generated:

    Mask 3D Coordinates (ProtMaskCoordinates) — User Manual
        Overview

        The Mask 3D Coordinates protocol filters a set of 3D coordinates
        according to one or more segmentation masks.

        Its main purpose is to retain only those coordinates that fall
        inside selected segmented regions of a tomogram.

        In cryo-electron tomography workflows, this protocol is useful
        when particle coordinates must be restricted to biologically
        meaningful regions, such as membranes, organelles, protein
        assemblies, or other segmented cellular structures.

        From a biological perspective, this protocol helps remove
        coordinates located in irrelevant background regions and
        preserves only those consistent with prior structural
        annotation.

        Inputs and General Workflow

        The protocol requires two inputs:

            - A SetOfCoordinates3D containing the coordinates to be
              filtered.
            - A SetOfTomoMasks containing segmentation masks associated
              with the corresponding tomograms.

        Each tomogram is processed independently.

        During execution, the protocol creates an output coordinate set
        and keeps it open while each tomogram is processed.

        For every tomogram, the protocol searches for a segmentation
        sharing the same tilt-series identifier. If a matching
        segmentation exists, a binary mask is computed and every
        coordinate is evaluated against that mask.

        Coordinates that satisfy the masking condition are copied into
        the output set.

        Segmentation Labels

        The protocol allows filtering by segmentation label.

        If the segmentation label is negative, all non-zero voxels of
        the segmentation are considered valid.

        This is useful when the segmentation contains only foreground
        and background information, or when all segmented structures are
        biologically relevant.

        If a non-negative label is provided, only voxels matching that
        exact label are considered valid.

        This is especially useful in multi-label segmentations where
        different cellular structures have been assigned distinct
        identifiers.

        Coordinate Evaluation

        For every coordinate, the protocol reads its 3D position in the
        bottom-left-corner reference system.

        The coordinate position is rounded to the nearest voxel index,
        and the corresponding voxel value is checked in the binary mask.

        A coordinate is accepted only if the corresponding voxel belongs
        to the selected segmented region.

        In practical terms, this means that only coordinates physically
        located inside the segmentation are preserved.

        Handling Missing Segmentations

        A tomogram may not always have an associated segmentation.

        The protocol provides two possible behaviors in this situation.

        When excludeUnsegmented is enabled, all coordinates belonging to
        tomograms without a matching segmentation are discarded.

        When excludeUnsegmented is disabled, coordinates from those
        tomograms are copied directly into the output set without
        filtering.

        From a biological point of view, enabling exclusion is more
        conservative and ensures that all output coordinates are
        supported by segmentation information.

        Outputs and Their Interpretation

        The protocol produces a new SetOfCoordinates3D.

        This output preserves the metadata of the original coordinate
        set but contains only coordinates that satisfy the segmentation
        criterion.

        Biologically, the resulting coordinate set represents a subset
        spatially constrained to the segmented structures of interest.

        This filtered set can be directly used in downstream workflows
        such as subtomogram extraction, classification, or targeted
        analysis of specific cellular compartments.

        Validation and Consistency Checks

        Before execution, the protocol verifies that the coordinate set
        and the segmentation set are geometrically compatible.

        Two checks are performed:

            - The sampling rate of coordinates and segmentations must
              match.
            - The tomogram dimensions and segmentation dimensions must
              also match.

        These validations are essential because incorrect sampling or
        mismatched dimensions would lead to biologically meaningless
        spatial filtering.

        Practical Recommendations

        In routine tomography workflows, this protocol is most reliable
        when coordinates and segmentations originate from the same
        tomographic dataset and have not been independently resampled.

        When working with multi-label segmentations, it is good
        practice to carefully verify the biological meaning of the
        selected label before filtering.

        Visual inspection of the filtered coordinates is strongly
        recommended, especially near segmentation boundaries, where
        rounding effects may influence inclusion or exclusion.

        Final Perspective

        For cryo-ET users, this protocol provides a simple but highly
        practical way to incorporate prior segmentation knowledge into
        coordinate analysis.

        Although computationally straightforward, this filtering step is
        often biologically important because it restricts downstream
        analyses to spatially meaningful regions of the tomogram.
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
    