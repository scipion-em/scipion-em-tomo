# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
import enum
import numpy as np
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from .protocol_base import ProtTomoPicking
import tomo.constants as const
from ..objects import Coordinate3D, SetOfCoordinates3D, SetOfSubTomograms, SetOfTomograms, SubTomogram


class Output3dCoordExtraction(enum.Enum):
    coordinates3d = SetOfCoordinates3D


class ProtTomoExtractCoords(ProtTomoPicking):
    """
    Extract the coordinates information from a set of subtomograms.

    This protocol is useful when we want to re-extract the subtomograms
    (maybe resulting from classification) with the
    original dimensions. It can be also handy to visualize the resulting
    subtomograms in their location on the tomograms.
    """

    _label = 'extract 3D coordinates'
    _devStatus = BETA
    _possibleOutputs = Output3dCoordExtraction

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
        self.outputCoords = None
        self._tomoDict = None
        self._inputAreSubtomos = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSubTomos', params.PointerParam,
                      pointerClass=[SetOfSubTomograms, SetOfCoordinates3D],
                      label='Subtomograms or 3D coordinates', important=True,
                      help='Select the subtomograms from which you want\n'
                           'to extract the coordinates. The coordinate belonging to '
                           'each subtomogram should be already associated to an initial '
                           'tomogram.')

        form.addParam('inputTomos', params.PointerParam,
                      pointerClass=SetOfTomograms,
                      label='Tomograms', important=True,
                      help='Select the tomograms to which you want to\n'
                           'associate the coordinates from the subtomograms.')

        form.addParam('boxSize', params.IntParam,
                      allowsNull=True, expertLevel=params.LEVEL_ADVANCED, label='Box Size',
                      help='Determine the box size of the extracted coordinates. By default, '
                           'the program assigns the box size directly from the coordinates '
                           'associated to the subtomograms.')

        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.extractCoordinatesStep)
        self._insertFunctionStep(self.createOutputStep)

    def getTomogramFromItem(self, item):
        """ Returns the tomogram associated with the item. Item could be either a subtomogram or a 3D coordinate."""

        tomoDict = self.getTomogramDictionary()

        # get the coordinate
        coord = self.getCoordFromItem(item)

        # If no coordinate, we are dealing with imported subtomograms
        if coord is None:

            # From file name and subtomo:
            file = pwutils.removeBaseExt(item.getVolName())

            tomo = tomoDict[file]

        # There are coordinates
        else:
            tomo = tomoDict.get(coord.getTomoId(), None)

        # last resource: vol identifier.
        if tomo is None:
            # Try by vol id
            tomo = tomoDict[item.getVolId()]

        return tomo

    def getTomogramDictionary(self):
        """ Returns a dictionary of tomogram where the key is any of the
        possible tomogram identifiers to do the matching"""

        if not self._tomoDict:
            self._tomoDict = dict()
            # iterate over the SetOfTomograms
            for tomo in self.getInputTomos().iterItems():
                # Clone the tomogram
                tomoClone = tomo.clone()

                # Add the tomogram based on the filename
                filename = pwutils.removeBaseExt(tomo.getFileName())
                self._tomoDict[filename] = tomoClone

                # Add it by tilt series id
                self._tomoDict[tomoClone.getTsId()] = tomoClone

                # Add it by row identifier (the weakest due to join sets renumbering identifiers)
                self._tomoDict[tomoClone.getObjId()] = tomoClone

        return self._tomoDict

    def extractCoordFromItem(self, item, boxSize, scaleCoords, scaleShifts):
        coord = self.getCoordFromItem(item)
        tomo = self.getTomogramFromItem(item)

        if tomo is None:
            self.warning("Tomogram not found for %s" % item)
            return None
        else:
            newCoord = Coordinate3D()

            coord.setVolume(tomo)
            newCoord.copyObjId(coord)
            x, y, z = coord.getPosition(const.SCIPION)
            newCoord.setVolume(tomo)
            newCoord.setPosition(x * scaleCoords, y * scaleCoords, z * scaleCoords, const.SCIPION)

            newCoord.setBoxSize(boxSize)
            transformation = self.checkMatrix(item)
            transformation[0, 3] *= scaleShifts
            transformation[1, 3] *= scaleShifts
            transformation[2, 3] *= scaleShifts
            newCoord.setMatrix(transformation)
            if coord.hasGroupId():
                newCoord.setGroupId(coord.getGroupId())

            return newCoord

    def areInputSubtomos(self):
        """
        Returns true if input re Subtomograms. (Lazy loaded)
        :param item: an item of the inputset
        :return:
        """
        if self._inputAreSubtomos is None:
            self._inputAreSubtomos = isinstance(self.getInputSubTomos(), SetOfSubTomograms)

        return self._inputAreSubtomos

    def getCoordFromItem(self, item):
        """ Returns the Coordinate 3D from the item"""
        if self.areInputSubtomos():
            coord = item.getCoordinate3D()
        else:
            coord = item
        return coord

    @classmethod
    def checkMatrix(cls, item):
        """ Returns the matrix of item (subtomo or coordinate3D"""
        transform = item.getTransform().getMatrix() if isinstance(item, SubTomogram) else item.getMatrix()
        if transform is not None:
            return transform
        else:
            return np.eye(transform.shape[0])

    def extractCoordinatesStep(self):
        """What must be considered here:

            1. If the input is a set of subtomograms:
                1.1. The coordinates associated will be at the scale of the tomograms they were picked from.
                1.2. The shifts of each subtomogram transformation matrix will be scaled properly.
            2. If the input is a set of 3d coordinates: both the coordinates and the shifts will be at the scale
               of the tomograms they were picked from.

        Thus, two scale factors will be necessary to considerate all the possible cases: one for the shifts and
        another for the coordinates. They will be equal in the case of introducing a set of 3d coordinates."""
        inTomos = self.getInputTomos()
        inSubTomos = self.getInputSubTomos()
        inCoords = self.getCoordinates()

        inTomosSRate = inTomos.getSamplingRate()
        scaleCoords = inCoords.getSamplingRate() / inTomosSRate
        scaleShifts = inSubTomos.getSamplingRate() / inTomosSRate

        # Create the output set of coordinates
        self.outputCoords = self._createSetOfCoordinates3D(inTomos)
        self.outputCoords.setSamplingRate(inTomosSRate)

        if self.boxSize.get() is None:

            if self.areInputSubtomos():
                boxSize = inSubTomos.getXDim() * scaleCoords
            else:
                boxSize = inSubTomos.getBoxSize() * scaleCoords
        else:
            boxSize = self.boxSize.get()

        # For each item (Subtomo or coordinate) in the input
        for item in inSubTomos:
            newCoord = self.extractCoordFromItem(item, boxSize, scaleCoords, scaleShifts)
            if newCoord:
                self.outputCoords.append(newCoord)

        self.outputCoords.setBoxSize(boxSize)

    def createOutputStep(self):
        if self.outputCoords.getSize() > 0:
            self._defineOutputs(**{Output3dCoordExtraction.coordinates3d.name: self.outputCoords})
            self._defineSourceRelation(self.inputSubTomos, self.outputCoords)
            self._defineSourceRelation(self.inputTomos, self.outputCoords)
        else:
            raise Exception("No coordinates were extracted from the input subtomograms probably "
                            "due to an issue during the association with the new tomograms. In case "
                            "the association was done by the filename, please, check that the tomograms where "
                            "subtomograms were extracted and new tomograms have the same file names "
                            "and try again.")

    # ------------- UTILS functions ----------------
    def getSuffix(self, suffix):
        return "_tmp%s" % suffix

    def getTmpOutputPath(self, suffix):
        return self._getPath("coordinates%s.sqlite" % self.getSuffix(suffix))

    def getInputTomos(self):
        return self.inputTomos.get()

    def getInputSubTomos(self):
        return self.inputSubTomos.get()

    def getCoordinates(self):
        if self.areInputSubtomos():
            return self.getInputSubTomos().getCoordinates3D()
        else:
            return self.getInputSubTomos()

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        # ps1 = self.getCoordinates().getSamplingRate()
        # ps2 = self.getInputTomos().getSamplingRate()
        # summary.append(u'Input subtomograms pixel size: *%0.3f* (Å/px)' % ps1)
        # summary.append(u'Input tomograms pixel size: *%0.3f* (Å/px)' % ps2)
        # summary.append('Scaling coordinates by a factor of *%0.3f*' % (ps1 / ps2))

        if hasattr(self, 'outputCoordinates3D'):
            summary.append('Output coordinates: *%d*'
                           % self.outputCoordinates3D.getSize())

        return summary

    def _methods(self):
        return self._summary()

    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of errors.
        """
        errors = []
        inputSubTomos = self.getInputSubTomos()
        first = inputSubTomos.getFirstItem()
        if self.getCoordFromItem(first) is None:
            errors.append('The input particles do not have coordinates!!!')

        return errors
