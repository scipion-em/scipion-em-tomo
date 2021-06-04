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

import numpy as np
import os

from pyworkflow import BETA
import pyworkflow.protocol.params as params

import pyworkflow.utils as pwutils
from pyworkflow.object import Integer

from .protocol_base import ProtTomoPicking

import tomo.constants as const
from ..objects import Coordinate3D


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

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSubTomos', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label='Subtomograms', important=True,
                      help='Select the subtomograms from which you want\n'
                           'to extract the coordinates. The coordinate belonging to '
                           'each subtomogram should be already associated to an initial '
                           'tomogram.')

        form.addParam('inputTomos', params.PointerParam,
                      pointerClass='SetOfTomograms',
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
        self._insertFunctionStep('extractCoordinatesStep')
        self._insertFunctionStep('createOutputStep')

    def extractCoordinatesStep(self):
        inTomos = self.getInputTomos()
        inSubTomos = self.getInputSubTomos()
        scale = inSubTomos.getSamplingRate() / inTomos.getSamplingRate()
        print("Scaling coordinates by a factor *%0.2f*" % scale)
        filesTomo = [pwutils.removeBaseExt(tomo.getFileName()) for tomo in inTomos.iterItems()]

        suffix = ''
        self.outputCoords = self._createSetOfCoordinates3D(inTomos, suffix=suffix)
        self.outputCoords.setSamplingRate(inTomos.getSamplingRate())

        def appendCoordFromSubTomo(subTomo, boxSize):
            coord = subTomo.getCoordinate3D()
            tomoKey = coord.getVolId()
            tomo = inTomos[tomoKey]

            if tomo is None:
                print("Key %s not found, trying to associate tomogram using filename" % tomoKey)
                try:
                    idx = filesTomo.index(pwutils.removeBaseExt(subTomo.getVolName()))
                except:
                    idx = None
                if idx is not None:
                    if len(filesTomo) == 1:
                        newCoord.setVolume(inTomos.getFirstItem())
                        coord.setVolume(inTomos.getFirstItem())
                    else:
                        newCoord.setVolume(inTomos[idx+1])
                        coord.setVolume(inTomos[idx+1])
                    x, y, z = coord.getPosition(const.SCIPION)
                    newCoord.copyObjId(subTomo)

                    newCoord.setPosition(x * scale, y * scale, z * scale, const.SCIPION)
                    newCoord.setBoxSize(boxSize)
                    newCoord.setMatrix(checkMatrix(subTomo, coord))
                    if coord.hasGroupId():
                        newCoord.setGroupId(coord.getGroupId())
                    self.outputCoords.append(newCoord)
            else:
                coord.setVolume(tomo)
                newCoord.copyObjId(subTomo)
                x, y, z = coord.getPosition(const.SCIPION)
                newCoord.setVolume(tomo)
                newCoord.setPosition(x * scale, y * scale, z * scale, const.SCIPION)

                newCoord.setBoxSize(boxSize)
                newCoord.setMatrix(checkMatrix(subTomo, coord))
                self.outputCoords.append(newCoord)

        def checkMatrix(subTomo, coord):
            transform_subTomo = subTomo.getTransform().getMatrix()
            transform_coordinate = coord.getMatrix()
            if not np.allclose(transform_coordinate, np.eye(transform_coordinate.shape[0])):
                return transform_coordinate
            elif not np.allclose(transform_subTomo, np.eye(transform_subTomo.shape[0])):
                return transform_subTomo
            else:
                return np.eye(transform_subTomo.shape[0])

        newCoord = Coordinate3D()
        if self.boxSize.get() is None:
            boxSize = inSubTomos.getXDim() * scale
        else:
            boxSize = self.boxSize.get()
        for subTomo in inSubTomos:
            appendCoordFromSubTomo(subTomo, boxSize)

        self.outputCoords.setBoxSize(boxSize)

    def createOutputStep(self):
        if self.outputCoords.getSize() > 0:
            self._defineOutputs(outputCoordinates3D=self.outputCoords)
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

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        ps1 = self.getInputSubTomos().getSamplingRate()
        ps2 = self.getInputTomos().getSamplingRate()
        summary.append(u'Input subtomograms pixel size: *%0.3f* (Å/px)' % ps1)
        summary.append(u'Input tomograms pixel size: *%0.3f* (Å/px)' % ps2)
        summary.append('Scaling coordinates by a factor of *%0.3f*' % (ps1 / ps2))

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
        if first.getCoordinate3D() is None:
            errors.append('The input particles do not have coordinates!!!')

        return errors