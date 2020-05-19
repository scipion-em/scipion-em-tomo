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

import pyworkflow.protocol.params as params

from .protocol_base import ProtTomoPicking

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

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSubTomos', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label='Input Subtomograms', important=True,
                      help='Select the subtomograms from which you want\n'
                           'to extract the coordinates.')

        form.addParam('inputTomos', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input Tomograms', important=True,
                      help='Select the tomograms to which you want to\n'
                           'associate the coordinates from the subtomograms.')

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

        suffix = ''
        self.outputCoords = self._createSetOfCoordinates3D(inTomos, suffix=suffix)

        def appendCoordFromSubTomo(subTomo, boxSize):
            coord = subTomo.getCoordinate3D()
            tomoKey = coord.getVolId()
            tomo = inTomos[tomoKey]

            if tomo is None:
                print("Skipping subtomogram, key %s not found" % tomoKey)
            else:
                newCoord.copyObjId(subTomo)
                x, y, z = coord.getPosition()
                newCoord.setPosition(x * scale, y * scale, z * scale)

                newCoord.setVolume(tomo)
                newCoord.setBoxSize(boxSize)
                newCoord.setMatrix(coord.getMatrix())
                self.outputCoords.append(newCoord)

        newCoord = Coordinate3D()
        boxSize = inSubTomos.getXDim() * scale
        for subTomo in inSubTomos:
            appendCoordFromSubTomo(subTomo, boxSize)

        self.outputCoords.setBoxSize(boxSize)

    def createOutputStep(self):
        self._defineOutputs(outputCoordinates=self.outputCoords)
        self._defineSourceRelation(self.inputSubTomos, self.outputCoords)
        self._defineSourceRelation(self.inputTomos, self.outputCoords)

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

        if hasattr(self, 'outputCoordinates'):
            summary.append('Output coordinates: *%d*'
                           % self.outputCoordinates.getSize())

        return summary

    def _methods(self):
        # No much to add to summary information
        return self._summary()

    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of errors.
        If the list is empty the protocol can be executed.
        """
        errors = []
        inputSubTomos = self.getInputSubTomos()
        first = inputSubTomos.getFirstItem()
        if first.getCoordinate3D() is None:
            errors.append('The input particles do not have coordinates!!!')

        return errors