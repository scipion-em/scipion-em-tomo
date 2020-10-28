# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             Tomas Majtner (tmajtner@cnb.csic.es)  -- streaming version
# *             David Herreros Calero (dherreros@cnb.csic.es) -- Tomo version
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************

import numpy as np

import pyworkflow.protocol.params as params

import pyworkflow.utils as pwutils

from tomo.protocols import ProtTomoPicking
from tomo.utils import delaunayTriangulation, normalFromMatrix

class ProtTomoScoreCoordinates(ProtTomoPicking):
    '''Scoring and (optional) filtering of coordinates based on different scoring
    functions (normals angle, carbon distance, neighbour distance)'''

    _label = 'score/filter coordinates'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input 3D coordinates", important=True,
                      help='Select the set of 3D coordinates to compare')
        form.addParam('outliers', params.BooleanParam, default=True,
                      label="Score outluiers?")
        form.addParam('outliersThreshold', params.FloatParam, default=10,
                      label="Outliers distance threshold", condition='outliers == True',
                      help='Distance expressed in pixels')
        form.addParam('angle', params.BooleanParam, default=True,
                      label="Score normals?")
        form.addParam('angleThreshold', params.FloatParam, default=0,
                      label="Angle threshold", condition='angle == True',
                      help='Angle expressed in ***')
        form.addParam('carbon', params.BooleanParam, default=True,
                      label="Score carbon closeness?")
        form.addParam('carbonThreshold', params.FloatParam, default=10,
                      label="Carbon distance threshold", condition='outliers == True',
                      help='Distance expressed in pixels')

    def _insertAllSteps(self):
        coordinates = self.inputCoordinates.get()
        self._insertFunctionStep('computeParams', coordinates)

    def computeParams(self, coordinates):
        import time
        time.sleep(10)
        tomos = coordinates.getPrecedents()
        tomoNames = [pwutils.removeBaseExt(tomo.getFileName()) for tomo in tomos]
        vesicleIds = set([coord._vesicleId.get() for coord in coordinates.iterCoordinates()])
        axis = np.array([0, 0, 1])
        tomo_vesicles = {tomoField: {'vesicles': [], 'normals': []}
                         for tomoField in tomoNames}

        for idt, tomo in enumerate(tomos.iterItems()):
            for idv in vesicleIds:
                vesicle = []
                normals = []
                for coord in coordinates.iterCoordinates(volume=tomo):
                    if coord._vesicleId == idv:
                        vesicle.append(coord.getPosition())
                        trMat = coord.getMatrix()
                        rotMat = np.linalg.inv(trMat[:3, :3])
                        normals.append(rotMat.dot(axis))
                tomo_vesicles[tomoNames[idt]]['vesicles'].append(np.asarray(vesicle))
                tomo_vesicles[tomoNames[idt]]['normals'].append(np.asarray(normals))
        print('Done')