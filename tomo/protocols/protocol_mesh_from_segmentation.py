# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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

from pwem.protocols import EMProtocol

import numpy as np
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam
from tomo.objects import SetOfMeshes, MeshPoint, SetOfTomoMasks
import tomo.constants as const

from tomo.protocols import ProtTomoBase
import mrcfile


class ProtMeshFromSegmentation(EMProtocol, ProtTomoBase):
    """
    Creates meshes based on segmentations (TomoMasks) voxels values.
    """
    _label = 'meshes from tomoMask'
    _devStatus = BETA
    # Output name
    _OUTPUT_NAME = 'Meshes'
    _possibleOutputs = {_OUTPUT_NAME: SetOfMeshes}
    listOfMeshCoords = {}

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Parameters')

        form.addParam('tomoMasks',
                      PointerParam,
                      important=True,
                      pointerClass=SetOfTomoMasks,
                      label='Tomo Masks',
                      help='Set of tomo mask from which the meshes will be created')

        line = form.addLine('Threshold to keep',
                             help="Only the voxels between these two values will be considered to create the meshes.")

        line.addParam('lowLimit', FloatParam, default=0.1, label='Lowest')
        line.addParam('highLimit', FloatParam, default=1, label='Highest')

        form.addParam('density',
                      FloatParam,
                      label='Percentage of density ',
                      default=5.0,
                      help='This parameter goes from 0 - 100 and defines the percentage of voxel of the tomoMask that'
                           'will be considered as points of the mesh.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        for tm in self.tomoMasks.get():
            tsId = tm.getTsId()
            self._insertFunctionStep(self.createMeshesStep, tsId)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------

    def createMeshesStep(self, tsId):

        settm = self.tomoMasks.get()
        tomoMask = settm[{'_tsId': tsId}]

        if tomoMask is None:
            self.warning('TomoMask not found for tsId %s' % tsId)
            return

        tomoFn = tomoMask.getFileName()
        mytomo = mrcfile.read(tomoFn)

        candidates = (self.lowLimit.get() <= mytomo) & (mytomo <= self.highLimit.get())

        coordinates = np.argwhere(candidates)

        segmentedPoints = len(coordinates)
        idxArray = np.arange(0, segmentedPoints)

        # Randomize the order of the array and get the permutation indices
        idxArrayPermutations = np.random.permutation(idxArray)

        rho = self.density.get()/100.0
        self.listOfMeshCoords[tsId] = coordinates[idxArrayPermutations[0:int(np.floor(segmentedPoints*rho))]]


    def createOutputStep(self):
        inputTomos = self.tomoMasks.get()
        outSet = self._createSetOfMeshes(inputTomos)
        sampling = inputTomos.getSamplingRate()
        outSet.setSamplingRate(sampling)
        outSet.setBoxSize(10)
        grId = 1
        for tomo in inputTomos.iterItems():
            tsId = tomo.getTsId()
            listOfCoords = self.listOfMeshCoords[tsId]
            x = listOfCoords[:,2]
            y = listOfCoords[:,1]
            z = listOfCoords[:,0]

            for idx in range(0,len(x)):
                mesh = MeshPoint()
                mesh.setVolume(tomo)
                mesh.setPosition(x[idx], y[idx], z[idx], const.BOTTOM_LEFT_CORNER)
                mesh.setGroupId(grId)
                outSet.append(mesh)
            grId += 1
        outSet.setPrecedents(inputTomos)
        self._defineOutputs(**{self._OUTPUT_NAME:outSet})
        self._defineSourceRelation(inputTomos, outSet)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        message = 'This is a Scipiontomo method'
        return [message]

    def _validate(self):
        errors = []
        if self.density.get()>100:
            errors.append('The density must be lesser than zero')
        if self.density.get()<0.0:
            errors.append('The density must be greater than zero')
        return errors

    def _summary(self):
        summary = ['A set of meshes has been created']
        return summary
