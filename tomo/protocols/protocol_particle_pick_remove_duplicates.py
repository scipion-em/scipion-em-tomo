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
"""
Consensus picking protocol
"""

import numpy as np

import pyworkflow.protocol.params as params
from pyworkflow.utils import getFiles, removeBaseExt

from .protocol_particle_pick_consensus import (ProtTomoConsensusPicking,
                                               consensusWorker, getReadyTomos)


class ProtTomoPickingRemoveDuplicates(ProtTomoConsensusPicking):
    """
    This protocol removes coordinates that are closer than a given threshold.
    The remaining coordinate is the average of the previous ones.
    """

    _label = 'remove duplicates'
    outputName = 'outputCoordinates'
    FN_PREFIX = 'purgedCoords_'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input coordinates", important=True,
                      help='Select the set of 3D coordinates to compare')
        form.addParam('consensusRadius', params.IntParam, default=10,
                      label="Radius",
                      help="All coordinates within this radius (in pixels) "
                           "are presumed to correspond to the same particle")

        # FIXME: It's not using more than one since
        #         self.stepsExecutionMode = STEPS_SERIAL
        # form.addParallelSection(threads=4, mpi=0)

#--------------------------- INSERT steps functions ----------------------------
    def insertNewCoorsSteps(self, tomos):
        deps = []
        for tomo in tomos:
            stepId = self._insertFunctionStep("removeDuplicatesStep",
                                              tomo.getObjId(),
                                              tomo.getFileName(),
                                              prerequisites=[])
            deps.append(stepId)
        return deps

    def _checkNewInput(self):
        # If continue from an stopped run, don't repeat what is done
        if not self.checkedTomos:
            for fn in getFiles(self._getExtraPath()):
                fn = removeBaseExt(fn)
                if fn.startswith(self.FN_PREFIX):
                    self.checkedTomos.update([self.getTomoId(fn)])
                    self.processedTomos.update([self.getTomoId(fn)])

        readyTomos, self.streamClosed = getReadyTomos(self.inputCoordinates.get())

        newTomosIds = readyTomos.difference(self.checkedTomos)

        if newTomosIds:
            self.checkedTomos.update(newTomosIds)

            inTomos = self.getMainInput().getPrecedents()
            newTomos = [inTomos[tomoId].clone() for tomoId in newTomosIds]

            fDeps = self.insertNewCoorsSteps(newTomos)
            outputStep = self._getFirstJoinStep()
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def getMainInput(self):
        return self.inputCoordinates.get()

    def defineRelations(self, outputSet):
        self._defineTransformRelation(self.getMainInput(), outputSet)

    def removeDuplicatesStep(self, tomoId, tomoName):
        print("Removing duplicates for tomogram %d: '%s'"
              % (tomoId, tomoName))

        coordArray = np.asarray([x.getPosition() for x in
                                 self.getMainInput().iterCoordinates(tomoId)],
                                 dtype=int)

        consensusWorker([coordArray], 1, self.consensusRadius.get(),
                        self._getTmpPath('%s%s.txt' % (self.FN_PREFIX, tomoId)))

        self.processedTomos.update([tomoId])

    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        return ["Radius = %d" % self.consensusRadius]
