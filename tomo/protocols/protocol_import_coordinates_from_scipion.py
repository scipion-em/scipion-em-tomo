# *
# * Authors:     Scipion Team
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************

from os.path import basename, exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FileParam, IntParam, PointerParam
from pyworkflow.utils import copyFile, Message, removeBaseExt
from ..constants import SCIPION
from .protocol_base import ProtTomoBase
from ..objects import SetOfTomograms


class ProtImportCoordinates3DFromScipion(EMProtocol, ProtTomoBase):
    """Protocol to import a set of tomograms to the project"""

    _label = 'import set of coordinates 3D from scipion sqlite file'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('sqliteFile', FileParam,
                      label='Scipion sqlite file')
        form.addParam('inTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input tomograms',
                      help='Select the tomograms to which the coordinates should be referred to. '
                           'The matching between coordinates and tomograms is made checking the tsId/tomoId '
                           'attribute. If no matches are found, then it tries to do it comparing the filenames. '
                           '*IMPORTANT*: the coordinates will be assumed to be at the same sampling rate as the '
                           'introduced tomograms.')
        form.addParam('boxSize', IntParam,
                      label='Box Size [pix]',
                      default=20)

    def _insertAllSteps(self):
        # JORGE
        import os
        fname = "/home/jjimenez/Desktop/test_JJ.txt"
        if os.path.exists(fname):
            os.remove(fname)
        fjj = open(fname, "a+")
        fjj.write('JORGE--------->onDebugMode PID {}'.format(os.getpid()))
        fjj.close()
        print('JORGE--------->onDebugMode PID {}'.format(os.getpid()))
        import time
        time.sleep(10)
        # JORGE_END
        self._insertFunctionStep(self.importCoordinatesStep)

    # --------------------------- STEPS functions -----------------------------

    def importCoordinatesStep(self):
        # Make a copy of the introduced file in the current protocol's results folder
        sqliteFile = self.sqliteFile.get()
        newFileName = self._getPath(basename(sqliteFile))
        copyFile(sqliteFile, newFileName)

        # Generate a set of 3d coordinates and assign the mapper of the introduced sqlite file
        inTomoSet = self.importTomograms.get()
        coordsSet = self._createSetOfCoordinates3D(inTomoSet)
        coordsSet.setSamplingRate(inTomoSet.getSamplingRate())
        coordsSet.setBoxSize(self.boxSize.get())
        coordsSet._mapperPath.set('%s, %s' % (newFileName, ''))
        coordsSet.load()

        # Check if the coordinates and the tomograms can be related via the tomoId or the filename
        outTomoSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        outTomoSet.copyInfo(inTomoSet)
        self._getMatchingTomos(inTomoSet, outTomoSet, coordsSet)

        # Define the outputs
        self._defineOutputs(outputSetOfCoordinates=coordsSet)
        self._defineSourceRelation(inTomoSet, coordsSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errorList = []
        if not exists(self.sqliteFile.get()):
            errorList.append('Introduced file was not found:\n\t%s' % self.sqliteFile.get())

        return errorList

    def _summary(self):
        pass

    # --------------------------- UTILS functions ----------------------------

    @staticmethod
    def _getMatchingTomos(inTomoSet, outTomoSet, coordsSet):
        notFoundCoords = []
        inTomoSetMatchingIndices = []
        pattern = 'Row %i - (x, y, x) = (%.2f, %.2f, %.2f)'
        tomoTsIdList, tomoBaseNameList = zip(*[(tomo.getTsId(), removeBaseExt(tomo.getFileName()))
                                               for tomo in inTomoSet])
        for coord in coordsSet:
            coordTomoId = coord.getTomoId()
            if coordTomoId in tomoTsIdList:
                indByTomoId = tomoTsIdList.index(coordTomoId) + 1
                coord.setVolume(inTomoSet[indByTomoId])
                inTomoSetMatchingIndices.append(indByTomoId)
            else:
                indexByName = ProtImportCoordinates3DFromScipion._getMatchingIndexByFileName(coordTomoId,
                                                                                             tomoBaseNameList)
                if indexByName:
                    coord.setVolume(inTomoSet[indexByName])
                    inTomoSetMatchingIndices.append(indexByName)
                else:
                    notFoundCoords.append(pattern % (coord.getObjId(), *coord.getPosition(SCIPION)))

        # Build a precedents set with only the matching tomograms, in case there are not all the ones present in the
        # input set
        if inTomoSetMatchingIndices:
            for ind in set(inTomoSetMatchingIndices):
                outTomoSet.append(inTomoSet[ind])

        # TODO: register thisoutput if necessary

        return notFoundCoords

    @staticmethod
    def _getMatchingIndexByFileName(coordTomoId, tomoBaseNameList):
        matchingIndex = None
        matches = list(map(lambda x: coordTomoId in x, tomoBaseNameList))
        if any(matches):
            matchingIndex = matches.index(True) + 1

        return matchingIndex
