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
from pyworkflow.object import String
from pyworkflow.protocol import FileParam, IntParam, PointerParam
from pyworkflow.utils import copyFile, Message, removeBaseExt
from ..constants import SCIPION, ERR_COORDS_FROM_SQLITE_NO_MATCH
from .protocol_base import ProtTomoBase
from ..objects import SetOfTomograms, SetOfCoordinates3D


class ProtImportCoordinates3DFromScipion(EMProtocol, ProtTomoBase):
    """Protocol to import a set of tomograms to the project"""

    _label = 'import set of coordinates 3D from scipion sqlite file'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.notFoundCoordsMsg = None

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('sqliteFile', FileParam,
                      label='Scipion sqlite file')
        form.addParam('importTomograms', PointerParam,
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
        self._insertFunctionStep(self.importCoordinatesStep)

    # --------------------------- STEPS functions -----------------------------

    def importCoordinatesStep(self):
        inTomoSet = self.importTomograms.get()
        coordsSet = SetOfCoordinates3D()

        # Make a copy of the introduced file in the current protocol's results folder
        sqliteFile = self.sqliteFile.get()
        newFileName = self._getPath(basename(sqliteFile))
        copyFile(sqliteFile, newFileName)

        # Generate a set of 3d coordinates and assign the mapper of the introduced sqlite file
        coordsSet.setSamplingRate(inTomoSet.getSamplingRate())
        coordsSet.setBoxSize(self.boxSize.get())
        coordsSet._mapperPath.set('%s, %s' % (newFileName, ''))
        coordsSet.load()

        # Check if the coordinates and the tomograms can be related via the tomoId or the filename
        outTomoSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        outTomoSet.copyInfo(inTomoSet)
        self._checkCoordinatesMatching(inTomoSet, outTomoSet, coordsSet)
        if self.notFoundCoordsMsg:
            self._store()

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
        summaryMsg = []
        if self.isFinished():
            statusMsg = getattr(self, 'notFoundCoordsMsg', None)
            if statusMsg:
                summaryMsg.append(statusMsg.get())

        return summaryMsg

    # --------------------------- UTILS functions ----------------------------

    def _checkCoordinatesMatching(self, inTomoSet, outTomoSet, coordsSet):
        notFoundCoords = []
        notFoundCoordsMsg = ''
        notFoundTomosMsg = ''
        inTomoSetMatchingIndices = []
        pattern = 'Row %i  -  tomoId = %s  -  (x, y, x) = (%.2f, %.2f, %.2f)'
        tomoTsIdList, tomoBaseNameList = zip(*[(tomo.getTsId(), removeBaseExt(tomo.getFileName()))
                                               for tomo in inTomoSet])
        for coord in coordsSet:
            coordTomoId = coord.getTomoId()
            if coordTomoId in tomoTsIdList:
                indByTomoId = tomoTsIdList.index(coordTomoId) + 1
                coord.setVolume(inTomoSet[indByTomoId])
                inTomoSetMatchingIndices.append(indByTomoId)
            else:
                indexByName = self._getMatchingIndexByFileName(coordTomoId, tomoBaseNameList)
                if indexByName:
                    coord.setVolume(inTomoSet[indexByName])
                    inTomoSetMatchingIndices.append(indexByName)
                else:
                    coord.setVolume(inTomoSet[1])  # 3D coordinate must be referred to a volume to get its origin
                    notFoundCoords.append(pattern % (coord.getObjId(), coordTomoId, *coord.getPosition(SCIPION)))

        # Build a precedents set with only the matching tomograms, in case there are not all the ones present in the
        # input set
        pattern = '\t-{}\n'
        if inTomoSetMatchingIndices:
            inTomoSetMatchingIndices = set(inTomoSetMatchingIndices)
            if len(inTomoSetMatchingIndices) < inTomoSet.getSize():
                for ind in set(inTomoSetMatchingIndices):
                    outTomoSet.append(inTomoSet[ind])
                coordsSet.setPrecedents(outTomoSet)
                # Format the non-matching coordinates message and add the header
                notFoundCoordsMsg += 'The following coordinates have a tomoId which was not found in the tsId ' \
                                     'attribute of none of the tomograms introduced nor contained in their basename:' \
                                     '\n\n%s' % (pattern * len(notFoundCoords))
                notFoundCoordsMsg.format(*notFoundCoords)

                # Generate a message to report about the non-matching tomograms found
                notMatchingTomoFiles = self._getNotMatchingTomoFiles(inTomoSet, inTomoSetMatchingIndices)
                notFoundTomosMsg += 'The following tomograms were excluded from the set because no coordinates are ' \
                                    'referred to them:\n\n%s' % (pattern * len(notMatchingTomoFiles))
                notFoundTomosMsg.format(*notMatchingTomoFiles)
        else:
            raise Exception(ERR_COORDS_FROM_SQLITE_NO_MATCH)

        if notFoundCoords:
            # Format the non-matching coordinates message and add the header
            notFoundCoordsMsg += ('The following coordinates have a tomoId which was not found in the tsId '
                                  'attribute of none of the tomograms introduced nor contained in their basename:'
                                  '\n%s' % (pattern * len(notFoundCoords))).format(*notFoundCoords)

        self.notFoundCoordsMsg = String(notFoundTomosMsg + '\n\n' + notFoundCoordsMsg if
                                        notFoundTomosMsg else notFoundCoordsMsg)

    @staticmethod
    def _getMatchingIndexByFileName(coordTomoId, tomoBaseNameList):
        matchingIndex = None
        matches = list(map(lambda x: coordTomoId in x, tomoBaseNameList))
        if any(matches):
            matchingIndex = matches.index(True) + 1

        return matchingIndex

    @staticmethod
    def _getNotMatchingTomoFiles(inTomoSet, inTomoSetMatchingIndices):
        return [inTomoSet[ind + 1].getFileName() for ind in range(inTomoSetMatchingIndices)
                if ind not in inTomoSetMatchingIndices]


