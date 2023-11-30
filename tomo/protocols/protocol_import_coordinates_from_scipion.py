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
from enum import Enum
from os.path import exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol import FileParam, IntParam, PointerParam
from pyworkflow.utils import Message, removeBaseExt, yellowStr
from ..constants import SCIPION, ERR_COORDS_FROM_SQLITE_NO_MATCH
from .protocol_base import ProtTomoBase
from ..objects import SetOfCoordinates3D


class outputObjs(Enum):
    coordinates = SetOfCoordinates3D


class ProtImportCoordinates3DFromScipion(EMProtocol, ProtTomoBase):
    """Protocol to import a set of 3d coordinates from Scipion sqlite file"""

    _label = 'import 3D coordinates from scipion'
    _devStatus = BETA
    _possibleOutputs = outputObjs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.notMatchingMsg = None

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
        inCoordsSet = SetOfCoordinates3D()
        outCoordsSet = self._createSetOfCoordinates3D(self.importTomograms)

        # Generate a set of 3d coordinates and assign the mapper of the introduced sqlite file
        inCoordsSet.setSamplingRate(inTomoSet.getSamplingRate())
        inCoordsSet.setBoxSize(self.boxSize.get())
        inCoordsSet._mapperPath.set('%s, %s' % (self.sqliteFile.get(), ''))
        inCoordsSet.load()

        # Check if the coordinates and the tomograms can be related via the tomoId or the filename
        self._checkCoordinatesMatching(inTomoSet, inCoordsSet, outCoordsSet)
        if self.notMatchingMsg:
            self._store()

        # Set some set attributes
        outCoordsSet.setSamplingRate(inTomoSet.getSamplingRate())
        outCoordsSet.setBoxSize(self.boxSize.get())
        # Set volume pointers (required by getCoordinates().getX, Y and Z
        tomoIdsDict = {tomo.getTsId(): tomo.clone() for tomo in inTomoSet}
        for coord in outCoordsSet.iterCoordinates():
            coord.setVolume(tomoIdsDict[coord.getTomoId()])

        # Define outputs and relations
        self._defineOutputs(**{outputObjs.coordinates.name: outCoordsSet})
        self._defineSourceRelation(inTomoSet, outCoordsSet)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        errorList = []
        if not exists(self.sqliteFile.get()):
            errorList.append('Introduced file was not found:\n\t%s' % self.sqliteFile.get())

        return errorList

    def _summary(self):
        summaryMsg = []
        if self.isFinished():
            if getattr(self, 'outputTomograms', None):
                summaryMsg.append('A *set of tomograms was generated* containing only the ones which there are\n'
                                  'at least one coordinate referred to.\n')
            statusMsg = getattr(self, 'notMatchingMsg', None)
            if statusMsg:
                summaryMsg.append(statusMsg.get())

        return summaryMsg

    # --------------------------- UTILS functions ----------------------------

    def _checkCoordinatesMatching(self, inTomoSet, inCoordsSet, outCoordsSet):
        notFoundCoords = []
        notFoundCoordsMsg = ''
        notFoundTomosMsg = ''
        inTomoSetMatchingIndices = []
        pattern = 'Row %i  -  tomoId = %s  -  (x, y, x) = (%.2f, %.2f, %.2f)'
        tomoTsIdList, tomoBaseNameList = zip(*[(tomo.getTsId(), removeBaseExt(tomo.getFileName()))
                                               for tomo in inTomoSet])
        for coord in inCoordsSet:
            coordTomoId = coord.getTomoId()
            if coordTomoId:
                if coordTomoId in tomoTsIdList:
                    indByTomoId = tomoTsIdList.index(coordTomoId) + 1
                    coord.setVolume(inTomoSet[indByTomoId])
                    inTomoSetMatchingIndices.append(indByTomoId)
                    # Add it to the output set of coordinates
                    outCoordsSet.append(coord)
                else:
                    indexByName = self._getMatchingIndexByFileName(coordTomoId, tomoBaseNameList)
                    if indexByName:
                        coord.setVolume(inTomoSet[indexByName])
                        inTomoSetMatchingIndices.append(indexByName)
                        # Add it to the output set of coordinates
                        outCoordsSet.append(coord)
                    else:
                        self._appendBaddCoordMsgToList(coord, notFoundCoords, inTomoSet, coordTomoId, pattern)

            else:
                self._appendBaddCoordMsgToList(coord, notFoundCoords, inTomoSet, 'NoTomoId', pattern)

        # Build a precedents set with only the matching tomograms, in case there are not all the ones present in the
        # input set
        pattern = '\t-{}\n'
        if not inTomoSetMatchingIndices:
            raise Exception(ERR_COORDS_FROM_SQLITE_NO_MATCH)

        if notFoundCoords:
            nOfNonMatchingCoords = len(notFoundCoords)
            # Format the non-matching coordinates message and add the header
            notFoundCoordsMsg += '*[%i] coordinates were excluded*.\nThey have a tomoId which was not found in the ' \
                                 'tsId attribute of none of the tomograms introduced nor contained in their basename.' \
                                 '\nThe details can be checked in the output log.' % nOfNonMatchingCoords

            # Print the detailed information in the output log
            print(yellowStr(('EXCLUDED COORDINATES [%i]:\n%s' %
                            (nOfNonMatchingCoords, pattern * nOfNonMatchingCoords)).format(*notFoundCoords)))

        self.notMatchingMsg = String(notFoundTomosMsg + '\n\n' + notFoundCoordsMsg if
                                     notFoundTomosMsg else notFoundCoordsMsg)

    @staticmethod
    def _getMatchingIndexByFileName(coordTomoId, tomoBaseNameList):
        matchingIndex = None
        matches = list(map(lambda x: coordTomoId in x, tomoBaseNameList))
        if any(matches):
            matchingIndex = matches.index(True) + 1

        return matchingIndex

    @staticmethod
    def _appendBaddCoordMsgToList(coord, notFoundCoordsList, inTomoSet, coordTomoId, pattern):
        coord.setVolume(inTomoSet[1])  # 3D coordinate must be referred to a volume to get its origin
        notFoundCoordsList.append(pattern % (coord.getObjId(), coordTomoId, *coord.getPosition(SCIPION)))


