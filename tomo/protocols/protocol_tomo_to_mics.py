# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
from os.path import basename
import mrcfile
import numpy as np
from pwem.objects import SetOfMicrographs, Micrograph, Acquisition
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, GT
from pwem.protocols import EMProtocol
from pyworkflow.utils import replaceExt, removeExt
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfCoordinates3D, Coordinate3D, SetOfTomograms


class ProtTomoToMicsOutput(enum.Enum):
    outputMicrographs = SetOfMicrographs()


class ProtTomoToMics(EMProtocol):
    """ Turns tomograms into set of micrographs to apply SPA picking methods."""
    _label = 'tomograms to micrographs'
    _devStatus = BETA
    _possibleOutputs = ProtTomoToMicsOutput

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      important=True,
                      help='Select the tomograms to be turned into micrographs')
        form.addParam('slicesGap', IntParam,
                      label="Slices gap",
                      default=10,
                      validators=[GT(0)],
                      help='Number of slices to skip when turning the tomogram slices into micrographs, starting from '
                           'the central slice. For example, if a tomogram has a thickness of 300 pixels and the gap '
                           'is set to 20, then the list of slices that will be used as reference for each micrograph '
                           'will be, considering that the slices are numbered from 0 to the thickness value - 1:'
                           '\n\t- Central slice = 149'
                           '\n\t- Going up = 169, 189, 209, 229, 249, 269, 289'
                           '\n\t- Going down = 129, 109, 89, 69, 49, 29, 9.'
                           '\nSo 15 micrographs will be generated.')
        form.addParam('noSlicesToAvg', IntParam,
                      label='No. slices to sum',
                      default=5,
                      validators=[GT(0)],
                      help='It represents the number of adjacent slices that will be sum for each reference slice to '
                           'get the corresponding micrographs (see help of the parameter "Slices gap"). For example, '
                           'if the number is 5, 2 slices above and 2 slices below will be sum for each reference '
                           'slice, so 5 slices will be sum in total to get each micrograph. If set to 1, no sum will '
                           'be performed.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def getInputTomograms(self) -> SetOfTomograms:
        return self.input.get()

    def createOutputStep(self):
        input = self.input.get()
        output = SetOfMicrographs.create(self._getPath())
        output.setSamplingRate(input.getSamplingRate())

        # Acquisition
        micAcq = Acquisition()
        tomoAcq = input.getAcquisition()
        micAcq.setMagnification(tomoAcq.getMagnification())
        micAcq.setSphericalAberration(tomoAcq.getSphericalAberration())
        micAcq.setVoltage(tomoAcq.getVoltage())
        micAcq.setAmplitudeContrast(tomoAcq.getAmplitudeContrast())
        micAcq.setDoseInitial(tomoAcq.getDoseInitial())
        micAcq.setDosePerFrame(tomoAcq.getDosePerFrame())

        output.setAcquisition(micAcq)

        # For each tomogram
        for tomo in input:
            self.appendMicsFromTomogram(output, tomo)

        self._defineOutputs(**{ProtTomoToMicsOutput.outputMicrographs.name: output})
        self._defineSourceRelation(self.input, output)

    # --------------------------- UTILS functions --------------------------------------------
    def appendMicsFromTomogram(self, output, tomo):
        nSlicesAvg = self.noSlicesToAvg.get()
        gap = self.slicesGap.get()
        sRate = tomo.getSamplingRate()
        self.info("Creating micrographs for %s with a %s gap and adding %s contiguous slices." %
                  (tomo.getTsId(), gap, nSlicesAvg))

        # Load the tomogram image to get the dimensions
        with mrcfile.mmap(tomo.getFileName(), mode='r+') as tomoMrc:
            tomogramData = tomoMrc.data
            thk, width, height = tomogramData.shape
        shape = (1, width, height)

        # For each slice
        for indexLims in self.genCenterList(thk, gap, nSlicesAvg):
            downLim, center, upLim = indexLims[:]
            micName = tomoSliceToMicName(tomo, center)
            self.info("Creating micrograph (%s) from %s to %s slices of %s" % (micName, downLim, upLim-1,tomo.getTsId()))
            outputMicName = replaceExt(self._getExtraPath(micName), "mrc")
            with mrcfile.new_mmap(outputMicName, shape, overwrite=True) as mic:
                micData = np.empty(shape, dtype=tomogramData.dtype)
                micData[0] = np.sum(tomogramData[downLim:upLim], axis=0)
                mic.set_data(micData)
                mic.update_header_from_data()
                mic.update_header_stats()
                mic.voxel_size = sRate

            # Create the micrograph metadata object
            newMic = Micrograph()
            newMic.setFileName(outputMicName)
            newMic.setMicName(micName)
            newMic.setSamplingRate(tomo.getSamplingRate())

            # Append it
            output.append(newMic)

    @staticmethod
    def genCenterList(dataLen, gap, nSlicesAvg):
        """Generates a list of lists in which each element is composed by the lower limit, the center and the upper
        limit of the intervals of indices with a given step (gap) between them and calculating the adjacent indices
        given the total number of slices that will include the interval (nSlicesAvg)."""
        center = dataLen // 2 - 1  # Indices begin in 0, not in 1
        slicesUp = (nSlicesAvg - 1) // 2  # Minus the gap index
        slicesDown = nSlicesAvg - 1 - slicesUp  # The remaining ones, to cover even and odd cases of nSlicesAvg
        indexListOfLists = [[center - slicesDown, center, center + slicesUp + 1]]
        for i in range(gap, center + 1, gap):
            if center - i - slicesDown >= 0:
                index = center - i
                indexListOfLists.append([index - slicesDown, index, index + slicesUp + 1])
            if center + i + slicesUp <= dataLen:
                index = center + i
                indexListOfLists.append([index - slicesDown, index,  index + slicesUp + 1])
        # Sort the list of lists based on the center element of each subList
        return sorted(indexListOfLists, key=lambda l: l[1])

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('Slices gap = *%i*\nSlices for averaging = *%i*' %
                           (self.slicesGap.get(), self.noSlicesToAvg.get()))
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validateMsgs = []
        return validateMsgs


class Prot2DcoordsTo3DCoordsOutput(enum.Enum):
    outputCoordinates = SetOfCoordinates3D()


class Prot2DcoordsTo3DCoords(EMProtocol):
    """ Turns 2d coordinates into set of 3d coordinates. Works in coordination with
    'tomograms to micrographs' protocol"""
    _label = '2d coordinates to 3d coordinates'
    _devStatus = BETA
    _possibleOutputs = Prot2DcoordsTo3DCoordsOutput

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('tomograms', PointerParam, pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Select the tomograms to be associated to the 3D coordinates')
        form.addParam('coordinates', PointerParam, pointerClass='SetOfCoordinates', label="2D Coordinates",
                      help='Set of 2d coordinates picked on tomogram slices.')

    # --------------------------- INSERT steps functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        tomograms = self.tomograms.get()
        coordinates = self.coordinates.get()

        output = SetOfCoordinates3D.create(self._getPath())
        # output.copyInfo(input)
        output.setPrecedents(tomograms)
        output.setBoxSize(coordinates.getBoxSize())
        output.setSamplingRate(coordinates.getMicrographs().getSamplingRate())

        # Get a dictionary of tomograms by tsId
        tomoDict = dict()
        for tomogram in tomograms:
            tomo = tomogram.clone()
            tomoDict[removeExt(basename(tomo.getFileName()))] = tomo

        # For each 2d coordinate
        for cord2d in coordinates:
            # Extract Z
            micName = coordinates.getMicrographs()[cord2d.getMicId()].getFileName()
            z, fileName = sliceAndNameFromMicName(micName)
            newCoord = Coordinate3D()
            newCoord.setVolume(tomoDict[fileName])
            newCoord.setX(cord2d.getX(), BOTTOM_LEFT_CORNER)
            newCoord.setY(cord2d.getY(), BOTTOM_LEFT_CORNER)
            newCoord.setZ(z, BOTTOM_LEFT_CORNER)
            newCoord.setBoxSize(coordinates.getBoxSize())
            newCoord.setGroupId(1)
            output.append(newCoord)

        self._defineOutputs(**{Prot2DcoordsTo3DCoordsOutput.outputCoordinates.name: output})
        self._defineSourceRelation(self.coordinates, output)

    # --------------------------- UTILS functions --------------------------------------------

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validateMsgs = []
        return validateMsgs


def tomoSliceToMicName(tomo, slice):
    return f"S{slice:03}_{basename(tomo.getFileName())}"


def sliceAndNameFromMicName(micname):
    """ Extracts z and tomo name from a micname composed with tomoSliceToMicName"""

    parts = removeExt(basename(micname)).split("_")
    slice = parts[0]
    slice = int(slice.replace("S", ""))
    fileName = "_".join(parts[1:])
    return slice, fileName
