# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from os.path import join, exists
from pathlib import PureWindowsPath

import numpy as np

from pwem.emlib.metadata import (MetaData, MDL_XCOOR, MDL_YCOOR, MDL_ZCOOR)

import pyworkflow.utils as pwutils
from pyworkflow.utils import removeBaseExt, getParentFolder


class TomoImport:

    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = protocol.getCopyOrLink()

    def importCoordinates3D(self, fileName, addCoordinate):
        from tomo.objects import Coordinate3D
        if pwutils.exists(fileName):
            ext = pwutils.getExt(fileName)

        if ext == ".txt":
            md = MetaData()
            md.readPlain(fileName, "xcoor ycoor zcoor")
            for objId in md:
                x = md.getValue(MDL_XCOOR, objId)
                y = md.getValue(MDL_YCOOR, objId)
                z = md.getValue(MDL_ZCOOR, objId)
                coord = Coordinate3D()
                coord.setPosition(x, y, z)
                addCoordinate(coord)

        else:
            raise Exception('Unknown extension "%s" to import Eman coordinates' % ext)


def getMeshVolFileName(volId):
    return 'Meshes_Vol%d.txt' % volId

def setOfMeshes2Files(meshes, path):

    def writeFile():
        fnInputCoor = getMeshVolFileName(currentVolId)
        pathInputCoor = pwutils.join(path, fnInputCoor)
        np.savetxt(pathInputCoor, np.asarray(coords), fmt='%d', delimiter=",")

    currentVolId = None
    coords = []
    for coor in meshes.iterCoordinates(orderBy="_volId"):
        if coor.getVolId() != currentVolId:
            if currentVolId != None:
                writeFile()
            currentVolId = coor.getVolId()
            coords = []
        coords.append([coor.getX(), coor.getY(), coor.getZ(), coor.getGroupId()])
    if coords:
        writeFile()


class MDoc:

    def __init__(self, fileName, isImportingTsMovies=True):
        self._mdocFileName = fileName
        self.isImportingTsMovies = isImportingTsMovies
        self._tsId = None
        # Acquisition general attributes
        self._voltage = None
        self._magnification = None
        self._samplingRate = None
        # Acquisition specific attributes (per angle)
        self._angles = []
        self._angleMovieFiles = []
        self._accumulatedDoses = []

    def read(self):
        validateTSFromMdocErrMsgList = ''
        tsFile = None
        mdoc = self._mdocFileName
        headerDict, zSlices = self._parseMdoc()

        # Get acquisition general info
        self._getAcquisitionInfoFromMdoc(headerDict, zSlices[0])
        self._tsId = removeBaseExt(mdoc)
        parentFolder = getParentFolder(mdoc)
        if not self.isImportingTsMovies:
            tsFile = join(parentFolder, self._tsId + '.mrcs')
            validateTSFromMdocErrMsgList = self._validateTSFromMdoc(mdoc, tsFile)

        # Get acquisition specific (per angle) info
        self._getSlicesData(zSlices, tsFile)

        # Check Mdoc info read
        validateMdocContentsErrorMsgList = self._validateMdocInfoRead()

        # Check all the possible errors found
        exceptionMsg = ''
        if validateTSFromMdocErrMsgList:
            exceptionMsg += ' '.join(validateTSFromMdocErrMsgList)
        if validateMdocContentsErrorMsgList:
            exceptionMsg += ' '.join(validateMdocContentsErrorMsgList)

        return exceptionMsg

    def _parseMdoc(self):
        """
        Parse the mdoc file and return a list with a dict key=value for each
        of the [Zvalue = X] sections and a dictionary for the first lines global variables.
        :return: dictionary (header), list of dictonaries (Z slices)
        """
        headerDict = {}
        headerParsed = False
        zvalueList = []

        with open(self._mdocFileName) as f:
            for line in f:
                if line.startswith('[ZValue'):
                    # We have found a new Zvalue
                    headerParsed = True
                    zvalue = int(line.split(']')[0].split('=')[1])
                    if zvalue != len(zvalueList):
                        raise Exception("Unexpected ZValue = %d" % zvalue)
                    zvalueDict = {}
                    zvalueList.append(zvalueDict)
                else:
                    if not line.startswith('[T') and line.strip():
                        key, value = line.split('=')
                        if not headerParsed:
                            headerDict[key.strip()] = value.strip()
                        if zvalueList:
                            zvalueDict[key.strip()] = value.strip()

        return headerDict, zvalueList

    def _getAcquisitionInfoFromMdoc(self, headerDict, firstSlice):
        """Acquisition data is read from to data sources (from higher to lower priority):
            - From the first ZSlice data.
            - From the file header data."""
        VOLTAGE = 'Voltage'
        MAGNIFICATION = 'Magnification'
        PIXEL_SPACING = 'PixelSpacing'
        self._voltage = firstSlice.get(VOLTAGE, headerDict.get(VOLTAGE, None))
        self._magnification = firstSlice.get(MAGNIFICATION, headerDict.get(MAGNIFICATION, None))
        self._samplingRate = firstSlice.get(PIXEL_SPACING, headerDict.get(PIXEL_SPACING, None))

    def _getSlicesData(self, zSlices, tsFile):
        parentFolder = getParentFolder(self._mdocFileName)
        accumulatedDose = 0
        for zSlice in zSlices:
            # Angles
            self._angles.append(zSlice.get('TiltAngle', None))
            # Files
            self._angleMovieFiles.append(self._getAngleMovieFileName(parentFolder, zSlice, tsFile))
            # Doses
            accumulatedDose = self._getDoseFromMdoc(zSlice, accumulatedDose)
            self._accumulatedDoses.append(accumulatedDose)

    @staticmethod
    def _getAngleMovieFileName(parentFolder, zSlice, tsFile):
        if tsFile:
            return tsFile
        else:
            # PureWindowsPath pathlib is ised to make possible deal with different path separators, like \\
            return join(parentFolder, PureWindowsPath(zSlice['SubFramePath']).parts[-1])

    @staticmethod
    def _getDoseFromMdoc(zSlice, accumulatedDose):
        """It calculates the accumulated dose on the frames represented by zSlice, and add it to the
        previous accumulated dose"""
        EXPOSURE_DOSE = 'ExposureDose'  # Dose on specimen during camera exposure in electrons/sq. Angstrom
        FRAME_DOSES_AND_NUMBERS = 'FrameDosesAndNumbers'  # Dose per frame in electrons per square Angstrom followed
        # by number of frames at that dose
        DOSE_RATE = 'DoseRate'  # Dose rate to the camera, in electrons per unbinned pixel per second
        EXPOSURE_TIME = 'ExposureTime'  # Image exposure time
        MIN_MAX_MEAN = 'MinMaxMean'  # Minimum, maximum, and mean value for this image
        PIXEL_SIZE = 'PixelSpacing'  # Pixel spacing in Angstroms for individual image
        COUNTS_PER_ELECTRON = 'CountsPerElectron'

        def _keysInDict(listOfKeys):
            return all([key in zSlice.keys() for key in listOfKeys])

        # Different ways of calculating the dose, ordered by priority considering the possible variability between
        # different mdoc files
        newDose = 0

        # Directly from field ExposureDose
        if EXPOSURE_DOSE in zSlice:
            expDoseVal = zSlice[EXPOSURE_DOSE]
            if expDoseVal:
                newDose = float(expDoseVal)

        # Directly from field FrameDosesAndNumbers
        if not newDose and FRAME_DOSES_AND_NUMBERS in zSlice:
            frameDoseAndNums = zSlice[FRAME_DOSES_AND_NUMBERS]
            if frameDoseAndNums:
                newDose = float(list(frameDoseAndNums.strip())[-1])  # Get the mean from a string like '0 6'

        # Calculated from fields DoseRate and ExposureTime
        if not newDose and _keysInDict([DOSE_RATE, EXPOSURE_TIME]):
            doseRate = zSlice[DOSE_RATE]
            expTime = zSlice[EXPOSURE_TIME]
            if doseRate and expTime:
                newDose = float(doseRate) * float(expTime)

        # Calculated from fields MinMaxMean, PixelSpacing and CountsPerElectron
        if not newDose and _keysInDict([MIN_MAX_MEAN, PIXEL_SIZE, COUNTS_PER_ELECTRON]):
            minMaxMean = zSlice[MIN_MAX_MEAN]
            pixelSize = zSlice[PIXEL_SIZE]
            counts = zSlice[COUNTS_PER_ELECTRON]
            if all([minMaxMean, pixelSize, counts]):
                meanVal = list(minMaxMean.strip())[-1]  # Get the mean from a string like '-42 2441 51.7968'
                newDose = (float(meanVal) / int(counts)) / float(pixelSize) ** 2

        return newDose + accumulatedDose

    @staticmethod
    def _validateTSFromMdoc(mdoc, tsFile):
        errMsg = ''
        if not exists(tsFile):
            errMsg = '\nMdoc --> %s\nExpected tilt series file not found \n%s' % (mdoc, tsFile)

        return errMsg

    def _validateMdocInfoRead(self):
        validateMdocContentsErrorMsgList = []
        msg = ['\n*Data not found in file*\n%s:\n' % self._mdocFileName]
        missingFiles = [file for file in self._angleMovieFiles if not exists(file)]
        if not self._voltage:
            msg.append('*Voltage*')
        if not self._magnification:
            msg.append('*Magnification*')
        if not self._samplingRate:
            msg.append('*PixelSpacing*')
        if None in self._angles:
            indices = [str(i) for i, v in enumerate(self._angles) if v is None]
            msg.append('*TiltAngle*: %s' % ' '.join(indices))
        if missingFiles:
            msg.append('*Missing files*:\n%s' % ' '.join(missingFiles))
        if len(msg) > 1:
            validateMdocContentsErrorMsgList.append(' '.join(msg))

        return validateMdocContentsErrorMsgList

    def getFileName(self):
        return self._mdocFileName

    def getTsId(self):
        return self._tsId

    def getVoltage(self):
        return self._voltage

    def getMagnification(self):
        return self._magnification

    def getSamplingRate(self):
        return self._samplingRate

    def getTiltAngles(self):
        return self._angles

    def getAngleMovieFiles(self):
        return self._angleMovieFiles

    def getAccumDoses(self):
        return self._accumulatedDoses




