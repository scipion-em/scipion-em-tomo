import logging
import os
from datetime import datetime as dt
from os.path import join, exists
from pyworkflow.utils import getParentFolder, removeBaseExt
from pathlib import PureWindowsPath


logger = logging.getLogger(__name__)
SUB_FRAME_PATH = 'SubFramePath'
TS_PREFIX = "TS_"
# Dose on specimen during camera exposure in electrons/sq. Angstrom
EXPOSURE_DOSE = 'ExposureDose'
# Dose per frame in electrons per square Angstrom followed
# by number of frames at that dose
FRAME_DOSES_AND_NUMBERS = 'FrameDosesAndNumber'
# Dose rate to the camera, in electrons per unbinned pixel per second
DOSE_RATE = 'DoseRate'
# Image exposure time
EXPOSURE_TIME = 'ExposureTime'
# Minimum, maximum, and mean value for this image
MIN_MAX_MEAN = 'MinMaxMean'
COUNTS_PER_ELECTRON = 'CountsPerElectron'
DIVIDED_BY_TWO = 'DividedBy2'


class MDoc:
    """class to interact with the IMOD "autodoc" files used by serialEM.
    This format consists of keyword-value pairs organized into blocks
    called sections.
    A section begins with a bracketed key-value pair:
      [sectionType = name]
    where the section "name" or value will typically be unique.
    Lines below a section header of the form
      key = value
    provide data associated with that section.

    In addition, key-value pairs can occur at the beginning of the file,
    before any section header, and these are referred to as global values.
    Files with extension ".mdoc" provides data about an MRC file and has
    the same name as the image file, with the additional extension ".mdoc".
    In these files, the main section type is "ZValue" and the name
    for each section is the Z value of the image in the file, numbered from 0.
    A description of each key is available at URL:
    https://bio3d.colorado.edu/SerialEM/hlp/html/about_formats.htm

    Additional information may be stored in section headers of the type "T"
    (i.e. [T  a = b ]). In theory these information in also stored in the
    "titles" of the MRC files.

    Example of mdoc file

    DataMode = 6   # MRC type
    ImageSize = 4096 4096 # x, y dim
    ImageFile = 20211130_PKV_tomo_01.mrc # image with tilt series
    PixelSpacing = 3.64  # A/px
    Voltage = 200.00  # kvolt

    [T = Tomography: TALOS-D3558    21-Nov-30  17:42:06]
    [T =   TiltAxisAngle = -91.81  Binning = 1  SpotSize = 7]

    [ZValue = 0] # first movie
    TiltAngle = -0.01
    ...
    [ZValue = 1] # second movie
    TiltAngle = -10.00
    ...
    """

    def __init__(self, fileName, voltage=None, magnification=None,
                 samplingRate=None, doseProvidedByUser=None,
                 tiltAngleProvidedByUser=None):

        self._mdocFileName = fileName
        self._tsId = None

        # Dose related attributes
        self.doseProvidedByUser = doseProvidedByUser
        self.mdocHasDose = False

        # Acquisition general attributes
        self._voltage = voltage
        self._magnification = magnification
        self._samplingRate = samplingRate
        self._tiltAxisAngle = tiltAngleProvidedByUser

        # Acquisition specific attributes (per angle)
        self._tiltsMetadata = []

    @staticmethod
    def normalizeTSId(mdocFileName):
        mdocBaseName = removeBaseExt(mdocFileName)
        tsPrefixAdded = True if mdocBaseName[0].isdigit() else False
        normTsId = normalizeTSId(mdocFileName)
        return normTsId, tsPrefixAdded

    def read(self, isImportingTsMovies=True, ignoreFilesValidation=False):
        """ Parse mdoc file. There is a collection of getters that may be
        used to access this information:
        """
        validateTSFromMdocErrMsg = ''
        tsFile = None
        mdoc = self._mdocFileName

        try:
            headerDict, zSlices = self._parseMdoc()
            zSlices = self._sortByTimestamp(zSlices)

            # Get acquisition general info
            self._getAcquisitionInfoFromMdoc(headerDict, zSlices[0])
            self._tsId, tsPrefixAdded = self.normalizeTSId(mdoc)
            parentFolder = getParentFolder(mdoc)
            if not isImportingTsMovies:
                # Some mdoc files point to an .st file stored in the ImageFile
                # header line
                tsFile = join(parentFolder, headerDict.get("ImageFile", ''))
                expectedBaseName = self._tsId.replace(TS_PREFIX, '') if tsPrefixAdded else self._tsId
                if not os.path.exists(tsFile):
                    tsFile = join(parentFolder, expectedBaseName + '.mrcs')
                if not os.path.exists(tsFile):
                    tsFile = join(parentFolder, expectedBaseName + '.mrc')
                if not os.path.exists(tsFile):
                    tsFile = join(parentFolder, expectedBaseName + '.st')
                validateTSFromMdocErrMsg = self._validateTSFromMdoc(mdoc, tsFile)

            # Get acquisition specific (per angle) info
            self._getSlicesData(zSlices, tsFile)

            # Check Mdoc info read
            validateMdocContentsErrorMsgList = self._validateMdocInfoRead(
                ignoreFilesValidation=ignoreFilesValidation or not
                isImportingTsMovies)

            # Check all the possible errors found
            exceptionMsg = ''
            if validateTSFromMdocErrMsg:
                exceptionMsg += validateTSFromMdocErrMsg
            if validateMdocContentsErrorMsgList:
                exceptionMsg += ' '.join(validateMdocContentsErrorMsgList)

            # If ts images we are assuming slices follow the angle order
            if not isImportingTsMovies and not exceptionMsg:
                self._tiltsMetadata.sort(key=lambda x: float(x.getTiltAngle()),
                                         reverse=False)
                for index, tiMd in enumerate(self._tiltsMetadata):
                    tiMd.setAngleMovieFile((index + 1, tiMd.getAngleMovieFile()))

            return exceptionMsg
        except Exception as e:

            return "\n*CRITICAL mdoc parsing error: %s can't be parsed.*\n %s\n" % (mdoc, str(e))

    def _parseMdoc(self):
        """
        Parse the mdoc file and return a list with a dict key=value for each
        of the [Zvalue = X] sections and a dictionary for the first lines
        global variables.

        :return: dictionary (header), list of dictionaries (Z slices)
        """
        headerDict = {}
        headerParsed = False
        zvalueList = []  # list of dictionaries with
        with open(self._mdocFileName) as f:
            for line in f:
                if line.startswith('[ZValue'):  # each tilt movie
                    # We have found a new z value
                    headerParsed = True
                    zvalue = int(line.split(']')[0].split('=')[1])
                    if zvalue != len(zvalueList):
                        raise Exception("Unexpected Z value = %d" % zvalue)
                    zvalueDict = {}
                    zvalueList.append(zvalueDict)
                elif line.startswith('[T'):  # auxiliary global information
                    if self.getTiltAxisAngle():
                        # It's in the mdoc, but the user has specified
                        # it manually
                        continue
                    else:
                        strLine = line.strip().replace(' ', ''). \
                            replace(',', '').lower()
                        pattern = 'tiltaxisangle='
                        if pattern in strLine:
                            # Example of the most common syntax
                            # [T =     Tilt axis angle = 90.1, binning = 1
                            #                   spot = 9  camera = 0]
                            # [T =     TiltAxisAngle = -91.81  Binning = 1
                            #                                   SpotSize = 7]
                            tiltAxisAngle = \
                                strLine.split('=')[2].split('binning')[0]
                            # Check if it's a string which
                            # represents a float or not
                            if tiltAxisAngle.lstrip('-+'). \
                                    replace('.', '', 1).isdigit():
                                self._tiltAxisAngle = float(tiltAxisAngle)
                elif line.strip():  # global variables no in [T sections]
                    key, value = line.split('=')
                    if not headerParsed:
                        headerDict[key.strip()] = value.strip()
                    if zvalueList:
                        zvalueDict[key.strip()] = value.strip()

        return headerDict, zvalueList

    def _getAcquisitionInfoFromMdoc(self, headerDict, firstSlice):
        """Acquisition data is read from to data sources
           (from higher to lower priority):
            - From the values introduced in the form by
              the user (introduced as attributes of the MDoc object)
            - From the first ZSlice data.
            - From the file header data."""
        if not self.getVoltage():
            VOLTAGE = 'Voltage'
            self._voltage = \
                firstSlice.get(VOLTAGE, headerDict.get(VOLTAGE, self._voltage))
        if not self.getMagnification():
            MAGNIFICATION = 'Magnification'
            self._magnification = \
                firstSlice.get(MAGNIFICATION,
                               headerDict.get(MAGNIFICATION,
                                              self._magnification))
        if not self.getSamplingRate():
            PIXEL_SPACING = 'PixelSpacing'
            self._samplingRate = \
                firstSlice.get(PIXEL_SPACING,
                               headerDict.get(PIXEL_SPACING,
                                              self._samplingRate))

    def _getSlicesData(self, zSlices, tsFile):
        parentFolder = getParentFolder(self._mdocFileName)
        accumulatedDose = 0
        for counter, zSlice in enumerate(zSlices):
            if self.doseProvidedByUser:
                incomingDose = self.doseProvidedByUser
            else:
                incomingDose = \
                    self._getDoseFromMdoc(zSlice, self.getSamplingRate())
            accumulatedDose += incomingDose
            self._tiltsMetadata.append(TiltMetadata(
                angle=zSlice.get('TiltAngle', None),
                angleFile=self._getAngleMovieFileName(
                    parentFolder, zSlice, tsFile),
                acqOrder=counter + 1,
                accumDose=accumulatedDose,
                incomingDose=incomingDose
            ))
        # round is used to make the condition more robust
        # for cases like 0.0000000001
        if round(accumulatedDose) > 0:
            self.mdocHasDose = True
        else:
            logger.debug("Dose not found or almost 0 (%s) in %s" %
                         (accumulatedDose, self._mdocFileName))

    @staticmethod
    def _sortByTimestamp(zSlices):
        """ MDOC file is not necessarily sorted by acquisition order,
            use TimeStamp key to sort Z-slices.
        """
        dtValue = zSlices[0].get('DateTime', None)
        if dtValue is None:
            return zSlices
        else:
            if len(dtValue.split()[0]) == 9:  # year is two digits
                fmt = '%d-%b-%y  %H:%M:%S'
            else:
                fmt = '%d-%b-%Y  %H:%M:%S'

        zSlices_sorted = sorted(zSlices,
                                key=lambda d: dt.strptime(d['DateTime'], fmt))

        return zSlices_sorted

    @staticmethod
    def _getAngleMovieFileName(parentFolder, zSlice, tsFile):
        if tsFile:
            return tsFile
        else:
            # PureWindowsPath pathlib is used to make
            # possible deal with different path separators, like \\
            try:
                return join(parentFolder, PureWindowsPath(
                    zSlice[SUB_FRAME_PATH]).parts[-1])
            except Exception as e:
                raise ValueError("Slice section does not have %s field." % SUB_FRAME_PATH)

    @staticmethod
    def _getDoseFromMdoc(zSlice, pixelSize):
        """It calculates the accumulated dose on the frames represented by
        zSlice, and add it to the previous accumulated dose"""
        # Different ways of calculating the dose, ordered by priority
        # considering the possible variability between different mdoc files
        newDose = 0

        # The case  pixelSize = 1 covers the possibility of no sampling rate in form
        # nor mdoc, avoiding the error execution before the whole
        # information is gathered and the exception is raised
        pixelSize = float(pixelSize) if pixelSize else 1

        # Directly from field ExposureDose
        expDoseVal = float(zSlice.get(EXPOSURE_DOSE, 0))  # It may be '0'
        if expDoseVal:
            return expDoseVal

        # Directly from field FrameDosesAndNumbers
        frameDoseAndNums = zSlice.get(FRAME_DOSES_AND_NUMBERS, 0)
        if frameDoseAndNums:
            # Get the mean from a string like '0 6'
            data = frameDoseAndNums.split()
            dosePerFrame = float(data[0])
            nFrames = float(data[1])
            newDose = dosePerFrame * nFrames  # Any of them may be 0
            if newDose:
                return newDose

        # Calculated from fields DoseRate and ExposureTime
        doseRate = float(zSlice.get(DOSE_RATE, 0))
        expTime = float(zSlice.get(EXPOSURE_DOSE, 0))
        if doseRate and expTime:
            return doseRate * expTime / pixelSize ** 2

        # Calculated from fields MinMaxMean, PixelSpacing and CountsPerElectron
        minMaxMean = zSlice.get(MIN_MAX_MEAN, 0)
        counts = float(zSlice.get(COUNTS_PER_ELECTRON, 0))
        if minMaxMean and counts:
            # Get the mean from a string like '-42 2441 51.7968'
            meanVal = float(minMaxMean.split()[-1])
            newDose = (meanVal / counts) / pixelSize ** 2

            # Calculated as in Grigorieff paper --> https://dx.doi.org/10.7554/eLife.06980.001
            divByTwoFactor = 2 if zSlice.get(DIVIDED_BY_TWO, None) else 1
            return newDose * divByTwoFactor

        return newDose

    @staticmethod
    def _validateTSFromMdoc(mdoc, tsFile):
        errMsg = ''
        if not exists(tsFile):
            errMsg = '\nMdoc --> %s\nExpected tilt series file not found \n%s' \
                     % (mdoc, tsFile)

        return errMsg

    def _validateMdocInfoRead(self, ignoreFilesValidation=False):
        validateMdocContentsErrorMsgList = []
        msg = [f'\n{self._mdocFileName} is missing:\n']
        missingFiles = []
        missingAnglesIndices = []
        for i, tiltMetadata in enumerate(self._tiltsMetadata):
            # Check the angles
            if not tiltMetadata.getTiltAngle():
                missingAnglesIndices.append(str(i))
            # Check the angle stack files read
            if not ignoreFilesValidation:
                # Ignore the files validation is sometimes used for
                # test purposes
                file = tiltMetadata.getAngleMovieFile()
                if not exists(file):
                    missingFiles.append(file)

        if not self._voltage:
            msg.append(' - Voltage\n')
        if not self._magnification:
            msg.append(' - Magnification\n')
        if not self._samplingRate:
            msg.append(' - PixelSpacing\n')
        if not self.mdocHasDose:
            msg.append(' - Dose values. Related mdoc labels are: '
                       'ExposureDose or FrameDosesAndNumber or '
                       '(DoseRate and ExposureTime) or '
                       '(MinMaxMean and CountsPerElectron)\n')
        if not self.getTiltAxisAngle():
            msg.append(' - RotationAngle (tilt axis angle)\n')
        if missingAnglesIndices:
            msg.append(' - TiltAngle for Z values: %s\n' % ', '.join(missingAnglesIndices))
        if missingFiles:
            msg.append(' - Missing files: %s\n' % ', '.join(missingFiles))
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

    def getTiltsMetadata(self):
        return self._tiltsMetadata

    def getTiltAxisAngle(self):
        return self._tiltAxisAngle


class TiltMetadata:

    def __init__(self, angle=None, angleFile=None, acqOrder=None,
                 accumDose=None, incomingDose=None):
        self._angle = angle
        self._angleFile = angleFile
        self._acqOrder = acqOrder
        self._accumDose = accumDose
        self._incomingDose = incomingDose

    def setTiltAngle(self, tiltAngle):
        self._angle = tiltAngle

    def setAngleMovieFile(self, angleMovieFile):
        self._angleFile = angleMovieFile

    def setAcqOorder(self, order):
        self._acqOrder = order

    def setAccumDose(self, accumDose):
        self._accumDose = accumDose

    def setIncomingDose(self, incDose):
        self._incomingDose = incDose

    def getTiltAngle(self):
        return self._angle

    def getAngleMovieFile(self):
        return self._angleFile

    def getAcqOrder(self):
        return self._acqOrder

    def getAccumDose(self):
        return self._accumDose

    def getIncomingDose(self):
        return self._incomingDose


def normalizeTSId(rawTSId):
    """ Normalizes the name of a TS to prevent sqlite errors,
    it ends up as a table in a set"""
    # remove paths and extension
    normTSID = removeBaseExt(rawTSId)

    # Avoid dots case: TS_234.mrc.mdoc
    normTSID = normTSID.split(".")[0]

    if normTSID[0].isdigit():
        normTSID = TS_PREFIX + normTSID

    # Replace/remove from the TsId the special characters that may be problematic when registering the data in
    # the sqlite
    specialCharDict = {
        '-': '_',
        '.': '',
        '[': '',
        ']': '',
        '__': '_',
    }

    for specChar, okChar in specialCharDict.items():
        normTSID = normTSID.replace(specChar, okChar)

    return normTSID
