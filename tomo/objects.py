# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
import logging
from sqlite3 import OperationalError
from typing import Dict, Tuple, Any, Optional, Union

from pwem import ALIGN_NONE

logger = logging.getLogger(__name__)

import csv
import math
import os
import threading
from collections import OrderedDict
from datetime import datetime

import numpy as np
import pwem.objects.data as data
import pyworkflow.utils.path as path
import tomo.constants as const
from pwem.convert.transformations import euler_matrix
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.object import Integer, Float, String, Pointer, Boolean, CsvList


class MATRIX_CONVERSION:
    RELION = const.TR_RELION
    XMIPP = "xmipp"
    EMAN = const.TR_EMAN
    DYNAMO = const.TR_DYNAMO


def convertMatrix(M, convention=None, direction=None):
    """
            Parameters:
                - M --> Transformation matrix
                - convention --> One of the valid conventions to convert M. It can be:
                    * None         : Return the matrix stored in the metadata (Scipion convention).
                    * relion (str) : Relion matrix convention
                    * eman (str)   : Eman matrix convention
                - direction --> Determine how to perform the conversion (not considered if convention is None):
                    * 'get' (str)   : Convert the matrix stored in metadata (Scipion definition) to the given 'convention'
                    * 'set' (str)   : Convert the matrix from the given 'convention' to Scipion definition

            Scipion transformation matrix definition is described in detailed below:

               Notation:
                    - r      --> A 3D position vector
                    - f'(r)  --> Moved map
                    - f(r)   --> Reference
                    - M      --> Matrix to be stored by Scipion
                    - @      --> Matrix product

               Definition:

                   f'(r) = f(M@r)

               Example:

                   If r = (0,0,0,1) = o (origin vector) and M is a translation only transformation matrix
                   of the form:

                        M = [[1,0,0,0],[0,1,0,0],[0,0,0,1],[-dx,-dy,-dz,1]]

                   Being dx, dy, and dz a infinitesimally small displacement, then our transformation
                   verifies that:

                        f'(o) = f(M@o) = [-dx,-dy,-dz,1]

            We include conversions to the following matrix conventions:

                - Eman: Same as Scipion convention

                - Relion:
                    Notation:
                        - X'  -->  inverse of a matrix
                        - T   -->  translation matrix
                        - R   ---> rotation matrix
                        - M   -->  Scipion transformation matrix
                        - N   -->  Relion transformation matrix
                        - @   -->  Matrix multiplication

                    Conversion Scipion --> Relion
                        M = R@T' => T = M'@R
                        *** N = T@R = M'@R@R ***

                    Conversion Relion --> Scipion
                        N = M'@R@R => M' = N@R'@R' => *** M = R@R@N' ***
            """

    if convention is None or convention in [MATRIX_CONVERSION.EMAN, MATRIX_CONVERSION.DYNAMO]:
        return M
    elif direction == 'get' and convention in [MATRIX_CONVERSION.RELION, MATRIX_CONVERSION.XMIPP]:
        # Rotation matrix. Remove translation from the Scipion matrix
        R = np.eye(4)
        R[:3, :3] = M[:3, :3]
        Mi = np.linalg.inv(M)
        return Mi @ R @ R
    elif direction == 'set' and convention in [MATRIX_CONVERSION.RELION, MATRIX_CONVERSION.XMIPP]:
        # Rotation matrix. Remove translation from the Scipion matrix
        R = np.eye(4)
        R[:3, :3] = M[:3, :3]
        Mi = np.linalg.inv(M)
        return R @ R @ Mi


class TiltImageBase:
    """ Base class for TiltImageM and TiltImage. """
    TS_ID_FIELD = '_tsId'
    TILT_ANGLE_FIELD = '_tiltAngle'
    ACQ_ORDER_FIELD = '_acqOrder'
    INDEX_FIELD = '_index'

    def __init__(self, tsId=None, tiltAngle=None, acquisitionOrder=None, **kwargs):
        self._tiltAngle = Float(tiltAngle)
        self._tsId = String(tsId)
        self._acqOrder = Integer(acquisitionOrder)
        self._oddEvenFileNames = CsvList(pType=str)  # IMPORTANT: The odd is the first one and the even the second one

    def hasOddEven(self):
        return not self._oddEvenFileNames.isEmpty()

    def getOddEven(self):
        return self._oddEvenFileNames

    def getOdd(self):
        return self._oddEvenFileNames[0]

    def getEven(self):
        return self._oddEvenFileNames[1]

    def setOdd(self, fnOdd):
        self._oddEvenFileNames[0] = fnOdd

    def setEven(self, fnEven):
        self._oddEvenFileNames[1] = fnEven

    def setOddEven(self, listFileNames):
        self._oddEvenFileNames.set(listFileNames)

    def getTsId(self):
        """ Get unique TiltSerie ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def getTiltAngle(self):
        return self._tiltAngle.get()

    def setTiltAngle(self, value):
        self._tiltAngle.set(value)

    def getAcquisitionOrder(self):
        return self._acqOrder.get()

    def setAcquisitionOrder(self, value):
        self._acqOrder.set(value)

    def copyInfo(self, other, copyId=False, copyTM=True):
        self.copyAttributes(other, '_tiltAngle', '_tsId', '_acqOrder')
        if copyId:
            self.copyObjId(other)
        if copyTM and other.hasTransform():
            self.copyAttributes(other, '_transform')
        else:
            self.setTransform(None)

        if other.hasOddEven():
            self.copyAttributes(other, '_oddEvenFileNames')


class TiltImage(data.Image, TiltImageBase):
    """ Tilt image """

    def __init__(self, location=None, **kwargs):
        data.Image.__init__(self, location, **kwargs)
        TiltImageBase.__init__(self, **kwargs)

    def clone(self, copyEnable=True):
        return super().clone(copyEnable=copyEnable)

    def copyInfo(self, other, copyId=False, copyTM=True, copyStatus=True):
        data.Image.copyInfo(self, other)
        TiltImageBase.copyInfo(self, other, copyId=copyId, copyTM=copyTM)
        if copyStatus:
            self.setEnabled(other.isEnabled())

    def parseFileName(self, suffix="", extension=None):
        """
        This method returns the filename of the Tilt-Image adding a specified suffix and changing its extension.
        :param suffix: String to be added at the end of the location path (before extension).
        :param extension: String containing the new extension of the filename.
        :return: String containing the parsed filename with the specified suffix and extension.
        """

        fileName = os.path.basename(self.getFileName())
        fileName, fileExtension = os.path.splitext(fileName)

        if extension is not None:
            fileExtension = extension

        return fileName + suffix + fileExtension


TS_IGNORE_ATTRS = ['_mapperPath', '_size', '_hasAlignment', '_hasOddEven']


class TiltSeriesBase(data.SetOfImages):
    TS_ID_FIELD = '_tsId'
    ACQ_ORDER_FIELD = '_acqOrder'

    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)
        self._tsId = String(kwargs.get('tsId', None))
        # TiltSeries will always be used inside a SetOfTiltSeries
        # so, let's do not store the mapper path by default
        self._mapperPath.setStore(False)
        self._acquisition = TomoAcquisition()
        self._origin = Transform()
        self._anglesCount = Integer()
        self._hasAlignment = Boolean(False)
        self._hasOddEven = Boolean(False)
        self._interpolated = Boolean(False)
        self._ctfCorrected = Boolean(False)

    def getAnglesCount(self):
        return self._anglesCount

    def hasOddEven(self):
        return self._hasOddEven.get()

    def extractFileName(self, inputStr):
        return inputStr.split('@')[-1]

    def getOddFileName(self):
        firstItem = self.getFirstItem()
        return self.extractFileName(firstItem.getOdd())

    def getEvenFileName(self):
        firstItem = self.getFirstItem()
        return self.extractFileName(firstItem.getEven())

    def setAnglesCount(self, value):

        if isinstance(value, int):
            self._anglesCount.set(value)
        else:
            self._anglesCount = value

    def hasAlignment(self):

        return self._hasAlignment.get()

    def ctfCorrected(self):
        """ Returns true if ctf has been corrected"""
        return self._ctfCorrected.get()

    def setCtfCorrected(self, corrected):
        """ Sets the ctf correction status"""
        self._ctfCorrected.set(corrected)

    def interpolated(self):
        """ Returns true if tilt series has been interpolated"""
        return self._interpolated.get()

    def setInterpolated(self, interpolated):
        """ Sets the interpolation status of the tilt series"""
        self._interpolated.set(interpolated)

    def getTsId(self):
        """ Get unique TiltSeries ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def copyInfo(self, other, copyId=False):
        """ Copy basic information (id and other properties) but
        not _mapperPath or _size from other set of tilt series to current one.
        """
        self.copy(other, copyId=copyId, ignoreAttrs=TS_IGNORE_ATTRS)
        # self.copyAttributes(other, '_tsId', '_anglesCount')

    def write(self, properties=True):
        """ Do not save properties for this "Second level object"""

        super().write(properties=False)

    def append(self, tiltImage: TiltImageBase):
        tiltImage.setTsId(self.getTsId())
        data.SetOfImages.append(self, tiltImage)

        # TODO: Do it only once? Size =1?
        self._hasAlignment.set(tiltImage.hasTransform())
        self._hasOddEven.set(tiltImage.hasOddEven())

    def clone(self, ignoreAttrs=TS_IGNORE_ATTRS):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=ignoreAttrs)
        return clone

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfTiltSeries
        pass

    def getScannedPixelSize(self):
        mag = self._acquisition.getMagnification()
        return self._samplingRate.get() * 1e-4 * mag

    def generateTltFile(self, tltFilePath, reverse=False, excludeViews=False, presentAcqOrders=None):
        """ Generates an angle file in .tlt format in the specified location. If reverse is set to true the angles in
        file are sorted in the opposite order.
        :param tltFilePath: String containing the path where the file is created.
        :param reverse: Boolean indicating if the angle list must be reversed.
        :param excludeViews: boolean used to indicate if the tlt file should contain only the data concerning
        the non-excluded views (True) or all of them (False).
        :param presentAcqOrders: set containing the present acq orders in both the given TS and CTFTomoSeries. Used to
        filter the tilt angles that will be written in the tlt file generated. The parameter excludedViews is ignored
        if presentAcqOrders is provided, as the excluded views info may have been used to generate the presentAcqOrders
        (see tomo > utils > getCommonTsAndCtfElements)
        """

        if presentAcqOrders:
            angleList = [ti.getTiltAngle() for ti in self.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD) if
                         ti.getAcquisitionOrder() in presentAcqOrders]
        else:
            angleList = []
            for ti in self.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD):
                if excludeViews and not ti.isEnabled():
                    continue
                angleList.append(ti.getTiltAngle())
        if reverse:
            angleList.reverse()

        with open(tltFilePath, 'w') as f:
            f.writelines("%.3f\n" % angle for angle in angleList)


    def hasOrigin(self):
        """ Method indicating if the TiltSeries object has a defined origin. """

        return self._origin is not None

    def setOrigin(self, newOrigin):
        """ Method to set the origin of the TiltSeries object.
        :param newOrigin: Scipion Transform object indicating the origin to be set to the TiltSeries.
        """

        self._origin = newOrigin

    def getOrigin(self, force=False):
        """ Method to get the origin associated to the TiltSeries. If there is no origin associated to the object
        it may create a default one.
        :param force: Boolean indicating if the method must return a default origin in case the object has no one
        associated.
        """

        if self.hasOrigin():
            return self._origin
        else:
            if force:
                return self._getDefaultOrigin()
            else:
                return None

    def _getDefaultOrigin(self):
        sampling = self.getSamplingRate()
        t = Transform()
        x, y, z = self.getDim()
        if z > 1:
            z /= -2.

        t.setShifts(x / -2. * sampling, y / -2. * sampling, z * sampling)
        return t  # The identity matrix

    def getShiftsFromOrigin(self):
        """ Method to return the origin shift from the Scipion Transform object. """

        origin = self.getOrigin(force=True).getShifts()
        x = origin[0]
        y = origin[1]
        z = origin[2]
        return x, y, z
        # x, y, z are floats in Angstroms

    def updateOriginWithResize(self, resizeFactor):
        """ Method to update the origin after resizing the TiltSeries. """

        origin = self.getOrigin()

        xOri, yOri, zOri = self.getShiftsFromOrigin()

        origin.setShifts(xOri * resizeFactor,
                         yOri * resizeFactor,
                         zOri * resizeFactor)

        self.setOrigin(origin)
        # x, y, z are floats in Angstroms


def tiltSeriesToString(tiltSeries):
    s = []
    # Matrix info
    if tiltSeries.hasAlignment():
        s.append('+ali')

    # Interpolated
    if tiltSeries.interpolated():
        s.append('! interp')

    # CTF status
    if tiltSeries.ctfCorrected():
        s.append('+ctf')

    # Odd even associated
    if tiltSeries.hasOddEven():
        s.append('+oe')

    return (", " + ", ".join(s)) if len(s) else ""


class TiltSeries(TiltSeriesBase):
    ITEM_TYPE = TiltImage

    def __str__(self):

        s = super().__str__()

        # Matrix info
        s += tiltSeriesToString(self)

        return s

    def applyTransform(self, outputFilePath, swapXY=False):
        ih = ImageHandler()
        inputFilePath = self.getFirstItem().getFileName()
        newStack = True
        # TODO: Handle output tilt-series datatype format
        if self.getFirstItem().hasTransform():
            for index, ti in enumerate(self):
                if ti.hasTransform():
                    if newStack:
                        if swapXY:
                            ih.createEmptyImage(fnOut=outputFilePath,
                                                xDim=ti.getYDim(),
                                                yDim=ti.getXDim(),
                                                nDim=self.getSize())
                            newStack = False
                        else:
                            ih.createEmptyImage(fnOut=outputFilePath,
                                                xDim=ti.getXDim(),
                                                yDim=ti.getYDim(),
                                                nDim=self.getSize())
                            newStack = False
                    transform = ti.getTransform().getMatrix()
                    transformArray = np.array(transform)
                    if swapXY:
                        ih.applyTransform(inputFile=str(index + 1) + ':mrcs@' + inputFilePath,
                                          outputFile=str(index + 1) + '@' + outputFilePath,
                                          transformMatrix=transformArray,
                                          shape=(ti.getXDim(), ti.getYDim()))
                    else:
                        ih.applyTransform(inputFile=str(index + 1) + ':mrcs@' + inputFilePath,
                                          outputFile=str(index + 1) + '@' + outputFilePath,
                                          transformMatrix=transformArray,
                                          shape=(ti.getYDim(), ti.getXDim()))
                else:
                    raise Exception('ERROR: Some tilt-image is missing from transform object associated.')
        else:
            path.createAbsLink(os.path.abspath(inputFilePath), outputFilePath)

    def _dimStr(self):
        """ Return the string representing the dimensions. """

        return '%s x %s' % (self._firstDim[0],
                            self._firstDim[1])

    def getExcludedViewsIndex(self, caster=int, indexOffset=0):
        """Return a list with a list of the excluded views.

         :param caster: casting method to cast each index
         :param indexOffset: Value to add to the index. If you want to start the count in 0 pass -1"""
        excludeViewsList = []
        for ti in self.iterItems():
            if not ti.isEnabled():
                excludeViewsList.append(caster(ti.getIndex() + indexOffset))
        return excludeViewsList

    def _getExcludedViewsIndex(self):

        return self.getExcludedViewsIndex()

    def writeNewstcomFile(self, ts_folder, **kwargs):
        """Writes an artificial newst.com file"""
        newstcomPath = ts_folder + '/newst.com'
        pathi = self.getTsId()
        taperAtFill = kwargs.get('taperAtFill', (1, 0))
        offsetsInXandY = kwargs.get('offsetsInXandY', (0.0, 0.0))
        imagesAreBinned = kwargs.get('imagesAreBinned', 1.0)
        binByFactor = kwargs.get('binByFactor', 1)

        with open(newstcomPath, 'w') as f:
            f.write('$newstack -StandardInput\n\
InputFile {}.st\n\
OutputFile {}.ali\n\
TransformFile {}.xf\n\
TaperAtFill	{},{}\n\
AdjustOrigin\n\
OffsetsInXandY {},{}\n\
#DistortionField	.idf\n\
ImagesAreBinned	{}\n\
BinByFactor	{}\n\
#GradientFile {}.maggrad\n\
$if (-e ./savework) ./savework'.format(pathi, pathi, pathi,
                                       taperAtFill[0],
                                       taperAtFill[1], offsetsInXandY[0],
                                       offsetsInXandY[1], imagesAreBinned,
                                       binByFactor, pathi))

        return newstcomPath

    def writeTiltcomFile(self, ts_folder, **kwargs):
        """Writes an artificial tilt.com file"""
        tiltcomPath = ts_folder + '/tilt.com'
        pathi = self.getTsId()
        thickness = kwargs.get('thickness', 500)
        binned = kwargs.get('binned', 1)
        offset = kwargs.get('offset', 0.0)
        shift = kwargs.get('shift', (0.0, 0.0))
        radial = kwargs.get('radial', (0.35, 0.035))
        xAxisTill = kwargs.get('xAxisTill', 0.0)
        log = kwargs.get('log', 0.0)
        scale = kwargs.get('scale', (0.0, 1000.0))
        mode = kwargs.get('mode', 2)
        subsetStart = kwargs.get('subsetStart', (0, 0))
        actionIfGPUFails = kwargs.get('actionIfGPUFails', (1, 2))
        excludedViewsList = self.getExcludedViewsIndex(caster=str)
        excludedViewsIndexes = ''
        if excludedViewsList:
            excludedViewsIndexes = 'EXCLUDELIST %s \n' % ",".join(excludedViewsList)

        # The dimensions considered will be read, by default, from the corresponding tilt series. However, they
        # can be specified via th kwarg dims, as can be the case of a resized tomogram, in which the X and Y dimensions
        # considered in the tilt.com should be the ones corresponding to the tomogram
        intorducedDims = kwargs.get('dims', None)  #
        if intorducedDims:
            dims = intorducedDims
        else:
            dims = (self.getDim()[0], self.getDim()[1])
        # Swap dimensions case
        if kwargs.get('swapDims', False):
            dims = (dims[1], dims[0])

        with open(tiltcomPath, 'w') as f:
            f.write('$tilt -StandardInput\n\
InputProjections {}.ali\n\
OutputFile {}.rec\n\
IMAGEBINNED {} \n\
TILTFILE {}.tlt\n\
THICKNESS {}\n\
RADIAL {} {}\n\
FalloffIsTrueSigma 1\n\
XAXISTILT {}\n\
LOG	{}\n\
SCALE {} {}\n\
PERPENDICULAR\n\
MODE {}\n\
FULLIMAGE {} {}\n\
SUBSETSTART {} {}\n\
AdjustOrigin\n\
ActionIfGPUFails {},{}\n\
XTILTFILE {}.xtilt\n\
OFFSET {}\n\
SHIFT {} {}\n\
{}\
$if (-e ./savework) ./savework'.format(pathi, pathi, binned, pathi, thickness,
                                       radial[0], radial[1], xAxisTill, log,
                                       scale[0], scale[1], mode,
                                       dims[0], dims[1],
                                       subsetStart[0], subsetStart[1],
                                       actionIfGPUFails[0], actionIfGPUFails[1],
                                       pathi, offset, shift[0], shift[1],
                                       excludedViewsIndexes))

        return tiltcomPath

    def writeTltFile(self, ts_folder, excludeViews=False):
        """Writes a tlt file.
        :param ts_folder: path of the directory in which the tlt file will be generated.
        :param excludeViews: boolean used to indicate if the tlt file should contain only the data concerning
        the non-excluded views (True) or all of them (False).
        """
        xtiltPath = ts_folder + '/%s.tlt' % self.getTsId()
        with open(xtiltPath, 'w') as f:
            for ti in self.iterItems():
                if excludeViews and not ti.isEnabled():
                    continue
                f.write(str(ti.getTiltAngle()) + '\n')

    def writeXtiltFile(self, ts_folder):
        xtiltPath = ts_folder + '/%s.xtilt' % self.getTsId()
        with open(xtiltPath, 'w') as f:
            for ti in self:
                f.write('0.00\n')

    def writeXfFile(self, transformFilePath):
        """ This method takes a tilt series and the output transformation file
        path and creates an IMOD-based transform
        file in the location indicated. """

        tsMatrixTransformList = []

        for ti in self:
            if ti.getTransform() is not None:
                transform = ti.getTransform().getMatrix().flatten()
                transformIMOD = ['%.7f' % transform[0],
                                 '%.7f' % transform[1],
                                 '%.7f' % transform[3],
                                 '%.7f' % transform[4],
                                 '%.3f' % transform[2],
                                 '%.3f' % transform[5]]
            else:
                from pyworkflow.utils import yellowStr
                logging.info(
                    yellowStr('WARNING: The Tilt series lacks of alignment information (transformation matrices). '
                              'The identity transformation will be written in the .xf file'))

                #  This is the identity matrix
                transformIMOD = ['1.0000000',
                                 '0.0000000',
                                 '0.0000000',
                                 '1.0000000',
                                 '0.000',
                                 '0.000']
            tsMatrixTransformList.append(transformIMOD)

        with open(transformFilePath, 'w') as f:
            csvW = csv.writer(f, delimiter='\t')
            csvW.writerows(tsMatrixTransformList)

    def writeImodFiles(self, folderName, **kwargs):
        """Writes the following IMOD files:
        - newst.com
        - tilt.com
        - tlt
        - xf
        - xtilt
        :param folderName: path of the directory in which the files will be generated.
        :keyword tltIgnoresExcluded: boolean used to indicate if the tlt file should contain only the data concerning
        the non-excluded views (True) or all of them (False).
        """
        # Create a newst.com file
        self.writeNewstcomFile(folderName, **kwargs)
        # Create a tilt.com file
        self.writeTiltcomFile(folderName, **kwargs)
        # Create a .tlt file
        self.writeTltFile(folderName, excludeViews=kwargs.get('tltIgnoresExcluded', False))
        # Create a .xtilt file
        self.writeXtiltFile(folderName)
        # Create a .xf file
        transformFilePath = folderName + '/%s.xf' % self.getTsId()
        self.writeXfFile(transformFilePath)


class SetOfTiltSeriesBase(data.SetOfImages):
    EXPOSE_ITEMS = True
    USE_CREATE_COPY_FOR_SUBSET = True

    """ Base class for SetOfTiltImages and SetOfTiltImagesM.
    """

    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)
        self._anglesCount = Integer()
        self._acquisition = TomoAcquisition()
        self._hasAlignment = Boolean(False)
        self._hasOddEven = Boolean(False)
        self._ctfCorrected = Boolean(False)
        self._interpolated = Boolean(False)

    def getAcquisition(self):
        return self._acquisition

    def hasOddEven(self):
        return self._hasOddEven.get()

    def getAnglesCount(self):
        return self._anglesCount.get()

    def setAnglesCount(self, value):
        self._anglesCount.set(value)

    def ctfCorrected(self):
        """ Returns true if ctf has been corrected"""
        return self._ctfCorrected.get()

    def interpolated(self):
        """ Returns true if tilt series has been interpolated"""
        return self._interpolated.get()

    def hasAlignment(self):
        """ Returns true if at least one of its items has alignment information"""
        return self._hasAlignment.get()

    def copyInfo(self, other):
        """ Copy information (sampling rate and ctf)
        from other set of images to current one"""
        super().copyInfo(other)
        self.copyAttributes(other, '_anglesCount', '_hasAlignment', '_ctfCorrected', '_interpolated')

    def iterClassItems(self, iterDisabled=False):
        """ Iterate over the images of a class.
        Params:
            iterDisabled: If True, also include the disabled items. """
        for cls in self.iterItems():
            if iterDisabled or cls.isEnabled():
                for img in cls:
                    if iterDisabled or img.isEnabled():
                        yield img

    def _setItemMapperPath(self, item):
        """ Set the mapper path of this class according to the mapper
        path of the SetOfClasses and also the prefix according to class id
        """
        item._mapperPath.set('%s,%s' % (self.getFileName(), item.getTsId()))
        item.load()

    def _insertItem(self, item):
        """ Create the SetOfImages assigned to a class.
        If the file exists, it will load the Set.
        """
        self._setItemMapperPath(item)
        data.EMSet._insertItem(self, item)
        item.write(properties=False)  # Set.write(self)

    def __getitem__(self, itemId):
        """ Set the mapper of the TiltSerie (item) to point to the right table. The one with its own tilt images. """
        tiltSerie = data.SetOfImages.__getitem__(self, itemId)
        self._setItemMapperPath(tiltSerie)
        return tiltSerie

    def getFirstItem(self) -> TiltSeriesBase:
        classItem = data.EMSet.getFirstItem(self)
        self._setItemMapperPath(classItem)
        return classItem

    def iterItems(self, **kwargs) -> TiltSeriesBase:
        for item in data.EMSet.iterItems(self, **kwargs):
            self._setItemMapperPath(item)
            yield item

    def copyItems(self, inputTs,
                  orderByTs='id', updateTsCallback=None,
                  orderByTi='id', updateTiCallback=None,
                  itemSelectedCallback=None):
        """ Copy items (TiltSeries and TiltImages) from the input Set.
         Params:
            inputTs: input TiltSeries (or movies) from where to copy elements.
            orderByTs: optional orderBy value for iterating over TiltSeries
            updateTsCallback: optional callback after TiltSeries is created
            orderByTi: optional orderBy value for iterating over TiltImages
            updateTiCallback: optional callback after TiltImage is created
            itemSelectedCallback: Optional, callback receiving an item and
                returning true if it has to be copied
        """

        if itemSelectedCallback is None:
            itemSelectedCallback = data.SetOfImages.isItemEnabled

        for i, ts in enumerate(inputTs.iterItems(orderBy=orderByTs)):
            if itemSelectedCallback(ts):
                tsOut = self.ITEM_TYPE()
                tsOut.copyInfo(ts)
                tsOut.copyObjId(ts)
                if updateTsCallback:
                    updateTsCallback(i, ts, tsOut)
                self.append(tsOut)
                for j, ti in enumerate(ts.iterItems(orderBy=orderByTi)):
                    tiOut = tsOut.ITEM_TYPE()
                    tiOut.copyInfo(ti)
                    tiOut.setAcquisition(ti.getAcquisition())
                    tiOut.copyObjId(ti)
                    tiOut.setLocation(ti.getLocation())
                    tiOut.setEnabled(ti.isEnabled())  # Clone disabled tilt images
                    if updateTiCallback:
                        updateTiCallback(j, ts, ti, tsOut, tiOut)
                    tsOut.append(tiOut)

                self.update(tsOut)

    def update(self, item: TiltSeriesBase):

        self.setDim(item.getDim())
        self._anglesCount.set(item.getSize())
        self._hasAlignment.set(item.hasAlignment())
        self._interpolated.set(item.interpolated())
        self._ctfCorrected.set(item.ctfCorrected())
        self._hasOddEven.set(item.hasOddEven())

        super().update(item)

    def updateDim(self):
        """ Update dimensions of this set base on the first element. """

        logger.warning("TO DEVELOPERS: update is called always before this. This call to updateDim could be removed.")

        # firstItem = self.getFirstItem()
        # self.update(firstItem)

    def getScannedPixelSize(self):
        mag = self._acquisition.getMagnification()
        return self._samplingRate.get() * 1e-4 * mag

    def getTiltSeriesFromTsId(self, tsId):
        return self[{"_tsId": tsId}]

    def getTSIds(self):
        """ Returns al the Tilt series ids involved in the set."""
        return self.getUniqueValues(TiltSeries.TS_ID_FIELD)


class SetOfTiltSeries(SetOfTiltSeriesBase):
    ITEM_TYPE = TiltSeries

    def _dimStr(self):
        """ Return the string representing the dimensions. """

        s = '%s x %s x %s' % (self._anglesCount,
                              self._firstDim[0],
                              self._firstDim[1])
        s += tiltSeriesToString(self)

        return s


class TiltImageM(data.Movie, TiltImageBase):
    """ Tilt movie. """

    def __init__(self, location=None, **kwargs):
        data.Movie.__init__(self, location, **kwargs)
        TiltImageBase.__init__(self, **kwargs)

    def copyInfo(self, other, copyId=False, copyTM=True):
        data.Movie.copyInfo(self, other)
        TiltImageBase.copyInfo(self, other, copyId=copyId, copyTM=copyTM)


class TiltSeriesM(TiltSeriesBase):
    ITEM_TYPE = TiltImageM


class SetOfTiltSeriesM(SetOfTiltSeriesBase):
    ITEM_TYPE = TiltSeriesM

    def __init__(self, **kwargs):
        SetOfTiltSeriesBase.__init__(self, **kwargs)
        self._gainFile = String()
        self._darkFile = String()
        # Store the frames range to avoid loading the items
        self._firstFramesRange = data.FramesRange()

    def setGain(self, gain):
        self._gainFile.set(gain)

    def getGain(self):
        return self._gainFile.get()

    def setDark(self, dark):
        self._darkFile.set(dark)

    def getDark(self):
        return self._darkFile.get()

    def getFramesRange(self):
        return self._firstFramesRange

    def setFramesRange(self, value):
        self._firstFramesRange.set(value)

    def copyInfo(self, other):
        """ Copy SoM specific information plus inherited """
        SetOfTiltSeriesBase.copyInfo(self, other)
        self._gainFile.set(other.getGain())
        self._darkFile.set(other.getDark())
        self._firstFramesRange.set(other.getFramesRange())

    def _dimStr(self):
        """ Return the string representing the dimensions. """

        return '%s x %s' % (self._anglesCount, self._firstDim)


class TiltSeriesDict:
    """ Helper class that to store TiltSeries and TiltImage but
    using dictionaries for quick access.
    This class also contains some logic related to the streaming:
    - Check for new input items that needs to be processed
    - Check for items already done that needs to be saved.
    """

    def __init__(self, inputSet=None, outputSet=None,
                 newItemsCallback=None,
                 doneItemsCallback=None):
        """
        Initialize the dict.
        :param inputSet: The set with input items. It will be monitored
            for new items from streaming.
        :param newItemsCallback: When new items are discovered, this
            function will be called
        :param doneItemsCallback: When some items are done, this function
            will be called.
        """
        self.__dict = OrderedDict()
        self.__inputSet = inputSet
        if inputSet is not None:
            self.__inputClosed = inputSet.isStreamClosed()
        self.__lastCheck = None
        self.__finalCheck = False
        self.__newItemsCallback = newItemsCallback
        self.__doneItemsCallback = doneItemsCallback

        self.__new = set()
        self.__finished = set()  # Reported as finished tasks, but not saved
        self.__done = set()  # Finished and saved tasks
        self.__lock = threading.Lock()

        if outputSet is not None:
            for ts in outputSet:
                # We don't need tilt-images for done items
                self.addTs(ts, includeTi=False)
                self.__done.add(ts.getTsId())

    def addTs(self, ts, includeTi=False):
        """ Add a clone of the tiltseries. """
        self.__dict[ts.getTsId()] = (ts.clone(), OrderedDict())
        if includeTi:
            for ti in ts:
                self.addTi(ti)

    def hasTs(self, tsId):
        return tsId in self.__dict

    def getTs(self, tsId):
        return self.__dict[tsId][0]

    def addTi(self, ti):
        self.getTiDict(ti.getTsId())[ti.getObjId()] = ti.clone()

    def getTi(self, tsId, tiObjId):
        return self.getTiDict(tsId)[tiObjId]

    def getTiDict(self, tsId):
        return self.__dict[tsId][1]

    def getTiList(self, tsId):
        return list(self.getTiDict(tsId).values())

    def __iter__(self):
        for ts, d in self.__dict.values():
            yield ts

    # ---- Streaming related methods -------------
    def update(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        logging.debug("TiltSeriesDict._checkNewInput called.")

        inputSetFn = self.__inputSet.getFileName()
        mTime = datetime.fromtimestamp(os.path.getmtime(inputSetFn))
        # if self.__lastCheck:
        # print('Last check: %s, modification: %s'
        #       % (pwutils.prettyTime(self.__lastCheck),
        #          pwutils.prettyTime(mTime)))

        if self.__lastCheck is None or self.__lastCheck <= mTime:
            updatedSet = self.__inputSet.getClass()(filename=inputSetFn)
            updatedSet.loadAllProperties()
            newItems = []
            for ts in updatedSet:
                if not self.hasTs(ts.getTsId()):
                    self.addTs(ts, includeTi=True)
                    newItems.append(ts.getTsId())
            self.__inputClosed = updatedSet.isStreamClosed()
            updatedSet.close()
            if newItems:
                self.__newItemsCallback(newItems)
        self.__lastCheck = datetime.now()

    def _checkNewOutput(self):
        logger.debug("TiltSeriesDict._checkNewOutput")

        # First check that we have some items in the finished
        self.__lock.acquire()
        doneItems = list(self.__finished)
        self.__finished.clear()
        self.__lock.release()

        if doneItems or (self.allDone() and not self.__finalCheck):
            self.__done.update(doneItems)
            self.__doneItemsCallback(doneItems)
            if self.allDone():
                self.__finalCheck = True

    def setFinished(self, *tsIdList):
        """ Notify that all TiltSeries in the list of ids are finished. """
        self.__lock.acquire()
        self.__finished.update(tsIdList)
        self.__lock.release()

    def allDone(self):
        """ Return True if input stream is closed and all task are done. """
        # print(">>> DEBUG: allDone\n"
        #       "    inputClosed: %s\n"
        #       "    len(dict):   %s\n"
        #       "    len(done):   %s" % (self.__inputClosed, len(self.__dict),
        #                                len(self.__done)))
        return self.__inputClosed and len(self.__dict) == len(self.__done)


class TomoAcquisition(data.Acquisition):
    """ Tomography acquisition metadata object"""

    def __init__(self, angleMin=None, angleMax=None, step=None,
                 accumDose=None, tiltAxisAngle=None, **kwargs):
        data.Acquisition.__init__(self, **kwargs)
        self._angleMin = Float(angleMin)
        self._angleMax = Float(angleMax)
        self._step = Float(step)
        self._accumDose = Float(accumDose)
        self._tiltAxisAngle = Float(tiltAxisAngle)

    def getTiltAxisAngle(self):
        return self._tiltAxisAngle.get()

    def setTiltAxisAngle(self, value):
        self._tiltAxisAngle.set(value)

    def getAngleMax(self):
        return self._angleMax.get()

    def setAngleMax(self, value):
        self._angleMax.set(value)

    def getAngleMin(self):
        return self._angleMin.get()

    def setAngleMin(self, value):
        self._angleMin.set(value)

    def getStep(self):
        return self._step.get()

    def setStep(self, value):
        return self._step.set(value)

    def getAccumDose(self):
        return self._accumDose.get()

    def setAccumDose(self, value):
        self._accumDose.set(value)


class Tomogram(data.Volume):
    """ Class to hold the tomogram abstraction inside Scipion. The origin (self._origin) of the volume is set as the
    location of the first coordinate loaded from the binary file. The volume may be displaced by setting a different
    origin using the methods implemented in the inherited class data.Image in scipion-em plugin.
    """
    TS_ID_FIELD = '_tsId'
    ORIGIN_MATRIX_FIELD = '_origin._matrix'

    def __init__(self, **kwargs):
        data.Volume.__init__(self, **kwargs)
        self._acquisition = None
        self._tsId = String(kwargs.get('tsId', None))
        self._dim = None
        self._ctfCorrected = Boolean(False)

    def getTsId(self):
        """ Get unique TiltSeries ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def getAcquisition(self):
        return self._acquisition

    def setAcquisition(self, acquisition):
        self._acquisition = acquisition

    def hasAcquisition(self):
        return (self._acquisition is not None
                and self._acquisition.getAngleMin() is not None
                and self._acquisition.getAngleMax() is not None)

    def copyInfo(self, other):
        """ Copy basic information """
        super().copyInfo(other)
        self.copyAttributes(other, '_acquisition', self.TS_ID_FIELD)
        if other.hasOrigin():
            self.copyAttributes(other, '_origin')

    def ctfCorrected(self):
        """ Returns true if ctf has been corrected"""
        return self._ctfCorrected.get()

    def setCtfCorrected(self, corrected):
        """ Sets the ctf correction status"""
        self._ctfCorrected.set(corrected)


class SetOfTomograms(data.SetOfVolumes):
    ITEM_TYPE = Tomogram
    EXPOSE_ITEMS = False

    def __init__(self, *args, **kwargs):
        data.SetOfVolumes.__init__(self, **kwargs)
        self._acquisition = TomoAcquisition()
        self._hasOddEven = Boolean(False)
        self._ctfCorrected = Boolean(False)

    def hasOddEven(self):
        return self._hasOddEven.get()

    def updateDim(self):
        """ Update dimensions of this set base on the first element. """
        self.setDim(self.getFirstItem().getDim())

    def __str__(self):
        sampling = self.getSamplingRate()
        tomoStr = "+oe," if self.hasOddEven() else ''
        # CTF status
        if self.ctfCorrected():
            tomoStr += '+ctf,'

        if not sampling:
            logger.error("FATAL ERROR: Object %s has no sampling rate!!!"
                         % self.getName())
            sampling = -999.0

        s = "%s (%d items, %s, %s %0.2f Å/px%s)" % \
            (self.getClassName(), self.getSize(),
             self._dimStr(), tomoStr, sampling, self._appendStreamState())
        return s

    def update(self, item: Tomogram):
        self._hasOddEven.set(item.hasHalfMaps())
        self.setCtfCorrected(item.ctfCorrected())
        super().update(item)

    def getTSIds(self):
        """ Returns al the Tilt series ids involved in the set."""
        return self.getUniqueValues(Tomogram.TS_ID_FIELD)

    def ctfCorrected(self):
        """ Returns true if ctf has been corrected"""
        return self._ctfCorrected.get()

    def setCtfCorrected(self, corrected):
        """ Sets the ctf correction status"""
        self._ctfCorrected.set(corrected)


class TomoMask(Tomogram):
    """ Object used to represent segmented tomograms
    """

    def __init__(self, **kwargs):
        Tomogram.__init__(self, **kwargs)
        self._volName = String()

    def getVolName(self):
        """ Get the reference tomogram file for the current tomoMask.
        """
        return self._volName.get()

    def setVolName(self, tomoName):
        """ Set the reference tomogram file for the current tomoMask.
        """
        self._volName.set(tomoName)

    def getTomogram(self):
        """ Generate the reference tomogram object for the current tomoMask.
        """
        tomo = Tomogram()
        tomo.setLocation(self.getVolName())
        tomo.setSamplingRate(self.getSamplingRate())
        tomo.setAcquisition(self.getAcquisition())
        return tomo


class SetOfTomoMasks(SetOfTomograms):
    ITEM_TYPE = TomoMask
    EXPOSE_ITEMS = True


class Coordinate3D(data.EMObject):
    """This class holds the (x,y,z) position and other information
    associated with a coordinate"""

    TOMO_ID_ATTR = "_tomoId"
    GROUP_ID_ATTR = "_groupId"
    SCORE_ATTR = "_score"

    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._boxSize = 0
        self._volumePointer = Pointer(objDoStore=False)
        self._x = Float()
        self._y = Float()
        self._z = Float()
        self._volId = Integer()
        self._eulerMatrix = data.Transform()
        self._groupId = Integer(0)  # This may refer to a mesh, ROI, vesicle or any group of coordinates
        self._tomoId = String(kwargs.get('tomoId', None))  # Used to access to the corresponding tomogram from each
        self._score = Float(0)
        # coord (it's the tsId)

    def _getOffset(self, dim, originFunction=const.SCIPION):
        """ Returns the offset to apply to a one of the coordinates
        :param dim integer to get the dimension from (X=0, Y=1, Z=2)
        :param originFunction function to call to do the conversion
        """
        if originFunction == const.SCIPION:
            return 0

        origin_Scipion = self.getVolumeOrigin()[dim]
        aux = originFunction(self.getVolume().getDim())
        origin = aux[dim] if aux is not None else -origin_Scipion
        return origin + origin_Scipion

    def getX(self, originFunction):
        """ See getPosition method for a full description of how "originFunction"
        works"""
        return self._x.get() - self._getOffset(0, originFunction)

    def setX(self, x, originFunction):
        """ See setPosition method for a full description of how "originFunction"
        works"""
        self._x.set(x + self._getOffset(0, originFunction))

    def shiftX(self, shiftX):
        self._x.sum(shiftX)

    def getY(self, originFunction):
        """ See getPosition method for a full description of how "originFunction"
        works"""

        return self._y.get() - self._getOffset(1, originFunction)

    def setY(self, y, originFunction):
        """ See setPosition method for a full description of how "originFunction"
        works"""

        self._y.set(y + self._getOffset(1, originFunction))

    def shiftY(self, shiftY):
        self._y.sum(shiftY)

    def getZ(self, originFunction):
        """ See getPosition method for a full description of how "originFunction"
        works"""
        return self._z.get() - self._getOffset(2, originFunction)

    def setZ(self, z, originFunction):
        """ See setPosition method for a full description of how "originFunction"
        works"""
        self._z.set(z + self._getOffset(2, originFunction))

    def shiftZ(self, shiftZ):
        self._z.sum(shiftZ)

    def setMatrix(self, matrix, convention=None):
        self._eulerMatrix.setMatrix(convertMatrix(matrix, direction=const.SET, convention=convention))

    def getMatrix(self, convention=None):
        return convertMatrix(self._eulerMatrix.getMatrix(), direction=const.GET, convention=convention)

    def hasTransform(self):
        return self._eulerMatrix is not None

    def euler2Matrix(self, r, p, y):
        # FIXME: Queremos mantener esta conversion? Puede ser muy lioso
        self._eulerMatrix.setMatrix(euler_matrix(r, p, y))

    def eulerAngles(self):
        R = self.getMatrix()
        sy = math.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
        singular = sy < 1e-6
        if not singular:
            x = math.atan2(R[2, 1], R[2, 2])
            y = math.atan2(-R[2, 0], sy)
            z = math.atan2(R[1, 0], R[0, 0])

        else:
            x = math.atan2(-R[1, 2], R[1, 1])
            y = math.atan2(-R[2, 0], sy)
            z = 0

        return np.array([x, y, z])

    def scale(self, factor):
        """ Scale x, y and z coordinates by a given factor.
        """
        self._x.multiply(factor)
        self._y.multiply(factor)
        self._z.multiply(factor)

    def getPosition(self, originFunction):
        """Get the position a Coordinate3D refered to a given origin defined by originFunction.
        The input of the method is a funtion (originFunction) which moves the coordinate
        position refered to the bottom left corner to other origin (retrieved by originFunction) in the grid.

        Parameters:

            :param function originFunction: Function to return a Vector to refer a coordinate to the bottom left corner
                                            from a given convention.

        Example:

            >>> origin = originFunction((Lx, Ly, Lz))
            >>> (vx, vy, vz)  # Vector to refer (x,y,z) coordinate to an origin from the bottom left corner

        Firstly, the Scipion origin vector stored in the Tomogram associated to the Coordinate3D
        will be applied to refer the current coordinate to the bottom left coordinate of the Tomogram.
        """
        return self.getX(originFunction), self.getY(originFunction), self.getZ(originFunction)

    def setPosition(self, x, y, z, originFunction):
        """Set the position of the coordinate to be saved in the Coordinate3D object.
        The inputs of the method are the (x,y,z) position of the coordinate and a
        funtion (originFunction) which moves the current position to the bottom left
        corner of the Tomogram with dimensions Lx, Ly, Lz.

        Parameters:

            :param float x: Position of the coordinate in the X axis
            :param float y: Position of the coordinate in the Y axis
            :param float z: Position of the coordinate in the Z axis
            :param function originFunction: Function to return a Vector to refer a coordinate to the bottom left corner from a
                                            given convention.

        Example:

            >>> origin = originFunction((Lx, Ly, Lz))
            >>> (vx, vy, vz)  # Vector to refer (x,y,z) coordinate to the bottom left corner

        In this way, it is possible to apply the Scipion origin vector stored in the
        Tomogram associated to the Coordinate3D which moves the positions referred to
        the bottom left corner of a grid to the center of gravity of the grid (or any
        other origin specified by the user).

        IMPORTANT NOTE: For this method to work properly, it is needed to associate the Tomogram
        before doing a call to this method.

        Example:

            >>> coord = Coordinate3D()
            >>> coord.setPosition(x, y, z, originFunction)
            >>> Error: Tomogram is still NoneType
            >>> coord.setVolume(Tomogram)
            >>> coord.setPosition(x, y, z, originFunction)
            >>> Exit: Everything runs normally

        This requirement is only needed for "setPostion" method. The remaining attributes of the object
        can be set either before or after calling "setVolume" method.
        """
        self.setX(x, originFunction)
        self.setY(y, originFunction)
        self.setZ(z, originFunction)

    def getVolume(self):
        """ Return the tomogram object to which
        this coordinate is associated.
        """
        return self._volumePointer.get()

    def setVolume(self, volume):
        """ Set the micrograph to which this coordinate belongs. """
        self._volumePointer.set(volume)
        self._volId.set(volume.getObjId())
        if volume.getTsId():  # See getCoordinate3D() --> as a tomo is necessary to be created, the tomoId (tsId),
            # which may have been previously stored is deleted when calling setVolume
            self.setTomoId(volume.getTsId())

    def setBoxSize(self, boxSize):
        logger.info('Deprecated, use SetOfCoordinates3D box size instead.')
        self._boxSize = boxSize

    def getBoxSize(self):
        logger.info('Deprecated, use SetOfCoordinates3D box size instead.')
        return self._boxSize

    def getVolId(self):
        return self._volId.get()

    def setVolId(self, volId):
        self._volId.set(volId)

    def invertY(self):
        if not self.getVolume() is None:
            dims = self.getVolume().getDim()
            height = dims[1]
            self.setY(height - self.getY(const.SCIPION), const.SCIPION)
        # else: error TODO

    def getVolName(self):
        return self.getVolume().getFileName()

    def getGroupId(self):
        return self._groupId.get()

    def setGroupId(self, groupId):
        self._groupId.set(groupId)

    def hasGroupId(self):
        return self._groupId is not None

    def getVolumeOrigin(self, angstrom=False):
        """Return the vector that can be used to move the position of the Coordinate3D
        (referred to the center of the Tomogram or other origin specified by the user)
        to the bottom left corner of the Tomogram
        """
        vol = self.getVolume()
        if not vol:
            raise Exception("3D coordinate must be referred to a volume to get its origin.")

        # Tomogram origin
        origin = vol.getShiftsFromOrigin()

        if angstrom:
            return origin
        else:
            sr = vol.getSamplingRate()

            return origin[0] / sr, origin[1] / sr, origin[2] / sr

    def getTomoId(self):
        return self._tomoId.get()

    def setTomoId(self, tomoId):
        self._tomoId.set(tomoId)

    def getScore(self):
        return self._score.get()

    def setScore(self, val):
        self._score.set(val)

    def composeCoordId(self, sampligRate):
        return "%s,%s,%s,%s" % (self.getTomoId(),
                                int(sampligRate * self._x.get()),
                                int(sampligRate * self._y.get()),
                                int(sampligRate * self._z.get()))

    def __str__(self):
        return "%s, G%s" % (self.composeCoordId(1), self.getGroupId())


class SetOfCoordinates3D(data.EMSet):
    """ Encapsulate the logic of a set of volumes coordinates.
    Each coordinate has a (x,y,z) position and is related to a Volume
    The SetOfCoordinates3D can also have information about TiltPairs.
    """
    ITEM_TYPE = Coordinate3D

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
        self._boxSize = Integer()
        self._samplingRate = Float()
        self._precedentsPointer = Pointer()
        self._tomos = None

    def getBoxSize(self):
        """ Return the box size of the particles.
        """
        return self._boxSize.get()

    def setBoxSize(self, boxSize):
        """ Set the box size of the particles. """
        self._boxSize.set(boxSize)

    def getSamplingRate(self):
        """ Return the sampling rate of the particles. """
        return self._samplingRate.get()

    def setSamplingRate(self, sampling):
        """ Set the sampling rate of the particles. """
        self._samplingRate.set(sampling)

    def iterVolumes(self):
        """ Iterate over the objects set associated with this
        set of coordinates.
        """
        return self.getPrecedents()

    def iterVolumeCoordinates(self, volume):
        """ Iterates over the set of coordinates belonging to that micrograph.
        """
        pass

    def iterCoordinates(self, volume: Tomogram = None, orderBy='id') -> Coordinate3D:
        """ Iterate over the coordinates associated with a tomogram.
        If volume=None, the iteration is performed over the whole
        set of coordinates.

        IMPORTANT NOTE: During the storing process in the database, Coordinates3D will lose their
        pointer to ther associated Tomogram. This method overcomes this problem by retrieving and
        relinking the Tomogram as if nothing would ever happened.

        It is recommended to use this method when working with Coordinates3D, being the common
        "iterItems" deprecated for this set.

        Example:

            >>> for coord in coordSet.iterItems()
            >>>     print(coord.getVolName())
            >>>     Error: Tomogram associated to Coordinate3D is NoneType (pointer lost)
            >>> for coord in coordSet.iterCoordinates()
            >>>     print(coord.getVolName())
            >>>     '/path/to/Tomo.file' retrieved correctly

        """

        # Iterate over all coordinates if tomoId is None,
        # otherwise use tomoId to filter the where selection
        if volume is None:
            coordWhere = '1'
        elif isinstance(volume, int):
            logger.warning("FOR DEVELOPERS: Do not use volId, use volName")
            coordWhere = '_volId=%d' % volume
        elif isinstance(volume, Tomogram):
            coordWhere = '%s="%s"' % (Coordinate3D.TOMO_ID_ATTR, volume.getTsId())
        else:
            raise Exception('Invalid input tomogram of type %s'
                            % type(volume))

        # Iterate over all coordinates if tomoId is None,
        # otherwise use tomoId to filter the where selection
        for coord in self.iterItems(where=coordWhere, orderBy=orderBy):
            # Associate the tomogram
            self._associateVolume(coord)
            yield coord

    def _getTomogram(self, tsId):
        """ Returns  the tomogram from a tsId"""
        tomos = self._getTomograms()

        if tsId not in tomos.keys():
            tomo = self.getPrecedents()[{Tomogram.TS_ID_FIELD: tsId}]
            self._tomos[tsId] = tomo
            return tomo
        else:
            return self._tomos[tsId]

    def _getTomograms(self):
        if self._tomos is None:
            self.initTomos()

        return self._tomos

    def getPrecedents(self):
        """ Returns the SetOfTomograms associated with
                this SetOfCoordinates"""
        return self._precedentsPointer.get()

    def getPrecedent(self, tomoId):
        return self.getPrecedentsInvolved()[tomoId]

    def setPrecedents(self, precedents):
        """ Set the tomograms  or Tilt Series associated with this set of coordinates.
                Params:
                    tomograms: Either a SetOfTomograms or Tilt Series object or a pointer to it.
                """
        if precedents.isPointer():
            self._precedentsPointer.copy(precedents)
        else:
            self._precedentsPointer.set(precedents)

    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths

    def getSummary(self):
        summary = []
        summary.append("Number of particles picked: %s" % self.getSize())
        summary.append("Particle size: %s" % self.getBoxSize())
        return "\n".join(summary)

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not _mapperPath or _size
        from other set of objects to current one.
        """
        self.setBoxSize(other.getBoxSize())
        self.setSamplingRate(other.getSamplingRate())
        self.setPrecedents(other.getPrecedents())

    def __str__(self):
        """ String representation of a set of coordinates. """
        if self._boxSize.hasValue():
            boxSize = self._boxSize.get()
            boxStr = ' %d x %d x %d' % (boxSize, boxSize, boxSize)
        else:
            boxStr = 'No-Box'
        sRate = '%.2f' % self.getSamplingRate() if self.getSamplingRate() is not None else 'None!!!'
        s = "%s (%d items, %s, %s Å/px%s)" % (self.getClassName(), self.getSize(), boxStr,
                                              sRate, self._appendStreamState())

        return s

    def getFirstItem(self):
        coord = data.EMSet.getFirstItem(self)
        self._associateVolume(coord)
        return coord

    def _associateVolume(self, coord):
        coord.setVolume(self._getTomogram(coord.getTomoId()))

    def __getitem__(self, itemId):
        """Add a pointer to a Tomogram before returning the Coordinate3D"""
        coord = data.EMSet.__getitem__(self, itemId)
        self._associateVolume(coord)
        return coord

    def initTomos(self):
        """ Initialize internal _tomos to a dictionary if not done already"""
        if self._tomos is None:
            self._tomos = dict()

    def getPrecedentsInvolved(self):
        """ Returns a list  with only the tomograms involved in the subtomograms. May differ when
        subsets are done."""

        tomoId_attr = Coordinate3D.TOMO_ID_ATTR

        uniqueTomos = self.aggregate(['count'], tomoId_attr, [tomoId_attr])

        for row in uniqueTomos:
            tsId = row[tomoId_attr]
            # This should register the tomogram in the internal _tomos
            self._getTomogram(tsId)

        return self._tomos

    def append(self, item: Coordinate3D):
        if self.getBoxSize() is None and item._boxSize:
            self.setBoxSize(item._boxSize)
        super().append(item)

    def getTSIds(self):
        """ Returns all the TS ID (tomoId) present in this set"""
        return self.getUniqueValues(Coordinate3D.TOMO_ID_ATTR)


class SubTomogram(data.Volume):
    """The coordinate associated to each subtomogram is not scaled. To do that, the coordinates and the subtomograms
    sampling rates should be compared (because of how the extraction protocol works). But when shifts are applied to
    the coordinates, it has to be considered that if we're operating with coordinates coming from subtomogrmas, those
    shifts will be scaled, but if the coordinates come from coordinates, they won't be."""

    VOL_NAME_FIELD = "_volName"
    COORD_VOL_NAME_FIELD = "_coordinate.%s" % Coordinate3D.TOMO_ID_ATTR

    def __init__(self, **kwargs):
        data.Volume.__init__(self, **kwargs)
        self._acquisition = None
        self._volId = Integer()
        # This coordinate is NOT SCALED. To do that, the coordinates and subtomograms sampling rates
        # should be compared (because of how the extraction protocol works)
        self._coordinate = None
        self._volName = String()

    def hasCoordinate3D(self):
        return self._coordinate is not None

    def setCoordinate3D(self, coordinate):
        self._coordinate = coordinate
        self.setVolId(coordinate.getVolId())

    def getCoordinate3D(self) -> Coordinate3D:
        """Since the object Coordinate3D needs a volume, use the information stored in the
        SubTomogram to reconstruct the corresponding Tomogram associated to its Coordinate3D"""
        # We do not do this here but in the set iterator tha will "plug" the volume (tomogram) is exists
        # tomo = Tomogram()
        # subtomoOrigin = self.getOrigin()
        # if subtomoOrigin:
        #     tomo.setOrigin(subtomoOrigin)
        # tomo.setLocation(self.getVolName())
        # tomo.setSamplingRate(self.getSamplingRate())
        # coord = self._coordinate
        # coord.setVolume(tomo)
        # return coord
        return self._coordinate

    def getAcquisition(self):
        return self._acquisition

    def setAcquisition(self, acquisition):
        self._acquisition = acquisition

    def hasAcquisition(self):
        return self._acquisition is not None and \
            self._acquisition.getAngleMin() is not None and \
            self._acquisition.getAngleMax() is not None

    def getVolId(self):
        """ Return the tomogram id if the coordinate is not None.
        or have set the _volId property.
        """
        if self._volId.hasValue():
            return self._volId.get()
        if self.hasCoordinate3D():
            return self.getCoordinate3D().getVolId()

        return None

    def setVolId(self, volId):
        self._volId.set(volId)

    def getVolName(self):
        """ Return the tomogram filename if the coordinate is not None.
        or have set the _volName property.
        """
        if self._volName.hasValue():
            return self._volName.get()

        if self.hasCoordinate3D():
            return self.getCoordinate3D().getVolName()

        return "Missing"

    def setVolName(self, volName):
        self._volName.set(volName)

    def getVolumeOrigin(self, angstrom=False):
        """Return the vector that can be used to move the position of the Coordinate3D
        associated to the SubTomogram (referred to the center of the Tomogram or other
        origin specified by the user) to the bottom left corner of the Tomogram
        """
        if angstrom:
            return self.getShiftsFromOrigin()
        else:
            sr = self.getSamplingRate()
            origin = self.getShiftsFromOrigin()
            return int(origin[0] / sr), int(origin[1] / sr), int(origin[2] / sr)

    def setTransform(self, newTransform, convention=None):
        if newTransform is None:
            newTransform = Transform()
            matrix = np.eye(4)
        else:
            matrix = newTransform.getMatrix()
        newTransform.setMatrix(convertMatrix(matrix, direction=const.SET, convention=convention))
        self._transform = newTransform

    def getTransform(self, convention=None):

        if convention is not None:
            matrix = self._transform.getMatrix()
            return Transform(convertMatrix(matrix, direction=const.GET, convention=convention))
        else:
            return self._transform


class SetOfSubTomograms(data.SetOfVolumes):
    ITEM_TYPE = SubTomogram
    REP_TYPE = SubTomogram
    EXPOSE_ITEMS = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._acquisition = TomoAcquisition()
        self._coordsPointer = Pointer()
        self._tomos = None

    def copyInfo(self, other):
        """ Copy basic information (sampling rate and ctf)
        from other set of images to current one"""
        super().copyInfo(other)
        if hasattr(other, '_coordsPointer'):  # Like the vesicles in pyseg
            self.copyAttributes(other, '_coordsPointer')

    def append(self, subtomo: SubTomogram):
        # Set the alignment attribute value when adding the first element to the set
        if self.isEmpty():
            if subtomo.hasTransform():
                trMatrix = subtomo.getTransform().getMatrix()
                hasNonEyeTransform = not np.allclose(trMatrix, np.eye(4))
                if hasNonEyeTransform:
                    self.setAlignment3D()
            else:
                self.setAlignment(ALIGN_NONE)

        super().append(subtomo)

    def hasCoordinates3D(self):
        return self._coordsPointer.hasValue()

    def getCoordinates3D(self, asPointer=False):
        """ Returns the SetOfCoordinates associated with
        this SetOfSubTomograms"""

        return self._coordsPointer if asPointer else self._coordsPointer.get()

    def setCoordinates3D(self, coordinates):
        """ Set the SetOfCoordinates associated with
        this set of particles.
         """
        if isinstance(coordinates, Pointer):
            self._coordsPointer = coordinates
        else:
            self._coordsPointer.set(coordinates)

    def iterCoordinates(self, volume: Tomogram = None, orderBy='id') -> Coordinate3D:
        """ Mimics SetOfCoordinates.iterCoordinates so can be passed to viewers or protocols transparently"""
        if self.hasCoordinates3D():
            for subtomo in self.iterSubtomos(volume, orderBy=orderBy):
                coord = subtomo.getCoordinate3D()
                coord.setObjId(subtomo.getObjId())
                yield coord
        else:
            yield

    def iterSubtomos(self, volume: Tomogram = None, orderBy='id') -> SubTomogram:
        """ Iterates over the sutomograms, enriching them with the related tomogram if apply so coordinate getters and setters will work
        If volume=None, the iteration is performed over the whole
        set of subtomograms.

        IMPORTANT NOTE: During the storing process in the database, Coordinates3D will lose their
        pointer to the associated Tomogram. This method overcomes this problem by retrieving and
        relinking the Tomogram as if nothing would ever happend.

        It is recommended to use this method when working with subtomograms, anytime you want to properly use
        its coordinate3D attached object.

        Example:

            >>> for subtomo in subtomos.iterItems()
            >>>     print(subtomo.getCoordinate3D().getX(SCIPION))
            >>>     Error: Tomogram associated to Coordinate3D is NoneType (pointer lost)
            >>> for subtomo in subtomos.iterSubtomos()
            >>>     print(subtomo.getCoordinate3D().getX(SCIPION))
            >>>     330 retrieved correctly

        """
        # Iterate over all Subtomograms if tomoId is None,
        # otherwise use tomoId to filter the where selection
        if volume is None:
            subtomoWhere = '1'
        elif isinstance(volume, int):
            logger.warning("FOR DEVELOPERS: Do not use volId, use volName or tsId")
            subtomoWhere = '_volId=%d' % volume
        elif isinstance(volume, Tomogram):
            subtomoWhere = '%s="%s"' % (SubTomogram.VOL_NAME_FIELD, volume.getTsId())
        else:
            raise Exception('Invalid input tomogram of type %s'
                            % type(volume))

        for subtomo in self.iterItems(where=subtomoWhere, orderBy=orderBy):
            if subtomo.hasCoordinate3D():
                subtomo.getCoordinate3D().setVolume(self.getTomogram(subtomo))
            yield subtomo

    def getTomogram(self, subtomo):
        """ returns and caches the tomogram related with a subtomogram.
        If the subtomograms were imported and not associated to any tomogram returns None."""

        # Tomogram is stored with the coordinate data
        coord = subtomo.getCoordinate3D()

        # If there is no coordinate associated
        if coord is None:
            return None

        # Else, there are coordinates
        volId = coord.getVolId()
        tsId = coord.getTomoId()

        self.initTomos()

        # If tsId is not cached, save both identifiers.
        if tsId not in self._tomos:
            tomo = self.getCoordinates3D().getPrecedents()[{Tomogram.TS_ID_FIELD: tsId}]
            self._tomos[volId] = tomo
            self._tomos[tsId] = tomo
            return tomo
        else:
            return self._tomos[tsId]

    def initTomos(self):
        """ Initialize internal _tomos to a dictionary if not done already"""
        if self._tomos is None:
            self._tomos = dict()

    def getTomograms(self):
        """ Returns a list  with only the tomograms involved in the subtomograms. May differ when
        subsets are done."""

        tomoId_attr = SubTomogram.COORD_VOL_NAME_FIELD
        if self._tomos is None:

            self.initTomos()

            uniqueTomos = self.aggregate(['count'], tomoId_attr, [tomoId_attr])

            for row in uniqueTomos:
                tsId = row[tomoId_attr]
                tomo = self.getCoordinates3D().getPrecedents()[{Tomogram.TS_ID_FIELD: tsId}]
                self._tomos[tsId] = tomo

        return self._tomos


class AverageSubTomogram(SubTomogram):
    """Represents a Average SubTomogram.
        It is a SubTomogram but it is useful to differentiate outputs."""

    def __init__(self, **kwargs):
        SubTomogram.__init__(self, **kwargs)


class SetOfAverageSubTomograms(SetOfSubTomograms):
    """Represents a set of Averages.
    It is a SetOfSubTomograms but it is useful to differentiate outputs."""
    ITEM_TYPE = AverageSubTomogram
    REP_TYPE = AverageSubTomogram
    EXPOSE_ITEMS = True

    def __init__(self, **kwargs):
        SetOfSubTomograms.__init__(self, **kwargs)


class ClassSubTomogram(SetOfSubTomograms):
    """ Represent a Class that groups SubTomogram objects.
    The representative of the class is an AverageSubTomogram.
    """
    REP_TYPE = AverageSubTomogram

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not
        _mapperPath or _size from other set of SubTomograms to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])

    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfClasses
        pass


class SetOfClassesSubTomograms(data.SetOfClasses):
    """ Store results from a subtomogram averaging method. """
    ITEM_TYPE = ClassSubTomogram
    REP_TYPE = AverageSubTomogram
    REP_SET_TYPE = SetOfAverageSubTomograms

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._coordsPointer = Pointer()

    def copyInfo(self, other):
        """ Copy properties from other set of images to current one"""
        super().copyInfo(other)
        if other._coordsPointer.hasValue():
            self.copyAttributes(other, '_coordsPointer')
        else:
            logger.warning("The source %s seems an old execution and does not have coordinates associated."
                           " This may set may fail with some protocols treating contained classes "
                           "as SetOfSubtomograms with coordinates.")

    def setCoordinates3D(self, coordinates):
        """ Set the SetOfCoordinates associated with
        this set.
         """
        if isinstance(coordinates, Pointer):
            self._coordsPointer = coordinates
        else:
            self._coordsPointer.set(coordinates)

    def _setItemMapperPath(self, item: ClassSubTomogram):
        """ This will happen when retrieving any item from this set. We take this chance to 'inject' the coordinates."""
        super()._setItemMapperPath(item)
        item.setCoordinates3D(self._coordsPointer)


class LandmarkModel(data.EMObject):
    """Represents the set of landmarks belonging to a specific tilt-series."""

    def __init__(self,
                 tsId=None,
                 fileName=None,
                 modelName=None,
                 size=5,
                 applyTSTransformation=True,
                 hasResidualInfo=False,
                 **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._tsId = String(tsId)
        self._fileName = String(fileName)
        self._modelName = String(modelName)
        self._size = Integer(size)  # Diameter in Angstroms
        self._applyTSTransformation = Boolean(applyTSTransformation)
        self._tiltSeries = Pointer(objDoStore=False)
        self._count = Integer(0)
        self._chains = None
        self._hasResidualInfo = Boolean(hasResidualInfo)

    def getTiltSeries(self):
        """ Return the tilt-series associated with this landmark model. """

        return self._tiltSeries.get()

    def setTiltSeries(self, tiltSeries):
        """ Set the tilt-series from which this landmark model were calculated.
        :param tiltSeries: Either a TiltSeries object or a pointer to it.
        """

        if tiltSeries.isPointer():
            self._tiltSeries.copy(tiltSeries)

        else:
            self._tiltSeries.set(tiltSeries)

    def getTsId(self):
        return str(self._tsId)

    def setTsId(self, tsId):
        self._tsId.set(tsId)

    def getSize(self):
        return self._size.get()

    def setSize(self, size):
        self._size.set(size)

    def getCount(self):
        return self._count.get()

    def setCount(self, count):
        self._count.set(count)

    def applyTSTransformation(self):
        return self._applyTSTransformation.get()

    def setApplyTSTransformation(self, apply):
        self._applyTSTransformation.set(apply)

    def getFileName(self):
        return self._fileName.get()

    def setFileName(self, fileName):
        self._fileName.set(fileName)

    def getModelName(self):
        return self._modelName.get()

    def setModelName(self, modelName):
        self._modelName.set(modelName)

    def hasResidualInfo(self):
        return self._hasResidualInfo

    def setHasResidualInfo(self, hasResidualInfo):
        self._hasResidualInfo.set(hasResidualInfo)

    def addLandmark(self, xCoor, yCoor, tiltIm, chainId, xResid, yResid):
        fieldNames = ['xCoor', 'yCoor', 'tiltIm', 'chainId', 'xResid', 'yResid']

        mode = "a" if os.path.exists(self.getFileName()) else "w"

        with open(self.getFileName(), mode) as f:
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldNames)
            if mode == "w":
                writer.writeheader()
            writer.writerow({'xCoor': xCoor,
                             'yCoor': yCoor,
                             'tiltIm': tiltIm,
                             'chainId': chainId,
                             'xResid': xResid,
                             'yResid': yResid})

            self._registerChain(chainId)

    def _registerChain(self, chainId):
        """ registers new chainId in a dictionary to later on store the chain count"""

        if self._chains is None:
            self._chains = dict()

        if chainId not in self._chains:
            self._chains[chainId] = None
            self.setCount(len(self._chains))

    def retrieveInfoTable(self):
        """ This method returns a table containing the information of the landkmark model. One landmark per line
        specifying in order: xCoor, YCoor, tiltIm, chainId, xResid, yResid"""

        fileName = self.getFileName()

        outputInfo = []

        with open(fileName) as f:
            reader = csv.reader(f)

            # Ignore header
            next(reader)

            for line in reader:
                vector = line[0].split()
                outputInfo.append(vector)

        return outputInfo

    def __str__(self):
        return "%s landmarks of %s Å %s to %s" \
            % (self.getCount(), self.getSize(),
               "to apply" if self.applyTSTransformation() else "applied",
               self.getTsId())


class SetOfLandmarkModels(data.EMSet):
    """Represents a class that groups a set of landmark models."""
    ITEM_TYPE = LandmarkModel

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._setOfTiltSeriesPointer = Pointer()
        self._hasResidualInfo = Boolean(False)

    def __getitem__(self, itemId):
        """Add a pointer to a tilt-series before returning the landmark model"""

        lm = super().__getitem__(itemId)

        return self.completeLandmarkModel(lm)

    def completeLandmarkModel(self, lm):
        """This method completes a landmark model object setting in execution time the tilt-series associated to it,
        since it is not possible to save pointers in the item classes of the set.

        IMPORTANT: this method must be used every time is necessary to retrieve information from the tilt-series
        associated to the landmark models that compose the set."""

        tsId = lm.getTsId()

        # Check for tilt series in set with coincident tsId
        for ts in self.getSetOfTiltSeries().iterItems(where="_tsId=='%s'" % tsId):
            lm.setTiltSeries(ts)

        return lm

    def getLandmarkModelFromTsId(self, tsId):
        """ This method return the landmark model belonging to the set that has a coincident input tsId.

        :param tsId: tilt-series ID to search the landmark model into the set."""

        for lm in self.iterItems(where="_tsId=='%s'" % tsId):
            return lm

    def getSetOfTiltSeries(self, pointer=False) -> SetOfTiltSeries:
        """ Return the set of tilt-series associated with this set of landmark models. """

        if pointer:
            return self._setOfTiltSeriesPointer

        else:
            return self._setOfTiltSeriesPointer.get()

    def setSetOfTiltSeries(self, setOfTiltSeries):
        """ Set the set of tilt-series from which this set of landmark models were calculated.
        :param tiltSeries: Either a TiltSeries object or a pointer to it.
        """

        if setOfTiltSeries.isPointer():
            self._setOfTiltSeriesPointer.copy(setOfTiltSeries, copyId=False)
        else:
            self._setOfTiltSeriesPointer.set(setOfTiltSeries)

    def hasResidualInfo(self):
        return self._hasResidualInfo

    def setHasResidualInfo(self, hasResidualInfo):
        self._hasResidualInfo.set(hasResidualInfo)


class MeshPoint(Coordinate3D):
    """Mesh object: it stores the coordinates of the points (specified by the user) needed to define
    the triangulation of a volume.
    A Mesh object can be considered as a point cloud in 3D containing the coordinates needed to divide a given region of
    space into planar triangles interconnected that will result in a closed surface."""

    def __init__(self, **kwargs):
        Coordinate3D.__init__(self, **kwargs)
        self._volumeName = String()
        self._description = None  # Algebraic description of fitted mesh

    def getVolumeName(self):
        return self._volumeName.get()

    def setVolumeName(self, volName):
        self._volumeName.set(volName)

    def getDescription(self):
        return self._description

    def setDescription(self, description):
        self._description = description

    def hasDescription(self):
        return self._description is not None


class SetOfMeshes(SetOfCoordinates3D):
    """ Store a series of meshes. """
    ITEM_TYPE = MeshPoint

    def __init__(self, **kwargs):
        SetOfCoordinates3D.__init__(self, **kwargs)
        self._numberOfMeshes = Integer()  # Indicates how many meshes are in the set

    def getNumberOfMeshes(self):
        return self._numberOfMeshes.get()

    def setNumberOfMeshes(self, n):
        self._numberOfMeshes.set(n)


class Ellipsoid(data.EMObject):
    """This class represent an ellipsoid. This is an instance class of description attribute of object MeshPoint"""

    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._center = String()
        self._radii = String()
        self._algebraicDesc = String()

    def getCenter(self):
        return self._center.get()

    def setCenter(self, center):
        self._center.set(center)

    def getRadii(self):
        return self._radii.get()

    def setRadii(self, radii):
        self._radii.set(radii)

    def getAlgebraicDesc(self):
        return self._center.get()

    def setAlgebraicDesc(self, algebraicDesc):
        self._algebraicDesc.set(algebraicDesc)

    def hasAlgebraicDesc(self):
        return self._algebraicDesc is not None


class CTFTomo(data.CTFModel):
    """ Represents a generic CTF model for a tilt-image. """
    ACQ_ORDER_FIELD = '_acqOrder'
    INDEX_FIELD = '_index'

    def __init__(self, index=None, acqOrder=None, **kwargs):
        super().__init__(**kwargs)
        self._index = Integer(index)
        self._acqOrder = Integer(acqOrder)

    def copyInfo(self, other, copyId=False):
        """ Copy info is similar to clone, but is often used to tranfer data from other type of objects with same attributes"""
        self.copy(other, copyId=copyId)

    def clone(self, copyEnable=True):
        return super().clone(copyEnable=copyEnable)

    def copy(self, other, copyId=True, ignoreAttrs=[], copyEnable=True):
        self.copyAttributes(other, '_defocusU', '_defocusV', '_defocusAngle', '_defocusRatio', '_psdFile',
                            '_resolution', '_fitQuality', self.INDEX_FIELD, self.ACQ_ORDER_FIELD)

        if copyEnable:
            self.setEnabled(other.isEnabled())

        if other.hasPhaseShift():
            self.setPhaseShift(other.getPhaseShift())

        if other.hasEstimationInfoAsList():
            if other.hasAstigmatismInfoAsList():
                self._defocusUList = CsvList(pType=float)
                self._defocusVList = CsvList(pType=float)
                self._defocusAngleList = CsvList(pType=float)
                self.setDefocusUList(other.getDefocusUList())
                self.setDefocusVList(other.getDefocusVList())
                self.setDefocusAngleList(other.getDefocusAngleList())

            else:
                self._defocusUList = CsvList(pType=float)
                self.setDefocusUList(other.getDefocusUList())

        if other.hasPhaseShiftInfoAsList():
            self._phaseShiftList = CsvList(pType=float)
            self.setPhaseShiftList(other.getPhaseShiftList())

        if other.hasCutOnFrequncyInfoAsList():
            self._cutOnFreqList = CsvList(pType=float)
            self.setCutOnFreqList(other.getCutOnFreqList())
            self.setCutOnFreq(other.getCutOnFreq())

        if copyId:
            self.copyObjId(other)

    @staticmethod
    def ctfModelToCtfTomo(ctfModel):
        newCTFTomo = CTFTomo()
        newCTFTomo.copyAttributes(ctfModel, '_defocusU', '_defocusV',
                                  '_defocusAngle', '_defocusRatio', '_psdFile',
                                  '_resolution', '_fitQuality')
        return newCTFTomo

    def getIndex(self):
        return self._index.get()

    def setIndex(self, value):
        self._index.set(value)

    def getAcquisitionOrder(self):
        return self._acqOrder.get()

    def setAcquisitionOrder(self, value):
        self._acqOrder.set(value)

    def getCutOnFreq(self):
        return self._cutOnFreq

    def setCutOnFreq(self, value):
        self._cutOnFreq = Float(value)

    " List data methods allow compatibility with IMOD metadata. "

    def getDefocusUList(self):
        return self._defocusUList.get()

    def setDefocusUList(self, defList):
        self._defocusUList.set(defList)

    def appendDefocusUList(self, value):
        self._defocusUList.append(value)

    def getDefocusVList(self):
        return self._defocusVList.get()

    def setDefocusVList(self, defList):
        self._defocusVList.set(defList)

    def appendDefocusVList(self, value):
        self._defocusVList.append(value)

    def getDefocusAngleList(self):
        return self._defocusAngleList.get()

    def setDefocusAngleList(self, defList):
        self._defocusAngleList.set(defList)

    def appendDefocusAngleList(self, value):
        self._defocusAngleList.append(value)

    def getPhaseShiftList(self):
        return self._phaseShiftList.get()

    def setPhaseShiftList(self, defList):
        self._phaseShiftList.set(defList)

    def appendPhaseShiftList(self, value):
        self._phaseShiftList.append(value)

    def getCutOnFreqList(self):
        return self._cutOnFreqList.get()

    def setCutOnFreqList(self, cutOnFreqList):
        self._cutOnFreqList.set(cutOnFreqList)

    def appendCutOnFreqList(self, value):
        self._cutOnFreqList.append(value)

    def hasEstimationInfoAsList(self):
        """ This method checks if the CTFTomo object contains estimation information in the form of a list. """

        if hasattr(self, "_defocusUList") or hasattr(self, "_defocusUList"):
            return True
        else:
            return False

    def hasAstigmatismInfoAsList(self):
        """ This method checks if the CTFTomo object contains astigmatism information in the form of a list. """

        if hasattr(self, "_defocusUList") and hasattr(self, "_defocusVList"):
            return True
        else:
            return False

    def hasPhaseShiftInfoAsList(self):
        """ This method checks if the CTFTomo object contains phase shift information in the form of a list. """

        if hasattr(self, "_phaseShiftList"):
            return True
        else:
            return False

    def hasCutOnFrequncyInfoAsList(self):
        """ This method checks if the CTFTomo object contains cut-on frequency information in the form of a list. """

        if hasattr(self, "_cutOnFreqList"):
            return True
        else:
            return False

    def completeInfoFromList(self):
        """ This method will set the _defocusU, _defocusV and _defocusAngle attributes from the provided CTF estimation
        information lists.

        Based on the IMOD program ctfphaseflip: "The program  will assign that defocus value to the midpoint of the
        range of views.  For a view at a given tilt angle, it will find the defocus either by interpolating between
        two surrounding midpoint angles, if there are such angles, or by taking the nearest defocus value, if the
        angle is beyond the range of the available midpoint angles. "
        - From IMOD documentation https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html

        This method will assign as the defocus value and angle the median of the estimation list. """

        " DEFOCUS INFORMATION -----------------------------------------------------------------------------------------"

        " Check that at least one list is provided "
        if not self.hasEstimationInfoAsList():
            raise Exception("CTFTomo object has no _defocusUList neither _defocusUList argument initialized. No "
                            "list information available.")

        " Get the number of provided list (1 or 2) "
        numberOfProvidedList = 2 if (hasattr(self, "_defocusUList") and hasattr(self, "_defocusVList")) else 1

        " No astigmatism is estimated (only one list provided) "
        if numberOfProvidedList == 1:
            providedDefocusUList = self.getDefocusUList() if hasattr(self, "_defocusUList") else self.getDefocusVList()
            providedDefocusUList = providedDefocusUList.split(",")

            " DefocusAngle is set to 0 degrees "
            self.setDefocusAngle(0)

            " DefocusU and DefocusV are set at the same value, equal to the middle estimation of the list "
            middlePoint = math.trunc(len(providedDefocusUList) / 2)

            " If the size of the defocus list is even, mean the 2 centre values "
            if len(providedDefocusUList) % 2 == 0:
                value = (float(providedDefocusUList[middlePoint]) + float(providedDefocusUList[middlePoint - 1])) / 2

                self.setDefocusU(value)
                self.setDefocusV(value)

            else:
                " If the size of defocus estimation is odd, get the centre value "

                value = providedDefocusUList[middlePoint]

                self.setDefocusU(value)
                self.setDefocusV(value)

        else:
            " Astigmatism is estimated (two lists are provided) "

            providedDefocusUList = self.getDefocusUList()
            providedDefocusUList = providedDefocusUList.split(",")

            providedDefocusVList = self.getDefocusVList()
            providedDefocusVList = providedDefocusVList.split(",")

            providedDefocusAngleList = self.getDefocusAngleList()
            providedDefocusAngleList = providedDefocusAngleList.split(",")

            " Check that the three list are equally long "
            if len(providedDefocusUList) != len(providedDefocusVList) or \
                    len(providedDefocusUList) != len(providedDefocusAngleList) or \
                    len(providedDefocusVList) != len(providedDefocusAngleList):
                raise Exception("DefocusUList, DefocusVList and DefocusAngleList lengths must be equal.")

            " DefocusU, DefocusV and DefocusAngle are set equal to the middle estimation of the list "
            middlePoint = math.trunc(len(providedDefocusUList) / 2)

            " If the size of the defocus list is even, mean the 2 centre values "
            if len(providedDefocusUList) % 2 == 0:
                defocusU = (float(providedDefocusUList[middlePoint]) +
                            float(providedDefocusUList[middlePoint - 1])) / 2
                defocusV = (float(providedDefocusVList[middlePoint]) +
                            float(providedDefocusVList[middlePoint - 1])) / 2
                defocusAngle = (float(providedDefocusAngleList[middlePoint]) +
                                float(providedDefocusAngleList[middlePoint - 1])) / 2

                self.setDefocusU(defocusU)
                self.setDefocusV(defocusV)
                self.setDefocusAngle(defocusAngle)

            else:
                " If the size of defocus estimation list is odd, get the centre value "

                defocusU = providedDefocusUList[middlePoint]
                defocusV = providedDefocusVList[middlePoint]
                defocusAngle = providedDefocusAngleList[middlePoint]

                self.setDefocusU(defocusU)
                self.setDefocusV(defocusV)
                self.setDefocusAngle(defocusAngle)

        " PHASE SHIFT INFORMATION -------------------------------------------------------------------------------------"

        " Check if phase shift information is also available "
        if hasattr(self, "_phaseShiftList"):
            providedPhaseShiftList = self.getPhaseShiftList()
            providedPhaseShiftList = providedPhaseShiftList.split(",")

            " Check that all the lists are equally long "
            if len(providedDefocusUList) != len(providedPhaseShiftList):
                raise Exception("PhaseShiftList length must be equal to DefocusUList, DefocusVList and "
                                "DefocusAngleList lengths.")

            " PhaseShift is set equal to the middle estimation of the list "
            middlePoint = math.trunc(len(providedPhaseShiftList) / 2)

            " If the size of the phase shift list is even, mean the 2 centre values "
            if len(providedPhaseShiftList) % 2 == 0:
                phaseShift = (float(providedPhaseShiftList[middlePoint]) +
                              float(providedPhaseShiftList[middlePoint - 1])) / 2

                self.setPhaseShift(phaseShift)

            else:
                " If the size of phase shift list estimation is odd, get the centre value "

                phaseShift = providedPhaseShiftList[middlePoint]

                self.setPhaseShift(phaseShift)

        " CUT-ON FREQUENCY INFORMATION --------------------------------------------------------------------------------"

        " Check if cut-on frequency information is also available "
        if hasattr(self, "_cutOnFreqList"):
            providedCutOnFreqList = self.getCutOnFreqList()
            providedCutOnFreqList = providedCutOnFreqList.split(",")

            " Check that all the lists are equally long "
            if len(providedPhaseShiftList) != len(providedCutOnFreqList):
                raise Exception("CutOnFreqList length must be equal to PhaseShiftList, DefocusUList, DefocusVList and "
                                "DefocusAngleList lengths.")

            " Cut-on frequency is set equal to the middle estimation of the list "
            middlePoint = math.trunc(len(providedCutOnFreqList) / 2)

            " If the size of the cut-on frequency shift list is even, mean the 2 centre values "
            if len(providedCutOnFreqList) % 2 == 0:
                cutOnFreq = (float(providedCutOnFreqList[middlePoint]) +
                             float(providedCutOnFreqList[middlePoint - 1])) / 2

                self.setCutOnFreq(cutOnFreq)

            else:
                " If the size of the cut-on frequency list estimation is odd, get the centre value "

                cutOnFreq = providedCutOnFreqList[middlePoint]

                self.setCutOnFreq(cutOnFreq)

        " Standardize the input values "
        self.standardize()


class CTFTomoSeries(data.EMSet):
    """ Represents a set of CTF models belonging to the same tilt-series. """
    ITEM_TYPE = CTFTomo
    TS_ID_FIELD = '_tsId'

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
        self._tiltSeriesPointer = Pointer(kwargs.get('tiltSeriesPointer', None))
        self._tsId = String(kwargs.get('tsId', None))
        self._isDefocusUDeviationInRange = Boolean(True)
        self._isDefocusVDeviationInRange = Boolean(True)

        # CtfModels will always be used inside a SetOfTiltSeries
        # so, let's do not store the mapper path by default
        self._mapperPath.setStore(False)

    def clone(self, ignoreAttrs=('_mapperPath', '_size')):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=ignoreAttrs)
        clone.setEnabled(self.isEnabled())
        return clone

    def __del__(self):
        # Cancel closing the mapper since this class is an item of a set and shares the mapper with its parent set.
        pass

    def getTiltSeries(self):
        """ Return the tilt-series associated with this CTF model series. """
        return self._tiltSeriesPointer.get()

    def setTiltSeries(self, tiltSeries):
        """ Set the tilt-series from which this CTF model series were estimated.
        :param tiltSeries: Either a TiltSeries object or a pointer to it.
        """
        if tiltSeries.isPointer():
            self._tiltSeriesPointer.copy(tiltSeries)
        else:
            self._tiltSeriesPointer.set(tiltSeries)

    def getTsId(self):
        """ Get unique TiltSeries ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def getNumberOfEstimationsInRange(self):
        """ Return the tilt-images range size used for estimation. """
        return self._estimationsInRange.get()

    def setNumberOfEstimationsInRange(self, estimationRange):
        """ Set the tilt-images range size used for estimation.
        :param estimationRange: Integer of the range size. """

        self._estimationsInRange = Integer(estimationRange)

    def getIMODDefocusFileFlag(self):
        """ Return the format file from which the CTF estimation information has been acquired. This parameter is
        useful for posterior information and format conversions between IMOD and Scipion. The flag value "is the sum of:

          1 if the file has astigmatism values
          2 if the astigmatism axis angle is in radians, not degrees
          4 if the file has phase shifts
          8 if the phase shifts are in radians, not degrees
         16 if tilt angles need to be inverted to match what the
             program expects (what Ctfplotter would produce)
             with the -invert option
         32 if the file has cut-on frequencies attenuating the phase
             at low frequencies"

             from https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html """

        return self._IMODDefocusFileFlag.get()

    def setIMODDefocusFileFlag(self, flag):
        """ Set the format file from which the CTF estimation information has been acquired.
        :param flag: Integer of the range size.

        This parameter is
        useful for posterior information and format conversions between IMOD and Scipion. The flag value "is the sum of:

          1 if the file has astigmatism values
          2 if the astigmatism axis angle is in radians, not degrees
          4 if the file has phase shifts
          8 if the phase shifts are in radians, not degrees
         16 if tilt angles need to be inverted to match what the
             program expects (what Ctfplotter would produce)
             with the -invert option
         32 if the file has cut-on frequencies attenuating the phase
             at low frequencies"

             from https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html """

        self._IMODDefocusFileFlag = Integer(flag)

    def setNumberOfEstimationsInRangeFromDefocusList(self):
        """ Set the tilt-images estimation range size used for estimation from the defocus info list size. """

        estimationRange = 0

        for ctfEstimation in self:
            # Check that at least one list is provided
            if not (hasattr(ctfEstimation, "_defocusUList") or hasattr(
                    ctfEstimation, "_defocusUList")):
                raise Exception(
                    "CTFTomo object has no _defocusUList neither _defocusUList argument initialized. No "
                    "list information available.")

            providedList = ctfEstimation.getDefocusUList() if hasattr(
                ctfEstimation, "_defocusUList") \
                else ctfEstimation.getDefocusVList()
            providedList = providedList.split(",")

            listLength = len(providedList) - 1

            if listLength > estimationRange:
                estimationRange = listLength

        self.setNumberOfEstimationsInRange(estimationRange)

    def getIsDefocusUDeviationInRange(self):
        return True

    def setIsDefocusUDeviationInRange(self, value):
        pass

    def getIsDefocusVDeviationInRange(self):
        return True

    def setIsDefocusVDeviationInRange(self, value):
        pass

    def calculateDefocusUDeviation(self, defocusUTolerance=20):
        pass

    def calculateDefocusVDeviation(self, defocusVTolerance=20):
        pass

    def getCtfTomoFromTi(self, ti: TiltImage, onlyEnabled: bool = True) -> Optional[CTFTomo]:
        """Get the corresponding CTFModel from a given tilt-image. If there's no match, it returns None.
        :param ti: Tilt-image.
        :param onlyEnabled: boolean used to indicate the matching behaves in terms of the value the attribute
        _objEnabled from both ti and the matching CTFModel attribute:
            - If True (default), the CTFModel found is returned only if both the ti and the CTFModel are enabled.
            - If False, the CTFModel found is returned no matter the value of _objEnabled.
        """
        if onlyEnabled and not ti.isEnabled():
            logger.debug('The introduced tilt-image is not enabled and working with onlyEnabled = True')
            return None
        try:
            # The method getItem raises an exception of type:
            #   - sqlite3.OperationalError if the key is not found
            #   - UnboundLocalError if the value is not found
            ctfTomo = self.getItem(CTFTomo.ACQ_ORDER_FIELD, ti.getAcquisitionOrder())
        except UnboundLocalError:
            return None
        except OperationalError:
            try:
                ctfTomo = self.getItem(CTFTomo.INDEX_FIELD, ti.getIndex())
                logger.warning('WARNING! The current CTF series does not have the attribute "acquisition order" '
                               '(_acqOrder). The matching between the CTF and the tilt-image is carried out '
                               'using the index --> LESS RELIABLE. CHECK THE RESULTS CAREFULLY')
            except (OperationalError, UnboundLocalError):
                logger.warning(f'No CTF found in the current CTF series {self.getTsId()} that matches the '
                               f'given tilt-image of tsId = {ti.getTsId()}.')
                return None

        if ctfTomo.isEnabled() or not onlyEnabled:
            return ctfTomo
        else:
            return None


class SetOfCTFTomoSeries(data.EMSet):
    """ Represents a set of CTF model series belonging to the same set of tilt-series. """
    ITEM_TYPE = CTFTomoSeries
    USE_CREATE_COPY_FOR_SUBSET = True

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
        self._setOfTiltSeriesPointer = Pointer(kwargs.get('tiltSeriesPointer', None))
        self._idDict = {}

    def copyInfo(self, other):
        data.EMSet.copyInfo(self, other)
        self.setSetOfTiltSeries(other.getSetOfTiltSeries(pointer=True))

    def copyItems(self, other, itemSelectedCallback=None):
        """ Copy items (CTFTomoSeries and CTFTomo) from the other Set.
         Params:
            other:  SetOfCTFTomoSeries from where to copy elements.
            
            itemSelectedCallback: Optional, callback receiving an item and
                returning true if it has to be copied
        """
        for i, ctfSerie in enumerate(other.iterItems()):
            if itemSelectedCallback(ctfSerie):
                ctfSerieOut = ctfSerie.clone()
                self.append(ctfSerieOut)

                for j, ctf in enumerate(ctfSerie.iterItems()):
                    ctfOut = ctf.clone()
                    ctfSerieOut.append(ctfOut)

                self.update(ctfSerieOut)

    def getSetOfTiltSeries(self, pointer=False):
        """ Return the tilt-series associated with this CTF model series. """
        return self._setOfTiltSeriesPointer.get() if not pointer else self._setOfTiltSeriesPointer

    def setSetOfTiltSeries(self, setOfTiltSeries):
        """ Set the tilt-series from which this CTF model series were estimated.
        :param setOfTiltSeries: Either a TiltSeries object or a pointer to it.
        """
        if setOfTiltSeries.isPointer():
            self._setOfTiltSeriesPointer.copy(setOfTiltSeries)
        else:
            self._setOfTiltSeriesPointer.set(setOfTiltSeries)

    def iterClassItems(self, iterDisabled=False):
        """ Iterate over the images of a class.
        Params:
            iterDisabled: If True, also include the disabled items. """
        for cls in self.iterItems():
            if iterDisabled or cls.isEnabled():
                for img in cls:
                    if iterDisabled or img.isEnabled():
                        yield img

    def _setItemMapperPath(self, item):
        """ Set the mapper path of this class according to the mapper
        path of the SetOfClasses and also the prefix according to class id
        """
        item._mapperPath.set('%s,id%s' % (self.getFileName(), item.getObjId()))
        item.load()

    def _insertItem(self, item):
        """ Create the SetOfImages assigned to a class.
        If the file exists, it will load the Set.
        """
        self._setItemMapperPath(item)
        data.EMSet._insertItem(self, item)
        item.write(properties=False)  # Set.write(self)

    def __getitem__(self, itemId):
        """ Setup the mapper classes before returning the item. """
        classItem = data.EMSet.__getitem__(self, itemId)

        objId = None
        for tiltSeries in self.getSetOfTiltSeries().iterItems(iterate=False):
            if tiltSeries.getTsId() == classItem.getTsId():
                objId = tiltSeries.getObjId()

        if objId is None:
            raise ("Could not find tilt-series with tsId = %s" % classItem.getTsId())

        classItem.setTiltSeries(self.getSetOfTiltSeries()[objId])

        self._setItemMapperPath(classItem)
        return classItem

    def getFirstItem(self):
        classItem = data.EMSet.getFirstItem(self)
        self._setItemMapperPath(classItem)
        return classItem

    def iterItems(self, orderBy='id', direction='ASC', **kwargs):
        for item in super().iterItems(orderBy=orderBy,
                                      direction=direction, **kwargs):

            ts = self._getTiltSeriesFromTsId(item.getTsId())
            if ts is None:
                raise Exception("Could not find tilt-series with tsId = %s" % item.getTsId())

            item.setTiltSeries(ts)
            self._setItemMapperPath(item)

            yield item

    def _getTiltSeriesFromTsId(self, tsId):
        if self._idDict:
            return self._idDict.get(tsId, None)
        else:
            self._idDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.getSetOfTiltSeries()}
            return self._idDict.get(tsId, None)

    def getTSIds(self):
        """ Returns al the Tilt series ids involved in the set."""
        return self.getUniqueValues(CTFTomoSeries.TS_ID_FIELD)


class TiltSeriesCoordinate(data.EMObject):
    """This class holds the (x,y,z) positions, in angstroms, and other information
    associated with a coordinate related to a tilt series"""

    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._x = Float()
        self._y = Float()
        self._z = Float()
        # Used to access to the corresponding tilt series from each coord (it's the tsId)
        self._tsId = String(kwargs.get('tsId', None))

    def copyInfo(self, coord):
        """ Copy information from other coordinate. """
        self.setPosition(*coord.getPosition())
        self.setTsId(coord.getTsId())
        self.setObjId(coord.getObjId())

    def getX(self):
        """ Returns the X dimension (Å) of the coordinate"""
        return self._x.get()

    def setX(self, x):
        """ Sets the x dimension (Å) of the coordinate """
        self._x.set(x)

    def getY(self):
        """ Returns the Y dimension (Å) of the coordinate"""
        return self._y.get()

    def setY(self, y):
        """ Sets the Y dimension of (Å) the coordinate"""
        self._y.set(y)

    def getZ(self):
        """ Returns the Z dimension (Å) of the coordinate"""
        return self._z.get()

    def setZ(self, z):
        """ Sets the Z dimension (Å) of the coordinate"""
        self._z.set(z)

    def getPosition(self, sampling_rate=1):
        """Returns the position a TiltSeriesCoordinate in a tuple at a specific sampling rate (optional)"""
        return self.getX() / sampling_rate, self.getY() / sampling_rate, self.getZ() / sampling_rate

    def setPosition(self, x, y, z, sampling_rate):
        """Set the position of the coordinate
            :param int x: Position of the coordinate in the X axis
            :param int y: Position of the coordinate in the Y axis
            :param int z: Position of the coordinate in the Z axis
            :param flat sampling_rate: sampling rate in which x,y,z are measured. Default 1 = Å
        """
        self.setX(x * sampling_rate)
        self.setY(y * sampling_rate)
        self.setZ(z * sampling_rate)

    def getTsId(self):
        return self._tsId.get()

    def setTsId(self, tsId):
        self._tsId.set(tsId)


class SetOfTiltSeriesCoordinates(data.EMSet):
    """ Encapsulate the logic of a set of tilt series coordinates.
    Each coordinate has a (x,y,z) position in scipion's convention.
    Scipion's convention is the center of a theoretical volume when applying
    the tilt series transformation matrix:
    X0 --> half x dimension
    Y0 --> half y dimension
    Z0 --> half z of a theoretical volume?
    """
    ITEM_TYPE = TiltSeriesCoordinate

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
        self._SetOfTiltSeriesPointer = Pointer()

    def getSetOfTiltSeries(self):
        """ Returns the Tilt Series associated with
                this SetOfTiltSeriesCoordinates"""
        return self._SetOfTiltSeriesPointer.get()

    def setSetOfTiltSeries(self, setOfTiltSeries):
        """ Set the Tilt Series associated with this set of coordinates.

            Params:

            setOfTiltSeries: Tilt Series object or a pointer to it.
                """
        if setOfTiltSeries.isPointer():
            self._SetOfTiltSeriesPointer.copy(setOfTiltSeries)
        else:
            self._SetOfTiltSeriesPointer.set(setOfTiltSeries)

    def getSummary(self):
        summary = []
        summary.append("Number of tilt series coordinates: %s" % self.getSize())
        return "\n".join(summary)

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not _mapperPath or _size
        from other set of objects to current one.
        """
        self.setSetOfTiltSeries(other.getSetOfTiltSeries())
