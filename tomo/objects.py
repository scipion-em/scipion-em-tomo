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

import os
from datetime import datetime
import threading
from collections import OrderedDict
import numpy as np
import math
import csv

import pyworkflow.object as pwobj
import pyworkflow.utils.path as path
import pwem.objects.data as data
from pwem.convert.transformations import euler_matrix
import pwem.emlib as emlib
from pwem.emlib.image import ImageHandler


class TiltImageBase:
    """ Base class for TiltImageM and TiltImage. """
    def __init__(self, **kwargs):
        self._tiltAngle = pwobj.Float(kwargs.get('tiltAngle', None))
        self._tsId = pwobj.String(kwargs.get('tsId', None))

        # Use the acquisition order as objId
        if 'acquisitionOrder' in kwargs:
            self.setObjId(int(kwargs['acquisitionOrder']))

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
        return self.getObjId()

    def copyInfo(self, other, copyId=False):
        self.copyAttributes(other, '_tiltAngle', '_tsId')
        if copyId:
            self.copyObjId(other)


class TiltImage(data.Image, TiltImageBase):
    """ Tilt image """
    def __init__(self, location=None, **kwargs):
        data.Image.__init__(self, location, **kwargs)
        TiltImageBase.__init__(self, **kwargs)

    def copyInfo(self, other, copyId=False):
        data.Image.copyInfo(self, other)
        TiltImageBase.copyInfo(self, other, copyId=copyId)


class TiltSeriesBase(data.SetOfImages):
    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)
        self._tsId = pwobj.String(kwargs.get('tsId', None))
        # TiltSeries will always be used inside a SetOfTiltSeries
        # so, let's do no store the mapper path by default
        self._mapperPath.setStore(False)

    def getTsId(self):
        """ Get unique TiltSerie ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def copyInfo(self, other, copyId=False):
        """ Copy basic information (id and other properties) but
        not _mapperPath or _size from other set of tomograms to current one.
        """
        self.copy(other, copyId=copyId, ignoreAttrs=['_mapperPath', '_size'])

    def append(self, tiltImage):
        tiltImage.setTsId(self.getTsId())
        data.SetOfImages.append(self, tiltImage)

    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfTiltSeries
        pass

    def getScannedPixelSize(self):
        mag = self._acquisition.getMagnification()
        return self._samplingRate.get() * 1e-4 * mag

    def generateTltFile(self, tltFilePath, reverse=False):
        """Generates an angle file in .tlt format in the specified location. If reverse is set to true the angles in
        file are sorted in the opposite order"""
        angleList = []
        for ti in self:
            angleList.append(ti.getTiltAngle())
        if reverse:
            angleList.reverse()
        with open(tltFilePath, 'w') as f:
            f.writelines("%s\n" % angle for angle in angleList)


class TiltSeries(TiltSeriesBase):
    ITEM_TYPE = TiltImage

    def applyTransform(self, outputFilePath):
        ih = ImageHandler()
        inputFilePath = self.getFirstItem().getFileName()
        newStack = True
        # TODO: Handle output tilt-series datatype format
        if self.getFirstItem().hasTransform():
            for index, ti in enumerate(self):
                if ti.hasTransform():
                    if newStack:
                        ih.createEmptyImage(fnOut=outputFilePath,
                                            xDim=ti.getXDim(),
                                            yDim=ti.getYDim(),
                                            nDim=self.getSize())
                        newStack = False
                    transform = ti.getTransform().getMatrix()
                    transformArray = np.array(transform)
                    ih.applyTransform(inputFile=str(index + 1) + ':mrcs@' + inputFilePath,
                                      outputFile=str(index + 1) + '@' + outputFilePath,
                                      transformMatrix=transformArray,
                                      shape=(ti.getXDim(), ti.getYDim()))
                else:
                    raise Exception('ERROR: Some tilt-image is missing from transform object associated.')
        else:
            path.createLink(inputFilePath, outputFilePath)


class SetOfTiltSeriesBase(data.SetOfImages):
    EXPOSE_ITEMS = True

    """ Base class for SetOfTiltImages and SetOfTiltImagesM.
    """
    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)

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
        """ Setup the mapper classes before returning the item. """
        classItem = data.SetOfImages.__getitem__(self, itemId)
        self._setItemMapperPath(classItem)
        return classItem

    def getFirstItem(self):
        classItem = data.EMSet.getFirstItem(self)
        self._setItemMapperPath(classItem)
        return classItem

    def iterItems(self, orderBy='id', direction='ASC'):
        for item in data.EMSet.iterItems(self, orderBy=orderBy,
                                         direction=direction):
            self._setItemMapperPath(item)
            yield item

    def copyItems(self, inputTs,
                  orderByTs='id', updateTsCallback=None,
                  orderByTi='id', updateTiCallback=None):
        """ Copy items (TiltSeries and TiltImages) from the input Set.
         Params:
            inputTs: input TiltSeries (or movies) from where to copy elements.
            orderByTs: optional orderBy value for iterating over TiltSeries
            updateTsCallback: optional callback after TiltSeries is created
            orderByTi: optional orderBy value for iterating over TiltImages
            updateTiCallback: optional callback after TiltImage is created
        """
        for i, ts in enumerate(inputTs.iterItems(orderBy=orderByTs)):
            tsOut = self.ITEM_TYPE()
            tsOut.copyInfo(ts)
            tsOut.copyObjId(ts)
            if updateTsCallback:
                updateTsCallback(i, ts, tsOut)
            self.append(tsOut)
            for j, ti in enumerate(ts.iterItems(orderBy=orderByTi)):
                tiOut = tsOut.ITEM_TYPE()
                tiOut.copyInfo(ti)
                tiOut.copyObjId(ti)
                tiOut.setLocation(ti.getLocation())
                if updateTiCallback:
                    updateTiCallback(j, ts, ti, tsOut, tiOut)
                tsOut.append(tiOut)

            self.update(tsOut)

    def updateDim(self):
        """ Update dimensions of this set base on the first element. """
        self.setDim(self.getFirstItem().getDim())

    def getScannedPixelSize(self):
        mag = self._acquisition.getMagnification()
        return self._samplingRate.get() * 1e-4 * mag


class SetOfTiltSeries(SetOfTiltSeriesBase):
    ITEM_TYPE = TiltSeries


class TiltImageM(data.Movie, TiltImageBase):
    """ Tilt movie. """
    def __init__(self, location=None, **kwargs):
        data.Movie.__init__(self, location, **kwargs)
        TiltImageBase.__init__(self, **kwargs)

    def copyInfo(self, other, copyId=False):
        data.Movie.copyInfo(self, other)
        TiltImageBase.copyInfo(self, other, copyId=copyId)


class TiltSeriesM(TiltSeriesBase):
    ITEM_TYPE = TiltImageM


class SetOfTiltSeriesM(SetOfTiltSeriesBase):
    ITEM_TYPE = TiltSeriesM

    def __init__(self, **kwargs):
        SetOfTiltSeriesBase.__init__(self, **kwargs)
        self._gainFile = pwobj.String()
        self._darkFile = pwobj.String()
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
        #self._firstFramesRange.set(other.getFramesRange())


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
        # print(">>> DEBUG: _checkNewInput ")

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
        # print(">>> DEBUG: _checkNewInput ")
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
    def __init__(self, **kwargs):
        data.Acquisition.__init__(self, **kwargs)
        self._angleMin = pwobj.Float(kwargs.get('angleMin', None))
        self._angleMax = pwobj.Float(kwargs.get('angleMax', None))
        self._step = pwobj.Integer(kwargs.get('step', None))
        self._angleAxis1 = pwobj.Float(kwargs.get('angleAxis1', None))
        self._angleAxis2 = pwobj.Float(kwargs.get('angleAxis2', None))

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

    def getAngleAxis1(self):
        return self._angleAxis1.get()

    def setAngleAxis1(self, value):
        self._angleAxis1.set(value)

    def getAngleAxis2(self):
        return self._angleAxis2.get()

    def setAngleAxis2(self, value):
        self._angleAxis2.set(value)


class Tomogram(data.Volume):
    def __init__(self, **kwargs):
        data.Volume.__init__(self, **kwargs)
        self._acquisition = None
        self._tsId = pwobj.String(kwargs.get('tsId', None))

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


class SetOfTomograms(data.SetOfVolumes):
    ITEM_TYPE = Tomogram
    EXPOSE_ITEMS = True

    def __init__(self, *args, **kwargs):
        data.SetOfVolumes.__init__(self, **kwargs)
        self._acquisition = TomoAcquisition()

    def updateDim(self):
        """ Update dimensions of this set base on the first element. """
        self.setDim(self.getFirstItem().getDim())


class Coordinate3D(data.EMObject):
    """This class holds the (x,y) position and other information
    associated with a coordinate"""
    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._volumePointer = pwobj.Pointer(objDoStore=False)
        self._x = pwobj.Integer(kwargs.get('x', None))
        self._y = pwobj.Integer(kwargs.get('y', None))
        self._z = pwobj.Integer(kwargs.get('z', None))
        self._volId = pwobj.Integer()
        self._volName = pwobj.String()
        self._eulerMatrix = data.Transform()

    def getX(self):
        return self._x.get()

    def setX(self, x):
        self._x.set(x)

    def shiftX(self, shiftX):
        self._x.sum(shiftX)

    def getY(self):
        return self._y.get()

    def setY(self, y):
        self._y.set(y)

    def shiftY(self, shiftY):
        self._y.sum(shiftY)

    def getZ(self):
        return self._z.get()

    def setZ(self, z):
        self._z.set(z)

    def setMatrix(self, matrix):
        self._eulerMatrix.setMatrix(matrix)

    def getMatrix(self):
        return self._eulerMatrix.getMatrix()

    def euler2Matrix(self, r, p, y):
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

    def getPosition(self):
        """ Return the position of the coordinate as a (x, y) tuple.
        mode: select if the position is the center of the box
        or in the top left corner.
        """
        return self.getX(), self.getY(), self.getZ()

    def setPosition(self, x, y, z):
        self.setX(x)
        self.setY(y)
        self.setZ(z)

    def getVolume(self):
        """ Return the micrograph object to which
        this coordinate is associated.
        """
        return self._volumePointer.get()

    def setVolume(self, volume):
        """ Set the micrograph to which this coordinate belongs. """
        self._volumePointer.set(volume)
        self._volId.set(volume.getObjId())
        self._volName.set(volume.getFileName())

    def copyInfo(self, coord):
        """ Copy information from other coordinate. """
        self.setPosition(*coord.getPosition())
        self.setObjId(coord.getObjId())
        self.setBoxSize(coord.getBoxSize())

    def getVolId(self):
        return self._volId.get()

    def setVolId(self, volId):
        self._volId.set(volId)

    def invertY(self):
        if not self.getVolume() is None:
            dims = self.getVolume().getDim()
            height = dims[1]
            self.setY(height - self.getY())
        # else: error TODO

    def setVolName(self, volName):
        self._volName.set(volName)

    def getVolName(self):
        return self._volName.get()


class SetOfCoordinates3D(data.EMSet):
    """ Encapsulate the logic of a set of volumes coordinates.
    Each coordinate has a (x,y,z) position and is related to a Volume
    The SetOfCoordinates3D can also have information about TiltPairs.
    """
    ITEM_TYPE = Coordinate3D

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
        self._boxSize = pwobj.Integer()
        self._samplingRate = pwobj.Float()
        self._precedentsPointer = pwobj.Pointer()

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

    def iterCoordinates(self, volume=None):
        """ Iterate over the coordinates associated with a tomogram.
        If tomogram=None, the iteration is performed over the whole
        set of coordinates.
        """
        if volume is None:
            volId = None
        elif isinstance(volume, int):
            volId = volume
        elif isinstance(volume, data.Volume):
            volId = volume.getObjId()
        else:
            raise Exception('Invalid input micrograph of type %s'
                            % type(volume))

        # Iterate over all coordinates if micId is None,
        # otherwise use micId to filter the where selection
        coordWhere = '1' if volId is None else '_volId=%d' % int(volId)

        for coord in self.iterItems(where=coordWhere):
            yield coord

    def getPrecedents(self):
        """ Returns the SetOfTomograms or Tilt Series associated with
                this SetOfCoordinates"""
        return self._precedentsPointer.get()

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

    def __str__(self):
        """ String representation of a set of coordinates. """
        if self._boxSize.hasValue():
            boxSize = self._boxSize.get()
            boxStr = ' %d x %d x %d' % (boxSize, boxSize, boxSize)
        else:
            boxStr = 'No-Box'
        s = "%s (%d items, %s%s)" % (self.getClassName(), self.getSize(),
                                     boxStr, self._appendStreamState())

        return s


class SubTomogram(data.Volume):
    def __init__(self, **kwargs):
        data.Volume.__init__(self, **kwargs)
        self._acquisition = None
        self._coordinate = None

    def hasCoordinate3D(self):
        return self._coordinate is not None

    def setCoordinate3D(self, coordinate):
        self._coordinate = coordinate

    def getCoordinate3D(self):
        return self._coordinate

    def getAcquisition(self):
        return self._acquisition

    def setAcquisition(self, acquisition):
        self._acquisition = acquisition

    def hasAcquisition(self):
        return self._acquisition is not None and \
               self._acquisition.getAngleMin() is not None and \
               self._acquisition.getAngleMax() is not None


class SetOfSubTomograms(data.SetOfVolumes):
    ITEM_TYPE = SubTomogram
    REP_TYPE = SubTomogram

    def __init__(self, **kwargs):
        data.SetOfVolumes.__init__(self, **kwargs)
        self._acquisition = TomoAcquisition()
        self._coordsPointer = pwobj.Pointer()

    def hasCoordinates3D(self):
        return self._coordsPointer.hasValue()

    def getCoordinates3D(self):
        """ Returns the SetOfCoordinates associated with
        this SetOfParticles"""
        return self._coordsPointer.get()

    def setCoordinates3D(self, coordinates):
        """ Set the SetOfCoordinates associates with
        this set of particles.
         """
        self._coordsPointer.set(coordinates)


class AverageSubTomogram(SubTomogram):
    """Represents a set of Averages.
    It is a SetOfParticles but it is useful to differentiate outputs."""
    def __init__(self, **kwargs):
        SubTomogram.__init__(self, **kwargs)


class SetOfAverageSubTomograms(SetOfSubTomograms):
    ITEM_TYPE = AverageSubTomogram
    REP_TYPE = AverageSubTomogram

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


class LandmarkModel(data.EMObject):
    """Represents the set of landmarks belonging to an specific Tilt-series."""
    def __init__(self, tsId=None, fileName=None, modelName=None, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._tsId = pwobj.String(tsId)
        self._fileName = pwobj.String(fileName)
        self._modelName = pwobj.String(modelName)

    def getTsId(self):
        return str(self._tsId)

    def getFileName(self):
        return str(self._fileName)

    def getModelName(self):
        return str(self._modelName)

    def setTsId(self, tsId):
        self._tsId = pwobj.String(tsId)

    def setFileName(self, fileName):
        self._fileName = pwobj.String(fileName)

    def setModelName(self, modelName):
        self._modelName = pwobj.String(modelName)

    def addLandmark(self, xCoor, yCoor,  tiltIm, chainId, xResid, yResid):
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


class SetOfLandmarkModels(data.EMSet):
    """Represents a class that groups a set of landmark models."""
    ITEM_TYPE = LandmarkModel

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)


class Mesh(data.EMObject):
    """Mesh object: it stores the coordinates of the points (specified by the user) needed to define
    the triangulation of a volume.
    A Mesh object can be consider as a point cloud in 3D containing the coordinates needed to divide a given region of
    space into planar triangles interconnected that will result in a closed surface."""
    def __init__(self, path=None, group=None, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self._path = pwobj.String(path)
        self._volume = None
        self._group = pwobj.Integer(group)

    def getPath(self):
        return self._path.get()

    def setPath(self, path):
        self._path.set(path)

    def getGroup(self):
        return self._group.get()

    def setGroup(self, group):
        self._group.set(group)

    def getMesh(self):
        filePath = self.getPath()
        mesh = []
        array = np.genfromtxt(filePath, delimiter=',')
        numMeshes = np.unique(array[:,3])
        for idm in numMeshes:
            mesh_group = array[array[:,3] == idm, :]
            mesh.append(mesh_group[:, 0:3])
        return mesh[self.getGroup() - 1]

    def getVolume(self):
        """ Return the tomogram object to which
        this mesh is associated.
        """
        return self._volume

    def setVolume(self, volume):
        """ Set the tomogram to which this mesh belongs. """
        self._volume = volume

    def __str__(self):
        return "Mesh (path=%s)" % self.getPath()


class SetOfMeshes(data.EMSet):
    """ Store a series of meshes. """
    ITEM_TYPE = Mesh

    def __init__(self, *args, **kwargs):
        data.EMSet.__init__(self, *args, **kwargs)
        self._volumesPointer = pwobj.Pointer()

    def setVolumes(self, volumes):
        """ Set the tomograms associated with this set of coordinates.
        Params:
            tomograms: Either a SetOfTomograms object or a pointer to it.
        """
        if volumes.isPointer():
            self._volumesPointer.copy(volumes)
        else:
            self._volumesPointer.set(volumes)

    def getVolumes(self):
        """ Returns the SetOfTomograms associated with
        this SetOfCoordinates"""
        return self._volumesPointer.get()

    def iterItems(self, orderBy='id', direction='ASC', where='1', limit=None):
        """ Redefine iteration to set the acquisition to images. """
        for mesh in data.EMSet.iterItems(self, orderBy=orderBy, direction=direction,
                                 where=where, limit=limit):
            yield mesh

