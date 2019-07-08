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

from collections import OrderedDict

import pyworkflow.object as pwobj
import pyworkflow.em.data as data


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


class TiltImage(data.Image, TiltImageBase):
    """ Tilt image """
    def __init__(self, location=None, **kwargs):
        data.Image.__init__(self, location, **kwargs)
        TiltImageBase.__init__(self, **kwargs)

    def copyInfo(self, other):
        data.Image.copyInfo(self, other)
        self.copyAttributes(other, '_tiltAngle', '_tsId')


class TiltSeriesBase(data.SetOfImages):
    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)
        self._tsId = pwobj.String(kwargs.get('tsId', None))

    def getTsId(self):
        """ Get unique TiltSerie ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but
        not _mapperPath or _size from other set of micrographs to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])

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


class TiltSeries(TiltSeriesBase):
    ITEM_TYPE = TiltImage


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
        item._mapperPath.setStore(False)
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

class TomoAcquisition:
    def __init__(self, **kwargs):
        self._acquisitionAngleMin = pwobj.Float(kwargs.get('acquisitionAngleMin', None))
        self._acquisitionAngleMax = pwobj.Float(kwargs.get('acquisitionAngleMax', None))
        self._step = pwobj.Integer(None)
        self._angleAxis1 = pwobj.Float(None)
        self._angleAxis2 = pwobj.Float(None)

    def getAcquisitionAngleMax(self):
        return self._acquisitionAngleMax.get()

    def getAcquisitionAngleMin(self):
        return self._acquisitionAngleMin.get()

    def setAcquisitionAngleMax(self, value):
        self._acquisitionAngleMax.set(value)

    def setAcquisitionAngleMin(self, value):
        self._acquisitionAngleMin.set(value)

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

class Tomogram(data.Volume, TomoAcquisition):
    def __init__(self, **kwargs):
        data.Volume.__init__(self, **kwargs)
        TomoAcquisition.__init__(self, **kwargs)
        self._tsId = pwobj.String(kwargs.get('tsId', None))

    def getTsId(self):
        """ Get unique TiltSeries ID, usually retrieved from the
        file pattern provided by the user at the import time.
        """
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

class SetOfTomograms(data.SetOfVolumes):
    ITEM_TYPE = Tomogram
    EXPOSE_ITEMS = True

class TiltSeriesDict:
    """ Helper class that to store TiltSeries and TiltImage but
    using dictionaries for quick access.
    """
    def __init__(self):
        self.__dict = OrderedDict()

    def addTs(self, tiltSeries, includeTi=False):
        """ Add a clone of the tiltseries. """
        self.__dict[tiltSeries.getTsId()] = (tiltSeries.clone(), OrderedDict())
        if includeTi:
            for ti in tiltSeries:
                self.addTi(ti)

    def getTs(self, tsId):
        return self.__dict[tsId][0]

    def addTi(self, ti):
        self.getTiDict(ti.getTsId())[ti.getObjId()] = ti.clone()

    def getTi(self, tsId, tiObjId):
        return self.getTiDict(tsId)[tiObjId]

    def getTiDict(self, tsId):
        return self.__dict[tsId][1]

    def getTiList(self, tsId):
        return self.getTiDict(tsId).values()

    def __iter__(self):
        for ts, d in self.__dict.values():
            yield ts


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
        return (self.getX(), self.getY(), self.getZ())

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
        self._volumesPointer = pwobj.Pointer()
        self._boxSize = pwobj.Integer()

    def getBoxSize(self):
        """ Return the box size of the particles.
        """
        return self._boxSize.get()

    def setBoxSize(self, boxSize):
        """ Set the box size of the particles. """
        self._boxSize.set(boxSize)

    def iterVolumes(self):
        """ Iterate over the micrographs set associated with this
        set of coordinates.
        """
        return self.getVolumes()

    def iterVolumeCoordinates(self, volume):
        """ Iterates over the set of coordinates belonging to that micrograph.
        """
        pass

    def iterCoordinates(self, volume=None):
        """ Iterate over the coordinates associated with a micrograph.
        If micrograph=None, the iteration is performed over the whole
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
        coordWhere = '1' if volId is None else '_volId=%d' % volId

        for coord in self.iterItems(where=coordWhere):
            yield coord

    def getVolumes(self):
        """ Returns the SetOfMicrographs associated with
        this SetOfCoordinates"""
        return self._volumesPointer.get()

    def setVolumes(self, volumes):
        """ Set the micrographs associated with this set of coordinates.
        Params:
            micrographs: Either a SetOfMicrographs object or a pointer to it.
        """
        if volumes.isPointer():
            self._volumesPointer.copy(volumes)
        else:
            self._volumesPointer.set(volumes)

    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths

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

class SubTomogram(data.Volume, TomoAcquisition):
    def __init__(self, **kwargs):
        data.Volume.__init__(self, **kwargs)
        TomoAcquisition.__init__(self, **kwargs)

class SetOfSubTomograms(data.SetOfVolumes):
    ITEM_TYPE = SubTomogram

    def __init__(self, **kwargs):
        data.SetOfVolumes.__init__(self, **kwargs)
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

