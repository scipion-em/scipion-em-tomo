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

import pyworkflow.object as pwobj
import pyworkflow.em.data as data


class TiltImageBase:
    """ Base class for TiltImageM and TiltImage. """
    def __init__(self, **kwargs):
        self._tiltAngle = pwobj.Float(kwargs.get('tiltAngle', None))
        self._acquisitionOrder = pwobj.Integer(kwargs.get('acquisitionOrder',
                                                          None))

    def getTiltAngle(self):
        return self._tiltAngle.get()

    def setTiltAngle(self, value):
        self._tiltAngle.set(value)

    def getAcquisitionOrder(self):
        return self._acquisitionOrder.get()

    def setAcquisitionOrder(self, value):
        self._acquisitionOrder.set(value)


class TiltImage(data.Image, TiltImageBase):
    """ Tilt image """
    def __init__(self, location=None, **kwargs):
        data.Image.__init__(self, location, **kwargs)
        TiltImageBase.__init__(self, **kwargs)


class TiltSeriesBase(data.SetOfImages):
    pass


class TiltSeries(TiltSeriesBase):
    ITEM_TYPE = TiltImage


class SetOfTiltSeriesBase(data.EMSet):
    """ Base class for SetOfTiltImages and SetOfTiltImagesM.
    """
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
        classPrefix = 'TS%03d' % item.getObjId()
        item._mapperPath.set('%s,%s' % (self.getFileName(), classPrefix))
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
        classItem = data.EMSet.__getitem__(self, itemId)
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


class Tomogram(data.Volume):
    pass


class SetOfTomograms(data.SetOfVolumes):
    ITEM_TYPE = Tomogram


