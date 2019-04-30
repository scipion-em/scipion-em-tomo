# **************************************************************************
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

import pyworkflow as pw
import pyworkflow.em as pwem
from pyworkflow.protocol.params import PointerParam
from pyworkflow.mapper.sqlite_db import SqliteDb

import tomo.objects


class ProtTomoBase:
    def __createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pw.utils.cleanPath(setFn)

        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj

    def _createSetOfTiltSeriesM(self, suffix=''):
            return self.__createSet(tomo.objects.SetOfTiltSeriesM,
                                    'tiltseriesM%s.sqlite', suffix)

    def _createSetOfTiltSeries(self, suffix=''):
            self._ouputSuffix = ''
            return self.__createSet(tomo.objects.SetOfTiltSeries,
                                    'tiltseries%s.sqlite', suffix)

    def _createSetOfCoordinates3D(self, volSet, suffix=''):
        coord3DSet = self.__createSet(tomo.objects.SetOfCoordinates3D,
                                      'coordinates%s.sqlite', suffix,
                                      indexes=['_volId'])
        coord3DSet.setVolumes(volSet)
        return coord3DSet

    def _createSetOfTomograms(self, suffix=''):
        return self.__createSet(tomo.objects.SetOfTomograms,
                                'tomograms%s.sqlite', suffix)


class ProtTomoReconstruct(pwem.EMProtocol, ProtTomoBase):
    """ Base class for Tomogram reconstruction protocols. """
    pass

class ProtTomoPicking(pwem.EMProtocol, ProtTomoBase):

    OUTPUT_PREFIX = 'output3DCoordinates'

    """ Base class for Tomogram boxing protocols. """
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputTomogram', PointerParam, label="Input Tomogram", important=True,
                      pointerClass='Tomogram',
                      help='Select the Tomogram to be used during picking.')

    def _getOutputSuffix(self):
        """ Get the name to be used for a new output.
        For example: output3DCoordinates7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(tomo.objects.SetOfCoordinates3D):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1 # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else '' # empty if not output



