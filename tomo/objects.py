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


class TiltImage(data.Micrograph):
    """ Tilt image """
    def __init__(self, location=None, **kwargs):
        data.Micrograph.__init__(self, location, **kwargs)
        self._tiltAngle = pwobj.Float()


class TiltSeries(data.SetOfMicrographs):
    ITEM_TYPE = TiltImage


class SetOfTiltSeries(data.EMSet):
    ITEM_TYPE = TiltSeries


class TiltImageM(data.Movie):
    """ Tilt movie. """
    def __init__(self, location=None, **kwargs):
        data.Movie.__init__(self, location, **kwargs)
        self._tiltAngle = pwobj.Float()


class TiltSeriesM(data.SetOfMovies):
    ITEM_TYPE = TiltImageM


class SetOfTiltSeriesM(data.EMSet):
    ITEM_TYPE = TiltSeriesM


class Tomogram(data.Volume):
    pass


class SetOfTomograms(data.SetOfVolumes):
    ITEM_TYPE = Tomogram


