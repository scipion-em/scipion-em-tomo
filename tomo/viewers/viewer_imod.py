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

import pyworkflow.viewer as pwviewer

import tomo.objects


class ImodViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.TiltSeries
    ]

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.TiltSeries):
            views.append(ImodTiltSerieView(obj))

        return views
    

class ImodTiltSerieView(pwviewer.CommandView):
    """ Wrapper to visualize different type of objects with the 3dmod.
    """
    def __init__(self, tiltSerie, **kwargs):
        # FIXME: Only will work for tiltSeries in a single stack
        pwviewer.CommandView.__init__(
            self, '3dmod "%s"' % tiltSerie.getFirstItem().getFileName())

