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
from pwem.viewers import ObjectView

from .views import ClassesSubTomogramsView
import tomo.objects


class TomoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfTiltSeriesM,
        tomo.objects.SetOfTiltSeries,
        tomo.objects.SetOfClassesSubTomograms,
        tomo.objects.SetOfMeshes
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        # For now handle both types of SetOfTiltSeries together
        if issubclass(cls, tomo.objects.SetOfTiltSeriesBase):
            # JMRT: Local import to avoid importing Tkinter stuff at top level
            from .views_tkinter_tree import TiltSeriesDialogView
            setTsView = TiltSeriesDialogView(self.getTkRoot(), self.protocol, obj)
            views.append(setTsView)

        elif issubclass(cls, tomo.objects.SetOfClassesSubTomograms):
            views.append(ClassesSubTomogramsView(self._project, obj.strId(),
                                                   obj.getFileName()))

        elif issubclass(cls, tomo.objects.SetOfMeshes):
            from .views_tkinter_tree import MeshesTreeProvider, TomogramsDialog
            meshList = [item.clone() for item in obj.iterItems()]

            tomoProvider = MeshesTreeProvider(meshList)

            path = self.protocol._getTmpPath()
            setView = TomogramsDialog(self._tkRoot, True, provider=tomoProvider, path=path)

        return views
