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

import os

import pyworkflow.utils as pwutils

import pyworkflow.viewer as pwviewer
from pwem.protocols import EMProtocol
from pwem.viewers import ObjectView
from pyworkflow.protocol import LabelParam

from .views import ClassesSubTomogramsView
from ..convert.convert import setOfMeshes2Files
import tomo.objects
from tomo.protocols import ProtTsCorrectMotion

SERIES_EVEN = "outputTiltSeriesEven"
SERIES_ODD = "outputTiltSeriesOdd"


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
            from .views_tkinter_tree import TomogramsTreeProvider, TomogramsDialog
            outputMeshes = obj
            tomoList = [item.clone() for item in outputMeshes.iterVolumes()]
            path = self.protocol._getExtraPath()
            tomoProvider = TomogramsTreeProvider(tomoList, path, 'txt', )
            path = os.path.join(path, '..')
            setOfMeshes2Files(outputMeshes, path)
            setView = TomogramsDialog(self._tkRoot, True, provider=tomoProvider, path=path)

        return views


class TSMotionCorrectionViewer(pwviewer.ProtocolViewer):
    """ Wrapper to visualize outputs of tilt series motion correction protocols
    """
    _label = 'Tilt series motion correction viewer'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [ProtTsCorrectMotion]

    def _defineParams(self, form):
        form.addSection(label='Visualization of tilt series')
        form.addParam('displayFullTiltSeries', LabelParam,
                      label='Display *full* frame aligned tilt series',
                      help='Shows full frames aligned set of tilt series')
        if self.hasEvenSet():
            form.addParam('displayEvenTiltSeries', LabelParam,
                          label='Display *even* frames aligned tilt series',
                          help='Shows even frames aligned set of tilt series')

            if self.hasOddSet():
                form.addParam('displayOddTiltSeries', LabelParam,
                              label='Display *odd* frames aligned tilt series',
                              help='Shows even frames aligned set of tilt series')

    def hasEvenSet(self):
        return hasattr(self.protocol, SERIES_EVEN)

    def hasOddSet(self):
        return hasattr(self.protocol, SERIES_ODD)

    def getEvenSet(self):
        return getattr(self.protocol, SERIES_EVEN)

    def getOddSet(self):
        return getattr(self.protocol, SERIES_ODD)

    def _displayEvenTiltSeries(self, param=None):
        return self._visualize(self.getEvenSet())

    def _displayOddTiltSeries(self, param=None):
        return self._visualize(self.getOddSet())

    def _displayFullTiltSeries(self, param=None):
        return self._visualize(self.protocol.outputTiltSeries)

    def _getVisualizeDict(self):
        return {
            'displayFullTiltSeries': self._displayFullTiltSeries,
            'displayEvenTiltSeries': self._displayEvenTiltSeries,
            'displayOddTiltSeries': self._displayOddTiltSeries,
        }

    def _visualize(self, setOfTiltSeries):
        from .views_tkinter_tree import TiltSeriesDialogView
        setTsView = TiltSeriesDialogView(self.getTkRoot(), self.protocol, setOfTiltSeries)

        return [setTsView]


class CtfEstimationTomoViewer(pwviewer.Viewer):
    """ This class implements a view using Tkinter CtfEstimationListDialog
    and the CtfEstimationTreeProvider.
    """
    _label = 'CTF estimation viewer'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [tomo.objects.SetOfCTFTomoSeries]

    def __init__(self, parent, protocol, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._tkParent = parent.root
        self._protocol = protocol
        self._title = 'CTF estimation viewer'

    def plot1D(self, ctfSet, ctfId):
        """ To be implemented in the viewers. """
        return None

    def plot2D(self, ctfSet, ctfId):
        """ To be implemented in the viewers. """
        return None

    def visualize(self, obj, windows=None, protocol=None):
        if not isinstance(obj, EMProtocol):
            self.visualizeSet(obj)
        else:
            for name, output in self._protocol._iterOutputsNew():
                if isinstance(output, tomo.objects.SetOfCTFTomoSeries):
                    self.visualizeSet(output)

    def visualizeSet(self, obj):
        # JMRT: Local import to avoid importing Tkinter stuff at top level
        from .views_tkinter_tree import (CtfEstimationTreeProvider,
                                         CtfEstimationListDialog)
        self._inputSetOfTiltSeries = obj.getSetOfTiltSeries()
        self._provider = CtfEstimationTreeProvider(self._tkParent,
                                                   self._protocol,
                                                   obj)

        CtfEstimationListDialog(self._tkParent, self._title, self._provider,
                                self._protocol, self._inputSetOfTiltSeries,
                                plot1Dfunc=self.getPlot1DCallback(),
                                plot2Dfunc=self.getPlot2DCallback())

    def getPlot1DCallback(self):
        if not pwutils.isSameFunction(self.plot1D,
                                      CtfEstimationTomoViewer.plot1D):
            return self.plot1D
        return None

    def getPlot2DCallback(self):
        if not pwutils.isSameFunction(self.plot2D,
                                      CtfEstimationTomoViewer.plot2D):
            return self.plot2D
        return None
