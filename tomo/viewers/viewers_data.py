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
from pwem.viewers import ObjectView, DataViewer, MODE, MODE_MD, VISIBLE, ORDER
from pyworkflow.protocol import LabelParam

from .views import ClassesSubTomogramsView
from ..convert.convert import setOfMeshes2Files
import tomo.objects

SERIES_EVEN = "TiltSeriesEven"
SERIES_ODD = "TiltSeriesOdd"
SERIES_DW = "TiltSeriesDW"


class TomoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    in a dedicated tree view
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfTiltSeriesM,
        tomo.objects.SetOfTomograms,
        tomo.objects.SetOfTiltSeries,
        tomo.objects.SetOfClassesSubTomograms,
        tomo.objects.SetOfMeshes
    ]
    _name = 'Tomo viewer'

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
        if issubclass(cls, tomo.objects.SetOfTiltSeriesBase) or issubclass(cls, tomo.objects.SetOfTomograms):
            # JMRT: Local import to avoid importing Tkinter stuff at top level
            from .views_tkinter_tree import TomoObjectDialogView
            setTsView = TomoObjectDialogView(self.getTkRoot(), self.protocol, obj)
            views.append(setTsView)

        elif issubclass(cls, tomo.objects.SetOfClassesSubTomograms):
            views.append(ClassesSubTomogramsView(self._project, obj.strId(),
                                                 obj.getFileName()))

        elif issubclass(cls, tomo.objects.SetOfMeshes):
            from .views_tkinter_tree import TomogramsTreeProvider, TomogramsDialog
            outputMeshes = obj
            tomos = outputMeshes.getPrecedents()
            volIds = outputMeshes.aggregate(["MAX", "COUNT"], "_volId", ["_volId"])
            volIds = [(d['_volId'], d["COUNT"]) for d in volIds]

            tomoList = []
            for objId in volIds:
                tomogram = tomos[objId[0]].clone()
                tomogram.count = objId[1]
                tomoList.append(tomogram)
            path = self.protocol._getExtraPath()
            tomoProvider = TomogramsTreeProvider(tomoList, path, 'txt', )
            path = os.path.join(path, '..')
            setOfMeshes2Files(outputMeshes, path)
            TomogramsDialog(self._tkRoot, True, provider=tomoProvider, path=path, lockGui=False)

        return views


class CtfEstimationTomoViewer(pwviewer.Viewer):
    """ This class implements a view using Tkinter CtfEstimationListDialog
    and the CtfEstimationTreeProvider.
    """
    _name = 'Tomo CTF'
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

    def plotExtra(self, ctfSet, ctfId):
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
                                plot2Dfunc=self.getPlot2DCallback(),
                                plotExtrafunc=self.getPlotExtraCallback())

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

    def getPlotExtraCallback(self):
        if not pwutils.isSameFunction(self.plotExtra,
                                      CtfEstimationTomoViewer.plotExtra):
            return self.plotExtra
        return None


# Register specific sets in pwem dataviewer.
DataViewer.registerConfig(tomo.objects.SetOfSubTomograms,
                          config={MODE: MODE_MD,
                                  VISIBLE: 'id _filename _volName _coordinate._x _coordinate._y _coordinate._z '
                                           '_transform._matrix '})

DataViewer.registerConfig(tomo.objects.SetOfCoordinates3D,
                          config={MODE: MODE_MD,
                                  VISIBLE: 'id _tomoId _x _y _z _groupId _eulerMatrix._matrix _score'})

DataViewer.registerConfig(tomo.objects.SetOfTiltSeriesCoordinates)
DataViewer.registerConfig(tomo.objects.SetOfMeshes)
DataViewer.registerConfig(tomo.objects.SetOfLandmarkModels)
DataViewer.registerConfig(tomo.objects.SetOfCTFTomoSeries)

labels = ('id %s _filename _samplingRate _ctfCorrected _halfMapFiles %s' %
          (tomo.objects.Tomogram.TS_ID_FIELD, tomo.objects.Tomogram.ORIGIN_MATRIX_FIELD))
DataViewer.registerConfig(tomo.objects.SetOfTomograms,
                          config={MODE:MODE_MD,
                                  ORDER: labels,
                                  VISIBLE: labels,})

# Register Tilt series sets
labels = ('id enabled _tsId _samplingRate _anglesCount _acquisition._tiltAxisAngle _hastCtf _hasOddEven _ctfCorrected'
          ' _hasAlignment _acquisition._angleMin _acquisition._angleMax'
          ' _acquisition._step _size _acquisition._accumDose'
)
DataViewer.registerConfig(tomo.objects.SetOfTiltSeries,
                          config={MODE: MODE_MD,
                                  ORDER: labels,
                                  VISIBLE: labels})

# TomoMasks
labels = 'id _tsId _filename _volName _samplingRate'
DataViewer.registerConfig(tomo.objects.SetOfTomoMasks,
                          config={MODE: MODE_MD,
                                  ORDER: labels,
                                  VISIBLE: labels})
