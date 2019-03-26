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

import pyworkflow.object as pwobj
import pyworkflow.em as pwem
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog

import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
from pyworkflow.viewer import View

import tomo.objects


class TiltSeriesTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_TS = 'Tilt Series'
    COL_TI = 'Info / Path'
    COL_TI_ID = 'Acq. Order'
    COL_TI_ANGLE = 'Tilt Angle'
    COL_TI_DEFOCUS_U = 'Defocus'
    ORDER_DICT = {COL_TI_ANGLE: '_tiltAngle',
                  COL_TI_DEFOCUS_U: '_ctfModel._defocusU'}

    def __init__(self, protocol, tiltSeries):
        self.protocol = protocol
        self.tiltseries = tiltSeries
        self._hasCtf = tiltSeries.getFirstItem().getFirstItem().hasCTF()
        TreeProvider.__init__(self, sortingColumnName=self.COL_TS)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200

    def getObjects(self):
        # Retrieve all objects of type className
        project = self.protocol.getProject()
        objects = []

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for ts in self.tiltseries:
            tsObj = ts.clone()
            tsObj._allowsSelection = True
            tsObj._parentObject = None
            objects.append(tsObj)
            for ti in ts.iterItems(orderBy=orderBy, direction=direction):
                tiObj = ti.clone()
                tiObj._allowsSelection = False
                tiObj._parentObject = tsObj
                objects.append(tiObj)

        return objects

    def _sortObjects(self, objects):
        pass

    def objectKey(self, pobj):
        pass

    def getColumns(self):
        cols = [
            (self.COL_TS, 100),
            (self.COL_TI_ID, 80),
            (self.COL_TI_ANGLE, 80),
            (self.COL_TI, 400),
        ]

        if self._hasCtf:
            cols.insert(3, (self.COL_TI_DEFOCUS_U, 80))

        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj):
        objId = obj.getObjId()
        tsId = obj.getTsId()

        if isinstance(obj, tomo.objects.TiltSeriesBase):
            key = obj.getObjId()
            text = tsId
            values = ['', '', '', str(obj)]
            opened = True
        else:  # TiltImageBase
            key = '%s.%s' % (tsId, obj.getObjId())
            text = ''
            values = [objId, obj.getTiltAngle(), str(obj.getLocation())]
            if self._hasCtf:
                values.insert(2, "%.03f" % obj.getCTF().getDefocusU())
            opened = False

        return {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': False,
            'parent': obj._parentObject
        }

    def getObjectActions(self, obj):
        actions = []

        if isinstance(obj, tomo.objects.TiltSeries):
            viewers = pwem.findViewers(obj.getClassName(),
                                       pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                def createViewer(viewerClass, obj):
                    proj = self.protocol.getProject()
                    item = self.tiltseries[obj.getObjId()]  # to load mapper
                    return lambda : viewerClass(project=proj).visualize(item)
                actions.append(('Open with %s' % viewerClass.__name__,
                                createViewer(viewerClass, obj)))

        return actions


class TiltSeriesDialogView(View):
    """ This class implements a view using Tkinter ListDialog
    and the TiltSeriesTreeProvider.
    """
    def __init__(self, parent, protocol, tiltSeries, **kwargs):
        """
         Params:
            parent: Tkinter parent widget


        From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback:
                    a callback function to validate selected items.
                allowSelect: if set to False, the 'Select' button will not
                    be shown.
                allowsEmptySelection: if set to True, it will not validate
                    that at least one element was selected.
        """
        self._tkParent = parent
        self._protocol = protocol
        self._provider = TiltSeriesTreeProvider(protocol, tiltSeries)

    def show(self):
        dlg = ListDialog(self._tkParent, 'TiltSeries display', self._provider)
