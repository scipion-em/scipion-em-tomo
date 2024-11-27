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
import ast
import glob
import os.path
import threading
from tkinter import messagebox, BOTH, RAISED

import numpy
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from pwem.convert.transformations import euler_from_matrix
from pwem.emlib.image.image_readers import ImageReadersRegistry, ImageStack
from pwem.viewers import showj
from pwem.viewers.showj import runJavaIJapp
from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog, ToolbarListDialog, showInfo
import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.plugin import Domain

import tomo.objects
from . import TomoDataViewer
from ..convert.convert import getMeshVolFileName
from ..objects import CTFTomo

# How many standard deviations to truncate above and below the mean when increasing contrast:
CONTRAST_STD = 2.0
afterId = None  # Global variable to store the reference of the auto navigate task
DIRECTION_DOWN = 0
DIRECTION_UP = 1
RESULT_CREATE_WITH = 0
RESULT_RESTACK = 1
THUMBNAIL_SIZE = 512


class TSTreeToolTip:
    def __init__(self, widget, text=""):
        self.widget = widget
        self.text = text
        self.tooltip = None

    def show(self, event):
        if self.tooltip is None:
            self.tooltip = tk.Toplevel(self.widget)
            self.tooltip.wm_overrideredirect(True)
            self.tooltip.wm_geometry(f"+{event.x_root+20}+{event.y_root+20}")
            label = tk.Label(self.tooltip, text=self.text, relief="solid", borderwidth=1)
            label.pack()

    def hide(self, event):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip = None


class TiltSerieState:
    EXCLUDED = 'excludedTs'
    INCLUDED = 'includedTs'


class TiltImageStates:
    EXCLUDED = 'excluded'
    INCLUDED = 'included'
    ODD = 'odd'
    EVEN = 'even'
    CHECK_MARK = "\u2611"    # ☑
    CHECK_UNMARK = "\u2610"  # ☐


class TiltSeriesTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_TS = 'Tilt series'
    COL_TI = 'Path'
    COL_TI_ANGLE = 'Tilt angle'
    COL_TI_ENABLED = 'Excluded'
    COL_TI_ACQ_ORDER = 'Order'
    COL_TI_DOSE = "Dose"
    COL_TI_TRANSFORM = "T. Matrix"
    COL_TI_ROT_ANGLE = "Rot"
    COL_TI_SHIFTX = 'ShiftX'
    COL_TI_SHIFTY = 'ShiftY'
    ORDER_DICT = {COL_TI_ANGLE: '_tiltAngle',
                  COL_TI_ACQ_ORDER: '_acqOrder',
                  COL_TI_DOSE: '_acquisition._accumDose'}

    def __init__(self, protocol, tiltSeries):
        self.protocol = protocol
        self.tiltSeries = tiltSeries
        TreeProvider.__init__(self, sortingColumnName=self.COL_TS)
        self.excludedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200
        self.updatedCount = 0
        self.changes = 0
        self.objects = []

    def configureTags(self, tree):
        self.tree = tree
        standardFont = getDefaultFont()
        self.tree.tag_configure(TiltSerieState.EXCLUDED, font=standardFont, foreground='red',
                                image=getImage(Icon.CHECKED))
        self.tree.tag_configure(TiltSerieState.INCLUDED, font=standardFont, foreground='black',
                                image=getImage(Icon.UNCHECKED))

        self.tree.tag_configure(TiltImageStates.EXCLUDED,  font=standardFont, foreground='red')
        self.tree.tag_configure(TiltImageStates.INCLUDED,  font=standardFont, foreground='black')
        self.tree.tag_configure(TiltImageStates.EVEN, background='#F2F2F2', foreground='black')
        self.tree.tag_configure(TiltImageStates.ODD, background='#E6E6E6', foreground='black')
        self.tree.bind('<space>', self.onSpace)
        tsTooltip = TSTreeToolTip(self.tree)
        self.tree.bind("<Motion>", lambda event: self.onMouseMotion(event, self.tree, tsTooltip))
        self.tree.bind("<Leave>", lambda event: tsTooltip.hide(event))

    def onMouseMotion(self, event, tsTree, tsTooltip):
        region = tsTree.identify_region(event.x, event.y)
        if region == "cell":
            column = tsTree.identify_column(event.x)
            if column == "#6" or column == "#7" or column == "#8":  # Display the matrix transformation
                tsTooltip.hide(event)
                item_id = tsTree.identify_row(event.y)
                if item_id and len(tsTree.item(item_id)['values']) > 5:
                    tMatrix = tsTree.item(item_id)['values'][8]
                    if tMatrix:
                        mlist = ast.literal_eval(tMatrix)
                        matriz = [mlist[i:i + 3] for i in range(0, len(mlist), 3)]
                        formattedMatrix = "\nTransformation matrix:\n\n[ " + "\n".join(
                            ["  ".join(f"{num:.4f}" for num in fila) for fila in matriz]) + " ]\n"
                        tsTooltip.text = f"{formattedMatrix}"
                        tsTooltip.show(event)
            else:
                tsTooltip.hide(event)
        else:
            tsTooltip.hide(event)

    def getTree(self):
        return self.tree

    def onSpace(self, event):
        selectedItem = self.tree.selection()[0]
        if selectedItem:
            obj = event.widget._objDict[selectedItem]
            self._itemSelected(obj)

    def _itemSelected(self, obj):
        _, y, _, _ = self.tree.bbox(self.tree.selection()[0])
        selectedItem = self.tree.identify_row(y + 1)
        itemValues = self.tree.item(selectedItem, 'values')

        if isinstance(obj, tomo.objects.TiltSeriesBase):
            excludedTS = TiltSerieState.INCLUDED if obj.isEnabled() else TiltSerieState.EXCLUDED
            obj.setEnabled(not obj.isEnabled())
            tags = TiltSerieState.INCLUDED
            if TiltSerieState.EXCLUDED in self.tree.item(selectedItem, 'tags'):
                tags = TiltSerieState.INCLUDED

            if excludedTS == TiltSerieState.INCLUDED:
                self.tree.item(selectedItem, tags=(TiltSerieState.EXCLUDED, tags,))
                excludedTi = TiltImageStates.CHECK_MARK
            else:
                self.tree.item(selectedItem, tags=(TiltSerieState.INCLUDED, tags,))
                excludedTi = TiltImageStates.CHECK_UNMARK

            children = self.tree.get_children(selectedItem)
            index, ts = self.getTiltSerie(obj)
            size = ts.getSize()
            for tiIndex in range(size):
                obj = self.objects[tiIndex + index + 1]
                if excludedTS == TiltSerieState.INCLUDED:
                    obj.setEnabled(False)
                else:
                    obj.setEnabled(True)
                item = children[tiIndex]
                itemValues = self.tree.item(item, 'values')
                self.excludeTiltImage(obj, item, itemValues, excludedTi)

        else:
            excludedTi = TiltImageStates.CHECK_MARK if obj.isEnabled() else TiltImageStates.CHECK_UNMARK
            obj.setEnabled(not obj.isEnabled())
            self.excludeTiltImage(obj, selectedItem, itemValues, excludedTi)
            # Verify if it is necessary to check or uncheck the tiltserie
            parentItem = self.tree.parent(selectedItem)
            parentObj = obj._parentObject

            if excludedTi == TiltImageStates.CHECK_UNMARK:
                self.tree.item(parentItem, tags=(TiltSerieState.INCLUDED, TiltSerieState.INCLUDED))
                parentObj.setEnabled(True)
            # Checking the tiltserie if all tiltimages are checked
            elif obj._parentObject.getSize() == len(self.excludedDict[self.tree.item(self.tree.parent(selectedItem))['text']]):
                self.tree.item(self.tree.parent(selectedItem), tags=(TiltSerieState.EXCLUDED, TiltSerieState.EXCLUDED))
                parentObj.setEnabled(False)

    def excludeTiltImage(self, obj,  selectedItem, itemValues, excluded):
        if isinstance(obj, tomo.objects.TiltImageM):
            newValues = (itemValues[0], itemValues[1], excluded, itemValues[3], itemValues[4])
        else:
            newValues = (itemValues[0], itemValues[1], excluded, itemValues[3], itemValues[4], itemValues[5],
                         itemValues[6], itemValues[7], itemValues[8])
        tags = TiltImageStates.EVEN
        if TiltImageStates.ODD in self.tree.item(selectedItem, 'tags'):
            tags = TiltImageStates.ODD

        if excluded == TiltImageStates.CHECK_MARK:
            self.tree.item(selectedItem, tags=(TiltImageStates.EXCLUDED, tags,))
            self.excludedDict[self.tree.item(self.tree.parent(selectedItem))['text']][obj.getObjId()] = True
            self.updatedCount += 1
            self.changes += 1
        else:
            self.tree.item(selectedItem, tags=(TiltImageStates.INCLUDED, tags,))
            self.excludedDict[self.tree.item(self.tree.parent(selectedItem))['text']].pop(obj.getObjId())
            self.updatedCount -= 1
            self.changes -= 1
        self.tree.item(selectedItem, values=newValues)

    def _toggleExclusion(self):
        itemSelected = self.tree.selection()
        self.tree.focus(itemSelected)
        self.tree.selectChild(itemSelected)
        itemId = self.tree.item(itemSelected)['text']
        itemParent = self.tree.parent(itemSelected)
        if not itemParent:
            _, obj = self.getTiltSerie(itemId)
        else:
            tsId = itemSelected[0].split('.')[0]
            _, obj = self.getTiltImage(tsId, itemId)
        self._itemSelected(obj)

    def getObjects(self):
        # Retrieve all objects of type className

        if self.objects:
            return self.objects

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for ts in self.tiltSeries:
            self.excludedDict[ts.getTsId()] = {}
            tsObj = ts.clone(ignoreAttrs=['_mapperPath'])
            tsObj._allowsSelection = True
            tsObj._parentObject = None
            tsObj.setEnabled(ts.isEnabled())
            self.objects.append(tsObj)
            for ti in ts.iterItems(orderBy=orderBy, direction=direction):
                tiObj = ti.clone()  # For some reason .clone() does not clone the enabled nor the creation time
                if not ti.isEnabled():
                    self.excludedDict[ts.getTsId()][ti.getObjId()] = True
                tiObj.setEnabled(ti.isEnabled())
                tiObj._allowsSelection = False
                tiObj._parentObject = tsObj
                self.objects.append(tiObj)

        return self.objects

    def getTiltSerieRepresentative(self, tiltSerie):
        """This method returns the central tiltImage of the set."""
        size = tiltSerie.getSize()
        objects = self.getObjects()
        for index, obj in enumerate(objects):
            if obj == tiltSerie:
                return self.objects[index + int(size/2)]

    def getTiltSerie(self, tiltSerie):
        objects = self.getObjects()
        is_tsId = isinstance(tiltSerie, str)

        for index, obj in enumerate(objects):
            if isinstance(obj, tomo.objects.TiltSeriesBase):
                if (is_tsId and obj.getTsId() == tiltSerie) or \
                        (not is_tsId and obj.getObjId() == tiltSerie.getObjId()):
                    return index, obj

        return None, None

    def getTiltImage(self, tsId, treeId):
        objects = self.getObjects()
        for index, obj in enumerate(objects):
            if isinstance(obj, tomo.objects.TiltImageBase):
                if obj.getTsId() == tsId and obj.getObjId() == treeId:
                    return index, obj

    def getExcludedViews(self):
        return self.excludedDict

    def getUpdatedCount(self):
        return self.updatedCount

    def getchanges(self):
        return self.changes

    def resetChanges(self):
        self.changes = 0

    def _sortObjects(self, objects):
        pass

    def objectKey(self, pobj):
        pass

    def getColumns(self):
        cols = [
            (self.COL_TS, 80),
            (self.COL_TI_ACQ_ORDER, 50),
            (self.COL_TI_ANGLE, 70),
            (self.COL_TI_ENABLED, 65),
            (self.COL_TI_DOSE, 55),
            (self.COL_TI, 400),
        ]
        if not isinstance(self.tiltSeries, tomo.objects.SetOfTiltSeriesM):
            cols.append((self.COL_TI_ROT_ANGLE, 100))
            cols.append((self.COL_TI_SHIFTX, 100))
            cols.append((self.COL_TI_SHIFTY, 100))
            # cols.append((self.COL_TI_TRANSFORM, 200))


        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj: tomo.objects.TiltImageBase):
        objId = obj.getObjId()
        tsId = obj.getTsId()

        if isinstance(obj, tomo.objects.TiltSeriesBase):
            key = objId
            text = tsId
            values = ['', '', '', '', str(obj)]
            opened = False
            tags = TiltSerieState.INCLUDED
            if not obj.isEnabled():
                obj.setEnabled(False)
                tags = TiltSerieState.EXCLUDED

        else:  # TiltImageBase
            key = '%s.%s' % (tsId, objId)
            text = objId

            dose = obj.getAcquisition().getAccumDose() if hasattr(obj.getAcquisition(), '_accumDose') else None
            adqOrder = obj.getAcquisitionOrder() if hasattr(obj, '_acqOrder') else None
            excluded = TiltImageStates.CHECK_UNMARK if obj.isEnabled() else TiltImageStates.CHECK_MARK

            values = [str("%d" % adqOrder) if adqOrder is not None else "",
                      str("%0.2f" % obj.getTiltAngle()),
                      excluded,
                      round(dose, 2) if dose is not None else "",
                      "%d@%s" % (obj.getIndex() or 1, obj.getFileName()),
                      ]

            if not isinstance(obj, tomo.objects.TiltImageM):
                matrix = ''
                angle = ''
                shiftX = ''
                shiftY = ''
                if obj.hasTransform():
                    transform = obj.getTransform()
                    # TODO These lines will be removed when transform.getEulerAngles is released
                    rotation = transform.getRotationMatrix()
                    _, _, rot = euler_from_matrix(rotation)
                    angle = str("%0.2f" % numpy.rad2deg(rot))
                    matrixList = transform.getMatrixAsList()
                    shiftX = str("%0.2f" % matrixList[2])
                    shiftY = str("%0.2f" % matrixList[5])
                    matrix = str([round(num, 4) for num in matrixList])
                values.append(angle)
                values.append(shiftX)
                values.append(shiftY)
                values.append(matrix)

            opened = False

            tags = TiltImageStates.INCLUDED
            if not obj.isEnabled():
                obj.setEnabled(False)
                tags = TiltImageStates.EXCLUDED

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': False,
            'parent': obj._parentObject
        }

        if obj.getObjId() % 2 == 0:
            item['tags'] = (tags, TiltImageStates.ODD)
        else:
            item['tags'] = (tags, TiltImageStates.EVEN)

        return item

    def getObjectActions(self, obj):
        actions = []

        if isinstance(obj, tomo.objects.TiltSeriesBase):
            viewers = Domain.findViewers(obj.getClassName(),
                                         pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                def createViewer(viewerClass, obj):
                    proj = self.protocol.getProject()
                    item = self.tiltSeries[obj.getObjId()]  # to load mapper
                    return lambda: viewerClass(project=proj, protocol=self.protocol).visualize(item)
                actions.append((viewerClass.getName(), createViewer(viewerClass, obj)))

        actionName = "Exclude" if obj.isEnabled() else "Include"
        actions.append((actionName, lambda: self._itemSelected(obj)))

        return actions


class TiltSeriesDialog(ToolbarListDialog):
    def __init__(self, parent, title, provider, tiltSeries, protocol, **kwargs):

        self._tiltSeries = tiltSeries
        self._protocol = protocol
        self._provider = provider
        self._applyContrastCallback = kwargs.get('applyContrastCallback', None)


        toolbarButtons = []

        if isinstance(self._tiltSeries, tomo.objects.SetOfTiltSeries):
            viewers = Domain.findViewers(tomo.objects.SetOfTiltSeries.getClassName(),
                                         pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                if viewerClass is not TomoDataViewer:
                    def launchViewer():
                        proj = self._protocol.getProject()
                        viewerInstance = viewerClass(project=proj, protocol=self._protocol)
                        return lambda event: self.launchViewer(viewerInstance)
                    toolbarButtons.append(dialog.ToolbarButton(viewerClass.getName(), launchViewer(), Icon.ACTION_RESULTS))

        toolbarButtons.append(dialog.ToolbarButton('|', None))
        toolbarButtons.append(dialog.ToolbarButton('Save', self._saveExcluded, Icon.ACTION_SAVE,
                              tooltip="Create a new output with excluded views marked", shortcut='<Control-s>'))
        toolbarButtons.append(dialog.ToolbarButton('Help', self._showHelp, Icon.ACTION_HELP))

        ToolbarListDialog.__init__(self, parent, title, provider, toolbarButtons=toolbarButtons, **kwargs)

    def body(self, bodyFrame):
        ToolbarListDialog.body(self, bodyFrame)
        firstTiltImage = self.tree.get_children(self.tree.get_children()[0])[0]
        if firstTiltImage:
            self.tree.selection_set(firstTiltImage)

    def validateClose(self):
        if self._provider.getUpdatedCount() and self._provider.getchanges():
            msg = "Do you want to exit and loose your changes ? "
            result = messagebox.askquestion("Loosing changes", msg, icon='warning', **{'parent': self})
            return result == messagebox.YES
        return True

    def launchViewer(self, viewerInstance):
        itemSelected = self.tree.selection() if not self.tree.parent(self.tree.selection()) else self.tree.parent(self.tree.selection())
        obj = self.tree._objects[self.tree.index(itemSelected) + 1]
        item = self._tiltSeries[obj.getObjId()]
        viewerInstance.visualize(item)

    def _showHelp(self, event=None):
        showInfo('TiltSeries viewer help',
                 'This viewer allows you to exclude or include TiltImages.\n\n'
                 '1. Toggle exclusion button to exclude or include a selected tiltimage or use click over the checkbox or space over the item\n'
                 '2. Increase contrast button to enhance the tiltimage contrast.\n'
                 '3. Save button to create a new set with excluded views marked.', self)

    def _saveExcluded(self, event=None):
        updatedCount = self._provider.getUpdatedCount()
        changes = self._provider.getchanges()

        msg = ("Are you going to create a new set of tiltseries with the excluded views?\n\n"
               "What do you do ? \n\n"
               "Yes: The set will be created with the excluded views. \n"
               "Re-stack: Delete excluded views and create a new TS stack") if updatedCount and changes \
            else ("This set of TiltSeries has already been saved previously or is the original one. \n\n"
                  "Are you still sure you want to generate a new set?")

        d = GenericDialog(self, "Create a new set of tiltseries", msg,
                          Icon.ALERT,
                          buttons=[('Yes', RESULT_CREATE_WITH),
                                   ('Re-stack', RESULT_RESTACK),
                                   ('Cancel', RESULT_CANCEL)],
                          default='Yes',
                          icons={RESULT_CANCEL: Icon.BUTTON_CANCEL,
                                 RESULT_CREATE_WITH: Icon.BUTTON_SELECT,
                                 RESULT_RESTACK: Icon.ACTION_EXECUTE})

        result = d.result
        if result != RESULT_CANCEL:

            self.update_idletasks()
            restack = result == RESULT_RESTACK
            outputSetOfTiltSeries = tomo.objects.SetOfTiltSeries.create(self._protocol.getPath(),
                                                                        suffix=str(self._protocol.getOutputsSize()))
            outputSetOfTiltSeries.copyInfo(self._tiltSeries)
            outputSetOfTiltSeries.setDim(self._tiltSeries.getDim())
            excludedViews = self._provider.getExcludedViews()
            outputName = self._protocol.getNextOutputName('TiltSeries_')
            outputPath = os.path.join(self._protocol._getExtraPath(), outputName)

            if restack:  # Creating a new directory to write a new stack
                if not os.path.exists(outputPath):
                    os.mkdir(outputPath)
            hasOddEven = self._tiltSeries.hasOddEven()
            count = len(self._tiltSeries.getTSIds())
            for tsIndex, ts in enumerate(self._tiltSeries):
                self.info('Processing %s of %s tiltseries ...' % (tsIndex + 1, count))
                self.update_idletasks()
                tsId = ts.getTsId()
                _, obj = self._provider.getTiltSerie(tsId)
                if obj.isEnabled():  # We will exclude the TS that are checked even if the restack option is not chosen.
                    newBinaryName = os.path.join(outputPath, tsId + '.mrcs')
                    if hasOddEven:
                        oddFileName = ts.getOddFileName()
                        newOddBinaryName = os.path.join(outputPath, tsId + '_odd.mrcs')
                        evenFileName = ts.getEvenFileName()
                        newEvenBinaryName = os.path.join(outputPath, tsId + '_even.mrcs')

                    index = 1
                    newTs = tomo.objects.TiltSeries()
                    newTs.copyInfo(ts)
                    outputSetOfTiltSeries.append(newTs)

                    # New Stacks
                    properties = {"sr": obj.getSamplingRate()}
                    newStack = ImageStack(properties=properties)
                    oddFileNames = ImageStack(properties=properties)
                    evenFileNames = ImageStack(properties=properties)

                    for ti in ts.iterItems():
                        included = False if ti.getObjId() in excludedViews[tsId] else True
                        if not restack or (included and restack):
                            newTi = ti.clone()
                            newTi.copyInfo(ti, copyId=False)
                            newTi.setObjId(None)
                            newTi.setAcquisition(ti.getAcquisition())
                            # For some reason .clone() does not clone the enabled nor the creation time
                            newTi.setEnabled(included)
                            if restack:
                                oldIndex = str(ti.getIndex())
                                newStack.append(ImageReadersRegistry.open(oldIndex + '@' + ti.getFileName()))
                                newTi.setLocation((index, newBinaryName))

                                if hasOddEven:
                                    oddFileNames.append(ImageReadersRegistry.open(oldIndex + '@' + oddFileName))
                                    evenFileNames.append(ImageReadersRegistry.open(oldIndex + '@' + evenFileName))
                                    newTi.setOddEven([newOddBinaryName, newEvenBinaryName])

                                index += 1
                            newTs.append(newTi)

                    if restack:
                        ImageReadersRegistry.write(newStack, newBinaryName, isStack=True)
                        if hasOddEven:
                            ImageReadersRegistry.write(oddFileNames, newOddBinaryName, isStack=True)
                            ImageReadersRegistry.write(evenFileNames, newEvenBinaryName, isStack=True)

                    if len(excludedViews[ts.getTsId()]) == ts.getSize():
                        newTs.setEnabled(False)
                    newTs.setDim(ts.getDim())
                    newTs.setAnglesCount(newTs.getSize())
                    newTs.write()
                    outputSetOfTiltSeries.update(newTs)
                    outputSetOfTiltSeries.write()

            self._protocol._defineOutputs(**{outputName: outputSetOfTiltSeries})
            self.info('The new set (%s) has been created successfully' % outputName)
            self._provider.resetChanges()


class TiltSeriesDialogView(pwviewer.View):
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
        self._tiltSeries = tiltSeries
        self._provider = TiltSeriesTreeProvider(self._protocol, self._tiltSeries)
        self._preview = None  # To store preview widget

    def show(self):
        previewCallback = self.previewTiltSeries
        if isinstance(self._tiltSeries, tomo.objects.SetOfTiltSeriesM):
            previewCallback = None

        TiltSeriesDialog(self._tkParent, 'Tilt series viewer', self._provider, self._tiltSeries, self._protocol,
                         lockGui=False, previewCallback=previewCallback,
                         itemOnClick=self.itemOnClick, allowSelect=False, cancelButton=True)

    def getPreviewWidget(self, frame):

        if not self._preview:

            actionBar = tk.Frame(frame, bd=1, relief=tk.SUNKEN)
            actionBar.grid(row=0, column=0, sticky='ns')

            self.startButton = tk.Button(actionBar, image=getImage(Icon.ACTION_CONTINUE), command=self.autoNavigate,
                                         width=25, height=25, relief=tk.RAISED)
            self.startButton.grid(row=0, column=0, sticky='ns')
            ToolTip(self.startButton, text='Auto navigate', delay=0)

            self.next = tk.Button(actionBar, image=getImage(Icon.ACTION_FIND_NEXT), command=lambda: self.navigate(DIRECTION_DOWN),
                                         width=25, height=25, relief=tk.RAISED)
            self.next.grid(row=0, column=1, sticky='ns')
            ToolTip(self.next, text='Next', delay=0)

            self.previous = tk.Button(actionBar, image=getImage(Icon.ACTION_FIND_PREVIOUS), command=lambda: self.navigate(DIRECTION_UP),
                                  width=25, height=25, relief=tk.RAISED)
            self.previous.grid(row=0, column=2, sticky='ns')
            ToolTip(self.previous, text='Previous', delay=0)

            self.stopButton = tk.Button(actionBar, image=getImage(Icon.ACTION_STOP), command=self.stopNavigate, width=25, height=25)
            self.stopButton.grid(row=0, column=3, sticky='ne')
            ToolTip(self.stopButton, text='Stop navigate', delay=0)

            self.autoContrast = tk.Button(actionBar, image=getImage(Icon.ACTION_CONTRAST),
                                          command=self.increaseContrast,
                                          width=25, height=25)
            self.autoContrast.grid(row=0, column=4, sticky='ns')
            ToolTip(self.autoContrast, text='Apply contrast to the selected tiltimage', delay=0)

            self.toogleExclusion = tk.Button(actionBar, image=getImage(Icon.ACTION_CLOSE),
                                             command=self._provider._toggleExclusion,
                                             width=25, height=25)
            self.toogleExclusion.grid(row=0, column=5, sticky='ns', command=None)
            ToolTip(self.toogleExclusion, text='Exclude or include the selection', delay=0)

            from pyworkflow.gui.matplotlib_image import ImagePreview
            self._preview = ImagePreview(frame, dim=500, label="Tilt series")
            self._preview.grid(row=2, column=0)
        return self._preview

    def stopNavigate(self):
        global afterId
        if afterId:
            tree = self._provider.getTree()
            tree.after_cancel(afterId)
            afterId = None
            self.startButton.configure(relief=tk.RAISED, state=tk.NORMAL)
            self.autoContrast.configure(state=tk.NORMAL)
            self.next.configure(state=tk.NORMAL)
            self.previous.configure(state=tk.NORMAL)

    def navigate(self, direction=DIRECTION_DOWN):
        tree = self._provider.getTree()
        itemSelected = tree.selection()
        if itemSelected:
            item = itemSelected[0]
            if not tree.parent(item):
                item = tree.get_children(item)[0]
            if direction == DIRECTION_DOWN:
                item = tree.next(item)
            else:
                item = tree.prev(item)
            tree.selection_set(item)
            tree.see(item)

    def autoNavigate(self, item=None, direction=DIRECTION_DOWN):

        # Making sure that we do not press the navigate button more than once.
        global afterId
        tree = self._provider.getTree()

        if item is None:
            item = tree.selection()[0]
            if not tree.parent(item):
                item = tree.get_children(item)[0]

        self.startButton.configure(relief=tk.SUNKEN, state=tk.DISABLED)
        self.autoContrast.configure(state=tk.DISABLED)
        self.next.configure(state=tk.DISABLED)
        self.previous.configure(state=tk.DISABLED)
        tree.selection_set(item)
        tree.see(item)
        # Obtain the following element
        if direction == DIRECTION_DOWN:
            nextItem = tree.next(item)
        else:
            nextItem = tree.prev(item)

        if nextItem:
            afterId = tree.after(100, self.autoNavigate, nextItem, direction)
        else:
            # Continue with the following item after a short delay
            children = tree.get_children(item)
            if not children:
                # Start navigating in the other direction
                if direction == DIRECTION_DOWN:
                    afterId = tree.after(100, self.autoNavigate, item, DIRECTION_UP)
                else:
                    afterId = tree.after(100, self.autoNavigate, item, DIRECTION_DOWN)

    def increaseContrast(self):

        # Try to increase the contrast. TODO: Move this to somewhere more resusable
        try:
            data = self._preview.figureimg.get_array()
            imgmean = data.mean()
            imgstd = data.std()
            low = imgmean - CONTRAST_STD * imgstd
            high = imgmean + CONTRAST_STD * imgstd
            data = np.clip(data, low, high)
            # Normalize data to the range [0, 255] if necessary
            data = 255 * (data - low) / (high - low)
            data = np.clip(data, 0, 255).astype(np.uint8)
            self._preview._update(data)

        except Exception as e:
            print(e)

    def itemOnClick(self, e=None):
        x, y, widget = e.x, e.y, e.widget
        elem = widget.identify("element", x, y)

        if not elem:
            return

        tree = self._provider.getTree()
        selectedItem = widget.identify_row(y)
        tree.focus(selectedItem)
        tree.selectChild(selectedItem)
        column = widget.identify_column(x)
        colNumber = int(column.replace('#', '')) - 1
        itemValue = tree.item(selectedItem, "values")[colNumber]

        if "image" in elem:  # click on the checkbox
            tsId = tree.item(selectedItem)['text']
            _, obj = self._provider.getTiltSerie(tsId)
        elif itemValue in [TiltImageStates.CHECK_UNMARK, TiltImageStates.CHECK_MARK]:
            tsId = tree.item(tree.parent(selectedItem))['text']
            tiId = tree.item(selectedItem)['text']
            _, obj = self._provider.getTiltImage(tsId, tiId)
        else:
            return

        self._provider._itemSelected(obj)

    def previewTiltSeries(self, obj, frame):

        preview = self.getPreviewWidget(frame)
        if isinstance(obj, tomo.objects.TiltSeriesBase):
            text = "Tilt Axis angle: %s" % obj.getAcquisition().getTiltAxisAngle()
            obj = self._provider.getTiltSerieRepresentative(obj)
        elif isinstance(obj, tomo.objects.TiltImageBase):
            text = "Tilt image at %sº" % obj.getTiltAngle()

        imageStk = ImageReadersRegistry.open(str(obj.getIndex()) + '@' + obj.getFileName())
        npImage = imageStk.getImage()  # Pass this when scipion-em is released with this functionality: pilImage=True)
        image = self._normalize(npImage)  # Delete this once pilImage is used
        # Get original size
        width, height = image.size

        # Calculate the new dimension keeping the proportion
        newWidth, newHeight = (THUMBNAIL_SIZE, int(THUMBNAIL_SIZE * height / width)) \
            if width > height else (int(THUMBNAIL_SIZE * width / height), THUMBNAIL_SIZE)

        # Resize the image creating a thumbnail
        newImage = image.copy()
        newImage.thumbnail((newWidth, newHeight))
        data = np.array(newImage)
        preview._update(data)
        preview.setLabel(text)
    
    def _normalize(cls, npImage):
        iMax = npImage.max()
        iMin = npImage.min()
        im255 = ((npImage - iMin) / (iMax - iMin) * 255).astype(numpy.uint8)
        return Image.fromarray(im255)


class TomogramsTreeProvider(TreeProvider):
    """ Populate Tree from SetOfTomograms. """

    def __init__(self, tomoList, path, mode):
        TreeProvider.__init__(self)
        self.tomoList = tomoList
        self._path = path
        self._mode = mode

    def getColumns(self):
        return [('Tomogram', 300), ("# coords", 100), ('status', 150)]

    def getObjectInfo(self, tomo):
        if self._mode == 'txt':
            tomogramName = os.path.basename(tomo.getFileName())
            tomogramName = os.path.splitext(tomogramName)[0]
            filePath = os.path.join(self._path, tomogramName + ".txt")
        elif self._mode == 'json':
            tomogramName = os.path.basename(tomo.getFileName())
            tomogramName = os.path.splitext(tomogramName)[0]

            outFile = '*%s_info.json' % pwutils.removeBaseExt(tomogramName.split("__")[0])
            pattern = os.path.join(self._path, outFile)
            files = glob.glob(pattern)

            filePath = ''
            if files:
                filePath = files[0]

        if not os.path.isfile(filePath):
            return {'key': tomogramName, 'parent': None,
                    'text': tomogramName, 'values': (tomo.count, "TODO"),
                    'tags': ("pending")}
        else:
            return {'key': tomogramName, 'parent': None,
                    'text': tomogramName, 'values': (tomo.count, "DONE"),
                    'tags': ("done")}

    def getObjectPreview(self, obj):
        return (None, None)

    def getObjectActions(self, obj):
        return []

    def _getObjectList(self):
        """Retrieve the object list"""
        return self.tomoList

    def getObjects(self):
        objList = self._getObjectList()
        return objList

    def configureTags(self, tree):
        tree.tag_configure("pending", foreground="red")
        tree.tag_configure("done", foreground="green")


class MeshesTreeProvider(TreeProvider):
    """Populate tree from SetOfMeshes"""

    def __init__(self, meshList):
        TreeProvider.__init__(self)
        self._parentDict = {}
        self.meshList = meshList
        self.tomoList = []
        self.tomo_names = set()
        id = 1
        for mesh in self.meshList:
            tomo = mesh.getVolume()
            if tomo.getFileName() not in self.tomo_names:
                tomo.setObjId(id)
                self.tomoList.append(tomo)
                self.tomo_names.add(tomo.getFileName())
                id += 1
        self.tomo_names = sorted(list(self.tomo_names), reverse=False)
        self.tomoList.sort(key=lambda x: x.getFileName(), reverse=False)

    def getColumns(self):
        return [('Tomogram', 300), ('Number of Meshes', 150)]

    def getObjectInfo(self, obj):
        if isinstance(obj, tomo.objects.MeshPoint):
            meshName = 'Mesh %d' % obj.getObjId()
            tomoName = pwutils.removeBaseExt(obj.getVolume().getFileName())
            return {'key': tomoName + '-' + str(obj.getObjId()),
                    'parent': self._parentDict.get(obj.getObjId(), None),
                    'text': meshName, 'values': ('')}
        elif isinstance(obj, tomo.objects.Tomogram):
            tomoName = pwutils.removeBaseExt(obj.getFileName())
            numMeshes = 0
            for mesh in self.meshList:
                if mesh.getVolume().getFileName() == obj.getFileName():
                    numMeshes += 1
            return {'key': tomoName, 'parent': None,
                    'text': tomoName, 'values': (numMeshes)}

    def getObjectActions(self, mesh):
        return []

    def _getObjectList(self):
        """Retrieve the object list"""
        return self.tomoList

    def getObjects(self):
        objList = self._getObjectList()
        self._parentDict = {}
        childs = []
        for obj in self.meshList:
            childs += self._getChilds(obj)
        objList += childs

        return objList

    def _getChilds(self, obj):
        childs = []
        childs.append(obj)
        for idx in range(len(self.tomo_names)):
            if obj.getVolume().getFileName() == self.tomo_names[idx]:
                self._parentDict[obj.getObjId()] = self.tomoList[idx]
        return childs


class TomogramsDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an ImageJ subprocess from a list of Tomograms.
    """

    def __init__(self, parent, viewer, lockGui=False, **kwargs):
        self.path = kwargs.get("path", None)
        self.provider = kwargs.get("provider", None)
        itemDoubleClick = self.doubleClickViewer if viewer else self.doubleClickOnTomogram
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=itemDoubleClick,
                                   allowSelect=False,
                                   lockGui=lockGui,
                                   cancelButton=True,
                                   **kwargs)

    def refresh_gui_viewer(self):
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui_viewer)
        else:
            meshFile = os.path.join(self.path, pwutils.removeBaseExt(self.tomo.getFileName()) + '.txt')
            self.tomo.count = np.loadtxt(meshFile, delimiter=',').shape[0]
            pwutils.cleanPath(self.macroPath)
            self.tree.update()

    def refresh_gui(self):
        self.tree.update()
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui)
        else:
            meshFile = os.path.join(self.path, pwutils.removeBaseExt(self.tomo.getFileName()) + '.txt')
            self.tomo.count = np.loadtxt(meshFile, delimiter=',').shape[0]
            pwutils.cleanPath(self.macroPath)
            self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchIJForTomogram, args=(self.path, self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def doubleClickViewer(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchIJForViewing, args=(self.path, self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui_viewer)

    def lanchIJForTomogram(self, path, tomogram):
        self.macroPath = os.path.join(self.path, "AutoSave_ROI.ijm")
        tomogramFile = tomogram.getFileName()
        tomogramName = os.path.basename(tomogramFile)

        macro = r"""path = "%s";
    file = "%s"
    
    // --------- Initialize Roi Manager ---------
    roiManager("Draw");
    setTool("polygon");
    
    newClass = "Yes";
    outPath = path + file + ".txt";
    
    // --------- Load SetOfMeshes ---------
    if (File.exists(outPath)){
    group = 0;
    groups = loadMeshFile(outPath);
    numMeshes = roiManager("count");
    emptyOutFile(outPath);
    aux = editMeshes(groups, numMeshes, outPath);
    group = group + aux + 1;
    }
    else{
    emptyOutFile(outPath);
    group = 1;
    }
    
    // --------- Draw new Meshes and save them ---------
    while (newClass == "Yes") {
    roiManager("Reset");
    //group = classDialog();
    waitForRoi();
    saveMeshes(group, outPath);
    newClass = newClassDialog();
    group = group + 1;
    }
    
    // --------- Close ImageJ ---------
    run("Quit");
    
    
    // --------- Functions Definition ---------
    function classDialog(){
    Dialog.create("Class Selection");
    Dialog.addMessage("Determine the group of the labels to be drawn");
    Dialog.addNumber("Class Group", 1);
    Dialog.show();
    return floor(Dialog.getNumber());
    }
    
    function newClassDialog(){
    choices = newArray("Yes", "No");
    Dialog.create("Create new class?");
    Dialog.addChoice("Choice", choices);
    Dialog.show();
    return Dialog.getChoice();
    }
    
    function waitForRoi(){
    waitForUser("Draw the desired ROIs\n\nThen click Ok");
    wait(50);
    while(roiManager("count")==0){
    waitForUser("Draw the desired ROIs\n\nThen click Ok");
    wait(50);
    }
    }
    
    function emptyOutFile(outPath){
    fid = File.open(outPath);
    File.close(fid);
    }
    
    function saveMeshes(class, outPath){
    string = "";
    meshes = roiManager("count");
    for (i=0; i<meshes; i++){
    roiManager("select", i);
    Stack.getPosition(channel, slice, frame);
    getSelectionCoordinates(xpoints, ypoints);
    for (j = 0; j < xpoints.length; j++) {
    string = string + "" + xpoints[j] + "," + ypoints[j] + "," + slice + "," + class + "\n";
    }
    }
    lastJump = lastIndexOf(string, "\n");
    File.append(substring(string, 0, lastJump), outPath);
    }
    
    function loadMeshFile(meshPath){
    c = "";
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    contents = split(File.openAsString(meshPath), "\n");
    xpoints = newArray();
    ypoints = newArray();
    groups = newArray();
    for (idx=0; idx < contents.length; idx++){
    values = split(contents[idx], ",");
    if (idx+1 < contents.length){
    valuesNext = split(contents[idx+1], ",");
    }
    else{
    valuesNext =  newArray(-1,-1,-1,-1);
    }
    xpoints = Array.concat(xpoints, values[0]);
    ypoints = Array.concat(ypoints, values[1]);
    if (values[2] != valuesNext[2]){
    xpoints = Array.concat(xpoints, xpoints[0]);
    ypoints = Array.concat(ypoints, ypoints[0]);
    groups = Array.concat(groups, values[3]);
    Stack.setSlice(values[2]);
    makeSelection("polyline", xpoints, ypoints);
    Roi.setName("Class " + values[3]);
    Roi.setStrokeWidth(5);
    roiManager("add");
    count = roiManager("count");
    roiManager("select", count-1);
    roiManager("Set Color", c);
    if (values[3] != valuesNext[3]){
    c = "";
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    }
    xpoints = newArray();
    ypoints = newArray();
    }
    }
    return groups;
    }
    
    function editMeshes(classVect, numMeshes, outPath){
    waitForUser("Edit the input ROIs if needed\n\nThen click Ok");
    string = "";
    for (i=0; i<numMeshes; i++){
    roiManager("select", i);
    Stack.getPosition(channel, slice, frame);
    getSelectionCoordinates(xpoints, ypoints);
    for (j = 0; j < xpoints.length; j++) {
    string = string + "" + xpoints[j] + "," + ypoints[j] + "," + slice + "," + classVect[i] + "\n";
    groupEnd = classVect[i];
    }
    }
    lastJump = lastIndexOf(string, "\n");
    File.append(substring(string, 0, lastJump), outPath);
    
    return groupEnd;
    }
    """ % (os.path.join(path, ''), os.path.splitext(tomogramName)[0])
        macroFid = open(self.macroPath, 'w')
        macroFid.write(macro)
        macroFid.close()

        args = "-i %s -macro %s" % (tomogramFile, self.macroPath)
        viewParams = {showj.ZOOM: 50}
        for key, value in viewParams.items():
            args = "%s --%s %s" % (args, key, value)

        app = "xmipp.ij.commons.XmippImageJ"

        runJavaIJapp(4, app, args).wait()

    def lanchIJForViewing(self, path, tomogram):
        self.macroPath = os.path.join(self.path, "View_ROI.ijm")
        tomogramFile = tomogram.getFileName()
        tomogramName = os.path.basename(tomogramFile)

        meshFile = getMeshVolFileName(self.tomo.getObjId())

        macro = r"""path = "%s";
    file = "%s"
    meshFile = "%s"
    
    outPath = path + meshFile;
    
    // --------- Load SetOfMeshes ---------
    if (File.exists(outPath)){
    loadMeshFile(outPath);
    }
    
    // --------- Functions Definition ---------    
    function loadMeshFile(meshPath){
    c = "";
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    contents = split(File.openAsString(meshPath), "\n");
    xpoints = newArray();
    ypoints = newArray();
    groups = newArray();
    for (idx=0; idx < contents.length; idx++){
    values = split(contents[idx], ",");
    if (idx+1 < contents.length){
    valuesNext = split(contents[idx+1], ",");
    }
    else{
    valuesNext =  newArray(-1,-1,-1,-1);
    }
    xpoints = Array.concat(xpoints, values[0]);
    ypoints = Array.concat(ypoints, values[1]);
    if (values[2] != valuesNext[2]){
    xpoints = Array.concat(xpoints, xpoints[0]);
    ypoints = Array.concat(ypoints, ypoints[0]);
    groups = Array.concat(groups, values[3]);
    Stack.setSlice(values[2]);
    makeSelection("polyline", xpoints, ypoints);
    Roi.setName("Class " + values[3]);
    Roi.setStrokeWidth(5);
    roiManager("add");
    count = roiManager("count");
    roiManager("select", count-1);
    roiManager("Set Color", c);
    if (values[3] != valuesNext[3]){
    c = "";
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    c = c + toHex(255*random);
    }
    xpoints = newArray();
    ypoints = newArray();
    }
    }
    return groups;
    }
    """ % (os.path.join(path, ''), os.path.splitext(tomogramName)[0], meshFile)
        macroFid = open(self.macroPath, 'w')
        macroFid.write(macro)
        macroFid.close()

        args = "-i %s -macro %s" % (tomogramFile, self.macroPath)
        viewParams = {showj.ZOOM: 50}
        for key, value in viewParams.items():
            args = "%s --%s %s" % (args, key, value)

        app = "xmipp.ij.commons.XmippImageJ"

        runJavaIJapp(4, app, args).wait()


class CTFSerieStates:
    UNCHECKED = 'unchecked'
    CHECKED = 'checked'
    ODD = 'odd'
    EVEN = 'even'
    EXCLUDED = "\u2611"  # ☑
    INCLUDED = "\u2610"  # ☐


class CtfEstimationTreeProvider(TreeProvider, ttk.Treeview):
    """ Model class that will retrieve the information from SetOfCTFTomoSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_CTF_SERIE = 'Tilt Series'
    COL_ACQ_ORDER = 'Acq. Order'
    COL_TILT_ANG = 'Tilt Angle'
    COL_CTF_EXCLUDED = 'Excluded'
    COL_CTF_EST_DEFOCUS_U = 'DefocusU (Å)'
    COL_CTF_EST_DEFOCUS_V = 'DefocusV (Å)'
    COL_CTF_EST_AST = 'Astigmatism (Å)'
    COL_CTF_EST_RES = 'Resolution (Å)'
    COL_CTF_EST_FIT = 'CC value'
    COL_CTF_EST_PHASE = 'Phase shift (deg)'

    ORDER_DICT = {COL_CTF_EST_DEFOCUS_U: '_defocusU',
                  COL_CTF_EST_DEFOCUS_V: '_defocusV',
                  COL_CTF_EST_AST: '_defocusRatio',
                  COL_CTF_EST_RES: '_resolution',
                  COL_CTF_EST_PHASE: '_phaseShift',
                  COL_CTF_EST_FIT: '_fitQuality'}

    def __init__(self, master, protocol, outputSetOfCTFTomoSeries, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        self.protocol = protocol
        self.ctfSeries = outputSetOfCTFTomoSeries
        self._hasPhaseShift = outputSetOfCTFTomoSeries.getFirstItem().getFirstItem().hasPhaseShift()
        TreeProvider.__init__(self, sortingColumnName=self.COL_CTF_SERIE)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200
        self._checkedItems = 0

    def getObjects(self):
        objects = []

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for ctfSerie in self.ctfSeries:
            ctfEstObj = ctfSerie.clone()
            ctfEstObj._allowsSelection = True
            ctfEstObj._parentObject = None
            objects.append(ctfEstObj)
            for item in ctfSerie.iterItems(orderBy=orderBy, direction=direction):
                ctfEstItem = item.clone()
                ctfEstItem._allowsSelection = False
                ctfEstItem._parentObject = ctfEstObj
                objects.append(ctfEstItem)

        return objects

    def getCTFSeries(self):
        return self.ctfSeries

    def _sortObjects(self, objects):
        # TODO
        pass

    def objectKey(self, pobj):
        pass

    def getColumns(self):
        cols = [
            (self.COL_CTF_SERIE, 100),
            (self.COL_ACQ_ORDER, 100),
            (self.COL_TILT_ANG, 100),
            (self.COL_CTF_EXCLUDED, 100),
            (self.COL_CTF_EST_DEFOCUS_U, 100),
            (self.COL_CTF_EST_DEFOCUS_V, 100),
            (self.COL_CTF_EST_AST, 150),
            (self.COL_CTF_EST_RES, 100),
            (self.COL_CTF_EST_FIT, 100)
        ]
        if self._hasPhaseShift:
            cols.insert(5, (self.COL_CTF_EST_PHASE, 150))

        return cols

    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        return False

    @staticmethod
    def _getParentObject(pobj, default=None):
        return getattr(pobj, '_parentObject', default)

    def getObjectInfo(self, obj):
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            key = obj.getTsId()
            text = obj.getTsId()
            # TODO: show avg defocus for TomoSeries
            values = ['',
                      #CTFSerieStates.INCLUDED if obj.isEnabled() else CTFSerieStates.EXCLUDED
                      ]
            opened = False
            selected = obj.isEnabled()
        else:  # CTFTomo
            key = "%s.%s" % (obj._parentObject.getTsId(), str(obj.getObjId()))
            text = obj.getObjId()
            acqOrder = obj.getAcquisitionOrder()
            ts = obj._parentObject.getTiltSeries()
            tiltAngle = ts.getItem(CTFTomo.ACQ_ORDER_FIELD, acqOrder).getTiltAngle()
            ast = obj.getDefocusU() - obj.getDefocusV()
            phSh = obj.getPhaseShift() if obj.hasPhaseShift() else 0

            values = [str(acqOrder),
                      str("%0.2f" % tiltAngle),
                      CTFSerieStates.INCLUDED if obj.isEnabled() else CTFSerieStates.EXCLUDED,
                      str("%d" % obj.getDefocusU()),
                      str("%d" % obj.getDefocusV()),
                      str("%d" % ast),
                      str("%0.1f" % obj.getResolution()) if obj.getResolution() else "-",
                      str("%0.3f" % obj.getFitQuality()) if obj.getFitQuality() else "-"]

            if self._hasPhaseShift:
                values.insert(4, str("%0.2f" % phSh))

            opened = False
            selected = False

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': selected,
            'parent': obj._parentObject
        }
        if isinstance(obj, tomo.objects.CTFTomoSeries):
            tags = CTFSerieStates.UNCHECKED
            if not obj.isEnabled():
                obj.setEnabled(True)
                tags = CTFSerieStates.CHECKED
                self._checkedItems += 1
        else:
            tags = CTFSerieStates.INCLUDED
            if not obj.isEnabled():
                tags = CTFSerieStates.EXCLUDED

        item['tags'] = (tags, CTFSerieStates.ODD) if obj.getObjId() % 2 == 0 else (tags, CTFSerieStates.EVEN)

        return item


class CTFEstimationTree(BoundTree):
    def __init__(self, master, provider,  **opts):
        BoundTree.__init__(self, master, provider, **opts)
        self.selectedItem = None
        self.provider = provider
        self._checkedItems = provider._checkedItems

    def toogleExclusion(self, item, isCTFTomoSerie):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        tags = CTFSerieStates.ODD if CTFSerieStates.ODD in self.item(item, 'tags') else CTFSerieStates.EVEN
        self.selectedItem = item
        if isCTFTomoSerie:

            if CTFSerieStates.UNCHECKED in self.item(item, 'tags'):
                self.item(item, tags=(CTFSerieStates.CHECKED, tags,))
                self._checkedItems += 1
                self.getSelectedObj().setEnabled(False)
                self.item(item)['selected'] = True
                excluded = True
            else:
                self.item(item, tags=(CTFSerieStates.UNCHECKED, tags,))
                self.getSelectedObj().setEnabled(True)
                self._checkedItems -= 1
                self.item(item)['selected'] = False
                excluded = False

            ctfTreeChildrenIds = list(self.get_children(item))
            for ctfChildrenId in ctfTreeChildrenIds:
                tags = CTFSerieStates.ODD if CTFSerieStates.ODD in self.item(ctfChildrenId, 'tags') else CTFSerieStates.EVEN
                self.toogleCtf(ctfChildrenId, tags, excluded)
        else:
            self.toogleCtf(item, tags, CTFSerieStates.INCLUDED in self.item(item, 'tags'))

    def toogleCtf(self, item, tags, excluded):
        itemValues = self.item(item)
        values = itemValues['values']
        if excluded:
            values[2] = CTFSerieStates.EXCLUDED
            self.item(item, tags=(CTFSerieStates.EXCLUDED, tags,))
            self.item(item, values=values)
            self._objDict[item].setEnabled(False)
            self.item(item)['selected'] = True
        else:
            values[2] = CTFSerieStates.INCLUDED
            self.item(item, tags=(CTFSerieStates.INCLUDED, tags,))
            self.item(item, values=values)
            self._objDict[item].setEnabled(True)
            self.item(item)['selected'] = False

    def _onItemClick(self, event=None):
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        self.selectedItem = self.identify_row(y)
        if "image" in elem:  # click on the CTFTomoSerie checkbox
            self.toogleExclusion(self.selectedItem, True)
        elif x > 308 and x < 317:  # Position of exclude/include checkbox(CTFTomo)
            self.toogleExclusion(self.selectedItem, False)

    def getSelectedItem(self):
        return self.selectedItem

    def getSelectedObj(self):
        obj = None
        if self.selectedItem:
            obj = self._objDict[self.getFirst()]
        return obj


class CtfEstimationListDialog(ListDialog):
    def __init__(self, parent, title, provider, protocol, inputTS, lockGui=False, **kwargs):
        self._project = protocol.getProject()
        self._protocol = protocol
        self._inputSetOfTiltSeries = inputTS
        self._checkedItems = provider._checkedItems
        # the funcs below should be implemented in viewers
        self._show1DPLot = kwargs.pop('plot1Dfunc', None)
        self._show2DPLot = kwargs.pop('plot2Dfunc', None)
        self._showExtraPlot = kwargs.pop('plotExtrafunc', None)
        ListDialog.__init__(self, parent, title, provider,
                            allowSelect=False,
                            cancelButton=True,
                            lockGui=lockGui,
                            **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config()
        self._col = 1
        self._fillCTFEstimationGUI(bodyFrame)

    def _addButton(self, frame, text, image, command, sticky='news', state=tk.NORMAL):
        btn = tk.Label(frame, text=text, image=self.getImage(image),
                       compound=tk.LEFT, cursor='hand2', state=state)
        btn.bind('<Button-1>', command)
        btn.grid(row=0, column=self._col, sticky=sticky,
                 padx=(0, 5), pady=5)
        self._col += 1
        return btn

    def _fillCTFEstimationGUI(self, bodyFrame):
        # Create a top panel to put the search box and buttons
        topPanel = tk.Frame(bodyFrame)
        topPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._createTopPanel(topPanel)

        # Create a bottom panel to put the tree and the plotter
        bottomPanel = tk.Frame(bodyFrame)
        bottomPanel.grid(row=1, column=0, padx=0, pady=0, sticky='news')
        self._createBottomPanel(bottomPanel)

    def _createTopPanel(self, topPanel):
        self._createFilterBox(topPanel)
        topRigthPanel = tk.Frame(topPanel)
        topRigthPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')
        self._createSubsetButton(topRigthPanel)
        if self._show1DPLot is not None:
            self._createShowFit(topRigthPanel)
        if self._showExtraPlot is not None:
            self._createShowExtra(topRigthPanel)
        self._createViewerHelp(topRigthPanel)

    def _createSubsetButton(self, topRigthPanel):
        state = tk.NORMAL
        if self._checkedItems or self._checkedItems == len(self.provider.getCTFSeries()):
            state = tk.DISABLED
        self.generateSubsetButton = self._addButton(topRigthPanel,
                                                    'Generate subsets',
                                                    pwutils.Icon.PROCESSING,
                                                    self._actionCreateSets,
                                                    sticky='ne',
                                                    state=state)

    def _createShowFit(self, topRigthPanel):
        self.createShowFitButton = self._addButton(topRigthPanel, '1D fit',
                                                   pwutils.Icon.ACTION_RESULTS, self._show1DFit,
                                                   state=tk.DISABLED,
                                                   sticky='ne')

    def _createViewerHelp(self, topRigthPanel):
        self._addButton(topRigthPanel, pwutils.Message.LABEL_HELP,
                        pwutils.Icon.ACTION_HELP, self._showHelp, sticky='ne')

    def _show1DFit(self, event=None):
        itemSelected = self.tree.getSelectedItem()
        obj = self.tree.getSelectedObj()
        if self.tree.parent(itemSelected):  # child item
            if obj is not None:
                for ctfSerie in self.provider.getCTFSeries():
                    if ctfSerie.getTsId() in itemSelected:
                        # TODO: sort ctfSerie by id
                        ctfId = int(itemSelected.split('.')[-1])
                        plot = self._show1DPLot(ctfSerie, ctfId)
                        plot.show()
                        break

    def _createShowExtra(self, topRigthPanel):
        self.createShowFitButton = self._addButton(topRigthPanel, 'Extra plots',
                                                   pwutils.Icon.ACTION_RESULTS, self._showExtra,
                                                   state=tk.DISABLED,
                                                   sticky='ne')

    def _showExtra(self, event=None):
        itemSelected = self.tree.getSelectedItem()
        obj = self.tree.getSelectedObj()
        if self.tree.parent(itemSelected):  # child item
            if obj is not None:
                for ctfSerie in self.provider.getCTFSeries():
                    if ctfSerie.getTsId() in itemSelected:
                        # TODO: sort ctfSerie by id
                        ctfId = int(itemSelected.split('.')[-1])
                        plot = self._showExtraPlot(ctfSerie, ctfId)
                        plot.show()
                        break

    def _actionCreateSets(self, event=None):
        if self.generateSubsetButton['state'] == tk.NORMAL:
            self.info('Processing...')
            protocol = self.provider.protocol
            ctfSeries = self.provider.getCTFSeries()
            suffix = self._getSuffix(protocol)
            goodCTFName = 'goodCtf%s' % suffix

            outputSetOfgoodCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                                prefix=goodCTFName,
                                                                copyInfo=True)
            outputSetOfbadCTFTomoSeries = String()
            outputSetOfbadCTFTomoSeries.set('')
            ctfSerieSeek = 0
            for ctfSerie in ctfSeries.iterItems():
                ctfSerieClon = ctfSerie.clone()
                ctfSerieClon.setEnabled(True)
                goodCTF = CTFSerieStates.UNCHECKED in self.tree.item(ctfSerie.getTsId(), 'tags')
                ctfSeek = 1
                if goodCTF:
                    # Adding the ctfSerie to the good set of ctfTomoSeries
                    outputSetOfgoodCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfgoodCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)
                    objList = list(self.tree._objDict)
                    for item in ctfSerie.iterItems():
                        ctfEstItem = item.clone()
                        obj = self.tree.getObjectFromId(objList[ctfSerieSeek + ctfSeek])
                        ctfEstItem.setEnabled(obj.isEnabled())
                        ctfSerieClon.append(ctfEstItem)
                        ctfSeek += 1

                    ctfSerieClon.write()
                    outputSetOfgoodCTFTomoSeries.update(ctfSerieClon)
                    outputSetOfgoodCTFTomoSeries.write()
                else:
                    # Adding the ctfSerie to the bad set of ctfTomoSeries
                    outputSetOfbadCTFTomoSeries.set(outputSetOfbadCTFTomoSeries.get() + ctfSerie.getTsId() + ' ')
                    ctfSeek += ctfSerie.getSize()

                ctfSerieSeek = ctfSeek

            outputgoodCTFSetName = 'goodSetOfCTFTomoSeries%s' % suffix
            outputbadCTFSetName = 'badSetOfCTFTomoSeries%s' % suffix

            if len(outputSetOfgoodCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputgoodCTFSetName: outputSetOfgoodCTFTomoSeries})

            if not outputSetOfbadCTFTomoSeries.empty() > 0:
                outputSetOfbadCTFTomoSeries.set(", ".join(outputSetOfbadCTFTomoSeries.get().split()))
                protocol._defineOutputs(**{outputbadCTFSetName: outputSetOfbadCTFTomoSeries})

            protocol._store()
            self.info('The output has been created successfully')

    def _showHelp(self, event=None):
        showInfo('CTFTomoSeries viewer help',
                 'This viewer allows you to create two '
                 'subsets of CTFTomoSeries which are called good '
                 'and bad respectively.\n\n'
                 'Note: The items that are excluded(checked) are the ones that '
                 'represent the bad CTFTomoSeries', self)

    def _getSuffix(self, protocol):
        """
        Return the index of the last output in order to complete the new
        output with a suffix
        """
        maxCounter = -1
        patterns = ['goodSetOfCTFTomoSeries', 'badSetOfCTFTomoSeries']
        for attrName, _ in protocol.iterOutputAttributes():
            for pattern in patterns:
                suffix = attrName.replace(pattern, '')
                try:
                    counter = int(suffix)
                except:
                    counter = 1  # when there is no number, assume 1
                maxCounter = max(counter, maxCounter)

        return str(maxCounter + 1) if maxCounter > 0 else ''

    def _createBottomPanel(self, bottomPanel):
        self._createCTFEstimationGUI(bottomPanel)
        self.initial_focus = self.tree

    def _createCTFEstimationGUI(self, bottomPanel):
        # Create a division Paned
        pw = tk.PanedWindow(bottomPanel, orient=tk.HORIZONTAL)
        # Create a left panel to put the tree
        bottomleftPanel = tk.Frame(pw)
        bottomleftPanel.grid(row=0, column=0, padx=0, pady=0, sticky='news')
        self._createTree(bottomleftPanel)
        pw.add(bottomleftPanel)
        # Create a right panel to put the plotter
        self.bottomRightPanel = ttk.Frame(pw)
        self.bottomRightPanel.grid(row=0, column=1, padx=0, pady=0, sticky='news')
        self._createPlotter()
        pw.add(self.bottomRightPanel)
        pw.pack(fill=BOTH, expand=True)
        # This method is used to show sash
        pw.configure(sashrelief=RAISED)

    def _createTree(self, parent):
        gui.configureWeigths(parent)
        self.tree = CTFEstimationTree(parent, self.provider,
                                      selectmode=self._selectmode)
        item = self.tree.identify_row(0)
        self.tree.selection_set(item)
        self.tree.focus_force()
        self.tree.selectedItem = item
        self.im_checked = gui.getImage(Icon.CHECKED)
        self.im_unchecked = gui.getImage(Icon.UNCHECKED)
        self.tree.tag_configure(CTFSerieStates.UNCHECKED,
                                image=self.im_unchecked)
        self.tree.tag_configure(CTFSerieStates.CHECKED,
                                image=self.im_checked)
        self.tree.tag_configure(CTFSerieStates.EVEN, background='#F2F2F2')
        self.tree.tag_configure(CTFSerieStates.ODD, background='#E6E6E6')
        self.tree.tag_configure(CTFSerieStates.CHECKED, foreground='red')
        self.tree.tag_configure(CTFSerieStates.UNCHECKED, foreground='black')
        self.tree.tag_configure(CTFSerieStates.INCLUDED, background='#F2F2F2', foreground='black')
        self.tree.tag_configure(CTFSerieStates.EXCLUDED, background='#E6E6E6', foreground='red')
        self.tree.bind("<ButtonRelease-1>", self._onItemClick)
        self.tree.bind("<Up>", self._onUp)
        self.tree.bind("<Down>", self._onDown)
        self.tree.bind("<space>", self._onSpace)

    def getTree(self):
        return self.tree

    def _onSpace(self, event):
        selectedItem = self.tree.selectedItem
        if selectedItem:
            obj = self.tree._objDict[selectedItem]
            isCTFTomoSerie = True if isinstance(obj, tomo.objects.CTFTomoSeries) else False
            self.tree.toogleExclusion(selectedItem, isCTFTomoSerie)

    def _onDown(self, event):
        selectedItem = self.tree.selectedItem
        if selectedItem:
            nextItem = self.tree.next(selectedItem)
            if nextItem:
                self.tree.selectedItem = nextItem
                self.tree.selectChild(nextItem)
                self.tree.selection_set(nextItem)
                self._createPlotter()

    def _onUp(self, event):
        selectedItem = self.tree.selectedItem
        if selectedItem:
            prevItem = self.tree.prev(selectedItem)
            if prevItem:
                self.tree.selectedItem = prevItem
                self.tree.selectChild(prevItem)
                self.tree.selection_set(prevItem)
                self._createPlotter()

    def _onItemClick(self, event):
        self.tree._onItemClick(event)
        self._createPlotter()

    def plotterChildItem(self, itemSelected):
        plotterPanel = tk.Frame(self.bottomRightPanel)
        if self._show2DPLot is None:
            self.plotterParentItem(self.tree.parent(itemSelected),
                                   int(itemSelected.split('.')[-1]))
        else:
            for ctfSerie in self.provider.getCTFSeries():
                if ctfSerie.getTsId() in itemSelected:
                    ctfId = int(itemSelected.split('.')[-1])
                    # TODO: sort ctfSerie by id
                    fig = self._show2DPLot(ctfSerie, ctfId)
                    if fig is None:
                        return
                    canvas = FigureCanvasTkAgg(fig, master=plotterPanel)
                    canvas.draw()
                    canvas.get_tk_widget().pack(fill=BOTH, expand=0)
                    canvas.get_tk_widget().bind("<Up>", self._onKeyPress)
                    canvas.get_tk_widget().bind("<Down>", self._onKeyPress)
                    canvas.get_tk_widget().bind("<space>", self._onKeyPress)
                    break
        plotterPanel.grid(row=0, column=1, sticky='news')

    def plotterParentItem(self, itemSelected, tiltSelected=None):
        plotterPanel = tk.Frame(self.bottomRightPanel)
        angDict = {}
        defocusUList = []
        defocusVList = []
        phShList = []
        resList = []

        for ts in self._inputSetOfTiltSeries:
            if ts.getTsId() == itemSelected:
                for item in ts:
                    # Related to excluded views:
                    # We will initially save all the angles
                    # and then, only exclude taking into account the defocus values.
                    angDict[item.getAcquisitionOrder()] = item.getTiltAngle()
                break

        for ctfSerie in self.provider.getCTFSeries():
            if ctfSerie.getTsId() == itemSelected:
                angList = []
                for item in ctfSerie.iterItems(orderBy='id'):
                    defocusU = item.getDefocusU()
                    # Related to excluded views:
                    # Add the angle corresponding to the defocus value
                    if item.getAcquisitionOrder() in angDict:
                        angList.append(angDict[item.getAcquisitionOrder()])
                        #   pwem method setWrongDefocus assigns:
                        #   ctfModel.setDefocusU(-999)
                        #   ctfModel.setDefocusV(-1)
                        #   ctfModel.setDefocusAngle(-999)
                        # If it's the case, the corresponding point won't be added to be plotted as
                        # it will widen the representation range, what would make the represented region
                        # of interest smaller
                        defocusUList.append(defocusU)
                        defocusVList.append(item.getDefocusV())
                        phShList.append(
                            item.getPhaseShift() if item.hasPhaseShift() else 0)
                        resList.append(item.getResolution())

                fig = Figure(figsize=(7, 7), dpi=100)
                defocusPlot = fig.add_subplot(111)
                defocusPlot.grid()
                defocusPlot.set_title(itemSelected)
                defocusPlot.set_xlabel('Tilt angle (deg)')
                defocusPlot.set_ylabel('Defocus (Å)')
                defocusPlot.plot(angList, defocusUList, marker='.',
                                 color='tab:red', label='DefocusU (Å)')
                defocusPlot.plot(angList, defocusVList, marker='.',
                                 color='tab:blue', label='DefocusV (Å)')

                if tiltSelected is not None:
                    tomoCtf = ctfSerie[tiltSelected]
                    angle = angDict[tomoCtf.getAcquisitionOrder()]
                    defocusU = tomoCtf.getDefocusU()
                    defocusV = tomoCtf.getDefocusV()
                    defocusPlot.scatter(angle, defocusU, marker='o', facecolors='none',
                                     edgecolors='tab:red', s=150)
                    defocusPlot.scatter(angle, defocusV, marker='o', facecolors='none',
                                     edgecolors='tab:blue', s=150)

                if item.hasPhaseShift():
                    phShPlot = defocusPlot.twinx()
                    phShPlot.set_ylim(0, 180)
                    phShPlot.set_ylabel('Phase shift', color='tab:green')
                    phShPlot.plot(angList, phShList, marker='.',
                                  color='tab:green', label='Phase shift (deg)')
                else:  # no phase shift, plot resolution instead
                    resPlot = defocusPlot.twinx()
                    resPlot.set_ylim(0, 30)
                    resPlot.set_ylabel('Resolution (A)', color='tab:green')
                    resPlot.plot(angList, resList, marker='.',
                                 color='tab:green', label='Resolution (Å)')

                fig.legend()
                canvas = FigureCanvasTkAgg(fig, master=plotterPanel)
                canvas.draw()
                canvas.get_tk_widget().pack(fill=BOTH, expand=0)
                plotterPanel.grid(row=0, column=1, sticky='news')
                canvas.get_tk_widget().bind("<Up>", self._onKeyPress)
                canvas.get_tk_widget().bind("<Down>", self._onKeyPress)
                canvas.get_tk_widget().bind("<space>", self._onKeyPress)
                break

    def _onKeyPress(self, event):
        if event.keysym == 'Up':
            self._onUp(event)
        elif event.keysym == 'Down':
            self._onDown(event)
        elif event.keysym == 'space':
            self._onSpace(event)

    def _createPlotter(self):
        itemSelected = self.tree.getSelectedItem()
        obj = self.tree.getSelectedObj()
        self._checkedItems = self.tree._checkedItems
        self.generateSubsetButton['state'] = tk.NORMAL

        if obj is not None:
            if self.tree.parent(itemSelected):  # child item
                if self._show1DPLot is not None:
                    self.createShowFitButton['state'] = tk.NORMAL
                if self._showExtraPlot is not None:
                    self.createShowFitButton['state'] = tk.NORMAL
                self.plotterChildItem(itemSelected)

            else:  # parent item
                if self._show1DPLot is not None:
                    self.createShowFitButton['state'] = tk.DISABLED
                if self._showExtraPlot is not None:
                    self.createShowFitButton['state'] = tk.DISABLED

                self.plotterParentItem(itemSelected)
