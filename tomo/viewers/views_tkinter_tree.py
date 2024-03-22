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

import glob
import threading
from tkinter import messagebox, BOTH, RAISED

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from pwem.viewers import showj
from pwem.viewers.showj import runJavaIJapp
from pyworkflow.gui import *
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog, ToolbarListDialog, showInfo
import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain

import tomo.objects
from ..convert.convert import getMeshVolFileName

# How many standard deviations to truncate above and below the mean when increasing contrast:
CONTRAST_STD = 0.5

class TiltImageStates:
    EXCLUDED = 'excluded'
    INCLUDED = 'included'
    ODD = 'odd'
    EVEN = 'even'


class TiltSeriesTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_TS = 'Tilt series'
    COL_TI = 'Path'
    COL_TI_ANGLE = 'Tilt angle'
    COL_TI_ENABLED = 'Excluded'
    COL_TI_ACQ_ORDER = 'Order'
    COL_TI_DOSE = "Accum. dose"
    COL_TI_TRANSFORM = "T. Matrix"
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
        self.tree.tag_configure(TiltImageStates.EXCLUDED,  font=standardFont, foreground='red')
        self.tree.tag_configure(TiltImageStates.INCLUDED,  font=standardFont, foreground='black')
        self.tree.tag_configure(TiltImageStates.EVEN, background='#F2F2F2', foreground='black')
        self.tree.tag_configure(TiltImageStates.ODD, background='#E6E6E6', foreground='black')

    def _itemSelected(self, obj):
        excluded = 'True' if obj.isEnabled() else 'False'
        obj.setEnabled(not obj.isEnabled())
        _, y, _, _ = self.tree.bbox(self.tree.selection()[0])
        selectedItem = self.tree.identify_row(y+1)
        itemValues = self.tree.item(selectedItem, 'values')
        newValues = (itemValues[0], itemValues[1], excluded, itemValues[3], itemValues[4], itemValues[5])
        tags = TiltImageStates.EVEN
        if TiltImageStates.ODD in self.tree.item(selectedItem, 'tags'):
            tags = TiltImageStates.ODD

        if excluded == 'True':
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
            (self.COL_TS, 100),
            (self.COL_TI_ACQ_ORDER, 100),
            (self.COL_TI_ANGLE, 100),
            (self.COL_TI_ENABLED, 100),
            (self.COL_TI_DOSE, 100),
            (self.COL_TI, 400),
        ]
        if not isinstance(self.tiltSeries, tomo.objects.SetOfTiltSeriesM):
            cols.append((self.COL_TI_TRANSFORM, 300))

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
        tags = None

        if isinstance(obj, tomo.objects.TiltSeriesBase):
            key = objId
            text = tsId
            values = [str(obj)]
            opened = False
        else:  # TiltImageBase
            key = '%s.%s' % (tsId, objId)
            text = objId

            dose = obj.getAcquisition().getAccumDose() if hasattr(obj.getAcquisition(), '_accumDose') else None
            adqOrder = obj.getAcquisitionOrder() if hasattr(obj, '_acqOrder') else None

            values = [str("%d" % adqOrder) if adqOrder is not None else "",
                      str("%0.2f" % obj.getTiltAngle()),
                      str(not obj.isEnabled()),
                      round(dose, 2) if dose is not None else "",
                      "%d@%s" % (obj.getIndex() or 1, obj.getFileName()),
                      ]

            if not isinstance(obj, tomo.objects.TiltImageM):
                matrix = "" if not obj.hasTransform() else obj.getTransform().getMatrixAsList()
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

        if tags is not None:
            if obj.getObjId() % 2 == 0:
                item['tags'] = (tags, TiltImageStates.ODD)
            else:
                item['tags'] = (tags, TiltImageStates.EVEN)

        return item

    def getObjectActions(self, obj):
        actions = []

        if isinstance(obj, tomo.objects.TiltSeries):
            viewers = Domain.findViewers(obj.getClassName(),
                                         pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                def createViewer(viewerClass, obj):
                    proj = self.protocol.getProject()
                    item = self.tiltSeries[obj.getObjId()]  # to load mapper
                    return lambda: viewerClass(project=proj, protocol=self.protocol).visualize(item)
                actions.append((viewerClass.getName(), createViewer(viewerClass, obj)))
        else:
            actionName = "Exclude" if obj.isEnabled() else "Include"
            actions.append((actionName, lambda: self._itemSelected(obj)))

        return actions


class TiltSeriesDialog(ToolbarListDialog):
    def __init__(self, parent, title, provider, tiltSeries, protocol, **kwargs):

        self._tiltSeries = tiltSeries
        self._protocol = protocol
        self._provider = provider
        self._applyContrastCallback = kwargs.get('applyContrastCallback', None)


        toolbarButtons = [
            dialog.ToolbarButton('Toggle exclusion', self._toggleExclusion, Icon.ACTION_CLOSE,
                                 tooltip="Exclude or include the selected tiltimage", shortcut='<space>'),
            dialog.ToolbarButton('Increase contrast', self._applyContrastCallback, Icon.ACTION_CONTRAST,
                                 tooltip="Apply contrast to the selected tiltimage", shortcut='<Control-c>'),
            dialog.ToolbarButton('Save', self._saveExcluded, Icon.ACTION_SAVE,
                                 tooltip="Create a new output with excluded views marked", shortcut='<Control-s>'),
            dialog.ToolbarButton('|', None)
        ]

        if isinstance(self._tiltSeries, tomo.objects.SetOfTiltSeries):
            viewers = Domain.findViewers(tomo.objects.TiltSeries.getClassName(),
                                         pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                def launchViewer():
                    proj = self._protocol.getProject()
                    viewerInstance = viewerClass(project=proj, protocol=self._protocol)
                    return lambda event: self.launchViewer(viewerInstance)
                toolbarButtons.append(dialog.ToolbarButton(viewerClass.getName(), launchViewer(), Icon.ACTION_RESULTS))

        toolbarButtons.append(dialog.ToolbarButton('|', None))
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
                 '1. Toggle exclusion button to exclude or include a selected tiltimage or use right-click.\n'
                 '2. Increase contrast button to enhance the tiltimage contrast.\n'
                 '3. Save button to create a new set with excluded views marked.', self)

    def _toggleExclusion(self, event=None):
        itemSelected = self.tree.selection()
        obj = self.tree._objects[self.tree.index(itemSelected) + 1]
        self._provider._itemSelected(obj)

    def _saveExcluded(self, event=None):
        updatedCount = self._provider.getUpdatedCount()
        changes = self._provider.getchanges()
        msg = "Are you sure you want to create a new set of TiltSeries?" if updatedCount and changes \
            else "This set of TiltSeries has already been saved previously or is the original one. Are you still sure you want to generate a new set?"
        result = messagebox.askquestion("Confirmation", msg, icon='info', **{'parent': self})
        if result == messagebox.YES:
            outputSetOfTiltSeries = tomo.objects.SetOfTiltSeries.create(self._protocol.getPath(),
                                                                        suffix=str(self._protocol.getOutputsSize()))
            outputSetOfTiltSeries.copyInfo(self._tiltSeries)
            outputSetOfTiltSeries.setDim(self._tiltSeries.getDim())
            excludedViews = self._provider.getExcludedViews()
            for ts in self._tiltSeries:
                newTs = ts.clone()
                newTs.copyInfo(ts)
                outputSetOfTiltSeries.append(newTs)
                for ti in ts.iterItems():
                    newTi = ti.clone()
                    newTi.copyInfo(ti, copyId=True)
                    newTi.setAcquisition(ti.getAcquisition())
                    newTi.setLocation(ti.getLocation())
                    # For some reason .clone() does not clone the enabled nor the creation time
                    included = False if ti.getObjId() in excludedViews[ts.getTsId()] else True
                    newTi.setEnabled(included)
                    newTs.append(newTi)

                newTs.setDim(ts.getDim())
                newTs.write()

                outputSetOfTiltSeries.update(newTs)
                outputSetOfTiltSeries.write()
            outputName = 'TiltSeries_' + str(self._protocol.getOutputsSize()+1)
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
        TiltSeriesDialog(self._tkParent, 'Tilt series viewer', self._provider, self._tiltSeries, self._protocol,
                         lockGui=True, previewCallback=self.previewTiltSeries, allowSelect=False, cancelButton=True,
                         applyContrastCallback=self.increaseContrast)

    def getPreviewWidget(self, frame):

        if not self._preview:
            from pyworkflow.gui.matplotlib_image import ImagePreview
            self._preview = ImagePreview(frame, dim=500, label="Tilt series", listenersDict={"<Button-1>": self.increaseContrast})
            self._preview.grid(row=0, column=0)
        return self._preview

    def increaseContrast(self, event):

        # Try to increase the contrast. TODO: Move this to somewhere more resusable
        try:
            
            data = self._preview.figureimg.get_array()
            imgmean = data.mean()
            imgstd = data.std()
            low = imgmean - CONTRAST_STD * imgstd
            high = imgmean + CONTRAST_STD * imgstd
            data[data < low] = low
            data[data > high] = high

            self._preview._update(data)
        except Exception as e:
            print(e)

    def previewTiltSeries(self, obj, frame):
        preview = self.getPreviewWidget(frame)
        if isinstance(obj, tomo.objects.TiltSeries):
            text = "Tilt Axis angle: %s" % obj.getAcquisition().getTiltAxisAngle()
            obj = self._provider.getTiltSerieRepresentative(obj)
        else:
            text = "Tilt image at %sº" % obj.getTiltAngle()

        image = obj.getImage()
        data = image.getData()
        preview._update(data)
        preview.setLabel(text)

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
    FAILED = 'Failed'
    OK = 'Ok'


class CtfEstimationTreeProvider(TreeProvider, ttk.Treeview):
    """ Model class that will retrieve the information from SetOfCTFTomoSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_CTF_SERIE = 'Tilt Series'
    COL_TILT_ANG = 'Tilt Angle'
    COL_CTF_EST_DEFOCUS_U = 'DefocusU (A)'
    COL_CTF_EST_DEFOCUS_V = 'DefocusV (A)'
    COL_CTF_EST_AST = 'Astigmatism (A)'
    COL_CTF_EST_RES = 'Resolution (A)'
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
            (self.COL_TILT_ANG, 100),
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
                      #CTFSerieStates.OK if obj.isEnabled() else CTFSerieStates.FAILED
                      ]
            opened = False
            selected = obj.isEnabled()
        else:  # CTFTomo
            key = "%s.%s" % (obj._parentObject.getTsId(), str(obj.getObjId()))
            text = obj.getIndex()
            ts = obj._parentObject.getTiltSeries()
            tiltAngle = ts[int(text)].getTiltAngle()
            ast = obj.getDefocusU() - obj.getDefocusV()
            phSh = obj.getPhaseShift() if obj.hasPhaseShift() else 0

            values = [str("%0.2f" % tiltAngle),
                      #CTFSerieStates.OK if obj.isEnabled() else CTFSerieStates.FAILED,
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

            if obj.getObjId() % 2 == 0:
                item['tags'] = (tags, CTFSerieStates.ODD)
            else:
                item['tags'] = (tags,  CTFSerieStates.EVEN)
        else:
            if obj.getObjId() % 2 == 0:
                item['tags'] = (CTFSerieStates.ODD,)
            else:
                item['tags'] = (CTFSerieStates.EVEN,)
        return item


class CTFEstimationTree(BoundTree):
    def __init__(self, master, provider,  **opts):
        BoundTree.__init__(self, master, provider, **opts)
        self.selectedItem = None
        self._checkedItems = provider._checkedItems

    def check_item(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        tags = CTFSerieStates.EVEN
        if CTFSerieStates.ODD in self.item(item, 'tags'):
            tags = CTFSerieStates.ODD

        if CTFSerieStates.UNCHECKED in self.item(item, 'tags'):
            self.item(item, tags=(CTFSerieStates.CHECKED, tags,))
            self._checkedItems += 1
            self.getSelectedObj().setEnabled(False)
            self.item(item)['selected'] = True
        else:
            self.item(item, tags=(CTFSerieStates.UNCHECKED, tags,))
            self.getSelectedObj().setEnabled(True)
            self._checkedItems -= 1
            self.item(item)['selected'] = False

    def _onClick(self, event=None):
        self._unpostMenu()
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        self.selectedItem = self.identify_row(y)
        self.focus(self.selectedItem)
        if "image" in elem:  # click on the checkbox
            self.check_item(self.selectedItem)

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
            protocol = self.provider.protocol
            ctfSeries = self.provider.getCTFSeries()
            suffix = self._getSuffix(protocol)
            goodCTFName = 'goodCtf%s' % suffix
            badCTFName = 'badCtf%s' % suffix

            outputSetOfgoodCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                                prefix=goodCTFName,
                                                                copyInfo=True)
            outputSetOfbadCTFTomoSeries = ctfSeries.createCopy(protocol._getPath(),
                                                               prefix=badCTFName,
                                                               copyInfo=True)
            for ctfSerie in ctfSeries:
                ctfSerieClon = ctfSerie.clone()
                if CTFSerieStates.UNCHECKED in self.tree.item(ctfSerie.getTsId(),
                                                              'tags'):
                    # Adding the ctfSerie to the good set of ctfTomoSeries
                    outputSetOfgoodCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfgoodCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)

                else:
                    # Adding the ctfSerie to the bad set of ctfTomoSeries
                    outputSetOfbadCTFTomoSeries.append(ctfSerieClon)
                    outputSetOfbadCTFTomoSeries.setSetOfTiltSeries(self._inputSetOfTiltSeries)

                for item in ctfSerie.iterItems():
                    ctfEstItem = item.clone()
                    ctfSerieClon.append(ctfEstItem)

            outputgoodCTFSetName = 'goodSetOfCTFTomoSeries%s' % suffix
            outputbadCTFSetName = 'badSetOfCTFTomoSeries%s' % suffix

            if len(outputSetOfgoodCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputgoodCTFSetName: outputSetOfgoodCTFTomoSeries})

            if len(outputSetOfbadCTFTomoSeries) > 0:
                protocol._defineOutputs(**{outputbadCTFSetName: outputSetOfbadCTFTomoSeries})

            protocol._store()
            self.cancel()

    def _showHelp(self, event=None):
        showInfo('CTFTomoSeries viewer help',
                 'This viewer allows you to create two '
                 'subsets of CTFTomoSeries which are called good '
                 'and bad respectively.\n\n'
                 'Note: The series that are checked are the ones that '
                 'represent the bad CTFTomoSeries', self.parent)

    def _getSuffix(self, protocol):
        """
        Return the number of the last output in order to complete the new
        output with a suffix
        """
        maxCounter = -1
        pattern = 'goodSetOfCTFTomoSeries'
        for attrName, _ in protocol.iterOutputAttributes():
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
        self._createPloter(self.bottomRightPanel)
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
        self.tree.focus(item)
        self.tree.selectedItem = item
        self.im_checked = gui.getImage(Icon.CHECKED)
        self.im_unchecked = gui.getImage(Icon.UNCHECKED)
        self.tree.tag_configure(CTFSerieStates.UNCHECKED,
                                image=self.im_unchecked)
        self.tree.tag_configure(CTFSerieStates.CHECKED,
                                image=self.im_checked)
        self.tree.tag_configure(CTFSerieStates.EVEN, background='#F2F2F2',
                                foreground='black')
        self.tree.tag_configure(CTFSerieStates.ODD, background='#E6E6E6',
                                foreground='black')
        self.tree.bind("<Button-1>", self._createPloter, True)

    def plotterChildItem(self, itemSelected):
        plotterPanel = tk.Frame(self.bottomRightPanel)
        if self._show2DPLot is None:
            self.plotterParentItem(self.tree.parent(itemSelected))
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
                    break
        plotterPanel.grid(row=0, column=1, sticky='news')

    def plotterParentItem(self, itemSelected):
        plotterPanel = tk.Frame(self.bottomRightPanel)
        angList = []
        defocusUList = []
        defocusVList = []
        phShList = []
        resList = []

        for ts in self._inputSetOfTiltSeries:
            if ts.getTsId() == itemSelected:
                for item in ts:
                    # Related to excluded views:
                    # Only represent the enabled tilt images
                    if item.isEnabled():
                        angList.append(item.getTiltAngle())
                break

        for ctfSerie in self.provider.getCTFSeries():
            if ctfSerie.getTsId() == itemSelected:
                for item in ctfSerie.iterItems(orderBy='id'):
                    defocusU = item.getDefocusU()
                    # Related to excluded views:
                    #   pwem method setWrongDefocus assigns:
                    #   ctfModel.setDefocusU(-999)
                    #   ctfModel.setDefocusV(-1)
                    #   ctfModel.setDefocusAngle(-999)
                    # If it's the case, the corresponding point won't be added to be plotted as
                    # it will widen the representation range, what would make the represented region
                    # of interest smaller
                    if item.isEnabled():
                        defocusUList.append(defocusU)
                        defocusVList.append(item.getDefocusV())
                        phShList.append(
                            item.getPhaseShift() if item.hasPhaseShift() else 0)
                        resList.append(item.getResolution())

                fig = Figure(figsize=(7, 7), dpi=100)
                defocusPlot = fig.add_subplot(111)
                defocusPlot.grid()
                defocusPlot.set_title(itemSelected)
                defocusPlot.set_xlabel('Tilt angle')
                defocusPlot.set_ylabel('DefocusU', color='tab:red')
                defocusPlot.plot(angList, defocusUList, marker='.',
                                 color='tab:red', label='DefocusU (A)')
                defocusPlot.set_ylabel('DefocusV', color='tab:blue')
                defocusPlot.plot(angList, defocusVList, marker='.',
                                 color='tab:blue', label='DefocusV (A)')

                if item.hasPhaseShift():
                    phShPlot = defocusPlot.twinx()
                    phShPlot.set_ylim(0, 180)
                    phShPlot.set_ylabel('Phase shift', color='tab:green')
                    phShPlot.plot(angList, phShList, marker='.',
                                  color='tab:green', label='Phase shift (deg)')
                else:  # no phase shift, plot resolution instead
                    resPlot = defocusPlot.twinx()
                    resPlot.set_ylim(0, 30)
                    resPlot.set_ylabel('Resolution', color='tab:green')
                    resPlot.plot(angList, resList, marker='.',
                                 color='tab:green', label='Resolution (A)')

                fig.legend()
                canvas = FigureCanvasTkAgg(fig, master=plotterPanel)
                canvas.draw()
                canvas.get_tk_widget().pack(fill=BOTH, expand=0)
                plotterPanel.grid(row=0, column=1, sticky='news')
                break

    def _createPloter(self, event):
        itemSelected = self.tree.getSelectedItem()
        obj = self.tree.getSelectedObj()
        self._checkedItems = self.tree._checkedItems
        if self._checkedItems and self._checkedItems != len(self.provider.getCTFSeries()):
            self.generateSubsetButton['state'] = tk.NORMAL
        else:
            self.generateSubsetButton['state'] = tk.DISABLED

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
