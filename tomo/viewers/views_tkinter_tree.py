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
import os
import threading
import numpy as np

from pwem.viewers import showj
from pwem.viewers.showj import runJavaIJapp
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.plugin import Domain
from pyworkflow.gui.dialog import ListDialog, ToolbarListDialog
import pyworkflow.config as conf

import pyworkflow.viewer as pwviewer

import pyworkflow.utils as pwutils

import tomo.objects
from tomo.convert.convert import getMeshVolFileName


class TiltSeriesTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from TiltSeries and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_TS = 'Tilt Series'
    COL_TI = 'Path'
    COL_TI_IX = 'Index'
    COL_TI_ID = 'Order'
    COL_TI_ANGLE = 'Angle'
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
            (self.COL_TI_ID, 50),
            (self.COL_TI_ANGLE, 50),
            (self.COL_TI_IX, 50),
            (self.COL_TI,400),
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
            values = [objId, obj.getTiltAngle(),
                      str(obj.getLocation()[0]),
                      str(obj.getLocation()[1])]
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
            viewers = Domain.findViewers(obj.getClassName(),
                                         pwviewer.DESKTOP_TKINTER)
            for viewerClass in viewers:
                def createViewer(viewerClass, obj):
                    proj = self.protocol.getProject()
                    item = self.tiltseries[obj.getObjId()]  # to load mapper
                    return lambda : viewerClass(project=proj).visualize(item)
                actions.append(('Open with %s' % viewerClass.__name__,
                                createViewer(viewerClass, obj)))

        return actions


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
        self._provider = TiltSeriesTreeProvider(protocol, tiltSeries)

    def show(self):
        dlg = ListDialog(self._tkParent, 'TiltSeries display', self._provider)


class TomogramsTreeProvider(TreeProvider):
    """ Populate Tree from SetOfTomograms. """

    def __init__(self, tomoList, path, mode):
        TreeProvider.__init__(self)
        self.tomoList = tomoList
        self._path = path
        self._mode = mode

    def getColumns(self):
        return [('Tomogram', 300), ('status', 150)]

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
                    'text': tomogramName, 'values': ("TODO"),
                    'tags': ("pending")}
        else:
            return {'key': tomogramName, 'parent': None,
                    'text': tomogramName, 'values': ("DONE"),
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
        if isinstance(obj, tomo.objects.Mesh):
            meshName = 'Mesh %d' % obj.getObjId()
            tomoName = pwutils.removeBaseExt(obj.getVolume().getFileName())
            return {'key': tomoName + '-' + str(obj.getObjId()), 'parent': self._parentDict.get(obj.getObjId(), None),
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

    def __init__(self, parent, viewer, **kwargs):
        self.path = kwargs.get("path", None)
        self.provider = kwargs.get("provider", None)
        if viewer:
            ToolbarListDialog.__init__(self, parent,
                                       "Tomogram List",
                                       allowsEmptySelection=False,
                                       itemDoubleClick=self.doubleClickViewer,
                                       **kwargs)
        else:
            ToolbarListDialog.__init__(self, parent,
                                       "Tomogram List",
                                       allowsEmptySelection=False,
                                       itemDoubleClick=self.doubleClickOnTomogram,
                                       **kwargs)

    def refresh_gui_viewer(self):
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui_viewer)
        else:
            pwutils.cleanPath(os.path.join(self.path, 'mesh.txt'))
            pwutils.cleanPath(self.macroPath)

    def refresh_gui(self):
        self.tree.update()
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui)
        else:
            pwutils.cleanPath(self.macroPath)

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
        self.macroPath = os.path.join(conf.Config.SCIPION_TMP, "AutoSave_ROI.ijm")
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
        self.macroPath = os.path.join(conf.Config.SCIPION_TMP, "View_ROI.ijm")
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
