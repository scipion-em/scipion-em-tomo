# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
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


"""
This module implements visualization program
for input tomograms.
"""

import os
from distutils.spawn import find_executable

import pyworkflow.protocol.params as params
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
import pwem.viewers as viewers
from pwem.emlib.image import ImageHandler


from tomo.protocols import ProtImportTomograms, ProtImportSubTomograms

TOMOGRAM_SLICES = 1
TOMOGRAM_CHIMERA = 0


class ViewerProtImportTomograms(ProtocolViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """

    _label = 'viewer input tomogram'
    _targets = [ProtImportTomograms, ProtImportSubTomograms]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization of input tomograms')
        form.addParam('displayTomo', params.EnumParam,
                      choices=['chimera', 'slices'],
                      default=TOMOGRAM_CHIMERA,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Display tomogram with',
                      help='*chimera*: display tomograms as surface with '
                           'Chimera.\n *slices*: display tomograms as 2D slices '
                           'along z axis.\n If number of tomograms == 1, '
                           'a system of coordinates is shown'
                      )

    def visualize(self, obj, **args):
        views = self._showTomogramsSlices()
        if views:
            for v in views:
                v.show()

    def _getVisualizeDict(self):
        return {
            'displayTomo': self._showTomograms,
        }

    def _validate(self):
        if (self.displayTomo == TOMOGRAM_CHIMERA
                and find_executable(viewers.viewer_chimera.Chimera.getProgram()) is None):
            return ["chimera is not available. "
                    "Either install it or choose option 'slices'. "]
        return []

    # =========================================================================
    # ShowTomograms
    # =========================================================================

    def _showTomograms(self, paramName=None):
        if self.displayTomo == TOMOGRAM_CHIMERA:
            return self._showTomogramsChimera()

        elif self.displayTomo == TOMOGRAM_SLICES:
            return self._showTomogramsSlices()

    def _createSetOfObjects(self):
        if hasattr(self.protocol, 'outputTomograms'):
            setOfObjects = self.protocol.outputTomograms
            sampling = self.protocol.outputTomograms.getSamplingRate()
        elif hasattr(self.protocol, 'outputSubTomograms'):
            setOfObjects = self.protocol.outputSubTomograms
            sampling = self.protocol.outputSubTomograms.getSamplingRate()
        return sampling, setOfObjects

    def _showTomogramsChimera(self):
        """ Create a chimera script to visualize selected tomograms. """
        tmpFileNameCMD = self.protocol._getTmpPath("chimera.cmd")
        f = open(tmpFileNameCMD, "w")
        sampling, _setOfObjects = self._createSetOfObjects()
        count = 0  # first model in chimera is a tomogram

        if len(_setOfObjects) == 1:
            count = 1  # first model in chimera is the bild file
            # if we have a single tomogram then create axis
            # as bild file. Chimera must read the bild file first
            # otherwise system of coordinates will not
            # be in the center of the window

            dim = _setOfObjects.getDim()[0]
            tmpFileNameBILD = os.path.abspath(self.protocol._getTmpPath(
                "axis.bild"))
            viewers.viewer_chimera.Chimera.createCoordinateAxisFile(dim,
                                                                    bildFileName=tmpFileNameBILD,
                                                                    sampling=sampling)
            f.write("open %s\n" % tmpFileNameBILD)
            f.write("cofr 0,0,0\n")  # set center of coordinates
            count = 1  # skip first model because is not a 3D map

        for tomo in _setOfObjects:
            localTomo = os.path.abspath(ImageHandler.removeFileType(
                tomo.getFileName()))
            if localTomo.endswith("stk"):
                self.showError("Extension .stk is not supported")
            f.write("open %s\n" % localTomo)
            f.write("volume#%d style surface voxelSize %f\n" %
                    (count, sampling))
            count += 1

        if len(_setOfObjects) > 1:
            f.write('tile\n')
        else:
            x, y, z = tomo.getShiftsFromOrigin()
            f.write("volume#1 origin %0.2f,%0.2f,%0.2f\n" % (x, y, z))
        f.close()
        return [viewers.viewer_chimera.ChimeraView(tmpFileNameCMD)]

    def _showTomogramsSlices(self):
        # Write an sqlite with all tomograms selected for visualization.
        sampling, setOfObjects = self._createSetOfObjects()

        viewParams= {viewers.showj.MODE: viewers.showj.MODE_MD}
        view = self.objectView(setOfObjects, viewParams=viewParams)
        view.setMemory(viewers.showj.getJvmMaxMemory() + 2)

        return [view]
