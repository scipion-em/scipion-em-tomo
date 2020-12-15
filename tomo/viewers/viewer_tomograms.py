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
import tempfile

import pyworkflow.protocol.params as params
from pwem.viewers import EmProtocolViewer, ChimeraView
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
import pwem.viewers as viewers
from pwem.emlib.image import ImageHandler


from tomo.protocols import ProtImportTomograms, \
    ProtImportSubTomograms
from tomo.objects import SetOfSubTomograms


TOMOGRAM_SLICES = 1
TOMOGRAM_CHIMERA = 0


class ViewerProtImportTomograms(EmProtocolViewer):
    """ Wrapper to visualize tomo objects
    with default viewers such as xmipp/showj and chimera. """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtImportTomograms, ProtImportSubTomograms, SetOfSubTomograms]
    _label = 'viewer input tomogram'

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

#    def visualize(self, obj, **args):
#        views = self._showTomogramsSlices()
#        if views:
#            for v in views:
#                v.show()
#
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

    def _sourceIsOutput(self):
        """If source is a SetOfSubtomograms"""
        return isinstance(self.protocol, SetOfSubTomograms)

    def _getOutput(self):
        if self._sourceIsOutput():
            # Protocol here is a SetOfSubtomograms
            return self.protocol
        elif hasattr(self.protocol, 'outputTomograms'):
            return self.protocol.outputTomograms
        elif hasattr(self.protocol, 'outputSubTomograms'):
            return self.protocol.outputSubTomograms

    def _getSetAndSampling(self):
        setOfObjects = self._getOutput()
        sampling = setOfObjects.getSamplingRate()
        return sampling, setOfObjects

    def _getExtraFolder(self, fileName):
        if self._sourceIsOutput():
            return os.path.join(tempfile.mkdtemp(), fileName)
        else:
            return self.protocol._getExtraPath(fileName)

    def _showTomogramsChimera(self):
        """ Create a chimera script to visualize selected tomograms. """
        tmpFileNameCMD = self._getExtraFolder("chimera.cxc")
        f = open(tmpFileNameCMD, "w")
        sampling, _setOfObjects = self._getSetAndSampling()
        count = 1  # first model in chimera is a tomogram

        if len(_setOfObjects) == 1:
            count = 2  # first model in chimera is the bild file
            # if we have a single tomogram then create axis
            # as bild file. Chimera must read the bild file first
            # otherwise system of coordinates will not
            # be in the center of the window

            dim = _setOfObjects.getDim()[0]
            tmpFileNameBILD = os.path.abspath(self._getExtraFolder("axis.bild"))
            viewers.viewer_chimera.Chimera.createCoordinateAxisFile(dim,
                                                                    bildFileName=tmpFileNameBILD,
                                                                    sampling=sampling)
            f.write("open %s\n" % tmpFileNameBILD)
            f.write("cofr 0,0,0\n")  # set center of coordinates
        oldFileName = ""
        for tomo in _setOfObjects:
            fileName = tomo.getFileName()
            if fileName == oldFileName:
                continue
            else:
                oldFileName = fileName
            localTomo = os.path.abspath(ImageHandler.removeFileType(
                tomo.getFileName()))
            if localTomo.endswith("stk"):
                self.showError("Extension .stk is not supported")
            f.write("open %s\n" % localTomo)
            f.write("volume #%d style surface voxelSize %f\n" %
                    (count, sampling))
            f.write("volume #%d level 1\n" % count)
            count += 1

        if len(_setOfObjects) > 1:
            f.write('tile\n')
        else:
            x, y, z = tomo.getShiftsFromOrigin()
            f.write("volume #2 origin %0.2f,%0.2f,%0.2f\n" % (x, y, z))
        f.close()
        view = ChimeraView(tmpFileNameCMD)
        return [view]

    def _showTomogramsSlices(self):
        # Write an sqlite with all tomograms selected for visualization.
        sampling, setOfObjects = self._getSetAndSampling()

        viewParams= {viewers.showj.MODE: viewers.showj.MODE_MD}
        view = viewers.views.ObjectView(self._project, setOfObjects.strId(),
                                        setOfObjects.getFileName(), viewParams=viewParams)
        view.setMemory(viewers.showj.getJvmMaxMemory() + 2)

        return [view]
