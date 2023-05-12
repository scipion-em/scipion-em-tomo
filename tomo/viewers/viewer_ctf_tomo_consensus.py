# **************************************************************************
# *
# * Authors: Daniel Marchan Torres    (da.marchan@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import matplotlib.pyplot as plt

from .viewers_data import CtfEstimationTomoViewer
import tomo.objects
import pyworkflow.viewer as pwviewer
from tomo.protocols import ProtCTFTomoSeriesConsensus
from pwem.viewers import ObjectView, DataViewer, MODE, MODE_MD, VISIBLE
from pyworkflow.protocol import LabelParam

class CtfConsensusTomoViewer(CtfEstimationTomoViewer):
    """ This class implements a view using Tkinter CtfEstimationListDialog.
    """

    _label = 'CTF consensus viewer'
    _targets = [ProtCTFTomoSeriesConsensus]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def plotExtra(self, ctfSet, ctfId):
        tsId =ctfSet.getTsId()
        fnAstigPlot = self._protocol.getAstigErrorPlot(str(tsId))
        fnDefPlot = self._protocol.getDefocusErrorPlot(str(tsId))
        fnResPlot = self._protocol.getResolutionErrorPlot(str(tsId))

        if os.path.exists(fnAstigPlot):
            fig = plt.figure()
            imageAstig = plt.imread(fnAstigPlot)
            figAstig = plt.imshow(imageAstig)
            figAstig.axes.get_xaxis().set_visible(False)
            figAstig.axes.get_yaxis().set_visible(False)

        if os.path.exists(fnDefPlot):
            fig = plt.figure()
            imageDef = plt.imread(fnDefPlot)
            figDef= plt.imshow(imageDef)
            figDef.axes.get_xaxis().set_visible(False)
            figDef.axes.get_yaxis().set_visible(False)

        if os.path.exists(fnResPlot):
            fig = plt.figure()
            imageRes = plt.imread(fnResPlot)
            figRes = plt.imshow(imageRes)
            figRes.axes.get_xaxis().set_visible(False)
            figRes.axes.get_yaxis().set_visible(False)

        plt.show()


