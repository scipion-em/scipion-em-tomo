# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.object import Set, List, String
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.constants import CENTER_GRAVITY
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils

METADATA_INPUT_COORDINATES = "fiducialCoordinates.xmd"


class XmippProtTsConvertCoords3d(EMProtocol, ProtTomoBase):
    """
    Scipion protocol to convert a set of tilt-series coordinates 3d to a set of coordinates 3d associated to a set of
    tomograms.
    """

    _label = 'detect misaligned TS'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfCoordinates',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeriesCoordinates',
                      important=True,
                      label='Input set of coordinates 3D',
                      help='Set of 3D coordinates indicating the position in space of the fiducials. This set should '
                           'be obtained from the previous alignment step of the tilt-series.')

        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.convertCoordinates)

    # --------------------------- STEPS functions ----------------------------

    def convertCoordinates(self):
        sotsc3d = self.inputSetOfCoordinates.get()

        self.getOutputSetOfCoordinates3Ds()

        for coor3d in sotsc3d:
            tsId = coor3d.getTsId()
            tomo = self.getTomoFromTsId(tsId)

            newCoord3D = tomoObj.Coordinate3D()
            newCoord3D.setVolume(tomo)
            newCoord3D.setX(coor3d.getX(), CENTER_GRAVITY)
            newCoord3D.setY(coor3d.getY(), CENTER_GRAVITY)
            newCoord3D.setZ(coor3d.getZ(), CENTER_GRAVITY)

            newCoord3D.setVolId(tsObjId)
            self.outputSetOfCoordinates3D.append(newCoord3D)
            self.outputSetOfCoordinates3D.update(newCoord3D)

            self._store()

    def getOutputSetOfCoordinates3Ds(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()

        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.inputSetOfTomograms.get(),
                                                                      suffix='Coords3d')

            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTomograms.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTomograms.get())
            outputSetOfCoordinates3D.setBoxSize(32)

            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTomograms.get(), outputSetOfCoordinates3D)

        return self.outputSetOfCoordinates3D

    def getTomoFromTsId(self, tsId):
        for tomo in self.inputSetOfTomograms.get():
            if tomo.getTsId() == tsId:
                return tomo
