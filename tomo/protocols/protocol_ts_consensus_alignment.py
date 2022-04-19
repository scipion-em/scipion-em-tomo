# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol

import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from tomo import utils


class ProtAssignTransformationMatrixTiltSeries(EMProtocol, ProtTomoBase):
    """
    Perform a consensus of two alignments for the same tilt series.
    """

    _label = 'Tilt-series consensus alignment'
    _devStatus = BETA

    _tsIdList = []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('setOfTiltSeries1',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='First set of tilt-series to be analyzed in the consensus alignment.',
                      label='First set of tilt-series')

        form.addParam('setOfTiltSeries2',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Second set of tilt-series to be analyzed in the consensus alignment.',
                      label='Second set of tilt-series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.generateTsIdList)

        for tsId in self._tsIdList:
            self._insertFunctionStep(self.consensusAlignment, tsId)

    # --------------------------- STEPS functions ----------------------------
    def generateTsIdList(self):
        for ts in self.setOfTiltSeries1.get():
            self._tsIdList.append(ts.getTsid())

        print(self._tsIdList)

        tmpTsIdList = []
        for ts in self.setOfTiltSeries2.get():
            tmpTsIdList.append(ts.getTsid())

        print(tmpTsIdList)

        for tsId in self._tsIdList:
            if tsId not in tmpTsIdList:
                self._tsIdList.remove(tsId)

        print(self._tsIdList)

        if len(self._tsIdList) == 0:
            raise Exception("None matching tilt-series between two sets.")

    def consensusAlignment(self, tsId):
        ts1 = self.setOfTiltSeries1.get().getTiltSeriesFromTsId(tsId)
        ts2 = self.setOfTiltSeries2.get().getTiltSeriesFromTsId(tsId)

        for ti1, ti2 in zip(ts1, ts2):
            ra1, sh1 = utils.getRotationAngleAndShiftFromTM(ti1)
            ra2, sh2 = utils.getRotationAngleAndShiftFromTM(ti2)

    # --------------------------- UTILS functions ----------------------------


