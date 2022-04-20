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
from pyworkflow.object import Set

from pwem.protocols import EMProtocol
from pwem.objects import Float

from tomo.objects import TiltSeries, TiltImage
from tomo.protocols import ProtTomoBase
from tomo import utils


class ProtConsensusAlignmentTS(EMProtocol, ProtTomoBase):
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
                      help='First set of tilt-series to be analyzed in the consensus alignment. The unrelated '
                           'alignment information will be taken form this set.',
                      label='First set of tilt-series')

        form.addParam('setOfTiltSeries2',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Second set of tilt-series to be analyzed in the consensus alignment.',
                      label='Second set of tilt-series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.consensusAlignment)

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def consensusAlignment(self):
        tsIdList = self.generateTsIdList()

        self.getOutputAlignmentConsensusSetOfTiltSeries()

        for tsId in tsIdList:
            ts1 = self.setOfTiltSeries1.get().getTiltSeriesFromTsId(tsId)
            ts2 = self.setOfTiltSeries2.get().getTiltSeriesFromTsId(tsId)

            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(ts1)
            self.outputAlignmentConsensusSetOfTiltSeries.append(newTs)

            for ti1, ti2 in zip(ts1, ts2):
                newTi = TiltImage()
                newTi.copyInfo(ti1, copyId=True)
                newTi.setLocation(ti1.getLocation())
                newTi.setTransform(ti1.getTransform())

                ra1, sh1 = utils.getRotationAngleAndShiftFromTM(ti1)
                ra2, sh2 = utils.getRotationAngleAndShiftFromTM(ti2)

                newTi._angleDiff = Float(abs(ra1-ra2))
                newTi._shiftDiff = Float((abs(sh1[0]-sh2[0]) + abs(sh1[1]-sh2[1]))/2)  # Shift average in both directions

                newTs.append(newTi)

            newTs.setDim(ts1.getDim())
            newTs.write()

            self.outputAlignmentConsensusSetOfTiltSeries.update(newTs)
            self.outputAlignmentConsensusSetOfTiltSeries.updateDim()
            self.outputAlignmentConsensusSetOfTiltSeries.write()

        self._store()

    def closeOutputSetsStep(self):
        self.outputAlignmentConsensusSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.outputAlignmentConsensusSetOfTiltSeries.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def generateTsIdList(self):
        tsIdList = []
        tmpTsIdList = []

        for ts in self.setOfTiltSeries1.get():
            tsIdList.append(ts.getTsId())

        for ts in self.setOfTiltSeries2.get():
            tmpTsIdList.append(ts.getTsId())

        for tsId in tsIdList:
            if tsId not in tmpTsIdList:
                tsIdList.remove(tsId)

        if len(tsIdList) == 0:
            raise Exception("None matching tilt-series between two sets.")

        print(tsIdList)
        return tsIdList

    def getOutputAlignmentConsensusSetOfTiltSeries(self):
        if not hasattr(self, "outputSetOfTiltSeries"):
            outputAlignmentConsensusSetOfTiltSeries = self._createSetOfTiltSeries(suffix='AliConsensus')
            outputAlignmentConsensusSetOfTiltSeries.copyInfo(self.setOfTiltSeries1.get())
            outputAlignmentConsensusSetOfTiltSeries.setDim(self.setOfTiltSeries1.get().getDim())

            self._defineOutputs(outputAlignmentConsensusSetOfTiltSeries=outputAlignmentConsensusSetOfTiltSeries)
            self._defineSourceRelation(self.setOfTiltSeries1, outputAlignmentConsensusSetOfTiltSeries)
        return self.outputAlignmentConsensusSetOfTiltSeries

