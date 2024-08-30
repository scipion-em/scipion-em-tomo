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
import logging
from enum import Enum

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.object import Set
from pwem.protocols import EMProtocol
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries
from tomo.protocols import ProtTomoBase

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtCorrectTiltOffset(EMProtocol, ProtTomoBase):
    """
    This protocol allows to correct the offset of a tilt series. This is particularly
    usefull for lamellas. The nominal tilt angle of the sample holder can be 0º
    however, the real tilt angle of the lamella is different, let's say 20º. In this
    it is neccesary to correct the tilt angle by 20º. This means that 20º will be added
    to all tilt angles. Example, a tilt series with nominal tilt angles in steps of 3º
    from -60º to 60º and a tilt offset of 15º will result in a tilt range from -45 to 75º
    """

    _label = 'correct tilt offset'
    _devStatus = BETA
    _possibleOutputs = outputObjects

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputTiltSeries',
                      PointerParam,
                      important=True,
                      label="Tilt series",
                      pointerClass='SetOfTiltSeries',
                      help='Select several sets of tilt-series where to evaluate the consensus in their alignment. '
                           'Output set will bring the information from the first selected set.')

        form.addParam('tiltOffset',
                      FloatParam,
                      label="Tilt Offset (degrees)",
                      default=10,
                      help='This is the offset of the tilt series. This magnitude will be added to all'
                           'tilt angles of the tilt series.')


    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.correctOffSetStep)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def correctOffSetStep(self):
        inputData = self.inputTiltSeries.get()
        print(inputData)

        offset = self.tiltOffset.get()

        suffix = '_offset_corr'
        newSetTs = self._createSet(SetOfTiltSeries, 'tiltseries%s.sqlite', suffix)
        newSetTs.copyInfo(inputData)
        newSetTs.copyAttributes(inputData)
        print(newSetTs)

        for ts in inputData:
            newTs = TiltSeries(tsId=ts.getTsId())
            newTs.copyInfo(ts)
            newSetTs.append(newTs)
            tomoAcq = ts.getAcquisition()
            tomoAcq.setAngleMax(tomoAcq.getAngleMax() + offset)
            tomoAcq.setAngleMax(tomoAcq.getAngleMin() + offset)
            newTs.setAcquisition(tomoAcq)
            for ti in ts:
                newTi = TiltImage()
                newTi.copyInfo(ti)
                newTi.setTiltAngle(ti.getTiltAngle() + offset)
                newTs.append(newTi)

            newTs.setDim(ts.getDim())
            newTs.write()
            newSetTs.update(ts)

        newSetTs.write()
        self._store()

        self._defineOutputs(**{self._possibleOutputs.tiltSeries.name: newSetTs})
        self._defineSourceRelation(inputData, newSetTs)

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary


