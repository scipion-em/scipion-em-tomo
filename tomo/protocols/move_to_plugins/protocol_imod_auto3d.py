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

import os

import pyworkflow as pw
import pyworkflow.protocol.params as params

from tomo.objects import TiltSeriesDict, TiltSeries
from tomo.protocols import ProtTomoReconstruct
from tomo.convert import writeTiStack


class ProtImodAuto3D(ProtTomoReconstruct):
    """
    Simple protocol to do a quick Tomogram reconstruction with IMOD.
    (Sample scripts provided by Javi Chichon)
    """

    _label = 'imod auto3d'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputTiltSeries', params.PointerParam,
                      pointerClass='TiltSeries,SetOfTiltSeries',
                      important=True,
                      label='Input Tilt-Series',
                      help='???')

        form.addParam('bin', params.BooleanParam, default=2,
                      label='Bin the input images?',
                      help='Binning of the input images.')

        form.addParam('zWidth', params.IntParam,
                      label='Z width of the tomograph',
                      help='???')

    # -------------------------- INSERT steps functions ---------------------
    def _loadInputTs(self):
        """ Load input TiltSeries (set or single item). """
        self._tsDict = TiltSeriesDict()
        inputTs = self.inputTiltSeries.get()

        if isinstance(inputTs, TiltSeries):
            self._tsDict.addTs(inputTs, includeTi=True)
        else:  # SetOfTiltSeries
            for ts in inputTs:
                print("Adding ts: %s" % ts)
                self._tsDict.addTs(ts, includeTi=True)

    def _insertAllSteps(self):
        self._loadInputTs()

        for ts in self._tsDict:  # Read if we input SetOfTiltSeries:
            self._insertFunctionStep('processTsStep', ts.getTsId())

        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ----------------------------
    def processTsStep(self, tsId):
        ts = self._tsDict.getTs(tsId)

        workingFolder = self._getTmpPath(tsId)
        prefix = os.path.join(workingFolder, tsId)

        pw.utils.makePath(workingFolder)

        # Write new stack discarding excluded tilts
        excludeList = [1]  # FIXME: Take this from input

        print("Full list:")
        for ti in self._tsDict.getTiList(tsId):
            print(str(ti))

        tiList = [ti for ti in self._tsDict.getTiList(tsId)
                  if ti.getObjId() not in excludeList]
        tiList.sort(key=lambda ti: ti.getTiltAngle())
        print("Used list: ")
        for ti in tiList:
            print(str(ti))
        writeTiStack(tiList,
                     outputStackFn=prefix + '.st',
                     outputTltFn=prefix + '.rawtlt')

        if not pw.utils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pw.utils.cleanPath(workingFolder)

    def createOutputStep(self):
        pass

