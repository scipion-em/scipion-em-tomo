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

from tomo.objects import TiltSeriesDict, TiltSeries, Tomogram
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
        form.addParam('excludeList', params.StringParam, default='',
                      label='Exclusion list',
                      help='Provide tilt images IDs (usually starting at 1) '
                           'that you want to exclude from the processing. ')
        form.addParam('bin', params.IntParam, default=2,
                      label='Bin the input images',
                      help='Binning of the input images.')
        form.addParam('zWidth', params.IntParam,
                      label='Z width of the tomograph',
                      help='???')

        group = form.addGroup('Fiducial markers')
        group.addParam('useRaptor', params.BooleanParam, default=True,
                       label='Use RAPTOR for automatic markers tracking?',
                       help='???')
        group.addParam('markersDiameter', params.IntParam, default=20,
                       label='Markers diameter (nm)',
                       help='Size of gold beads in nanometers.')
        group.addParam('markersNumber', params.IntParam, default=20,
                       label='Number of markers to track',
                       help='???')

    # -------------------------- INSERT steps functions ---------------------
    def _loadInputTs(self):
        """ Load input TiltSeries (set or single item). """
        inputTs = self.inputTiltSeries.get()
        self._tsDict = TiltSeriesDict()

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
        excludeList = map(int, self.excludeList.get().split())

        tiList = self._tsDict.getTiList(tsId)
        tiList.sort(key=lambda ti: ti.getTiltAngle())
        tsStack = prefix + '.st'

        writeTiStack(tiList,
                     outputStackFn=tsStack,
                     outputTltFn=prefix + '.rawtlt',
                     excludeList=excludeList)

        here = os.path.abspath(os.path.dirname(__file__))
        args = os.path.join(here, 'script_imod_auto3d.py')
        tomoName = self._getTomoName(tsId)
        args += ' --output %s ' % tomoName
        args += '%s --widthz %d --bin %d ' % (os.path.basename(tsStack),
                                              self.zWidth, self.bin)
        if self.useRaptor:
            args += ('--raptor --markers_diameter %d --markers_number %d'
                     % (self.markersDiameter, self.markersNumber))

        self.runJob('python', args, cwd=workingFolder)

        tomoPath = os.path.join(workingFolder, tomoName)

        if os.path.exists(tomoPath):
            pw.utils.moveFile(tomoPath, self._getPath(tomoName))
        else:
            print("ERROR: The expected tomogram for Tilt-Series %s "
                  "was not properly generated. " % tsId)

        if not pw.utils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pw.utils.cleanPath(workingFolder)

    def createOutputStep(self):
        inputTs = self.inputTiltSeries.get()
        outTomos = self._createSetOfTomograms()
        samplingRate = inputTs.getSamplingRate()

        if self.bin > 1:
            samplingRate *= self.bin.get()

        outTomos.setSamplingRate(samplingRate)

        if not hasattr(self, '_tsDict'):
            self._loadInputTs()

        for ts in self._tsDict:  # Read if we input SetOfTiltSeries:
            t = Tomogram(location=self._getPath(self._getTomoName(ts.getTsId())))
            outTomos.append(t)

        self._defineOutputs(outputTomograms=outTomos)

    # --------------------------- UTILS functions ----------------------------
    def _getTomoName(self, tsId):
        return '%s_tomo.mrc' % tsId

