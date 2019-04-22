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


class ProtImodEtomo(ProtTomoReconstruct):
    """
    Simple wrapper around etomo to manually reconstruct a Tomogram.

    More info:
        https://bio3d.colorado.edu/imod/doc/etomoTutorial.html
    """

    _label = 'imod etomo'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('action', params.EnumParam,
                      choices=['Build tomogram', 'Register tomogram'],
                      default=0,
                      label='Action', important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Choose *Register tomogram" option when you want to generate'
                           'the output Tomogram from your processing with '
                           'etomo. ')

        group = form.addGroup('Build Tomogram',
                              condition='action==0')
        group.addParam('inputTiltSeries', params.PointerParam,
                      pointerClass='TiltSeries',
                      important=True,
                      label='Input Tilt-Series',
                      help='???')

        group.addParam('excludeList', params.StringParam, default='',
                      label='Exclusion list',
                      help='Provide tilt images IDs (usually starting at 1) '
                           'that you want to exclude from the processing. ')

        group.addParam('binning', params.IntParam, default=2,
                      label='Bin the input images',
                      help='Binning of the input images.')

        group.addParam('markersDiameter', params.IntParam, default=20,
                      label='Fiducial markers diameter (nm)',
                      help='Size of gold beads in nanometers.')

        group.addParam('rotationAngle', params.FloatParam,
                      label='Tilt rotation angle (deg)',
                      help='Angle from the vertical to the tilt axis in raw '
                           'images.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        ts = self.inputTiltSeries.get()
        tsId = ts.getTsId()
        self._insertFunctionStep('convertInputStep', tsId)
        self._insertFunctionStep('runEtomoStep', tsId, interactive=True)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsId):
        ts = self.inputTiltSeries.get()
        workingFolder = self._getExtraPath(tsId)
        prefix = os.path.join(workingFolder, tsId)
        pw.utils.makePath(workingFolder)

        # Write new stack discarding excluded tilts
        excludeList = map(int, self.excludeList.get().split())
        tsStack = prefix + '.st'
        tiList = [ti.clone() for ti in ts]
        tiList.sort(key=lambda ti: ti.getTiltAngle())

        writeTiStack(tiList,
                     outputStackFn=tsStack,
                     outputTltFn=prefix + '.rawtlt',
                     excludeList=excludeList)

        args = '-name %s ' % tsId
        args += '-gold %0.3f ' % self.markersDiameter
        args += '-pixel %0.3f ' % ts.getSamplingRate()
        args += '-binning %d ' % self.binning
        args += '-rotation %0.3f ' % self.rotationAngle
        args += '-userawtlt'

        self.runJob('copytomocoms', args, cwd=workingFolder)

        edfFn = os.path.join(workingFolder, '%s.edf' % tsId)
        minTilt = min(ti.getTiltAngle() for ti in tiList)
        self._writeEtomoEdf(edfFn,
                            {
                                'date': pw.utils.prettyTime(),
                                'name': tsId,
                                'pixelSize': ts.getSamplingRate(),
                                'version': pw.__version__,
                                'minTilt': minTilt,
                                'binning': self.binning,
                                'markerDiameter': self.markersDiameter,
                                'rotationAngle': self.rotationAngle
                            })

    def runEtomoStep(self, tsId):
        workingFolder = self._getExtraPath(tsId)
        self.runJob('etomo', '%s.edf' % tsId, cwd=workingFolder)

        tomoFn = os.path.join(workingFolder, '%s_full.rec' % tsId)
        if os.path.exists(tomoFn):
            self._createOutput(tomoFn)

    def _createOutput(self, tomoFn):
        inputTs = self.inputTiltSeries.get()
        outTomos = self._createSetOfTomograms()
        samplingRate = inputTs.getSamplingRate()

        if self.binning > 1:
            samplingRate *= self.binning.get()

        outTomos.setSamplingRate(samplingRate)

        t = Tomogram(location=tomoFn+':mrc')
        t.setObjId(inputTs.getObjId())
        t.setTsId(inputTs.getTsId())
        outTomos.append(t)

        self._defineOutputs(outputTomograms=outTomos)
        self._defineSourceRelation(self.inputTiltSeries, outTomos)

    # --------------------------- UTILS functions ----------------------------
    def _writeEtomoEdf(self, fn, paramsDict):
        template = """
#%(date)s  // Generated by Scipion, version: %(version)s
Setup.DataSource=CCD
ReconstructionState.InvalidEdgeFunctionsA=no result
Setup.BackupDirectory=
ReconstructionState.InvalidEdgeFunctionsB=no result
ProcessTrack.PostProcessing=Not started
Setup.Combine.ManualCleanup=false
Setup.AxisA.TiltAngle.Type=File
Setup.FiducialessAlignmentB=false
Setup.FiducialessAlignmentA=false
ProcessTrack.FinalAlignedStack-A=Not started
ProcessTrack.FinalAlignedStack-B=Not started
Setup.tiltalign.TargetPatchSizeXandY=700,700
Setup.AxisB.TiltAngle.RangeMin=%(minTilt)f
Setup.Combine.TempDirectory=
Setup.Combine.UseList=
Setup.Combine.PatchBoundaryYMax=0
Setup.WholeTomogramSampleA=false
Setup.DatasetName=%(name)s
Setup.FiducialDiameter=%(markerDiameter)f
Setup.WholeTomogramSampleB=false
Setup.SetFEIPixelSize=true
ReconstructionState.A.AdjustOrigin=true
Setup.Setup.OrigImageStackExt=.st
ProcessTrack.RevisionNumber=2.0
Setup.ProjectLog.FrameSize.Width=683
Setup.Combine.PatchBoundaryYMin=0
Setup.Combine.MaxPatchBoundaryZMax=0
ProcessTrack.CoarseAlignment-A=Not started
Setup.ViewType=Single View
ProcessTrack.CoarseAlignment-B=Not started
ReconstructionState.TrimvolFlipped=no result
Setup.Track.B.TrackMethod=Seed
Setup.Track.A.TrackMethod=Seed
Setup.A.SizeToOutputInXandY=/
ProcessTrack.TomogramGeneration-A=Not started
ProcessTrack.TomogramGeneration-B=Not started
ProcessTrack.CleanUp=Not started
Setup.AxisA.TiltAngle.RangeMin=%(minTilt)f
Setup.AxisB.TiltAngle.TiltAngleFilename=
Setup.Pos.A.NewDialog=true
Setup.Squeezevol.LinearInterpolation=false
Setup.Stack.B.Is.Twodir=false
Setup.Combine.RevisionNumber=1.2
Setup.Stack.B.CTF.AutoFit.RangeAndStep=-Infinity,-Infinity
Setup.UseLocalAlignmentsB=true
Setup.UseLocalAlignmentsA=true
Setup.AxisB.TiltAngle.RangeStep=1.0
Setup.Combine.FiducialMatchListA=
Setup.Combine.FiducialMatchListB=
Setup.B.SizeToOutputInXandY=/
Setup.Binning=%(binning)s
Setup.Combine.ModelBased=false
ReconstructionState.MadeZFactorsB=no result
ReconstructionState.MadeZFactorsA=no result
ReconstructionState.SqueezevolFlipped=no result
ProcessTrack.FiducialModel-B=Not started
ProcessTrack.FiducialModel-A=Not started
Setup.Combine.PatchBoundaryXMax=0
ProcessTrack.Setup=Complete
ProcessTrack.PreProcessing-A=Not started
ProcessTrack.PreProcessing-B=Not started
ProcessTrack.FineAlignment-B=Not started
ProcessTrack.FineAlignment-A=Not started
Setup.AxisA.TiltAngle.TiltAngleFilename=
Setup.AxisA.TiltAngle.RangeStep=1.0
Setup=-Infinity,-Infinity
Setup.Combine.PatchBoundaryXMin=0
Setup.MagGradientFile=
Setup.RevisionNumber=1.12
Setup.Track.B.SeedModel.Transfer=true
Setup.Track.A.Raptor.UseRawStack=false
Setup.PixelSize=%(pixelSize)f
ReconstructionState.B.AdjustOrigin=true
ReconstructionState.Combine.ScriptsCreated=no result
Setup.Combine.Transfer=true
Setup.DistortionFile=
ReconstructionState.UsedLocalAlignmentsA=no result
ReconstructionState.UsedLocalAlignmentsB=no result
Setup.ProjectLog.FrameSize.Height=230
Setup.Stack.B.Twodir=0.0
Setup.Combine.FiducialMatch=BothSides
Setup.AxisType=Single Axis
Setup.ImageRotationA=%(rotationAngle)f
ProcessTrack.TomogramPositioning-A=Not started
Setup.ImageRotationB=
Setup.Stack.A.Is.Twodir=false
Setup.Pos.B.NewDialog=true
ProcessTrack.TomogramPositioning-B=Not started
Setup.Combine.PatchBoundaryZMax=0
Setup.DefaultGpuProcessing=false
Setup.Track.A.SeedModel.Auto=true
Setup.Combine.PatchSize=M
Setup.AxisB.TiltAngle.Type=Extract
Setup.Combine.PatchBoundaryZMin=0
Setup.ProjectLog.FrameLocation.Y=55
Setup.ProjectLog.FrameLocation.X=95
Setup.DefaultParallel=false
Setup.ProjectLog.Visible=true
Setup.tiltalign.NumberOfLocalPatchesXandY=5,5
Setup.Combine.PatchRegionModel=
ReconstructionState.NewstFiducialessAlignmentA=no result
ReconstructionState.NewstFiducialessAlignmentB=no result
Setup.Stack.A.Twodir=0.0
ProcessTrack.TomogramCombination=Not started
        """
        with open(fn, 'w') as f:
            f.write(template % paramsDict)