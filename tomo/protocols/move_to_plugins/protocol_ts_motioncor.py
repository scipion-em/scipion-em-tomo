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

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
from pyworkflow.em.data import MovieAlignment

from tomo.protocols import ProtTsCorrectMotion

import motioncorr
from motioncorr.convert import *
from motioncorr.constants import *


class ProtTsMotionCorr(ProtTsCorrectMotion):
    """
    This protocol wraps motioncor2 movie alignment program developed at UCSF.

    Motioncor2 performs anisotropic drift correction and dose weighting
        (written by Shawn Zheng @ David Agard lab)
    """

    _label = 'tiltseries motioncor'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtTsCorrectMotion._defineParams(self, form)

        form.addSection(label="Motioncor2 params")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=cons.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " Motioncor2 can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParam('doApplyDoseFilter', params.BooleanParam, default=True,
                      label='Apply dose filter',
                      help='Apply a dose-dependent filter to frames before '
                           'summing them. Pre-exposure and dose per frame '
                           'should  be specified during movies import.')

        line = form.addLine('Number of patches',
                            help='Number of patches to be used for patch based '
                                 'alignment. Set to *0 0* to do only global motion '
                                 'correction. \n')
        line.addParam('patchX', params.IntParam, default=5, label='X')
        line.addParam('patchY', params.IntParam, default=5, label='Y')

        form.addParam('patchOverlap', params.IntParam, default=0,
                      label='Patches overlap (%)',
                      help='In versions > 1.0.1 it is possible to specify'
                           'the overlapping between patches. '
                           '\nFor example, overlap=20 means that '
                           'each patch will have a 20% overlapping \n'
                           'with its neighboring patches in each dimension.')

        form.addParam('group', params.IntParam, default='1',
                      label='Group N frames',
                      help='Group every specified number of frames by adding '
                           'them together. The alignment is then performed on '
                           'the summed frames. By default, no grouping is '
                           'performed.')

        form.addParam('tol', params.FloatParam, default='0.5',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Tolerance (px)',
                      help='Tolerance for iterative alignment, default *0.5px*.')

        form.addParam('doSaveUnweightedMic', params.BooleanParam, default=True,
                      condition='doApplyDoseFilter',
                      label="Save unweighted micrographs?",
                      help="Yes by default, if you have selected to apply a "
                           "dose-dependent filter to the frames")

        group = form.addGroup('Magnification correction')
        group.addParam('doMagCor', params.BooleanParam, default=False,
                       label='Correct anisotropic magnification?',
                       help='Correct anisotropic magnification by '
                            'stretching image along the major axis, '
                            'the axis where the lower magnification is '
                            'detected.')
        group.addParam('useEst', params.BooleanParam, default=True,
                       label='Use previous estimation?',
                       condition='doMagCor',
                       help='Use previously calculated parameters of '
                            'magnification anisotropy (from magnification '
                            'distortion estimation protocol).')
        group.addParam('inputEst', params.PointerParam,
                       pointerClass='ProtMagDistEst',
                       condition='useEst and doMagCor',
                       label='Input protocol',
                       help='Select previously executed estimation protocol.')
        group.addParam('scaleMaj', params.FloatParam, default=1.0,
                       condition='not useEst and doMagCor',
                       label='Major scale factor',
                       help='Major scale factor.')
        group.addParam('scaleMin', params.FloatParam, default=1.0,
                       condition='not useEst and doMagCor',
                       label='Minor scale factor',
                       help='Minor scale factor.')
        group.addParam('angDist', params.FloatParam, default=0.0,
                       condition='not useEst and doMagCor',
                       label='Distortion angle (deg)',
                       help='Distortion angle, in degrees.')

        form.addParam('gainRot', params.EnumParam,
                      choices=['no rotation', '90 degrees',
                               '180 degrees', '270 degrees'],
                      label="Rotate gain reference:",
                      default=NO_ROTATION,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help="Rotate gain reference counter-clockwise.")

        form.addParam('gainFlip', params.EnumParam,
                      expertLevel=cons.LEVEL_ADVANCED,
                      choices=['no flip', 'upside down', 'left right'],
                      label="Flip gain reference:", default=NO_FLIP,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help="Flip gain reference after rotation.")

        form.addParam('defectFile', params.FileParam, allowsNull=True,
                      label='Camera defects file',
                      help='Defect file that stores entries of defects on camera.\n'
                           'Each entry corresponds to a rectangular region in image. '
                           'The pixels in such a region are replaced by '
                           'neighboring good pixel values. Each entry contains '
                           '4 integers x, y, w, h representing the x, y '
                           'coordinates, width, and height, respectively.')

        form.addParam('extraParams2', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""*Extra parameters:*\n
        -Bft\t\t100\t\tBFactor for alignment, in px^2.\n
        -Iter\t\t5\t\tMaximum iterations for iterative alignment.\n
        -MaskCent\t\t0 0\t\tCenter of subarea that will be used for alignment,
        \t\t\t\tdefault *0 0* corresponding to the frame center.\n
        -MaskSize\t\t1.0 1.0\t\tThe size of subarea that will be used for alignment,
        \t\t\t\tdefault *1.0 1.0* corresponding full size.\n
        -Align\t\t1\t\tGenerate aligned sum (1) or simple sum (0).\n
        -FmRef\t\t-1\t\tSpecify which frame to be the reference to which
        \t\t\t\tall other frames are aligned, by default (-1) the
        \t\t\t\tthe central frame is chosen. The central frame is
        \t\t\t\tat N/2 based upon zero indexing where N is the
        \t\t\t\tnumber of frames that will be summed, i.e., not
        \t\t\t\tincluding the frames thrown away.\n
        -Tilt\t\t0 0\t\tTilt angle range for a dose fractionated tomographic
        \t\t\t\ttilt series, e.g. *-60 60*\n

    Since version *1.1.0*:\n
        -GpuMemUsage\t\t0.5\t\tSpecify how much GPU memory is used to buffer movie frames.
        \t\t\t\tIt is recommended when running side by side processes in
        \t\t\t\tthe same card. By default is 50% (i. e 0.5).\n
        -InFmMotion\t\t1\t\tTakes into account of motion-induced blurring of
        \t\t\t\teach frame. It has shown resolution improvement
        \t\t\t\tin some test cases. By default this option is off.\n
        -Bft\t\t500 150\t\tSince version 1.1.0 this option can take two arguments.
        \t\t\t\tFirst one is used in global-motion measurement and the
        \t\t\t\tsecond one is for local-motion. (default 500 150).""")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def _processTiltImageM(self, workingFolder, tiltImageM,
                           initialDose, dosePerFrame, gain, dark, *args):
        outputFn, outputFnDW = self._getOutputTiltImagePaths(tiltImageM)

        def _getPath(path):
            """ shortcut to get relative path from workingFolder. """
            return os.path.relpath(path, workingFolder)

        self.info("workingFolder: %s" % workingFolder)
        self.info("outputFn: %s" % outputFn)
        self.info("outputFn (rel): %s" % _getPath(outputFn))

        logFile = self._getPath(self._getMovieLogFile(tiltImageM))

        a0, aN = self._getRange(tiltImageM, 'align')

        logFileBase = (logFile.replace('0-Full.log', '').replace('0-Patch-Full.log', ''))
        # default values for motioncor2 are (1, 1)
        cropDimX = self.cropDimX.get() or 1
        cropDimY = self.cropDimY.get() or 1

        numbOfFrames = self._getNumberOfFrames(tiltImageM)
        order = tiltImageM.getAcquisitionOrder()

        # reset values = 1 to 0 (motioncor2 does it automatically,
        # but we need to keep this for consistency)
        if self.patchX == 1:
            self.patchX.set(0)
        if self.patchY == 1:
            self.patchY.set(0)

        argsDict = {
            '-OutMrc': '"%s"' % _getPath(outputFn),
            '-Patch': '%d %d' % (self.patchX, self.patchY),
            '-MaskCent': '%d %d' % (self.cropOffsetX, self.cropOffsetY),
            '-MaskSize': '%d %d' % (cropDimX, cropDimY),
            '-FtBin': self.binFactor.get(),
            '-Tol': self.tol.get(),
            '-Group': self.group.get(),
            '-FmDose': dosePerFrame,
            '-Throw': '%d' % a0,
            '-Trunc': '%d' % (abs(aN - numbOfFrames + 1)),
            '-PixSize': tiltImageM.getSamplingRate(),
            '-kV': tiltImageM.getAcquisition().getVoltage(),
            '-LogFile': logFileBase,
            '-InitDose': initialDose + order * dosePerFrame,
            '-OutStack': 0
        }

        if self.defectFile.get():
            argsDict['-DefectFile'] = self.defectFile.get()

        patchOverlap = self.getAttributeValue('patchOverlap', None)
        if patchOverlap:  # 0 or None is False
            argsDict['-Patch'] += " %d" % patchOverlap

        if self.doMagCor:
            if self.useEst:
                inputEst = self.inputEst.get().getOutputLog()
                input_params = parseMagEstOutput(inputEst)
                argsDict['-Mag'] = '%0.3f %0.3f %0.3f' % (
                    input_params[1],
                    input_params[2],
                    input_params[0])
            else:
                argsDict['-Mag'] = '%0.3f %0.3f %0.3f' % (self.scaleMaj,
                                                          self.scaleMin,
                                                          self.angDist)

        tiFn = tiltImageM.getFileName()
        inputFn = os.path.abspath(tiFn)

        self.info("inputFn: %s" % tiFn)
        self.info("inputFn (rel): %s" % _getPath(inputFn))

        ext = pwutils.getExt(inputFn).lower()

        if ext in ['.mrc', '.mrcs']:
            args = ' -InMrc "%s" ' % inputFn
        elif ext in ['.tif', '.tiff']:
            args = ' -InTiff "%s" ' % inputFn
        else:
            raise Exception("Unsupported format: %s" % ext)

        args += ' '.join(['%s %s' % (k, v)
                          for k, v in argsDict.iteritems()])

        if gain:
            args += ' -Gain "%s" ' % gain
            args += ' -RotGain %d' % self.gainRot.get()
            args += ' -FlipGain %d' % self.gainFlip.get()

        if dark:
            args += ' -Dark "%s"' % dark

        args += ' -Gpu %(GPU)s'
        args += ' ' + self.extraParams2.get()

        self.runJob(motioncorr.Plugin.getProgram(), args,
                    cwd=workingFolder,
                    env=motioncorr.Plugin.getEnviron())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        # Check base validation before the specific ones
        errors = ProtTsCorrectMotion._validate(self)

        inputTs = self.inputTiltSeriesM.get()

        if self.doApplyDoseFilter:
            doseFrame = inputTs.getAcquisition().getDosePerFrame()

            if doseFrame == 0.0 or doseFrame is None:
                errors.append('Dose per frame for input movies is 0 or not '
                              'set. You cannot apply dose filter.')

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getMovieLogFile(self, tiltImageM):
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s_0%s-Full.log' % (self._getTiltImageMRoot(tiltImageM),
                                    '-Patch' if usePatches else '')

    def _getRange(self, movie, prefix):

        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, prefix)

        if iniFrame != indxFrame:
            first -= iniFrame
            last -= iniFrame
        else:
            first -= 1
            last -= 1

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
            else:
                return movie.getNumberOfFrames()
        else:
            return movie.getNumberOfFrames()

    # ############## TO-REVIEW ##########################
    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
        The shifts are in pixels irrespective of any binning.
        """
        logPath = self._getExtraPath(self._getMovieLogFile(movie))
        xShifts, yShifts = parseMovieAlignment2(logPath)

        return xShifts, yShifts

    def writeZeroShifts(self, movie):
        # TODO: find another way to do this
        shiftsMd = self._getTmpPath('zero_shifts.xmd')
        pwutils.cleanPath(shiftsMd)
        xshifts = [0] * movie.getNumberOfFrames()
        yshifts = xshifts
        alignment = MovieAlignment(first=1, last=movie.getNumberOfFrames(),
                                   xshifts=xshifts, yshifts=yshifts)
        roiList = [0, 0, 0, 0]
        alignment.setRoi(roiList)
        movie.setAlignment(alignment)
        writeShiftsMovieAlignment(movie, shiftsMd,
                                  1, movie.getNumberOfFrames())
        return shiftsMd

