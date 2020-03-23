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
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler

from ..objects import TiltSeries, TiltImage
from .protocol_ts_base import ProtTsProcess


class ProtTsCorrectMotion(ProtTsProcess):
    """
    Base class for movie alignment protocols such as:
    motioncorr, crosscorrelation and optical flow

    Alignment parameters are defined in common. For example,
    the frames range used for alignment and final sum, the binning factor
    or the cropping options (region of interest)
    """
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputTiltSeriesM', params.PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      important=True,
                      label='Input Tilt-Series (movies)',
                      help='Select input tiltseries-movies that you want'
                           'to correct for beam-induced motion. ')

        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN',
                             help='Frames range to ALIGN on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie.')
        line.addParam('alignFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', params.IntParam, default=0,
                      label='to')
        group.addParam('useAlignToSum', params.BooleanParam, default=True,
                       label='Use ALIGN frames range to SUM?',
                       help="If *Yes*, the same frame range will be used to "
                            "ALIGN and to SUM. If *No*, you can selected a "
                            "different range for SUM (must be a subset).")
        line = group.addLine('Frames to SUM', condition="not useAlignToSum",
                             help='Frames range to SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to sum, it means that you will sum '
                                  'until the last frame of the movie.')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')
        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        line = group.addLine('Crop offsets (px)',
                             expertLevel=params.LEVEL_ADVANCED)
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')

        line = group.addLine('Crop dimensions (px)',
                             expertLevel=params.LEVEL_ADVANCED,
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, inputId):
        inputTs = self.inputTiltSeriesM.get()
        ih = ImageHandler()

        def _convert(path, tmpFn):
            if path:
                ih.convert(path, self._getTmpPath(tmpFn))

        _convert(inputTs.getGain(), 'gain.mrc')
        _convert(inputTs.getDark(), 'dark.mrc')

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        tiltImageM = self._tsDict.getTi(tsId, tiltImageId)
        workingFolder = self.__getTiltImageMWorkingFolder(tiltImageM)
        pw.utils.makePath(workingFolder)
        self._processTiltImageM(workingFolder, tiltImageM, *args)

        tiFn, tiFnDW = self._getOutputTiltImagePaths(tiltImageM)
        if not os.path.exists(tiFn):
            raise Exception("Expected output file '%s' not produced!" % tiFn)

        if not pw.utils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pw.utils.cleanPath(workingFolder)

    def processTiltSeriesStep(self, tsId):
        """ Create a single stack with the tiltseries. """
        ts = self._tsDict.getTs(tsId)
        ts.setDim([])

        tiList = self._tsDict.getTiList(tsId)
        tiList.sort(key=lambda ti: ti.getTiltAngle())

        ih = ImageHandler()

        tsFn = self._getOutputTiltSeriesPath(ts)
        tsFnDW = self._getOutputTiltSeriesPath(ts, '_DW')

        for i, ti in enumerate(tiList):
            tiFn, tiFnDW = self._getOutputTiltImagePaths(ti)
            newLocation = (i+1, tsFn)
            ih.convert(tiFn, newLocation)
            ti.setLocation(newLocation)
            pw.utils.cleanPath(tiFn)
            if os.path.exists(tiFnDW):
                ih.convert(tiFnDW, (i+1, tsFnDW))
                pw.utils.cleanPath(tiFnDW)

        self._tsDict.setFinished(tsId)

    def _updateOutputSet(self, outputSet, tsIdList):
        """ Override this method to convert the TiltSeriesM into TiltSeries.
        """
        for tsId in tsIdList:
            ts = TiltSeries()
            ts.copyInfo(self._tsDict.getTs(tsId), copyId=True)
            outputSet.append(ts)
            for ti in self._tsDict.getTiList(tsId):
                tiOut = TiltImage(location=ti.getLocation())
                tiOut.copyInfo(ti, copyId=True)
                ts.append(tiOut)

            outputSet.update(ts)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _initialize(self):
        inputTs = self._getInputTs()
        acq = inputTs.getAcquisition()
        gain, dark = self.getGainAndDark()
        self.__basicArgs = [
            acq.getDoseInitial(), acq.getDosePerFrame(), gain, dark]

    def _getArgs(self):
        return self.__basicArgs

    def _getInputTsPointer(self):
        return self.inputTiltSeriesM

    def _processTiltImageM(self, workingFolder, tiltImageM, *args):
        """ This function should be implemented in subclasses to really provide
        the processing step for this TiltSeries Movie.
        Output corrected image (and DW one) should be copied to expected name.
        """
        pass

    def getGainAndDark(self):
        """ Return temporary paths of gain and dark if relevant. """
        inputTs = self.inputTiltSeriesM.get()
        gain = self._getTmpPath('gain.mrc') if inputTs.getGain() else None
        dark = self._getTmpPath('dark.mrc') if inputTs.getDark() else None
        return gain, dark

    def _getFrameRange(self, n, prefix):
        """
        Params:
        :param n: Number of frames of the movies
        :param prefix: what range we want to consider, either 'align' or 'sum'
        :return: (i, f) initial and last frame range
        """
        # In case that the user select the same range for ALIGN and SUM
        # we also use the 'align' prefix
        if self._useAlignToSum():
            prefix = 'align'

        first = self.getAttributeValue('%sFrame0' % prefix)
        last = self.getAttributeValue('%sFrameN' % prefix)

        if first <= 1:
            first = 1

        if last <= 0:
            last = n

        return first, last

    def _getBinFactor(self):
        return self.getAttributeValue('binFactor', 1.0)

    # ----- Some internal functions ---------
    def _getTiltImageMRoot(self, tim):
        return '%s_%02d' % (tim.getTsId(), tim.getObjId())

    def __getTiltImageMWorkingFolder(self, tiltImageM):
        return self._getTmpPath(self._getTiltImageMRoot(tiltImageM))

    def _getOutputTiltImagePaths(self, tiltImageM):
        """ Return expected output path for correct movie and DW one.
        """
        base = self._getExtraPath(self._getTiltImageMRoot(tiltImageM))
        return base + '.mrc', base + '_DW.mrc'

    def _getOutputTiltSeriesPath(self, ts, suffix=''):
        return self._getExtraPath('%s%s.mrcs' % (ts.getTsId(), suffix))

    def _useAlignToSum(self):
        return True


class ProtTsAverage(ProtTsCorrectMotion):
    """
    Simple protocol to average TiltSeries movies as basic
    motion correction. It is used mainly for testing purposes.
    """
    _label = 'average tiltseries'

    def _processTiltImageM(self, workingFolder, tiltImageM, *args):
        """ Simple add all frames and divide by its number. """
        ih = ImageHandler()
        sumImg = ih.createImage()
        img = ih.createImage()

        n = tiltImageM.getNumberOfFrames()
        fn = tiltImageM.getFileName()

        sumImg.read((1, fn))

        for frame in range(2, n + 1):
            img.read((frame, fn))
            sumImg.inplaceAdd(img)

        # sumImg.inplaceDivide(float(n))
        outputFn = self._getOutputTiltImagePaths(tiltImageM)[0]
        sumImg.write(outputFn)
