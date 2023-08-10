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
from os.path import abspath, relpath

import numpy as np

import pyworkflow as pw
import pyworkflow.protocol.params as params
from pyworkflow.object import String, CsvList, Set
from pyworkflow.utils import removeBaseExt
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler

from ..objects import TiltSeries, TiltImage, SetOfTiltSeries
from .protocol_ts_base import ProtTsProcess

OUTPUT_TILT_SERIES_ODD = 'TiltSeriesOdd'
OUTPUT_TILT_SERIES_EVEN = 'TiltSeriesEven'
EVEN = 'even'
ODD = 'odd'
OUTPUT_TILT_SERIES_DW = 'TiltSeriesDW'


class ProtTsCorrectMotion(ProtTsProcess):
    """
    Base class for movie alignment protocols such as:
    motioncorr, crosscorrelation and optical flow

    Alignment parameters are defined in common. For example,
    the frames range used for alignment and final sum, the binning factor
    or the cropping options (region of interest)
    """
    _possibleOutputs = {'outputTiltSeries': SetOfTiltSeries,
                        OUTPUT_TILT_SERIES_DW: SetOfTiltSeries,
                        OUTPUT_TILT_SERIES_EVEN: SetOfTiltSeries,
                        OUTPUT_TILT_SERIES_ODD: SetOfTiltSeries}

    # Even / odd functionality
    evenOddCapable = False
    TiltSeriesOdd = None
    TiltSeriesEven = None
    TiltSeriesDW = None
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form, addEvenOddParam = True):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputTiltSeriesM', params.PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      important=True,
                      label='Input Tilt-Series (movies)',
                      help='Select input tilt-series movies that you want'
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

        if self.evenOddCapable and addEvenOddParam:
            form.addParam('splitEvenOdd', params.BooleanParam,
                          default=False,
                          label='Split & sum odd/even frames?',
                          help='(Used for denoising data preparation). If set to Yes, 2 additional movies/tilt '
                               'series will be generated, one generated from the even frames and the other from the '
                               'odd ones using the same alignment for the whole stack of frames.')

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

    def _doSplitEvenOdd(self):
        """ Returns if even/odd stuff has to be done"""
        if not self.evenOddCapable:
            return False
        else:
            return self.splitEvenOdd.get()

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        tiltImageM = self._tsDict.getTi(tsId, tiltImageId)
        workingFolder = self.__getTiltImageMWorkingFolder(tiltImageM)
        pw.utils.makePath(workingFolder)
        self._processTiltImageM(workingFolder, tiltImageM, *args)

        if self._doSplitEvenOdd():
            baseName = removeBaseExt(tiltImageM.getFileName())
            evenName = abspath(self._getExtraPath(baseName + '_avg_' + EVEN))
            oddName = abspath(self._getExtraPath(baseName + '_avg_' + ODD))
            alignedFrameStack = self._getExtraPath(baseName + '_aligned_movie.mrcs')
            # Get even/odd xmd files
            args = '--img %s ' % abspath(alignedFrameStack)
            args += '--type frames '
            args += '-o %s ' % (evenName + '.xmd')
            args += '-e %s ' % (oddName + '.xmd')
            args += '--sum_frames '
            self.runJob('xmipp_image_odd_even', args)

            # Update even and odd average lists
            oddfn = oddName + '_aligned.mrc'
            evenfn = evenName + '_aligned.mrc'

            tiltImageM.setOddEven([oddfn, evenfn])
            pw.utils.cleanPath(alignedFrameStack)

        tiFn, tiFnDW = self._getOutputTiltImagePaths(tiltImageM)
        if not os.path.exists(tiFn):
            raise Exception("Expected output file '%s' not produced!" % tiFn)

        if not pw.utils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pw.utils.cleanPath(workingFolder)

    def createDWTs(self, ts):
        """ Dose weighting creation"""

        if self._createOutputWeightedTS():
            if self.TiltSeriesDW is None:
                self.TiltSeriesDW = self._createOutputSet(suffix='_dose-weighted')
                self.TiltSeriesDW.setSamplingRate(self._getOutputSampling())
            else:
                self.TiltSeriesDW.enableAppend()

            tsObjDW = TiltSeries()
            tsObjDW.copyInfo(ts, copyId=True)
            self.TiltSeriesDW.append(tsObjDW)

            tsFnDW = self._getDWTiltSeriesPath(ts)

            for i, ti in enumerate(ts):
                tiOut = TiltImage(location=(i+1, tsFnDW))
                tiOut.copyInfo(ti, copyId=True)
                tiOut.setAcquisition(ti.getAcquisition())
                tiOut.setSamplingRate(self._getOutputSampling())
                tiOut.setObjId(ti.getIndex())
                tsObjDW.append(tiOut)

            self.TiltSeriesDW.update(tsObjDW)

    def createOddEvenTs(self, ts, odd=True):
        def addTiltImage(tiLocation, tsObject, mainti, tsIde, samplingRate):
            """
            :param tiLocation: Location of the aligned tilt image in the stack
            :param tsObject: Tilt Series to which the new Ti Image will be added
            :param mainti: Tilt Series Movies object
            :param tsIde: Tilt series identifier
            :param samplingRate: current Tilt Series sampling rate
            """
            ta = mainti.getTiltAngle()
            to = mainti.getAcquisitionOrder()
            acq = mainti.getAcquisition()
            ti = TiltImage(tiltAngle=ta, tsId=tsIde, acquisitionOrder=to)
            ti.setSamplingRate(samplingRate)
            ti.setAcquisition(acq)
            index, fname = tiLocation.split("@")
            ti.setLocation(int(index), fname)
            tsObject.append(ti)
            pw.utils.cleanPath(tiLocation)

        if odd:
            suffix = ODD
            outputAttr = OUTPUT_TILT_SERIES_ODD
            oddEvenIndex = 0
        else:
            suffix = EVEN
            outputAttr = OUTPUT_TILT_SERIES_EVEN
            oddEvenIndex = 1

        template = 'tiltseries%s.sqlite'
        sRate = self._getOutputSampling()

        output = getattr(self, outputAttr, None)
        if output:
            output.enableAppend()
        else:
            output = SetOfTiltSeries.create(self._getPath(), template=template, suffix=suffix)
            setattr(self, outputAttr, output)
            output.setSamplingRate(sRate)

        tsObj = TiltSeries()
        tsObj.copyInfo(ts, copyId=True)
        output.append(tsObj)

        tsId = ts.getTsId()
        for ti in ts:
            fnImg = ti.getOddEven()[oddEvenIndex]
            addTiltImage(fnImg, tsObj, ti, tsId, sRate)

        # update items and size info
        output.update(tsObj)

    def processTiltSeriesStep(self, tsId):
        """ Create a single stack with the tiltseries. """
        def createStack(tiList, tsFn, fngetter, locationSetter=None):
            """ This function creates a stack from individual images """
            for i, ti in enumerate(tiList):
                tiFn = fngetter(ti)
                #newLocation = (i + 1, tsFn)
                #ih.convert(tiFn, newLocation)
                #pw.utils.cleanPath(tiFn)
                if os.path.exists(tiFn):
                    newLocation = (i + 1, tsFn)
                    ih.convert(tiFn, newLocation)
                    pw.utils.cleanPath(tiFn)
                if locationSetter:
                    locationSetter(newLocation, ti)

        ts = self._tsDict.getTs(tsId)
        ts.setDim([])

        tiList = self._tsDict.getTiList(tsId)
        tiList.sort(key=lambda ti: ti.getTiltAngle())

        ih = ImageHandler()

        tsFn = self._getOutputTiltSeriesPath(ts)

        # Merge all micrographs from the same tilt images in a single "mrcs" stack file
        createStack(tiList, tsFn, self._getOutputTiltImagePath, locationSetter=lambda newloc, ti: ti.setLocation(newloc))

        # Dose weighted
        if self._createOutputWeightedTS():

            createStack(tiList, self._getDWTiltSeriesPath(ts), self._getOutputTiltImageDWPath)

        if self._doSplitEvenOdd():
            tsFnOdd = self._getOutputTiltSeriesPath(ts, '_odd')
            tsFnEven = self._getOutputTiltSeriesPath(ts, '_even')
            createStack(tiList, tsFnOdd, self._getOutputTiltImageOddPath,
                        locationSetter=lambda newloc, ti: ti.setOdd(ih.locationToXmipp(newloc)))
            createStack(tiList, tsFnEven, self._getOutputTiltImageEvenPath,
                        locationSetter=lambda newloc, ti: ti.setEven(ih.locationToXmipp(newloc)))

        self._tsDict.setFinished(tsId)

    def _getDWTiltSeriesPath(self, ts:TiltSeries):

        return self._getOutputTiltSeriesPath(ts, '_DW')

    def _updateOutput(self, tsIdList):
        """ Update the output set with the finished Tilt-series.
        Params:
            :param tsIdList: list of ids of finished tasks.
        """

        def writeAndStore(obj):
            obj.write()
            self._store(obj)

        # Flag to check the first time we save output
        self._createOutput = getattr(self, '_createOutput', True)

        outputSet = self._getOutputSet()

        if outputSet is None:
            # Special case just to update the outputSet status
            # but it only makes sense when there is outputSet
            if not tsIdList:
                return
            outputSet = self._createOutputSet()
        else:
            outputSet.enableAppend()
            self._createOutput = False

        # Call the sub-class method to update the output
        outputSet.setSamplingRate(self._getOutputSampling())
        self._updateOutputSet(outputSet, tsIdList)
        outputSet.setStreamState(Set.STREAM_OPEN)

        if self._doSplitEvenOdd():
            self.TiltSeriesEven.setStreamState(Set.STREAM_OPEN)
            self.TiltSeriesOdd.setStreamState(Set.STREAM_OPEN)

        if self._createOutputWeightedTS():
            self.TiltSeriesDW.setStreamState(Set.STREAM_OPEN)

        if self._createOutput:
            outputSet.updateDim()
            outputs = {self._getOutputName(): outputSet}
            if self._createOutputWeightedTS():
                self.TiltSeriesDW.updateDim()
                outputs.update({OUTPUT_TILT_SERIES_DW: self.TiltSeriesDW})

            if self._doSplitEvenOdd():
                self.TiltSeriesEven.updateDim()
                self.TiltSeriesOdd.updateDim()
                outputs.update({OUTPUT_TILT_SERIES_EVEN: self.TiltSeriesEven,
                                OUTPUT_TILT_SERIES_ODD: self.TiltSeriesOdd
                                })
                self._defineOutputs(**outputs)
                self._defineSourceRelation(self._getInputTsPointer(), self.TiltSeriesEven)
                self._defineSourceRelation(self._getInputTsPointer(), self.TiltSeriesOdd)
            else:
                self._defineOutputs(**outputs)

            if self._createOutputWeightedTS():
                self._defineSourceRelation(self._getInputTsPointer(), self.TiltSeriesDW)

            self._defineSourceRelation(self._getInputTsPointer(), outputSet)
            self._createOutput = False
        else:
            writeAndStore(outputSet)
            if self._doSplitEvenOdd():
                writeAndStore(self.TiltSeriesEven)
                writeAndStore(self.TiltSeriesOdd)
            if self._createOutputWeightedTS():
                writeAndStore(self.TiltSeriesDW)

        outputSet.close()
        if self._doSplitEvenOdd():
            self.TiltSeriesEven.close()
            self.TiltSeriesOdd.close()

        if self._createOutputWeightedTS():
            self.TiltSeriesDW.close()

        if self._tsDict.allDone():
            self._coStep.setStatus(params.STATUS_NEW)

    def _updateOutputSet(self, outputSet, tsIdList):
        """ Override this method to convert the TiltSeriesM into TiltSeries.
        """
        for tsId in tsIdList:
            ts = TiltSeries()
            ts.copyInfo(self._tsDict.getTs(tsId), copyId=True)
            ts.setSamplingRate(self._getOutputSampling())
            outputSet.append(ts)
            tList = self._tsDict.getTiList(tsId)
            ind = np.argsort([ti.getTiltAngle() for ti in tList])
            counter = 1

            for i in ind:  # Make each row of the sqlite file be sorted by
                # index after having been sorted by angle previously, in order to avoid tilt image mismatching in
                # another operations, such as the fiducial alignment, which expects the sqlite to be sorted that way
                ti = tList[i]
                tiOut = TiltImage(location=(counter, ti.getFileName()))
                tiOut.copyInfo(ti, copyId=True)
                tiOut.setAcquisition(ti.getAcquisition())
                tiOut.setSamplingRate(self._getOutputSampling())
                tiOut.setObjId(ti.getIndex())
                ts.append(tiOut)
                counter += 1

            outputSet.update(ts)

            # Create dose weighted set
            self.createDWTs(ts)

            # Even and odd stuff
            if self._doSplitEvenOdd():
                self.createOddEvenTs(ts, True)
                self.createOddEvenTs(ts, False)

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

    def _getOutputSampling(self):
        return self.inputTiltSeriesM.get().getSamplingRate() * self._getBinFactor()

    def _processTiltImageM(self, workingFolder, tiltImageM, *args):
        """ This function should be implemented in subclasses to really provide
        the processing step for this TiltSeries Movie.
        Output corrected image (and DW one) should be copied to expected name.
        """
        pass

    def getGainAndDark(self):
        """ Return temporary paths of gain and dark if relevant. """
        inputTs = self.inputTiltSeriesM.get()
        gain = os.path.abspath(self._getTmpPath('gain.mrc')) if inputTs.getGain() else None
        dark = os.path.abspath(self._getTmpPath('dark.mrc')) if inputTs.getDark() else None
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
        return self._getOutputTiltImagePath(tiltImageM), self._getOutputTiltImageDWPath(tiltImageM)

    def _getOutputTiltImagePath(self, tiltImageM):
        """ Return the main path for the tilt image corrected movie.
        """
        return self._getExtraPath(self._getTiltImageMRoot(tiltImageM)) + '.mrc'

    def _getOutputTiltImageDWPath(self, tiltImageM):
        """ Return the main path for the tilt image dose weighted corrected movie.
        """
        return self._getExtraPath(self._getTiltImageMRoot(tiltImageM)) + '_DW.mrc'

    def _getOutputTiltImageOddPath(self, tiltImageM):
        """ Return the path for the odd tilt image corrected movie.
        """
        return tiltImageM.getOdd()

    def _getOutputTiltImageEvenPath(self, tiltImageM):
        """ Return the path for the even tilt image corrected movie.
        """
        return tiltImageM.getEven()

    def _getOutputTiltSeriesPath(self, ts, suffix=''):
        return self._getExtraPath('%s%s.mrcs' % (ts.getTsId(), suffix))

    def _useAlignToSum(self):
        return True

    def _createOutputWeightedTS(self):
        return False


class ProtTsAverage(ProtTsCorrectMotion):
    """
    Simple protocol to average TiltSeries movies as basic
    motion correction. It is used mainly for testing purposes.
    """
    _label = 'average tilt-series movies'
    _devStatus = pw.BETA

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
