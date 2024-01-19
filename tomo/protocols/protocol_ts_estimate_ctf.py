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
from pwem.convert.headers import getFileFormat, MRC
from pyworkflow.object import Set, Integer
from pyworkflow.protocol import STATUS_NEW
from pyworkflow.utils.properties import Message
from pwem import emlib

from .protocol_ts_base import ProtTsProcess
from ..objects import SetOfCTFTomoSeries, CTFTomoSeries


class ProtTsEstimateCTF(ProtTsProcess):
    """
    Base class for estimating the CTF on TiltSeries
    """

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('inputTiltSeries', params.PointerParam, important=True,
                      pointerClass='SetOfTiltSeries, SetOfCTFTomoSeries',
                      label='Input tilt series')

        self._defineProcessParams(form)
        self._defineStreamingParams(form)

    def _defineProcessParams(self, form):
        """ Should be implemented in subclasses. """
        pass

    # --------------------------- STEPS functions ----------------------------

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        ti = self._tsDict.getTi(tsId, tiltImageId)

        # Create working directory for a tilt-image
        workingDir = self.__getTiworkingDir(ti)
        tiFnMrc = os.path.join(workingDir, self.getTiPrefix(ti) + '.mrc')
        pw.utils.makePath(workingDir)
        self._convertInputTi(ti, tiFnMrc)

        # Call the current estimation of CTF that is implemented in subclasses
        self._estimateCtf(workingDir, tiFnMrc, ti, *args)

        if not pw.utils.envVarOn(pw.SCIPION_DEBUG_NOCLEAN):
            pw.utils.cleanPath(workingDir)

        ti.setCTF(self.getCtf(ti))

    def _convertInputTi(self, ti, tiFn):
        """ This function will convert the input tilt-image
        taking into account the downFactor.
        It can be overwritten in subclasses if another behaviour is required.
        """
        downFactor = self.ctfDownFactor.get()
        ih = emlib.image.ImageHandler()

        if not ih.existsLocation(ti):
            raise Exception("Missing input file: %s" % ti)

        tiFName = ti.getFileName()
        # Make xmipp considers the input object as TS to work as expected
        if getFileFormat(tiFName) == MRC:
            tiFName = tiFName.split(':')[0] + ':mrcs'
        tiFName = str(ti.getIndex()) + '@' + tiFName

        if downFactor != 1:
            ih.scaleFourier(tiFName, tiFn, downFactor)
        else:
            ih.convert(tiFName, tiFn, emlib.DT_FLOAT)

    def _estimateCtf(self, workingDir, tiFn, tiltImage, *args):
        raise Exception("_estimateCTF function should be implemented!")

    def processTiltSeriesStep(self, tsId):
        """ Step called after all CTF are estimated for a given tilt series. """
        self._tsDict.setFinished(tsId)

    def _updateOutput(self, tsIdList):
        """ Update the output set with the finished Tilt-series.
        Params:
            :param tsIdList: list of ids of finished tasks.
        """
        ts = self._getTiltSeries(tsIdList[0])
        tsId = ts.getTsId()
        objId = ts.getObjId()
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

        newCTFTomoSeries = CTFTomoSeries()
        newCTFTomoSeries.copyInfo(ts)
        newCTFTomoSeries.setTiltSeries(ts)
        newCTFTomoSeries.setTsId(tsId)
        newCTFTomoSeries.setObjId(objId)

        outputSet.append(newCTFTomoSeries)

        index = 1
        for ti in self._tsDict.getTiList(tsId):
            newCTFTomo = ti._ctfModel
            newCTFTomo.setIndex(Integer(index))
            index += 1
            newCTFTomoSeries.append(newCTFTomo)

        newCTFTomoSeries.write(properties=False)
        outputSet.update(newCTFTomoSeries)

        if self._createOutput:
            self._defineOutputs(**{self._getOutputName(): outputSet})
            self._defineSourceRelation(self._getInputTs(pointer=True),
                                       outputSet)
            self._createOutput = False
        else:
            outputSet.write()
            self._store(outputSet)

        outputSet.close()
        self._store()

        if self._tsDict.allDone():
            self._coStep.setStatus(STATUS_NEW)

    def createOutputStep(self):
        outputSet = self._getOutputSet()
        outputSet.setStreamState(outputSet.STREAM_CLOSED)
        outputSet.write()
        self._store(outputSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfCTFTomoSeries'):
            inputLabel = 'Input Tilt-Series'
            if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
                inputLabel = 'Input CTFTomoSeries'
            summary.append(
                inputLabel + ": %d.\nNumber of CTF estimated: %d.\n"
                % (self._getInputTs().getSize(),
                   self.outputSetOfCTFTomoSeries.getSize()))
        else:
            summary.append("Output CTFs are not ready yet.")
        return summary

    # --------------------------- UTILS functions ----------------------------
    def _createOutputSet(self, suffix=''):
        """ Method to create the output set.
        By default will a SetOfTiltSeries, but can be re-defined in subclasses.
        """
        outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                             template='CTFmodels%s.sqlite')
        outputSetOfCTFTomoSeries.setSetOfTiltSeries(self._getInputTs(pointer=True))
        outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
        return outputSetOfCTFTomoSeries

    def _getTiltSeries(self, itemId):
        obj = None
        inputSetOfTiltseries = self._getInputTs()
        for item in inputSetOfTiltseries.iterItems(iterate=False):
            if item.getTsId() == itemId:
                obj = item
                if isinstance(obj, CTFTomoSeries):
                    obj = item.getTiltSeries()
                break

        if obj is None:
            raise ("Could not find tilt-series with tsId = %s" % itemId)

        return obj

    def _getOutputName(self):
        """ Return the output name, by default 'outputTiltSeries'.
        This method can be re-implemented in subclasses that have a
        different output. (e.g outputTomograms).
        """
        return 'outputSetOfCTFTomoSeries'

    def _getOutputSet(self):
        return getattr(self, self._getOutputName(), None)

    def _getInputTsPointer(self):
        return self.inputTiltSeries

    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def _initialize(self):
        """ This function define a dictionary with parameters used
        for CTF estimation that are common for all micrographs. """
        # Get pointer to input micrographs
        inputTs = self._getInputTs()
        downFactor = self.getAttributeValue('ctfDownFactor', 1.0)
        self._params = createCtfParams(inputTs, self.windowSize.get(), self.lowRes.get(), self.highRes.get(),
                                       self.minDefocus.get(), self.maxDefocus.get(), downFactor=downFactor)

    def getCtfParamsDict(self):
        """ Return a copy of the global params dict,
        to avoid overwriting values. """
        return self._params

    def getCtf(self, ti):
        """ Should be implemented in subclasses. """
        pass

    # ----- Some internal functions ---------
    def getTiRoot(self, tim):
        return '%s_%02d' % (tim.getTsId(), tim.getObjId())

    def __getTiworkingDir(self, tiltImage):
        return self._getTmpPath(self.getTiRoot(tiltImage))

    def _getOutputTiPaths(self, tiltImageM):
        """ Return expected output path for correct movie and DW one.
        """
        base = self._getExtraPath(self.getTiRoot(tiltImageM))
        return base + '.mrc', base + '_DW.mrc'

    def _getOutputTiltSeriesPath(self, ts, suffix=''):
        return self._getExtraPath('%s%s.mrcs' % (ts.getTsId(), suffix))

    def _useAlignToSum(self):
        return True

    def getTiPrefix(self, ti):
        return '%s_%03d' % (ti.getTsId(), ti.getObjId())

    def allowsDelete(self, obj):
        return True


def createCtfParams(inputTs, windowSize, lowRes, highRes, minDefocus, maxDefocus, downFactor=1.0):
    """ This function define a dictionary with parameters used
    for CTF estimation that are common for all micrographs. """
    # Get pointer to input micrographs
    acq = inputTs.getAcquisition()
    sampling = inputTs.getSamplingRate()
    if downFactor != 1.0:
        sampling *= downFactor
    return {'voltage': acq.getVoltage(),
            'sphericalAberration': acq.getSphericalAberration(),
            'magnification': acq.getMagnification(),
            'ampContrast': acq.getAmplitudeContrast(),
            'samplingRate': sampling,
            'scannedPixelSize': inputTs.getScannedPixelSize(),
            'windowSize': windowSize,
            'lowRes': lowRes,
            'highRes': highRes,
            'minDefocus': minDefocus,
            'maxDefocus': maxDefocus
            }
