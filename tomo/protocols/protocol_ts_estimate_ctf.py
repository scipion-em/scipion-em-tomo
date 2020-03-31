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
from pwem import emlib

from .protocol_ts_base import ProtTsProcess


class ProtTsEstimateCTF(ProtTsProcess):
    """
    Base class for estimating the CTF on TiltSeries
    """
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('inputTiltSeries', params.PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
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

        if downFactor != 1:
            # Replace extension by 'mrc' because there are some formats
            # that cannot be written (such as dm3)
            ih.scaleFourier(ti, tiFn, downFactor)
        else:
            ih.convert(ti, tiFn, emlib.DT_FLOAT)

    def _estimateCtf(self, workingDir, tiFn, tiltImage, *args):
        raise Exception("_estimateCTF function should be implemented!")

    def processTiltSeriesStep(self, tsId):
        """ Step called after all CTF are estimated for a given tilt series. """
        self._tsDict.setFinished(tsId)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _getInputTsPointer(self):
        return self.inputTiltSeries

    def _initialize(self):
        """ This function define a dictionary with parameters used
        for CTF estimation that are common for all micrographs. """
        # Get pointer to input micrographs
        inputTs = self._getInputTs()
        acq = inputTs.getAcquisition()
        downFactor = self.getAttributeValue('ctfDownFactor', 1.0)
        sampling = inputTs.getSamplingRate()
        if downFactor != 1.0:
            sampling *= downFactor
        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': acq.getMagnification(),
                        'ampContrast': acq.getAmplitudeContrast(),
                        'samplingRate': sampling,
                        'scannedPixelSize': inputTs.getScannedPixelSize(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        'minDefocus': self.minDefocus.get(),
                        'maxDefocus': self.maxDefocus.get()
                        }

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
