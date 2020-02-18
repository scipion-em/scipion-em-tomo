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
from pwem.emlib.image import ImageHandler, DT_FLOAT

from .protocol_ts_base import ProtTsProcess


class ProtTsEstimateCTF(ProtTsProcess):
    """
    Base class for estimating the CTF on TiltSeries
    """
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        form.addParam('inputTiltSeries', params.PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label=Message.LABEL_INPUT_MIC)
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      #condition='not recalculate',
                      help='Set to 1 for no downsampling. Non-integer downsample '
                           'factors are possible. This downsampling is only used '
                           'for estimating the CTF and it does not affect any '
                           'further calculation. Ideally the estimation of the '
                           'CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) '
                           'and not occupying the whole power spectrum (since '
                           'this downsampling might entail aliasing).')

        self._defineProcessParams(form)

        line = form.addLine('Resolution',
                            #condition='not recalculate',
                            help='Give a value in digital frequency '
                                 '(i.e. between 0.0 and 0.5). These cut-offs '
                                 'prevent the typical peak at the center of the'
                                 ' PSD and high-resolution terms where only '
                                 'noise exists, to interfere with CTF '
                                 'estimation. The default lowest value is 0.05 '
                                 'but for micrographs with a very fine sampling '
                                 'this may be lowered towards 0. The default '
                                 'highest value is 0.35, but it should be '
                                 'increased for micrographs with signals '
                                 'extending beyond this value. However, if '
                                 'your micrographs extend further than 0.35, '
                                 'you should consider sampling them at a finer '
                                 'rate.')
        line.addParam('lowRes', params.FloatParam, default=0.05, label='Lowest')
        line.addParam('highRes', params.FloatParam, default=0.35, label='Highest')
        line = form.addLine('Defocus search range (microns)',
                            #condition='not recalculate',
                            expertLevel=params.LEVEL_ADVANCED,
                            help='Select _minimum_ and _maximum_ values for '
                                 'defocus search range (in microns). Underfocus'
                                 ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=0.25,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=4.,
                      label='Max')

        form.addParam('windowSize', params.IntParam, default=256,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Window size',
                      #condition='not recalculate',
                      help='The PSD is estimated from small patches of this '
                           'size. Bigger patches allow identifying more '
                           'details. However, since there are fewer windows, '
                           'estimations are noisier.')

        form.addParallelSection(threads=2, mpi=1)

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

        if not pw.utils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pw.utils.cleanPath(workingDir)

        ti.setCTF(self.getCtf(ti))

    def _convertInputTi(self, ti, tiFn):
        """ This function will convert the input tilt-image
        taking into account the downFactor.
        It can be overriden in subclasses if another behaviour is required.
        """
        downFactor = self.ctfDownFactor.get()

        ih = ImageHandler()

        if downFactor != 1:
            # Replace extension by 'mrc' because there are some formats
            # that cannot be written (such as dm3)
            ih.scaleFourier(ti, tiFn, downFactor)
        else:
            ih.convert(ti, tiFn, DT_FLOAT)

    def _estimateCtf(self, workingDir, tiFn, tiltImage, *args):
        raise Exception("_estimateCTF function should be implemented!")

    def processTiltSeriesStep(self, tsId):
        """ Step called after all CTF are estimated for a given tiltseries. """
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
                        # Convert from microns to Angstroms
                        'minDefocus': self.minDefocus.get() * 1e+4,
                        'maxDefocus': self.maxDefocus.get() * 1e+4
                        }

    def getCtfParamsDict(self):
        """ Return a copy of the global params dict,
        to avoid overwriting values. """
        return self._params

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
