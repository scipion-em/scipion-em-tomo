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
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params
from pyworkflow.protocol import STEPS_PARALLEL

from pyworkflow.utils.properties import Message

from tomo.objects import TiltSeriesDict
from .protocol_base import ProtTomoBase


class ProtTsEstimateCTF(pwem.EMProtocol, ProtTomoBase):
    """
    Base class for estimating the CTF on TiltSeries
    """
    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

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

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._defineCtfParamsDict()
        inputTs = self.inputTiltSeries.get()
        ciStepId = self._insertFunctionStep('convertInputStep',
                                            inputTs.getObjId())
        self._tsDict = TiltSeriesDict()
        allSteps = []

        for ts in inputTs:
            tsId = ts.getTsId()
            self._tsDict.addTs(ts)
            tsAllSteps = []
            for i, ti in enumerate(ts):
                self._tsDict.addTi(ti)
                tiStepId = self._insertFunctionStep('estimateCtfStep',
                                                    tsId, ti.getObjId(),
                                                    *self._getArgs(),
                                                    prerequisites=[ciStepId])
                tsAllSteps.append(tiStepId)

            tsStepId = self._insertFunctionStep('processTsStep', tsId,
                                                prerequisites=tsAllSteps)
            allSteps.append(tsStepId)

        self._insertFunctionStep('createOutputStep',
                                 prerequisites=allSteps)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, inputId):
        pass

    def estimateCtfStep(self, tsId, tiltImageId, *args):
        tiltImage = self._tsDict.getTi(tsId, tiltImageId)

        # Create working directory for a tilt-image
        workingFolder = self.__getTiWorkingFolder(tiltImage)
        pw.utils.makePath(workingFolder)

        # Call the current estimation of CTF that is implemented in subclasses
        self._estimateCtf(workingFolder, tiltImage, *args)

        if not pw.utils.envVarOn('SCIPION_DEBUG_NOCLEAN'):
            pw.utils.cleanPath(workingFolder)

    def _estimateCtf(self, workingFolder, tiltImage):
        raise Exception("_estimateCTF function should be implemented!")

    def processTsStep(self, tsId):
        """ Step called after all CTF are estimated for a given tiltseries. """
        pass

    def createOutputStep(self):
        inputTs = self.inputTiltSeries.get()
        outputTs = self._createSetOfTiltSeries()
        outputTs.copyInfo(inputTs)

        def _updateCtf(j, ts, ti, tsOut, tiOut):
            tiOut.setCTF(self.getCtf(ti))

        outputTs.copyItems(inputTs, updateTiCallback=_updateCtf)

        self._defineOutputs(outputTiltSeries=outputTs)
        self._defineSourceRelation(self.inputTiltSeries, outputTs)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _defineCtfParamsDict(self):
        """ This function define a dictionary with parameters used
        for CTF estimation that are common for all micrographs. """
        # Get pointer to input micrographs
        inputTs = self.inputTiltSeries.get()
        acq = inputTs.getAcquisition()
        sampling = inputTs.getSamplingRate()
        downFactor = self.getAttributeValue('ctfDownFactor', 1.0)
        if downFactor != 1.0:
            sampling *= downFactor
        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': acq.getMagnification(),
                        'ampContrast': acq.getAmplitudeContrast(),
                        'samplingRate': sampling,
                        #'scannedPixelSize': inputTs.getScannedPixelSize(),
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

    def _getArgs(self):
        """ Return a list with parameters that will be passed to the process
        TiltSeries step. It can be redefined by subclasses.
        """
        return []

    # ----- Some internal functions ---------
    def getTiRoot(self, tim):
        return '%s_%02d' % (tim.getTsId(), tim.getObjId())

    def __getTiWorkingFolder(self, tiltImage):
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

