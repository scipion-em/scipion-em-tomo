# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
import numpy as np
import math
import csv
import pwem.objects as data
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase

MISALIGNED_TS_NAME = "MisalignedTiltSeries"
INTERPOLATED_TS_NAME = "InterpolatedTiltSeries"


class ProtTomoMisalignTiltSeries(EMProtocol, ProtTomoBase):
    """
    Introduce misalignment in the transformation matrix of a tilt-series.
    NOTE: The Interpolated tilt series in this case resembles a not aligned tilt series
    or an aligned one in case you want to apply the inverse of the misalignment
    transformation matrix.
    """
    _label = 'misalign tilt-series '
    _devStatus = BETA
    _possibleOutputs = {MISALIGNED_TS_NAME: tomoObj.SetOfTiltSeries,
                        INTERPOLATED_TS_NAME: tomoObj.SetOfTiltSeries}

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.MisalignedTiltSeries = None
        self.InterpolatedTiltSeries = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        """ Options to introduce misalignment in the X axis shift"""
        form.addParam('shiftXNoiseToggle',
                      params.BooleanParam,
                      default=False,
                      label='Introduce misalignment in shift X?',
                      important=True,
                      help='Introduce noise in the shift alignment value in the X axis. Characterize the noise '
                           'behaviour through the parameters in the following formula:\n'
                           '\n'
                           'dx = a0 + a1 * i + a2 * sin((i + a3) / S * pi) + a4 * sin((i + a5) / S * 2 * pi) '
                           '+ N(0,a6)\n'
                           '\n'
                           'Being i the index position of the image inside the tilt-series, S the size of it and N a '
                           'normal distribution.'
                           'These parameters characterize the following behaviours:\n'
                           '- Constant (a0): an offset error (a0) is introduced in every image of the tilt-series.\n'
                           '- Incremental (a1): a constant incremental error (a1) is propagated through the '
                           'tilt-series.\n'
                           '- Sine lobe (a2, a3): the introduced error presents a half sine shape, characterized by '
                           'the error amplitude (a2) and the phase to displace the error function a given number of '
                           'images inside the tilt-series (a3).\n'
                           '- Sine cycle (a4, a5): the introduced error presents a full sine cycle shape, '
                           'characterized by the error amplitude (a4) and the phase to displace the error function a '
                           'given number of images inside the tilt-series (a5).\n'
                           '- Random (a6): a random error is introduced in every image of the tilt-series given a '
                           'sigma value (a6).\n')

        groupShiftX = form.addGroup('Misalignment parameters in shift X',
                                    condition='shiftXNoiseToggle==True')

        groupShiftX.addParam('a0param',
                             params.FloatParam,
                             default=0.0,
                             label='Offset error (a0)',
                             help='Offset shift error introduced in the X axis for every image of the tilt-series.')

        groupShiftX.addParam('a1param',
                             params.FloatParam,
                             default=0.0,
                             label='Incremental error (a1)',
                             help='Incremental shift error introduced in the X axis for every image of the '
                                  'tilt-series.')

        groupShiftX.addParam('a2param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine lobe error amplitude (a2)',
                             help='Maximum amplitude of the sine lobe error function introduced in the X axis.')

        groupShiftX.addParam('a3param',
                             params.IntParam,
                             default=0,
                             label='Sine lobe error phase (a3)',
                             help='Phase (displacement) of the sine lobe error function. The introduced number '
                                  'corresponds to the number of images from the tilt-series that the origin of the '
                                  'error function will be displaced.')

        groupShiftX.addParam('a4param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine error amplitude (a4)',
                             help='Maximum amplitude of the sine error function introduced in the X axis.')

        groupShiftX.addParam('a5param',
                             params.IntParam,
                             default=0,
                             label='Sine error phase (a5)',
                             help='Phase (displacement) of the sine error function. The introduced number corresponds '
                                  'to the number of images from the tilt-series that the origin of the error function '
                                  'will be displaced.')

        groupShiftX.addParam('a6param',
                             params.FloatParam,
                             default=0.0,
                             label='Random error sigma (a6)',
                             help='Sigma value for the random error introduced in the shift X.')

        """ Options to introduce misalignment in the Y axis shift"""
        form.addParam('shiftYNoiseToggle',
                      params.BooleanParam,
                      default=False,
                      label='Introduce misalignment in shift Y?',
                      important=True,
                      help='Introduce noise in the shift alignment value in the Y axis. Characterize the noise '
                           'behaviour through the parameters in the following formula:\n'
                           '\n'
                           'dY = b0 + b1 * i + b2 * sin((i + b3) / S * pi) + b4 * sin((i + b5) / S * 2 * pi) '
                           '+ N(0,b6)\n'
                           '\n'
                           'Being i the index position of the image inside the tilt-series, S the size of it and N a '
                           'normal distribution.'
                           'These parameters characterize the following behaviours:\n'
                           '- Constant (b0): an offset error (b0) is introduced in every image of the tilt-series.\n'
                           '- Incremental (b1): a constant incremental error (b1) is propagated through the '
                           'tilt-series.\n'
                           '- Sine lobe (b2, b3): the introduced error presents a half sine shape, characterized by '
                           'the error amplitude (b2) and the phase to displace the error function a given number of '
                           'images inside the tilt-series (b3).\n'
                           '- Sine cycle (b4, b5): the introduced error presents a full sine cycle shape, '
                           'characterized by the error amplitude (b4) and the phase to displace the error function a '
                           'given number of images inside the tilt-series (b5).\n'
                           '- Random (b6): a random error is introduced in every image of the tilt-series given a '
                           'sigma value (b6).\n')

        groupShiftY = form.addGroup('Misalignment parameters in shift Y',
                                    condition='shiftYNoiseToggle==True')

        groupShiftY.addParam('b0param',
                             params.FloatParam,
                             default=0.0,
                             label='Offset error (b0)',
                             help='Offset shift error introduced in the Y axis for every image of the tilt-series.')

        groupShiftY.addParam('b1param',
                             params.FloatParam,
                             default=0.0,
                             label='Incremental error (b1)',
                             help='Incremental shift error introduced in the Y axis for every image of the '
                                  'tilt-series.')

        groupShiftY.addParam('b2param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine lobe error amplitude (b2)',
                             help='Maximum amplitude of the sine lobe error function introduced in the Y axis.')

        groupShiftY.addParam('b3param',
                             params.IntParam,
                             default=0,
                             label='Sine lobe error phase (b3)',
                             help='Phase (displacement) of the sine lobe error function. The introduced number '
                                  'corresponds to the number of images from the tilt-series that the origin of the '
                                  'error function will be displaced.')

        groupShiftY.addParam('b4param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine error amplitude (b4)',
                             help='Maximum amplitude of the sine error function introduced in the Y axis.')

        groupShiftY.addParam('b5param',
                             params.IntParam,
                             default=0,
                             label='Sine error phase (b5)',
                             help='Phase (displacement) of the sine error function. The introduced number corresponds '
                                  'to the number of images from the tilt-series that the origin of the error function '
                                  'will be displaced.')

        groupShiftY.addParam('b6param',
                             params.FloatParam,
                             default=0.0,
                             label='Random error sigma (b6)',
                             help='Sigma value for the random error introduced in the shift Y.')

        """ Options to introduce misalignment in the angle"""
        form.addParam('angleNoiseToggle',
                      params.BooleanParam,
                      default=False,
                      label='Introduce misalignment in angle?',
                      important=True,
                      help='Introduce noise in the angle alignment value. Characterize the noise behaviour through the '
                           'parameters in the following formula:\n'
                           '\n'
                           'dA = c0 + c1 * i + c2 * sin((i + c3) / S * pi) + c4 * sin((i + c5) / S * 2 * pi) '
                           '+ N(0,c6)\n'
                           '\n'
                           'Being i the index position of the image inside the tilt-series, S the size of it and N a '
                           'normal distribution.'
                           'These parameters characterize the following behaviours:\n'
                           '- Constant (c0): an offset error (c0) is introduced in every image of the tilt-series.\n'
                           '- Incremental (c1): a constant incremental error (c1) is propagated through the '
                           'tilt-series.\n'
                           '- Sine lobe (c2, c3): the introduced error presents a half sine shape, characterized by '
                           'the error amplitude (c2) and the phase to displace the error function a given number of '
                           'images inside the tilt-series (c3).\n'
                           '- Sine cycle (c4, c5): the introduced error presents a full sine cycle shape, '
                           'characterized by the error amplitude (c4) and the phase to displace the error function a '
                           'given number of images inside the tilt-series (c5).\n'
                           '- Random (c6): a random error is introduced in every image of the tilt-series given a '
                           'sigma value (c6).\n')

        groupAngle = form.addGroup('Misalignment parameters in angle',
                                   condition='angleNoiseToggle==True')

        groupAngle.addParam('c0param',
                            params.FloatParam,
                            default=0.0,
                            label='Offset error (c0)',
                            help='Constant angle error to add for every image of the tilt-series. Angles are measured '
                                 'in degrees.')

        groupAngle.addParam('c1param',
                            params.FloatParam,
                            default=0.0,
                            label='Incremental error (c1)',
                            help='Initial angle error value for the first image (lowest angle) of the tilt-series. '
                                 'Angles are measured in degrees.')

        groupAngle.addParam('c2param',
                            params.FloatParam,
                            default=0.0,
                            label='Sine lobe error amplitude (c2)',
                            help='Maximum amplitude of the sine lobe error function introduced in the angle. Angles '
                                 'are measured in degrees.')

        groupAngle.addParam('c3param',
                            params.IntParam,
                            default=0,
                            label='Sine lobe error phase (c3)',
                            help='Phase (displacement) of the sine lobe error function. The introduced number '
                                 'corresponds to the number of images from the tilt-series that the origin of the '
                                 'error function will be displaced.')

        groupAngle.addParam('c4param',
                            params.FloatParam,
                            default=0.0,
                            label='Sine error amplitude (c4)',
                            help='Maximum amplitude of the sine error function introduced in the angle. Angles are '
                                 'measured in degrees.')

        groupAngle.addParam('c5param',
                            params.IntParam,
                            default=0,
                            label='Sine error phase (c5)',
                            help='Phase (displacement) of the sine error function. The introduced number corresponds '
                                 'to the number of images from the tilt-series that the origin of the error function '
                                 'will be displaced.')

        groupAngle.addParam('c6param',
                            params.FloatParam,
                            default=0.0,
                            label='Random error sigma (c6)',
                            help='Sigma value for random error introduced in the angle. Angles are measured in '
                                 'degrees.')

        """ Options for misalignment interpolation"""
        form.addParam('applyMatrix', params.BooleanParam,
                      default=False,
                      label='Generate interpolated tilt-series',
                      important=True,
                      help='Generate tilt-series applying the'
                           'obtained misalignment transformation matrix.')

        """ Options for aligned tilt series"""
        form.addParam('addInverseMatrix', params.BooleanParam,
                      default=False,
                      label='With inverted matrix',
                      important=False,
                      help='Save the inverse of the misalignment transformation matrix in the interpolated set.',
                      condition='applyMatrix==True')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.introduceRandomMisalignment, ts.getObjId())

            if self.applyMatrix.get():
                self._insertFunctionStep(self.interpolateTiltSeries, ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def introduceRandomMisalignment(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        outputMisalignedSetOfTiltSeries = self.getOutputMisalignedSetOfTiltSeries()
        missAliTs = tomoObj.TiltSeries(tsId=tsId)
        missAliTs.copyInfo(ts)
        outputMisalignedSetOfTiltSeries.append(missAliTs)

        for index, ti in enumerate(ts):
            missAliTi = tomoObj.TiltImage()
            missAliTi.copyInfo(ti, copyId=True)
            missAliTi.setLocation(ti.getLocation())

            if ti.hasTransform():
                transformMat = ti.getTransform().getMatrix()
            else:
                transformMat = np.identity(3)
            newTransformMat = self.modifyTransformMatrix(transformMat, index, ts.getSize(), tsId)

            newTransform = data.Transform()
            newTransform.setMatrix(newTransformMat)
            missAliTi.setTransform(newTransform)

            missAliTs.append(missAliTi)

        outputMisalignedSetOfTiltSeries.update(missAliTs)
        outputMisalignedSetOfTiltSeries.write()

        self._store()

    def interpolateTiltSeries(self, tsObjId):
        missAliTs = self.MisalignedTiltSeries[tsObjId]
        tsId = missAliTs.getTsId()

        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        extraPrefix = self._getExtraPath(tsId)
        outputTsFileName = os.path.join(extraPrefix, "%s_missAli.st" % tsId)

        """Apply the transformation form the input tilt-series"""
        missAliTs.applyTransform(outputTsFileName)

        missAliInterTs = tomoObj.TiltSeries(tsId=tsId)
        missAliInterTs.copyInfo(missAliTs)
        outputInterpolatedSetOfTiltSeries.append(missAliInterTs)

        saveMatrix = self.addInverseMatrix.get()

        for index, tiltImage in enumerate(missAliTs):
            missAliInterTi = tomoObj.TiltImage()

            if saveMatrix:
                # Calculate the inverse of the transformation matrix
                tiltImage.getTransform().invert()

            missAliInterTi.copyInfo(tiltImage, copyId=True, copyTM=saveMatrix)

            missAliInterTi.setLocation(index + 1, outputTsFileName)
            missAliInterTs.append(missAliInterTi)

        missAliInterTs.write()

        outputInterpolatedSetOfTiltSeries.update(missAliInterTs)
        outputInterpolatedSetOfTiltSeries.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def modifyTransformMatrix(self, transformMatrix, index, size, tsId):
        """Shift in X axis modifications"""
        if self.shiftXNoiseToggle.get():
            incrementShiftX = self.a0param.get() + \
                              self.a1param.get() * index + \
                              self.a2param.get() * abs(np.sin((index + self.a3param.get()) / size * np.pi)) + \
                              self.a4param.get() * np.sin((index + self.a5param.get()) / size * 2 * np.pi)

            if self.a6param.get() != 0:
                incrementShiftX += np.random.normal(transformMatrix[0, 2], self.a6param.get())

            transformMatrix[0, 2] += incrementShiftX

        """Shift in Y axis modifications"""
        if self.shiftYNoiseToggle.get():
            incrementShiftY = self.b0param.get() + \
                              self.b1param.get() * index + \
                              self.b2param.get() * abs(np.sin((index + self.b3param.get()) / size * np.pi)) + \
                              self.b4param.get() * np.sin((index + self.b5param.get()) / size * 2 * np.pi)

            if self.b6param.get() != 0:
                incrementShiftY += np.random.normal(transformMatrix[1, 2], self.b6param.get())

            transformMatrix[1, 2] += incrementShiftY

        """Angle modifications"""
        if self.angleNoiseToggle.get():
            oldAngle = np.arccos(transformMatrix[0, 0])

            incrementAngle = self.c0param.get() + \
                             self.c1param.get() * index + \
                             self.c2param.get() * abs(np.sin((index + self.c3param.get()) / size * np.pi)) + \
                             self.c4param.get() * np.sin((index + self.c5param.get()) / size * 2 * np.pi)

            if self.c6param.get() != 0:
                incrementAngle += np.random.normal(oldAngle, self.c6param.get())

            newAngle = oldAngle + math.radians(incrementAngle)

            transformMatrix[0, 0] = np.cos(newAngle)
            transformMatrix[0, 1] = - np.sin(newAngle)
            transformMatrix[1, 0] = np.sin(newAngle)
            transformMatrix[1, 1] = np.cos(newAngle)

        fileName = "TM_misalignment_" + tsId + ".xf"
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)
        filePath = os.path.join(extraPrefix, fileName)

        if 'incrementShiftX' not in locals():
            incrementShiftX = 0
        if 'incrementShiftY' not in locals():
            incrementShiftY = 0
        if 'incrementAngle' not in locals():
            incrementAngle = 0

        vector = [np.cos(math.radians(incrementAngle)),
                  np.sin(math.radians(incrementAngle)),
                  - np.sin(math.radians(incrementAngle)),
                  np.cos(math.radians(incrementAngle)),
                  incrementShiftX,
                  incrementShiftY]

        mode = "a" if os.path.exists(filePath) else "w"

        with open(filePath, mode) as f:
            writer = csv.writer(f, delimiter='\t')

            writer.writerow(vector)

        fileName = "TM_final_" + tsId + ".xf"
        filePath = os.path.join(extraPrefix, fileName)

        vector = [transformMatrix[0, 0],
                  transformMatrix[1, 0],
                  transformMatrix[0, 1],
                  transformMatrix[1, 1],
                  transformMatrix[0, 2],
                  transformMatrix[1, 2]]

        mode = "a" if os.path.exists(filePath) else "w"

        with open(filePath, mode) as f:
            writer = csv.writer(f, delimiter='\t')

            writer.writerow(vector)

        return transformMatrix

    def getOutputMisalignedSetOfTiltSeries(self):
        if not self.MisalignedTiltSeries:
            self.debug("Creating %s output." % MISALIGNED_TS_NAME)
            outputMisalignedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Misaligned')
            outputMisalignedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputMisalignedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(**{MISALIGNED_TS_NAME: outputMisalignedSetOfTiltSeries})
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputMisalignedSetOfTiltSeries)
        return self.MisalignedTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):

        if not self.InterpolatedTiltSeries:
            self.debug("Creating %s output." % INTERPOLATED_TS_NAME)
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(**{INTERPOLATED_TS_NAME: outputInterpolatedSetOfTiltSeries})
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)

        return self.InterpolatedTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if self.MisalignedTiltSeries:
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d."
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.MisalignedTiltSeries.getSize()))

        if self.InterpolatedTiltSeries:
            summary.append("Interpolated Tilt-Series: %d.\n"
                           % self.InterpolatedTiltSeries.getSize())

        if self.AlignedTiltSeries:
            summary.append("Aligned Tilt-Series: %d.\n"
                           % self.AlignedTiltSeries.getSize())

        if len(summary) == 0:
            summary.append("Output not ready yet.")

        return summary

    def _methods(self):
        methods = []
        if self.MisalignedTiltSeries:
            methods.append("New transformation matrices has been calculated for %d Tilt-series."
                           % (self.MisalignedTiltSeries.getSize()))

        if self.InterpolatedTiltSeries:
            methods.append("Also, interpolation has been completed for %d Tilt-series."
                           % self.InterpolatedTiltSeries.getSize())

        if len(methods) == 0:
            methods.append("Output not ready yet.")

        return methods
