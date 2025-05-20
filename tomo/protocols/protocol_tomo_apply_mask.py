# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * National Center of Biotechnology, CSIC, Spain
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
import logging
from enum import Enum
from os import remove
from typing import Union
import mrcfile
import numpy as np
from scipy.ndimage import gaussian_filter, binary_dilation
from pwem.emlib.image.image_readers import MRCImageReader, ImageReadersRegistry, ImageStack
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String, Set
from pyworkflow.protocol import PointerParam, BooleanParam, IntParam, STEPS_PARALLEL, GE
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTomograms, SetOfTomoMasks, Tomogram

logger = logging.getLogger(__name__)
IN_TOMO_SET = 'inTomoSet'
IN_MASK_SET = 'inMaskSet'


class TomoApplyMaskOutputs(Enum):
    maskedTomograms = SetOfTomograms


class ProtTomoApplyMask(EMProtocol):
    """This protocol applies a set of masks to a given set of tomograms. The protocol
    will try to match the tomograms and the masks by tsId. Once the mask/s are applied."""

    _label = 'apply tomomasks to tomograms'
    _devStatus = BETA
    _possibleOutputs = TomoApplyMaskOutputs
    stepsExecutionMode = STEPS_PARALLEL
    _sRateTol = 1e-3  # Angstrom/px

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomosDict = None
        self.tomoMaskDict = None
        self.doSmooth = False
        self.nonMatchingTsIdsMsg = String()
        self.failedApixTsIds = []
        self.failedDimsTsIds = []

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMO_SET, PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')
        form.addParam(IN_MASK_SET, PointerParam,
                      pointerClass="VolumeMask, SetOfTomoMasks",
                      important=True,
                      label='Masks',
                      help='The protocol will try to match the tomograms and the masks by tsId.')
        group = form.addGroup('Mask operations')
        group.addParam('invertMask', BooleanParam,
                       default=False,
                       label='Invert mask?')
        group.addParam('dilationPixels', IntParam,
                       default=0,
                       label='Number of pixels for dilation',
                       validators=[GE(0)],
                       haelp='The dilation will expands the shape by the given number of pixels '
                             'in all directions.')
        group.addParam('sigmaGaussian', IntParam,
                       default=3,
                       label='Std for gaussian smoothing',
                       validators=[GE(0)],
                       help='A gaussian filter can be applied by setting this parameter with a value greater than 0. '
                            'It can be used to smooth the borders, providing a more continuous transition between the '
                            'mask and the background. It helps to avoid undesired mathematical artifacts when '
                            'processing later the resulting masked tomograms.\n\nThis parameter (named sigma) '
                            'determines how much neighboring pixels influence each other during smoothing. '
                            'A small value of sigma means a sharper image, less smoothing, '
                            'while a larger value of sigma means a blurrier image, more smoothing.')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeStepDeps = []
        for tsId in self.tomosDict.keys():
            smoothId = self._insertFunctionStep(self.processMaskStep, tsId,
                                                prerequisites=[],
                                                needsGPU=False)
            aMId = self._insertFunctionStep(self.applyMaskStep, tsId,
                                            prerequisites=smoothId,
                                            needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=aMId,
                                              needsGPU=False)
            closeStepDeps.append(cOutId)
        self._insertFunctionStep(self.closeOutputSetStep,
                                 prerequisites=closeStepDeps,
                                 needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        inTomos = self._getInTomoSet()
        inTomoMasks = self._getInTomoMasks()
        tomosTsIds = set(inTomos.getTSIds())
        tomoMaskTsIds = set(inTomoMasks.getTSIds())
        # Check the common elements
        matchingTsIds = tomosTsIds & tomoMaskTsIds
        nonMatchingTsIds = tomosTsIds ^ tomoMaskTsIds
        if not matchingTsIds:
            raise Exception('No matching tsIds were found among the given sets of tomograms and tomo masks.')
        if nonMatchingTsIds:
            msg = f'Some non-matching tsIds were found: {nonMatchingTsIds}'
            self.nonMatchingTsIdsMsg.set(msg)
            logger.info(cyanStr(msg))
            self._store(self.nonMatchingTsIdsMsg)
        if self.sigmaGaussian.get() > 0:
            self.doSmooth = True
        self.tomosDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos
                          if tomo.getTsId() in matchingTsIds}
        self.tomoMaskDict = {tomoMask.getTsId(): tomoMask.clone() for tomoMask in inTomoMasks
                             if tomoMask.getTsId() in matchingTsIds}

    def processMaskStep(self, tsId: str):
        if self.doSmooth:
            logger.info(cyanStr(f'tsId = {tsId}: processing the mask...'))
            mask = self.tomoMaskDict[tsId]
            maskFileName = mask.getFileName()
            # Read the mask
            with mrcfile.mmap(maskFileName, mode='r', permissive=True) as mrc:
                # Invert if required
                data = 1 - mrc.data if self.invertMask.get() else mrc.data
            # Dilate the mask
            dilationPixels = self.dilationPixels.get()
            if dilationPixels > 0:
                data = np.array(data, dtype=bool)  # Required for the binary dilation
                data = binary_dilation(data, iterations=dilationPixels)
            # Smooth the mask
            data = np.array(data, dtype=float)  # Required to be cast from uint8 to float for the gaussian filtering
            smoothData = gaussian_filter(data, sigma=self.sigmaGaussian.get())
            # Write the result
            smoothMaskFn = self._getSmoothedMaskFn(tsId)
            with mrcfile.new_mmap(smoothMaskFn, overwrite=True, shape=smoothData.shape,
                                  mrc_mode=2) as mrc:  # Mode 2 is float32 (see new_mmap)
                for i in range(len(smoothData)):
                    mrc.data[i, :, :] = smoothData[i, :, :]
                mrc.update_header_from_data()
                mrc.voxel_size = mask.getSamplingRate()

    def applyMaskStep(self, tsId: str):
        logger.info(cyanStr(f'tsId = {tsId}: applying the mask...'))
        mask = self.tomoMaskDict[tsId]
        tomo = self.tomosDict[tsId]
        # Check tomo by tomo (to cover heterogeneous sets) both the sampling rate (checked out only at set level
        # in the _validate) and the dimensions
        maskSRate = mask.getSamplingRate()
        tomoSRate = tomo.getSamplingRate()
        if abs(tomoSRate - maskSRate) > self._sRateTol:
            self.failedApixTsIds.append(tsId)
        else:
            maskFileName = self._getSmoothedMaskFn(tsId) if self.doSmooth else mask.getFileName()
            tomoFileName = tomo.getFileName()
            maskDims = MRCImageReader.getDimensions(maskFileName)
            tomoDims = MRCImageReader.getDimensions(tomoFileName)
            if not np.allclose(np.array(maskDims), np.array(tomoDims)):
                self.failedDimsTsIds.append(tsId)
            else:
                maskStack = ImageReadersRegistry.open(maskFileName)
                tomoStack = ImageReadersRegistry.open(tomoFileName)
                resultingImgList = [np.multiply(maskSlice, tomoSlice) for
                                    maskSlice, tomoSlice in zip(maskStack, tomoStack)]
                resultingStack = ImageStack(resultingImgList)
                MRCImageReader.write(resultingStack, self._getResultFn(tsId), samplingRate=tomoSRate)
                # Remove the smoothed mask from the protocol tmp directory to avoid the storage of multiple
                # big temporal files in execution time at once
                if self.doSmooth:
                    remove(self._getSmoothedMaskFn(tsId))

    def createOutputStep(self, tsId: str):
        if not (tsId in self.failedApixTsIds or tsId in self.failedDimsTsIds):
            with self._lock:
                logger.info(cyanStr(f'tsId = {tsId}: registering the output...'))
                currentTomo = self.tomosDict[tsId]
                outputTomos = self._getOutTomos()
                tomo = Tomogram()
                tomo.copyInfo(currentTomo)
                tomo.setFileName(self._getResultFn(tsId))
                outputTomos.append(tomo)
                outputTomos.update(tomo)
                outputTomos.write()
                self._store(outputTomos)

    def closeOutputSetStep(self):
        super()._closeOutputSet()
        if self.failedApixTsIds:
            failedApixTsIdList = String(str(self.failedApixTsIds))
            self._store(failedApixTsIdList)
        if self.failedDimsTsIds:
            failedDimsTsIdList = String(str(self.failedDimsTsIds))
            self._store(failedDimsTsIdList)

    # --------------------------- UTILS functions -----------------------------
    def _getInTomoSet(self, returnPointer: bool = False) -> Union[SetOfTomograms, Pointer]:
        inTomoSetPointer = getattr(self, IN_TOMO_SET)
        return inTomoSetPointer if returnPointer else inTomoSetPointer.get()

    def _getInTomoMasks(self, returnPointer: bool = False) -> Union[SetOfTomoMasks, Pointer]:
        inTomoMasksPointer = getattr(self, IN_MASK_SET)
        return inTomoMasksPointer if returnPointer else inTomoMasksPointer.get()

    def _getResultFn(self, tsId: str):
        return self._getExtraPath(f'{tsId}.mrc')

    def _getSmoothedMaskFn(self, tsId: str) -> str:
        return self._getTmpPath(f'{tsId}_smooth_mask.mrc')

    def _getOutTomos(self) -> SetOfTomograms:
        outputName = self._possibleOutputs.maskedTomograms.name
        outTomograms = getattr(self, outputName, None)
        if outTomograms:
            outTomograms.enableAppend()
        else:
            inSetPointer = self._getInTomoSet(returnPointer=True)
            outTomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            outTomograms.copyInfo(inSetPointer.get())
            outTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputName: outTomograms})
            self._defineSourceRelation(inSetPointer, outTomograms)
            self._defineSourceRelation(self._getInTomoMasks(returnPointer=True), outTomograms)
        return outTomograms

    # --------------------------- INFO functions ------------------------------
    def _summary(self) -> list:
        msgList = []
        nonMatchingTsIdsMsg = self.nonMatchingTsIdsMsg.get()
        if nonMatchingTsIdsMsg:
            msgList.append(f'*{nonMatchingTsIdsMsg}*')
        failedTsIdApixList = getattr(self, 'failedApixTsIdList', None)
        if failedTsIdApixList:
            msgList.append(f'*WARNING*: Failed tsIds because of different pixel size '
                           f'in the tomo mask and the tomogram: {failedTsIdApixList}')
        failedTsIdDims = getattr(self, 'failedDimsTsIdList', None)
        if failedTsIdDims:
            msgList.append(f'*WARNING*: Failed tsIds because of different dimensions '
                           f'in the tomo mask and the tomogram: {failedTsIdDims}')
        return msgList

    def _validate(self) -> list:
        errorList = []
        inTomos = self._getInTomoSet()
        inMasks = self._getInTomoMasks()
        # Check the sampling rate
        inTomosSRate = inTomos.getSamplingRate()
        inMasksSRate = inMasks.getSamplingRate()
        if abs(inTomosSRate - inMasksSRate) > self._sRateTol:
            errorList.append(f'The sampling rate of the given tomograms and mask/s are different within tolerance: '
                             f'abs({inTomosSRate:.3f} - {inMasksSRate:.3f} > {self._sRateTol:3f}')
        return errorList
