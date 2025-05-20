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
from SCons.Subst import escape_list
from scipy.ndimage import gaussian_filter, binary_dilation, distance_transform_edt

from pwem.emlib.image.image_readers import MRCImageReader, ImageReadersRegistry, ImageStack
from pwem.objects import VolumeMask
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
from pyworkflow.protocol import PointerParam, BooleanParam, IntParam, STEPS_PARALLEL, GE
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfTomograms, SetOfTomoMasks, TomoMask

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
        form.addParam('sigmaGaussian', IntParam,
                      default=3,
                      label='Standard deviation for gaussian kernel (smoothing)',
                      validators=[GE(0)],
                      help='A gaussian filter can be applied by setting thi parameter with a value greater than 0. '
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
            smoothId = self._insertFunctionStep(self.smoothMaskStep, tsId,
                                                prerequisites=[],
                                                needsGPU=False)
            aMId = self._insertFunctionStep(self.applyMaskStep, tsId,
                                            prerequisites=smoothId,
                                            needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=aMId,
                                              needsGPU=False)
            closeStepDeps.append(cOutId)
        self._insertFunctionStep(self._closeOutputSet,
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

    def smoothMaskStep(self, tsId: str):
        if self.doSmooth:
            logger.info(cyanStr(f'tsId = {tsId}: smoothing the mask...'))
            mask = self._getCurrentMask(tsId)
            maskFileName = mask.getFileName()
            # Read and smooth
            with mrcfile.mmap(maskFileName, mode='r', permissive=True) as mrc:
                doInvert = True
                data = 1 - mrc.data if doInvert else mrc.data
                data = np.array(data, dtype=float)
                smoothData = gaussian_filter(data, sigma=3)  # TODO: CHANGE THIS
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
        mask = self._getCurrentMask(tsId)
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
                # if self.doSmooth:
                #     remove(self._getSmoothedMaskFn(tsId))

    def createOutputStep(self, tsId: str):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _getInTomoSet(self, returnPointer: bool = False) -> Union[SetOfTomograms, Pointer]:
        inTomoSetPointer = getattr(self, IN_TOMO_SET)
        return inTomoSetPointer if returnPointer else inTomoSetPointer.get()

    def _getInTomoMasks(self) -> Union[SetOfTomoMasks, VolumeMask]:
        return getattr(self, IN_MASK_SET).get()

    def _getResultFn(self, tsId: str):
        return self._getExtraPath(f'{tsId}.mrc')

    def _getSmoothedMaskFn(self, tsId: str) -> str:
        return self._getTmpPath(f'{tsId}_smooth_mask.mrc')

    def _getCurrentMask(self, tsId: str) -> Union[TomoMask, VolumeMask]:
        if self.tomoMaskDict:
            return self.tomoMaskDict[tsId]
        else:
            return self.mask

    def _getSigmaForGaussianSmooth(self) -> float:
        return self.smoothMaskNPix.get() / 3

    def _smoothMask(self,
                    maskData: ImageStack,
                    outFileName: str,
                    samplingRate: float) -> None:
        """It smooths a mask applying a gaussian filtering and save it as specified in the input argument
         named outFileName."""
        # sigma = self._getSigmaForGaussianSmooth()
        # smoothedImgList = [gaussian_filter(zSlice, sigma=sigma) for zSlice in maskData]
        nPix = 5#self.smoothMaskNPix.get()
        smoothedImgList = [self._dilateWithDecay(zSlice, nPix) for zSlice in maskData]
        smoothedStack = ImageStack(smoothedImgList)
        MRCImageReader.write(smoothedStack, outFileName, samplingRate=samplingRate)

    @staticmethod
    def _dilateWithDecay(zSliceBinary: np.array, nPix: int) -> np.array:
        # dilatedZSlice = binary_dilation(zSliceBinary, iterations=nPix)
        euclideanDistanceArray = distance_transform_edt(1 - zSliceBinary)
        expDecayArray = np.exp(-euclideanDistanceArray / nPix)
        return expDecayArray #* dilatedZSlice

    # --------------------------- INFO functions ------------------------------
    def _summary(self) -> list:
        msgList = []
        nonMatchingTsIdsMsg = self.nonMatchingTsIdsMsg.get()
        if nonMatchingTsIdsMsg:
            msgList.append(f'*{nonMatchingTsIdsMsg}*')
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
