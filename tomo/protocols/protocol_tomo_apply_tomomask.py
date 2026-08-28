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


class ApplyTomoMaskOutputs(Enum):
    maskedTomograms = SetOfTomograms


class ApplyTomoMaskFormParams(Enum):
    IN_TOMO_SET = 'inTomoSet'
    IN_MASK_SET = 'inMaskSet'
    INVERT_MASK = 'invertMask'
    DILATION_PX = 'dilationPixels'
    SIGMA_GAUSSIAN = 'sigmaGaussian'


class ProtTomoApplyTomoMask(EMProtocol):
    """
    Applies a set of tomographic masks to a set of tomograms by matching both inputs through
    their tilt-series identifiers (tsId). Each tomogram is multiplied voxel-wise by its
    corresponding mask, producing a new masked tomogram while preserving the original
    tomographic metadata.

    AI Generated:

    Apply Tomo Masks to Tomograms (ProtTomoApplyTomoMask) — User Manual
        Overview

        The Apply Tomo Masks to Tomograms protocol applies one mask to each tomogram
        in an input set. Matching between tomograms and masks is performed using the
        tilt-series identifier (tsId), ensuring that every tomogram receives only the
        mask that belongs to the same acquisition or reconstruction series.

        In practical cryo-electron tomography workflows, masking is commonly used to
        isolate biologically relevant regions, suppress surrounding solvent, remove
        noisy peripheral densities, or focus downstream analyses such as segmentation,
        subvolume extraction, or visualization.

        For biological users, this protocol is especially useful when the tomograms
        already contain reconstructed cellular or macromolecular information, but only
        specific spatial regions should be preserved for further interpretation.

        Inputs and Matching Strategy

        The protocol requires two inputs:

        - A set of tomograms.
        - A single tomographic mask or a set of tomographic masks.

        Matching is entirely based on tsId values. Only tomograms and masks sharing
        the same tsId will be processed together.

        If some tomograms or masks do not have a corresponding partner, they are not
        processed. The protocol reports these non-matching identifiers so the user can
        verify dataset consistency.

        If no common tsId exists between both input sets, execution stops with an error.

        This design is especially important in batch tomography workflows where many
        tomograms are processed in parallel and preserving one-to-one correspondence
        between reconstructions and masks is essential.

        Mask Preprocessing

        Before applying the mask, the protocol optionally allows several operations
        that modify the mask itself.

        Inversion

        The mask can be inverted before application.

        This is useful when the provided mask represents regions that should be
        suppressed instead of preserved. After inversion, masked and unmasked
        regions are exchanged.

        Dilation

        The mask can be expanded by a user-defined number of pixels in all directions.

        From a biological perspective, dilation is useful when the original mask is
        too conservative and may cut away peripheral structural information.
        For example, membrane-associated densities, flexible domains, or low-contrast
        boundaries often benefit from a small dilation.

        Excessive dilation should be avoided because it may reintroduce noise or
        irrelevant surrounding density.

        Gaussian Smoothing

        If the Gaussian sigma parameter is larger than zero, the mask is smoothed
        before application.

        Rather than producing a sharp binary edge, smoothing generates a gradual
        transition between masked and unmasked regions. This is often biologically
        advantageous because abrupt boundaries can create artificial discontinuities
        that later affect filtering, segmentation, visualization, or mathematical
        operations.

        Small sigma values preserve sharper boundaries, while larger sigma values
        create more diffuse transitions.

        In practice, smoothing is often one of the most useful options when the goal
        is to keep biologically meaningful densities while avoiding masking artifacts.

        Validation Before Mask Application

        Before multiplying a tomogram by its mask, the protocol performs per-tomogram
        validation checks.

        Sampling Rate Consistency

        The voxel size of the tomogram and the mask must agree within a small tolerance.

        This validation is essential because voxel-wise multiplication only makes
        biological sense when both volumes represent the same physical sampling.
        A mismatch in pixel size means that densities would not correspond spatially.

        If the sampling rates differ beyond tolerance, the tomogram is skipped and
        its tsId is stored as failed.

        Dimension Consistency

        The dimensions of the tomogram and mask must also be identical.

        If dimensions differ, the protocol cannot apply voxel-wise multiplication.
        Such tomograms are skipped and reported separately.

        These validations are especially important in heterogeneous tomography projects
        where masks may originate from different processing branches or intermediate
        reconstruction stages.

        Mask Application

        Once validation succeeds, the protocol applies the mask by multiplying each
        tomographic slice with the corresponding mask slice.

        This operation is performed slice by slice through the entire 3D volume.

        Biologically, the resulting tomogram preserves signal only in regions weighted
        by the mask. Binary masks preserve selected voxels and suppress the rest,
        while smoothed masks produce gradual attenuation.

        The resulting masked tomogram is written as a new MRC volume while preserving
        the sampling rate of the original tomogram.

        Temporary Files and Memory Handling

        When smoothing is enabled, the protocol creates temporary smoothed masks
        inside the temporary working directory.

        After the masked tomogram is generated, these temporary files are immediately
        removed.

        This behavior is particularly important in tomography workflows because
        tomograms and masks are often very large volumes. Removing temporary files
        avoids unnecessary disk consumption during parallel execution.

        Parallel Execution Strategy

        The protocol executes independently for each matched tsId.

        For every tomogram-mask pair, the workflow follows three sequential stages:

        1. Process the mask (optional dilation, smoothing, inversion).
        2. Apply the mask to the tomogram.
        3. Register the resulting masked tomogram in the output set.

        Different tomograms can be processed independently, which makes the protocol
        naturally suitable for parallel execution when large tomography datasets are used.

        Outputs

        The protocol produces a new set of tomograms called maskedTomograms.

        For every successfully processed tomogram:

        - The original tomographic metadata is preserved.
        - The file location is updated to point to the newly generated masked volume.

        Tomograms that fail validation are not included in the output set.

        The protocol also stores warning information about:

        - Non-matching tsIds between input tomograms and masks.
        - Failed tsIds caused by sampling rate mismatch.
        - Failed tsIds caused by dimension mismatch.

        Practical Recommendations

        In biological practice, this protocol is most useful when the user wants to
        focus analysis on a known region of interest.

        Typical examples include:

        - Isolating a cellular compartment from a crowded tomogram.
        - Preserving only membrane-proximal densities.
        - Removing solvent or empty reconstruction regions.
        - Preparing tomograms for segmentation, particle picking, or subvolume analysis.

        When the mask is already biologically well defined, using no dilation and
        little or no smoothing often provides the most faithful result.

        When mask boundaries are uncertain or too sharp, a small dilation together
        with moderate Gaussian smoothing often improves the continuity of biologically
        meaningful densities.

        Final Perspective

        Applying tomographic masks is not merely a computational filtering step.
        It is often a biologically meaningful operation that determines which parts
        of the reconstructed volume remain visible and interpretable.

        Careful control of mask matching, sampling consistency, and boundary smoothing
        is essential to avoid introducing artifacts or losing relevant structural
        information.

        In most tomography workflows, this protocol serves as a reliable preparation
        step before interpretation, segmentation, visualization, or downstream
        quantitative analysis.
    """
    _label = 'apply tomomasks to tomograms'
    _devStatus = BETA
    _possibleOutputs = ApplyTomoMaskOutputs
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
        form.addParam(ApplyTomoMaskFormParams.IN_TOMO_SET.value,
                      PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')
        form.addParam(ApplyTomoMaskFormParams.IN_MASK_SET.value,
                      PointerParam,
                      pointerClass="VolumeMask, SetOfTomoMasks",
                      important=True,
                      label='Masks',
                      help='The protocol will try to match the tomograms and the masks by tsId.')
        group = form.addGroup('Mask operations')
        group.addParam(ApplyTomoMaskFormParams.INVERT_MASK.value,
                       BooleanParam,
                       default=False,
                       label='Invert mask?')
        group.addParam(ApplyTomoMaskFormParams.DILATION_PX.value,
                       IntParam,
                       default=0,
                       label='Number of pixels for dilation',
                       validators=[GE(0)],
                       haelp='The dilation will expands the shape by the given number of pixels '
                             'in all directions.')
        group.addParam(ApplyTomoMaskFormParams.SIGMA_GAUSSIAN.value,
                       IntParam,
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
        if self._getFormAttrib(ApplyTomoMaskFormParams.SIGMA_GAUSSIAN.value) > 0:
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
                data = mrc.data
            # Dilate the mask
            dilationPixels = self._getFormAttrib(ApplyTomoMaskFormParams.DILATION_PX.value)
            if dilationPixels > 0:
                data = np.array(data, dtype=bool)  # Required for the binary dilation
                data = binary_dilation(data, iterations=dilationPixels)
            # Smooth the mask
            data = np.array(data, dtype=float)  # Required to be cast from uint8 to float for the gaussian filtering
            sigma = self._getFormAttrib(ApplyTomoMaskFormParams.SIGMA_GAUSSIAN.value)
            smoothData = gaussian_filter(data, sigma=sigma)
            # Invert if required
            invertMask = self._getFormAttrib(ApplyTomoMaskFormParams.INVERT_MASK.value)
            smoothData = 1 - smoothData if invertMask else smoothData
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
        # Validate the sampling rate
        if abs(tomoSRate - maskSRate) > self._sRateTol:
            self.failedApixTsIds.append(tsId)
        else:
            maskFileName = self._getSmoothedMaskFn(tsId) if self.doSmooth else mask.getFileName()
            tomoFileName = tomo.getFileName()
            maskDims = MRCImageReader.getDimensions(maskFileName)
            tomoDims = MRCImageReader.getDimensions(tomoFileName)
            # Validate the dimensions
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
    def _getFormAttrib(self, attribName: str):
        return getattr(self, attribName).get()

    def _getInTomoSet(self, returnPointer: bool = False) -> Union[SetOfTomograms, Pointer]:
        inTomoSetPointer = getattr(self, ApplyTomoMaskFormParams.IN_TOMO_SET.value)
        return inTomoSetPointer if returnPointer else inTomoSetPointer.get()

    def _getInTomoMasks(self, returnPointer: bool = False) -> Union[SetOfTomoMasks, Pointer]:
        inTomoMasksPointer = getattr(self, ApplyTomoMaskFormParams.IN_MASK_SET.value)
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
