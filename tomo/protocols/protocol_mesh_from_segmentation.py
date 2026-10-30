# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
# *              Oier Lauzirika Zarrabeita (oierlauzi@bizkaia.eu)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import traceback
from enum import Enum
import math
from os.path import basename

import numpy as np
import skimage.morphology
from pwem.emlib.image.image_readers import ImageReadersRegistry, MRCImageReader
from pwem.protocols import EMProtocol
from pwem.objects import Set, Integer
from pyworkflow import BETA
from pyworkflow.protocol import STEPS_PARALLEL, LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam, GE, LE)
from pyworkflow.utils import yellowStr, cyanStr, redStr, createLink
from tomo.objects import (SetOfMeshes, MeshPoint, SetOfTomoMasks, TomoMask,
                          SetOfTomograms, Tomogram)
import tomo.constants as const

logger = logging.getLogger(__name__)
MRC_EXT = '.mrc'


class OutputMeshesFromSegmentation(Enum):
    meshes = SetOfMeshes


class ProtMeshFromSegmentation(EMProtocol):
    """
    Creates meshes based on segmentations or voxels values (TomoMasks).
    """
    _label = 'meshes from tomo mask'
    _devStatus = BETA
    _possibleOutputs = OutputMeshesFromSegmentation
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.failedTsIds = []

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMasks', PointerParam, pointerClass=SetOfTomoMasks,
                      label='Tomo Masks', important=True,
                      help='Set of tomo mask from which the meshes will be created')
        form.addParam('inputTomograms', PointerParam, pointerClass=SetOfTomograms,
                      label='Tomograms',
                      help='Tomograms to which the meshes will be asotiated')
        form.addParam('boxSize', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=32,
                      label="Box size (px)",
                      help="Required by some visualization tools")

        form.addSection(label='Parameters')
        form.addParam('smoothMask', BooleanParam, label='Smooth mask',
                      important=True, default=False,
                      help='Wether the input is a segmentation or a smooth mask')
        form.addParam('backgroundLabel', IntParam, label='Background label',
                      condition='not smoothMask', default=0,
                      help='Label in the segmentation to be considered as '
                           'background')

        line = form.addLine('Threshold to keep',
                            condition='smoothMask',
                            help="Only the voxels between these two values "
                                 "will be considered to create the meshes.")
        line.addParam('lowLimit', FloatParam, default=0.1, label='Lowest')
        line.addParam('highLimit', FloatParam, default=1, label='Highest')

        group = form.addGroup('Morphological operations',
                              help='Operations are performed in the same order '
                                   'as in the form')
        group.addParam('applyDilation', IntParam, label='Dilation',
                       default=0,
                       help='When set to a positive number, a dilation '
                            'operation is applied with the provided pixel count.')
        group.addParam('applySkeletonization', BooleanParam, label='Skeletonization',
                       default=True,
                       help='When set to yes, a skeletonization operation is applied.')

        form.addParam('density',
                      FloatParam,
                      label='Percentage of density ',
                      default=5.0,
                      validators=[GE(0), LE(100)],
                      help='This parameter goes from 0 - 100 and defines the '
                           'percentage of voxel of the tomoMask that '
                           'will be considered as points of the mesh.')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        closeSetDeps = []
        self._initialize()
        for tsId in self.masks.keys():
            cInPId = self._insertFunctionStep(self.convertInputStep, tsId,
                                              prerequisites=[],
                                              needsGPU=False)
            pId = self._insertFunctionStep(self.processTomogramStep, tsId,
                                           prerequisites=cInPId,
                                           needsGPU=False)
            closeSetDeps.append(pId)
        self._insertFunctionStep(self.closeOutputSet,
                                 prerequisites=[],
                                 needsGPU=False)

    # --------------------------- STEPS functions ------------------------------
    def convertInputStep(self, tsId: str):
        try:
            # Mask
            mask = self._getInputMask(tsId)
            maskFName = mask.getFileName()
            convertedMaskFn = self._getMaskFn(tsId)
            self._convertOrLink(tsId, maskFName, convertedMaskFn)
            # Tomogram
            tomogram = self._getInputTomogram(tsId)
            tomoFName = tomogram.getFileName()
            convertedTomoFn = self._getTomoFn(tsId)
            self._convertOrLink(tsId, tomoFName, convertedTomoFn)
        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed '
                                f'with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def processTomogramStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId} -> generating the mesh...'))
            mask = self._getInputMask(tsId)
            tomogram = self._getInputTomogram(tsId)

            with self._lock:
                outputMeshes = self._getOutputMeshes()
                maskData = self._loadMask(mask)

                if self.smoothMask:
                    self._processSmoothMask(
                        mesh=outputMeshes,
                        tomogram=tomogram,
                        mask=maskData
                    )
                else:
                    self._processSegmentation(
                        mesh=outputMeshes,
                        tomogram=tomogram,
                        segmentation=maskData
                    )

                outputMeshes.write()
                self._store()
        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> process tomogram step failed '
                                f'with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def closeOutputSet(self):
        self._closeOutputSet()
        outputAttribName = self._possibleOutputs.meshes.name
        output = getattr(self, outputAttribName, None)
        if not output or (output and len(output) == 0):
            raise Exception(f'No output/s {outputAttribName} were generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')

    # --------------------------- UTILS functions ------------------------------
    def _initialize(self):
        inTomoMasks = self.inputMasks.get()
        inTomograms = self.inputTomograms.get()
        presentMaskTsIds = set(inTomoMasks.getTSIds())
        presentTomoTsIds = set(inTomoMasks.getTSIds())
        commonTsIds = presentMaskTsIds & presentTomoTsIds
        nonPresenTsIds = presentMaskTsIds ^ presentTomoTsIds
        if nonPresenTsIds:
            logger.info(yellowStr(f'Some tsIds are not common to both sets of '
                                  f'tomograms and tomomasks: {nonPresenTsIds}'))
        self.masks = {tsId: mask.clone() for mask in inTomoMasks if
                      (tsId := mask.getTsId()) in commonTsIds}
        self.tomos = {tsId: tomo.clone() for tomo in inTomograms if
                      (tsId := tomo.getTsId()) in commonTsIds}
        self.baseGroupId = Integer(1)

    def _getOutputMeshes(self) -> SetOfMeshes:
        output: SetOfMeshes = getattr(self, self._possibleOutputs.meshes.name, None)
        if output is not None:
            output.enableAppend()
        else:
            output = SetOfMeshes.create(self._getPath(), template='meshes%s.sqlite')
            output.setPrecedents(self.inputTomograms)
            output.setSamplingRate(self.inputMasks.get().getSamplingRate())
            output.setBoxSize(self.boxSize.get())
            output.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{self._possibleOutputs.meshes.name: output})
            self._defineSourceRelation(self.inputMasks, output)

        return output

    def _getInputMask(self, tsId: str) -> TomoMask:
        return self.masks[tsId]

    def _getInputTomogram(self, tsId: str) -> Tomogram:
        return self.tomos[tsId]

    @staticmethod
    def _loadMask(mask: TomoMask) -> np.ndarray:
        maskFName = mask.getFileName()
        return MRCImageReader.open(maskFName)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        sRateTol = 0.01
        masksSRate = self.inputMasks.get().getSamplingRate()
        tomosSRate = self.inputTomograms.get().getSamplingRate()
        if abs(masksSRate - tomosSRate) > sRateTol:
            errors.append(f'The sampling rate of the introduced sets are not equal within '
                          f'the tolerance {sRateTol}.\n'
                          f'MasksSRate -> {masksSRate} =! TomosSRate -> {tomosSRate}.')
        return errors

    def _processSegmentation(self,
                             mesh: SetOfMeshes,
                             tomogram: Tomogram,
                             segmentation: np.ndarray) -> None:
        labels = np.unique(segmentation)
        for label in labels:
            if label != self.backgroundLabel.get():
                mask = (segmentation == label)
                self._processBinaryMask(
                    mesh=mesh,
                    tomogram=tomogram,
                    mask=mask
                )

    def _processSmoothMask(self,
                           mesh: SetOfMeshes,
                           tomogram: Tomogram,
                           mask: np.ndarray) -> None:
        mask = (self.lowLimit.get() <= mask) & (mask <= self.highLimit.get())
        self._processBinaryMask(
            mesh=mesh,
            tomogram=tomogram,
            mask=mask
        )

    def _processBinaryMask(self,
                           mesh: SetOfMeshes,
                           tomogram: Tomogram,
                           mask: np.ndarray):
        if self.applyDilation.get() > 0:
            footprint = skimage.morphology.ball(self.applyDilation.get())
            mask = skimage.morphology.binary_dilation(mask, footprint)

        if self.applySkeletonization.get():
            mask = skimage.morphology.skeletonize_3d(mask)

        probability = self.density.get() / 100.0
        coordinates = np.argwhere(mask)
        indices = np.arange(0, len(coordinates))
        nPoints = math.floor(probability * len(indices))
        selection = np.random.choice(indices, size=nPoints, replace=False)

        for z, y, x in coordinates[selection, :]:
            point = MeshPoint()
            point.setVolume(tomogram)
            point.setGroupId(self.baseGroupId)
            point.setPosition(x, y, z, const.BOTTOM_LEFT_CORNER)
            mesh.append(point)

        self.baseGroupId.increment()

    def _getMaskFn(self, tsId: str) -> str:
        return self._getExtraPath(f'{tsId}_mask{MRC_EXT}')

    def _getTomoFn(self, tsId: str) -> str:
        return self._getExtraPath(f'{tsId}{MRC_EXT}')

    @staticmethod
    def _isMrc(filename: str) -> bool:
        return MRC_EXT in filename

    def _convertOrLink(self, tsId: str, inFileName: str, outFileName: str) -> None:
        inBaseName = basename(inFileName)
        if self._isMrc(inFileName):
            logger.info(cyanStr(f'tsId = {tsId} {inBaseName} -> Creating a link...'))
            createLink(inFileName, outFileName)
        else:
            logger.info(cyanStr(f'tsId = {tsId} {inBaseName} -> Converting into {MRC_EXT}...'))
            imgStack = ImageReadersRegistry.open(inFileName)
            ImageReadersRegistry.write(imgStack, outFileName)
