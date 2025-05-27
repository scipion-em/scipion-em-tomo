# ***************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
from typing import Tuple, Union
from imod.constants import OUTPUT_TOMOGRAMS_NAME
from imod.protocols import ProtImodTomoNormalization
from imod.protocols.protocol_base import IN_TOMO_SET, BINNING_FACTOR
from imod.protocols.protocol_base_preprocess import NO_ADJUST
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr
from tomo.objects import SetOfTomograms
from tomo.protocols import ProtImportTomograms, ProtImportTomomasks
from tomo.protocols.protocol_tomo_apply_tomomask import ApplyTomoMaskFormParams, ProtTomoApplyTomoMask
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from tomo.tests import EMD_10439, DataSetEmd10439


class TestApplyTomoMask(TestBaseCentralizedLayer):
    ds = None
    importedTomos = None
    binnedTomos = None
    importedTomoMasksBin2 = None
    binFactor = 2
    unbinnedSRate = DataSetEmd10439.unbinnedSRate.value

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(EMD_10439)
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        try:
            print(cyanStr(
                '--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
            cls._runPreviousProtocols()
            print(
                cyanStr(
                    '\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))
        except Exception as e:
            raise Exception(f'Something Failed when executing the previous protocols -> {e}')

    @classmethod
    def _runPreviousProtocols(cls):
        cls._importTomograms()
        cls._runBinTomograms()
        cls._importTomoMasks()

    @classmethod
    def _importTomograms(cls) -> None:
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.value),
                                             samplingRate=cls.unbinnedSRate)
        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, 'Tomograms', None)
        cls.importedTomos = outputTomos

    @classmethod
    def _runBinTomograms(cls) -> None:
        print(magentaStr('\n==> Binning the tomograms with IMOD:'))
        protArgsDict = {
            IN_TOMO_SET: cls.importedTomos,
            BINNING_FACTOR: cls.binFactor,
            'floatDensities': NO_ADJUST
        }
        protBinTomos = cls.newProtocol(ProtImodTomoNormalization, **protArgsDict)
        cls.launchProtocol(protBinTomos)
        binnedTomos = getattr(protBinTomos, OUTPUT_TOMOGRAMS_NAME, None)
        cls.binnedTomos = binnedTomos

    @classmethod
    def _importTomoMasks(cls) -> None:
        print(magentaStr("\n==> Importing data - tomoMasks:"))
        protImportTomomasks = cls.newProtocol(ProtImportTomomasks,
                                              filesPath=cls.ds.getFile(DataSetEmd10439.tomoMaskByTardisBin2.value),
                                              inputTomos=cls.binnedTomos)
        cls.launchProtocol(protImportTomomasks)
        cls.importedTomoMasksBin2 = getattr(protImportTomomasks, protImportTomomasks._possibleOutputs.tomomasks.name, None)

    @classmethod
    def _runApplyTomoMask(cls,
                          invertMask: bool = False,
                          dilationPixels: int = 0,
                          sigmaGaussian: float = 3.0) -> Union[SetOfTomograms, None]:

        infoStr, objLabel = cls._getInfoStrs(invertMask, dilationPixels, sigmaGaussian)
        print(magentaStr(infoStr))
        protocolInputDict = {
            ApplyTomoMaskFormParams.IN_TOMO_SET.value: cls.binnedTomos,
            ApplyTomoMaskFormParams.IN_MASK_SET.value: cls.importedTomoMasksBin2,
            ApplyTomoMaskFormParams.INVERT_MASK.value: invertMask,
            ApplyTomoMaskFormParams.DILATION_PX.value: dilationPixels,
            ApplyTomoMaskFormParams.SIGMA_GAUSSIAN.value: sigmaGaussian
        }
        protApplyTomoMask = cls.newProtocol(ProtTomoApplyTomoMask, **protocolInputDict)
        protApplyTomoMask.setObjLabel(objLabel)
        cls.launchProtocol(protApplyTomoMask)
        maskedTomos = getattr(protApplyTomoMask, protApplyTomoMask._possibleOutputs.maskedTomograms.name, None)
        return maskedTomos

    @staticmethod
    def _getInfoStrs(invertMask: bool,
                     dilationPixels: int,
                     sigmaGaussian: float) -> Tuple[str, str]:
        infoStr = (f'\n==> Applying the tomoMasks:'
                   f'\n\t- Invert Masks = {invertMask}'
                   f'\n\t- Dilation pixels = {dilationPixels}'
                   f'\n\t- Gaussian sigma = {sigmaGaussian}')
        objLabel = f'Apply TomoMask inv={invertMask}, dpx={dilationPixels}, sgm={sigmaGaussian}'
        return infoStr, objLabel

    def checkTomos(self, tomoSet: SetOfTomograms):
        self.checkTomograms(tomoSet,
                            expectedSetSize=len(self.importedTomos),
                            expectedSRate=DataSetEmd10439.unbinnedSRate.value * self.binFactor,
                            expectedDimensions=DataSetEmd10439.getBinnedDims(self.binFactor))

    def testApplyTomoMasks_01(self):
        maskedTomos = self._runApplyTomoMask(invertMask=False,
                                             dilationPixels=0,
                                             sigmaGaussian=3)
        self.checkTomos(maskedTomos)

    def testApplyTomoMasks_02(self):
        maskedTomos = self._runApplyTomoMask(invertMask=True,
                                             dilationPixels=0,
                                             sigmaGaussian=1)
        self.checkTomos(maskedTomos)

    def testApplyTomoMasks_03(self):
        maskedTomos = self._runApplyTomoMask(invertMask=False,
                                             dilationPixels=12,
                                             sigmaGaussian=4)
        self.checkTomos(maskedTomos)


