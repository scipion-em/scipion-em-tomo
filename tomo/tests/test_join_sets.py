# *
# * Authors:     Scipion Team
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
from pwem.protocols import ProtUnionSet
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTs, ProtImportTomograms
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestJoinTomoSets(TestBaseCentralizedLayer):
    ds = None
    nTiltSeries = 2

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        # 5 TS with no. tilt-images:
        #   - TS_01 = 40
        #   - TS_03 = 40
        #   - TS_43 = 41
        #   - TS_45 = 41
        #   - TS_54 = 41
        #
        # 5 tomograms with a thickness of (px):
        #   - TS_01 = 340
        #   - TS_03 = 280
        #   - TS_43 = 300
        #   - TS_45 = 300
        #   - TS_54 = 280

    @classmethod
    def _runImportTs(cls, exclusionWords=None):
        print(magentaStr("\n==> Importing the tilt series:"))
        protTsImport = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=exclusionWords,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=DataSetRe4STATuto.unbinnedPixSize.value,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protTsImport)
        tsImported = getattr(protTsImport, protTsImport.OUTPUT_NAME, None)
        return tsImported

    @classmethod
    def _runImportTomograms(cls, filesPattern=None):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=filesPattern,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    # TODO: add assertions

    def test_join_ts_heterogeneous_sets(self):
        print(magentaStr("\n==> Join sets of TS with different number of tilt-images:"))
        tsSet41imgs1 = self._runImportTs(exclusionWords='output 01 03 43 45')
        tsSet41imgs2 = self._runImportTs(exclusionWords='output 01 03 54')
        protUnion = self.newProtocol(ProtUnionSet)
        protUnion.inputSets.append(tsSet41imgs1)
        protUnion.inputSets.append(tsSet41imgs2)
        protUnion.setObjLabel('join het tsSets')
        self.launchProtocol(protUnion)

    def test_join_ts_homogeneous_sets(self):
        print(magentaStr("\n==> Join sets of TS with the same number of tilt-images:"))
        tsSet41imgs = self._runImportTs(exclusionWords='output 01 03')
        tsSet40imgs = self._runImportTs(exclusionWords='output 43 45 54')
        protUnion = self.newProtocol(ProtUnionSet)
        protUnion.inputSets.append(tsSet41imgs)
        protUnion.inputSets.append(tsSet40imgs)
        protUnion.setObjLabel('join homogen tsSets')
        self.launchProtocol(protUnion)

    def test_join_tomos_heterogeneous_sets(self):
        print(magentaStr("\n==> Join sets of tomograms with different thickness:"))
        tomoSetThk300 = self._runImportTomograms(filesPattern='TS_4*.mrc')
        tomoSetThk340 = self._runImportTomograms(filesPattern='TS_01.mrc')
        tomoSetThk280_1 = self._runImportTomograms(filesPattern='TS_03.mrc')
        tomoSetThk280_2 = self._runImportTomograms(filesPattern='TS_54.mrc')
        protUnion = self.newProtocol(ProtUnionSet)
        protUnion.inputSets.append(tomoSetThk300)
        protUnion.inputSets.append(tomoSetThk340)
        protUnion.inputSets.append(tomoSetThk280_1)
        protUnion.inputSets.append(tomoSetThk280_2)
        protUnion.setObjLabel('join het tomoSets')
        self.launchProtocol(protUnion)

    # def test_join_tomos_homogeneous_sets(self):
    #     print(magentaStr("\n==> Join sets of tomograms with different thickness:"))
    #     tomoSetThk300 = self._runImportTomograms(filesPattern='TS_4*.mrc')
    #     tomoSetThk340 = self._runImportTomograms(filesPattern='TS_01.mrc')
    #     tomoSetThk280_1 = self._runImportTomograms(filesPattern='TS_03.mrc')
    #     tomoSetThk280_2 = self._runImportTomograms(filesPattern='TS_54.mrc')
    #     protUnion = self.newProtocol(ProtUnionSet)
    #     protUnion.inputSets.append(tomoSetThk300)
    #     protUnion.inputSets.append(tomoSetThk340)
    #     protUnion.inputSets.append(tomoSetThk280_1)
    #     protUnion.inputSets.append(tomoSetThk280_2)
    #     self.launchProtocol(protUnion)
