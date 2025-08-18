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
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTs, ProtImportTomograms, ProtTsFromTomos
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto, TS_03, TS_54
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestTsFromTomograms(TestBaseCentralizedLayer):
    importedTs = None
    importedTomos = None
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    expectedTsSetSize = 2
    testAcqObjDict, expectedDimensionsDict, anglesCountDict = DataSetRe4STATuto.genTestTsDicts((TS_03, TS_54))

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.importedTomos = cls._runImportTomograms()

    @classmethod
    def _runImportTs(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords='output',
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=DataSetRe4STATuto.tomosPattern.value,
                                          exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    def testTsFromTomos(self):
        """Five tilt series are imported, with tsIds = [TS_01, TS_03, TS_43, TS_45, TS_54] and two
        tomograms are imported, with tsIds = [TS_03, TS_54]. The resulting test of tilt-series should
        be composed of [TS_03, TS_54]."""
        print(magentaStr("\n==> Getting the tilt-series from the tomograms...:"))
        protTsFromTomos = self.newProtocol(ProtTsFromTomos,
                                           inTomoSet=self.importedTomos,
                                           inTsSet=self.importedTs)
        self.launchProtocol(protTsFromTomos)
        ouTsSet = getattr(protTsFromTomos, protTsFromTomos._possibleOutputs.tiltSeries.name, None)
        # Check the tilt-series
        self.checkTiltSeries(ouTsSet,
                             expectedSetSize=self.expectedTsSetSize,
                             expectedSRate=self.unbinnedSRate,
                             imported=True,
                             expectedDimensions=self.expectedDimensionsDict,
                             testAcqObj=self.testAcqObjDict,
                             anglesCount=self.anglesCountDict,
                             isHeterogeneousSet=True,
                             presentTsIds=[TS_03, TS_54])

