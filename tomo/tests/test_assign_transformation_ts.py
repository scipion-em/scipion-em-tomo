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
from typing import Optional, Type

from pwem import ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr, weakImport
from tomo.objects import SetOfTiltSeries
from tomo.protocols import ProtImportTs
from tomo.protocols.protocol_assignTransformationTS import ProtAssignTransformationMatrixTiltSeries
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto, TS_03, TS_54
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

with weakImport("imod"):
    from imod.constants import OUTPUT_TILTSERIES_NAME
    from imod.protocols import ProtImodImportTransformationMatrix, ProtImodTsNormalization


class TestAssignTransformationTS(TestBaseCentralizedLayer):
    ds = None
    binFactor = 2
    exclusionWords = None
    tsSetNoAlignment = None
    tsSetWithAlignment = None
    tsSetBinned = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls._runPrevProtocols()

    @classmethod
    def _runPrevProtocols(cls):
        print(cyanStr('\n--------------------------------- RUNNING PREVIOUS PROTOCOLS '
                      '---------------------------------'))
        cls.exclusionWords = DataSetRe4STATuto.exclusionWordsTs03ts54.value
        # Import TS_03 + TS_54 once (no alignment) — reused across all tests
        cls.tsSetNoAlignment = cls._runImportTs(exclusionWords=cls.exclusionWords)
        # Derive a set with alignment (transform matrices imported via IMOD)
        cls.tsSetWithAlignment = cls._runImportTrMatrix(cls.tsSetNoAlignment)
        # Derive a binned set via IMOD preprocessing (bin 2)
        cls.tsSetBinned = cls._runBinTs(cls.tsSetNoAlignment, binning=cls.binFactor)
        print(cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED '
                      '---------------------------------'))

    @classmethod
    def _runImportTs(cls, exclusionWords: Optional[str] = None):
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
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)
        cls.launchProtocol(protTsImport)
        tsImported = getattr(protTsImport, protTsImport.OUTPUT_NAME, None)
        return tsImported

    @classmethod
    def _runImportTrMatrix(cls, inTsSet: SetOfTiltSeries) -> Optional[SetOfTiltSeries]:
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=DataSetRe4STATuto.transformPattern.value,
                                             inputSetOfTiltSeries=inTsSet)
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runBinTs(cls, inTsSet: SetOfTiltSeries, binning: int = 2) -> Optional[SetOfTiltSeries]:
        print(magentaStr(f"\n==> Binning the tilt-series with IMOD (bin={binning}):"))
        protTSNormalization = cls.newProtocol(ProtImodTsNormalization,
                                              inputSetOfTiltSeries=inTsSet,
                                              binning=binning)
        cls.launchProtocol(protTSNormalization)
        outTsSet = getattr(protTSNormalization, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runAssignTransformation(cls,
                                 fromTsSet: SetOfTiltSeries,
                                 toTsSet: SetOfTiltSeries,
                                 objLabel: str = 'assign transformation') \
            -> Type[ProtAssignTransformationMatrixTiltSeries]:

        print(magentaStr("\n==> Assigning transformation matrices:"))
        protAssign = cls.newProtocol(ProtAssignTransformationMatrixTiltSeries,
                                     getTMSetOfTiltSeries=fromTsSet,
                                     setTMSetOfTiltSeries=toTsSet)
        protAssign.setObjLabel(objLabel)
        cls.launchProtocol(protAssign)
        return protAssign

    def test_assignTransformation_sameSRate(self):
        """Test assigning transformation matrices between two sets of tilt series
        with the same sampling rate. The 'from' set has alignment (imported
        transformation matrices) and the 'to' set is a plain import without
        alignment. The output should have alignment assigned."""
        print(magentaStr("\n==> TEST: Assign transformation with the same sampling rate"))
        protAssign = self._runAssignTransformation(
            self.tsSetWithAlignment, self.tsSetNoAlignment,
            objLabel='assign TM same sRate')

        outputTsSet = getattr(protAssign, protAssign._possibleOutputs.tiltSeries.name, None)
        self.assertIsNotNone(outputTsSet, "No output tilt series set was generated.")

        tsIdList = (TS_03, TS_54)
        testAcqObjDict, expectedDimensionsDict, anglesCountDict = \
            DataSetRe4STATuto.genTestTsDicts(tsIdList=tsIdList)

        self.checkTiltSeries(outputTsSet,
                             expectedSetSize=2,
                             expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value,
                             expectedDimensions=expectedDimensionsDict,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             presentTsIds=[TS_03, TS_54])

    def test_assignTransformation_differentSRate(self):
        """Test assigning transformation matrices between two sets of tilt series
        with different sampling rates. The 'from' set has alignment at the unbinned
        sampling rate (1.35 A/pix), and the 'to' set is preprocessed with IMOD
        (bin 2 -> 2.7 A/pix), producing genuinely binned files with correct
        sampling rate and dimensions. Verifies that the protocol correctly scales
        the transformation matrix shifts by the sampling rate ratio
        (toSRate / fromSRate = 2.0)."""
        print(magentaStr("\n==> TEST: Assign transformation with different sampling rates"))
        unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value  # 1.35
        binnedSRate = unbinnedSRate * self.binFactor  # 2.7

        # Verify sampling rates are genuinely different before running the protocol
        self.assertAlmostEqual(self.tsSetWithAlignment.getSamplingRate(), unbinnedSRate,
                               delta=0.01,
                               msg="'from' set sampling rate is not the expected unbinned value.")
        self.assertAlmostEqual(self.tsSetBinned.getSamplingRate(), binnedSRate,
                               delta=0.01,
                               msg="'to' set sampling rate is not the expected binned value.")
        self.assertNotAlmostEqual(self.tsSetWithAlignment.getSamplingRate(),
                                  self.tsSetBinned.getSamplingRate(),
                                  delta=0.01,
                                  msg="Sampling rates should differ between 'from' and 'to' sets.")

        protAssign = self._runAssignTransformation(
            self.tsSetWithAlignment, self.tsSetBinned,
            objLabel='assign TM diff sRate')

        outputTsSet = getattr(protAssign, protAssign._possibleOutputs.tiltSeries.name, None)
        self.assertIsNotNone(outputTsSet, "No output tilt series set was generated.")

        tsIdList = (TS_03, TS_54)
        testAcqObjDict, expectedDimensionsDict, anglesCountDict = \
            DataSetRe4STATuto.genTestTsDicts(tsIdList=tsIdList,
                                             binFactor=self.binFactor,
                                             isImod=True)

        self.checkTiltSeries(outputTsSet,
                             expectedSetSize=2,
                             expectedSRate=binnedSRate,
                             expectedDimensions=expectedDimensionsDict,
                             testAcqObj=testAcqObjDict,
                             anglesCount=anglesCountDict,
                             hasAlignment=True,
                             alignment=ALIGN_2D,
                             presentTsIds=[TS_03, TS_54])

        # Check that the matrix shifts scale were properly scaled
        self.checkTrMatrixShiftsScale(self.tsSetWithAlignment, outputTsSet)

    def test_assignTransformation_validate_noAlignment(self):
        """Test that the protocol's validation correctly reports errors when
        the 'from' set has no alignment (no transformation matrices).
        The protocol requires the 'from' set to have alignment, so validate()
        must return meaningful error messages without launching the protocol."""
        print(magentaStr("\n==> TEST: Validation - 'from' TS without alignment"))
        protAssign = self.newProtocol(ProtAssignTransformationMatrixTiltSeries,
                                      getTMSetOfTiltSeries=self.tsSetNoAlignment,
                                      setTMSetOfTiltSeries=self.tsSetNoAlignment)
        protAssign.setObjLabel('assign TM no alignment (validate)')

        # Validate directly — no protocol launch, no exception
        errors = protAssign.validate()
        self.assertGreater(len(errors), 0,
                           "Protocol should report validation errors when 'from' set "
                           "has no alignment.")
        errorMsg = '\n'.join(errors)
        self.assertIn("transformation matrix", errorMsg,
                      f"Validation error should mention missing transformation "
                      f"matrix. Got: {errorMsg}")