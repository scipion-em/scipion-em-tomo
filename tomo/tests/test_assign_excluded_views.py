# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
"""
Tests for ProtAssignExcludedViews — transfer excluded-view annotations
(_enabled state) from a source set of tilt-series to a target set.

Dataset: RE4_STA_TUTO (5 tilt-series: TS_01, TS_03, TS_43, TS_45, TS_54)

Three scenarios:
  1) All 5 TS as target, exclusions on TS_03 (4 views) and TS_43 (5 views).
  2) Subset TS_03 + TS_54 as target, exclusions on both (3 and 5 views).
  3) Same as (2), followed by a physical restack using IMOD's excludeviews.
"""
import copy
from typing import Optional, Dict, List

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, weakImport
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from tomo.protocols import ProtImportTs, ProtImportTsBase
from tomo.protocols.protocol_assign_excluded_views import ProtAssignExcludedViews
from tomo.tests import (
    RE4_STA_TUTO, DataSetRe4STATuto,
    TS_01, TS_03, TS_43, TS_45, TS_54
)
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

with weakImport("imod"):
    from imod.protocols import ProtImodExcludeViews
    from imod.constants import OUTPUT_TILTSERIES_NAME


class TestAssignExcludedViews(TestBaseCentralizedLayer):
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value

    # Test 1: all 5 TS as target; exclude on TS_03 (4 views) and TS_43 (5 views)
    excludedViewsTs03Ts43 = {
        TS_01: [],
        TS_03: [0, 1, 38, 39],
        TS_43: [0, 1, 2, 39, 40],
        TS_45: [],
        TS_54: [],
    }

    # Tests 2 & 3: subset TS_03 + TS_54; exclude on both (3 and 5 views)
    excludedViewsTs03Ts54 = {
        TS_01: [],
        TS_03: [0, 38, 39],
        TS_43: [],
        TS_45: [],
        TS_54: [0, 1, 38, 39, 40],
    }

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)

        allTsIds = (TS_01, TS_03, TS_43, TS_45, TS_54)
        subsetTsIds = (TS_03, TS_54)

        cls.testAcqObjDict5, cls.expectedDimsDict5, cls.anglesCountDict5 = \
            DataSetRe4STATuto.genTestTsDicts(allTsIds)
        cls.testAcqObjDict2, cls.expectedDimsDict2, cls.anglesCountDict2 = \
            DataSetRe4STATuto.genTestTsDicts(subsetTsIds)

        cls.expectedSetSize5 = len(allTsIds)
        cls.expectedSetSize2 = len(subsetTsIds)

        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        print(magentaStr('\n--- Importing tilt-series for assign-excluded-views tests ---'))

        # Target for test 1: all 5 TS (clean)
        cls.tsSetAll = cls._runImportTs(objLabel='target: all 5 TS')

        # Source for test 1: all 5 TS with exclusions on TS_03 and TS_43
        cls.tsSetAllWithEV = cls._runImportTs(objLabel='source: all 5 TS (excl)')
        cls._excludeSetViews(cls.tsSetAllWithEV, cls.excludedViewsTs03Ts43)

        # Target for tests 2 & 3: subset TS_03 + TS_54 (clean)
        cls.tsSubset = cls._runImportTs(
            exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
            objLabel='target: TS_03+TS_54')

        # Source for tests 2 & 3: TS_03 + TS_54 with exclusions
        cls.tsSubsetWithEV = cls._runImportTs(
            exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
            objLabel='source: TS_03+TS_54 (excl)')
        cls._excludeSetViews(cls.tsSubsetWithEV, cls.excludedViewsTs03Ts54)

    @classmethod
    def _runImportTs(cls, exclusionWords: str = 'output', objLabel: str = 'Import TS') \
            -> Optional[SetOfTiltSeries]:
        print(magentaStr(f"\n==> Importing tilt-series: {objLabel}"))
        protImportTs = cls.newProtocol(
            ProtImportTs,
            filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
            filesPattern=DataSetRe4STATuto.tsPattern.value,
            exclusionWords=exclusionWords,
            anglesFrom=2,  # From tlt file
            voltage=DataSetRe4STATuto.voltage.value,
            magnification=DataSetRe4STATuto.magnification.value,
            sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
            amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
            samplingRate=cls.unbinnedSRate,
            doseInitial=DataSetRe4STATuto.initialDose.value,
            dosePerFrame=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value,
            tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)
        protImportTs.setObjLabel(objLabel)
        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, ProtImportTsBase.OUTPUT_NAME, None)
        return tsImported

    @classmethod
    def _excludeSetViews(cls, inSet: SetOfTiltSeries, excludedViewsDict: Dict[str, List[int]]) -> None:
        """Mark specific tilt images as disabled (_objEnabled=False) in the
        given SetOfTiltSeries, modifying it in-place. Only TS whose tsId
        appears in excludedViewsDict (with a non-empty list) are modified."""
        objList = [obj.clone(ignoreAttrs=[]) for obj in inSet]
        for obj in objList:
            tsId = obj.getTsId()
            if tsId in excludedViewsDict and excludedViewsDict[tsId]:
                cls._excIntermediateSetViews(inSet, obj, excludedViewsDict[tsId])

    @staticmethod
    def _excIntermediateSetViews(inSet: SetOfTiltSeries,
                                 obj: TiltSeries,
                                 excludedViewsList: List[int]) -> None:
        tiList = [ti.clone() for ti in obj.iterItems(orderBy=TiltImage.INDEX_FIELD)]
        for i, ti in enumerate(tiList):
            if i in excludedViewsList:
                ti._objEnabled = False
                obj.update(ti)
        obj.write()
        inSet.update(obj)
        inSet.write()

    def _runAssignExcludedViews(self,
                                sourceTsSet: SetOfTiltSeries,
                                targetTsSet: SetOfTiltSeries,
                                objLabel: str = 'Assign excluded views') -> Optional[SetOfTiltSeries]:
        print(magentaStr(f"\n==> Running assign excluded views: {objLabel}"))
        prot = self.newProtocol(
            ProtAssignExcludedViews,
            inputSourceTiltSeries=sourceTsSet,
            inputTargetTiltSeries=targetTsSet)
        prot.setObjLabel(objLabel)
        self.launchProtocol(prot)
        outTsSet = getattr(prot, prot._possibleOutputs.tiltSeries.name, None)
        return outTsSet

    # ------------------------------------------------------------------
    # Test 1:
    # Source: set of 5 ts, exclusions on TS_03 (4 views) and TS_43 (5 views)
    # Target: the same set of 5 ts, unmodified.
    # ------------------------------------------------------------------
    def test_excludeViews_01(self):
        outTsSet = self._runAssignExcludedViews(
            self.tsSetAllWithEV,
            self.tsSetAll,
            objLabel='test_01')
        self.assertIsNotNone(outTsSet, "No output tilt-series set produced")

        self.checkTiltSeries(
            outTsSet,
            expectedSetSize=self.expectedSetSize5,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=self.expectedDimsDict5,
            testAcqObj=self.testAcqObjDict5,
            anglesCount=self.anglesCountDict5,
            isHeterogeneousSet=True,
            excludedViewsDict=self.excludedViewsTs03Ts43,
            presentTsIds=[TS_01, TS_03, TS_43, TS_45, TS_54])

    # ------------------------------------------------------------------
    # Test 2:
    # Source: subset (TS_03 + TS_54), exclusions on both (TS_03: 3 views, TS_54: 5 views).
    # Target: set of 5 ts, unmodified.
    # ------------------------------------------------------------------
    def test_excludeViews_02(self):
        outTsSet = self._runAssignExcludedViews(
            self.tsSubsetWithEV, self.tsSetAll,
            objLabel='test_02')
        self.assertIsNotNone(outTsSet, "No output tilt-series set produced")

        self.checkTiltSeries(
            outTsSet,
            expectedSetSize=self.expectedSetSize5,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=self.expectedDimsDict5,
            testAcqObj=self.testAcqObjDict5,
            anglesCount=self.anglesCountDict5,
            isHeterogeneousSet=True,
            excludedViewsDict=self.excludedViewsTs03Ts54,
            presentTsIds=[TS_01, TS_03, TS_43, TS_45, TS_54])

    # ------------------------------------------------------------------
    # Test 3:
    # Source: set of 5 ts, exclusions on TS_03 (4 views) and TS_43 (5 views)
    # Target: subset (TS_03 + TS_54)
    # ------------------------------------------------------------------
    def test_excludeViews_03(self):
        outTsSet = self._runAssignExcludedViews(
            self.tsSetAllWithEV, self.tsSubset,
            objLabel='test_03')
        self.assertIsNotNone(outTsSet, "No output tilt-series set produced")

        self.checkTiltSeries(
            outTsSet,
            expectedSetSize=self.expectedSetSize2,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=self.expectedDimsDict2,
            testAcqObj=self.testAcqObjDict2,
            anglesCount=self.anglesCountDict2,
            isHeterogeneousSet=True,
            excludedViewsDict=self.excludedViewsTs03Ts43,
            presentTsIds=[TS_03, TS_54])

    # ------------------------------------------------------------------
    # Test 4:
    # Source: set of 5 ts, exclusions on TS_03 (4 views) and TS_43 (5 views)
    # Target: subset (TS_03 + TS_54), exclusions on both (TS_03: 3 views, TS_54: 5 views).
    # ------------------------------------------------------------------
    def test_excludeViews_04(self):
        outTsSet = self._runAssignExcludedViews(
            self.tsSetAllWithEV, self.tsSubsetWithEV,
            objLabel='test_04')
        self.assertIsNotNone(outTsSet, "No output tilt-series set produced")

        self.checkTiltSeries(
            outTsSet,
            expectedSetSize=self.expectedSetSize2,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=self.expectedDimsDict2,
            testAcqObj=self.testAcqObjDict2,
            anglesCount=self.anglesCountDict2,
            isHeterogeneousSet=True,
            excludedViewsDict=self.excludedViewsTs03Ts43,
            presentTsIds=[TS_03, TS_54])

    # ------------------------------------------------------------------
    # Test 5:
    # Source:  subset (TS_03 + TS_54), exclusions on both (TS_03: 3 views, TS_54: 5 views).
    # Target: set of 5 ts, exclusions on TS_03 (4 views) and TS_43 (5 views).
    # ------------------------------------------------------------------
    def test_excludeViews_05(self):
        outTsSet = self._runAssignExcludedViews(
            self.tsSubsetWithEV, self.tsSetAllWithEV,
            objLabel='test_05')
        self.assertIsNotNone(outTsSet, "No output tilt-series set produced")

        # Check the results
        testExcludedViews = copy.deepcopy(self.excludedViewsTs03Ts54)
        testExcludedViews[TS_43] = self.excludedViewsTs03Ts43[TS_43]
        self.checkTiltSeries(
            outTsSet,
            expectedSetSize=self.expectedSetSize5,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=self.expectedDimsDict5,
            testAcqObj=self.testAcqObjDict5,
            anglesCount=self.anglesCountDict5,
            isHeterogeneousSet=True,
            excludedViewsDict=testExcludedViews,
            presentTsIds=[TS_01, TS_03, TS_43, TS_45, TS_54])

    # ------------------------------------------------------------------
    # Test 6:
    # 6.1:
    # Source: subset target (TS_03 + TS_54), exclusions on both (TS_03: 3 views, TS_54: 5 views),
    # then physically restack with IMOD's ProtImodExcludeViews to remove disabled images.
    # Target: set of 5 ts, unmodified.
    #
    # 6.2:
    # Source: set of 5 ts, exclusions on TS_03 (4 views) and TS_43 (5 views)
    # Target: subset target (TS_03 + TS_54), exclusions on both (TS_03: 3 views, TS_54: 5 views),
    # then physically restack with IMOD's ProtImodExcludeViews to remove disabled images.
    # ------------------------------------------------------------------
    def test_excludeViews_06(self):
        # Restack using IMOD excludeviews
        print(magentaStr("\n==> Restacking with IMOD:"))
        protRestack = self.newProtocol(
            ProtImodExcludeViews,
            inputSetOfTiltSeries=self.tsSubsetWithEV)
        protRestack.setObjLabel('Restack TS_03+TS_54')
        self.launchProtocol(protRestack)
        restackedTsSubset = getattr(protRestack, OUTPUT_TILTSERIES_NAME, None)
        self.assertIsNotNone(restackedTsSubset, "No restacked output produced")
        # Expected image counts after restack (original - excluded) and tilt-series dimensions
        restackedAnglesCount = {
            TS_03: self.anglesCountDict2[TS_03] - len(self.excludedViewsTs03Ts54[TS_03]),
            TS_54: self.anglesCountDict2[TS_54] - len(self.excludedViewsTs03Ts54[TS_54]),
        }
        expectedDimsRestack2 = copy.deepcopy(self.expectedDimsDict2)
        expectedDimsRestack2[TS_03][-1] = restackedAnglesCount[TS_03]
        expectedDimsRestack2[TS_54][-1] = restackedAnglesCount[TS_54]

        # Test 6.1: assign excluded views
        assignedTsSet = self._runAssignExcludedViews(
            restackedTsSubset, self.tsSetAll,
            objLabel='test_06.1')
        self.assertIsNotNone(assignedTsSet, "No output from assign-excluded-views step")

        # Check the results
        self.checkTiltSeries(
            assignedTsSet,
            expectedSetSize=self.expectedSetSize5,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=self.expectedDimsDict5,
            testAcqObj=self.testAcqObjDict5,
            anglesCount=self.anglesCountDict5,
            isHeterogeneousSet=True,
            excludedViewsDict=self.excludedViewsTs03Ts54,
            presentTsIds=[TS_01, TS_03, TS_43, TS_45, TS_54])

        # Test 6.2: assign excluded views
        assignedTsSet = self._runAssignExcludedViews(
            self.tsSetAllWithEV, restackedTsSubset,
            objLabel='test_06.2')
        self.assertIsNotNone(assignedTsSet, "No output from assign-excluded-views step")

        # Check the results
        testExcludedViews = copy.deepcopy(self.excludedViewsTs03Ts54)
        for tsId in testExcludedViews.keys():
            testExcludedViews[tsId] = []
        testExcludedViews[TS_03] = [0]  # Because of re-stack there is a re-indexation and this is the only present
        # # Also the acquisition needs to be built adapted to the test because of the re-stack
        testAcqObjDictReStacked = {}
        acq_TS_03 = DataSetRe4STATuto.testAcq03.value.clone()
        acq_TS_03.setAccumDose(111.)
        acq_TS_03.setAngleMin(-54)
        acq_TS_03.setAngleMax(54)
        testAcqObjDictReStacked[TS_03] = acq_TS_03

        acq_TS_54 = DataSetRe4STATuto.testAcq54.value.clone()
        acq_TS_54.setAccumDose(108.)
        acq_TS_54.setAngleMin(-54)
        acq_TS_54.setAngleMax(51)
        testAcqObjDictReStacked[TS_54] = acq_TS_54

        self.checkTiltSeries(
            assignedTsSet,
            expectedSetSize=self.expectedSetSize2,
            expectedSRate=self.unbinnedSRate,
            imported=True,
            expectedDimensions=expectedDimsRestack2,
            testAcqObj=testAcqObjDictReStacked,
            anglesCount=restackedAnglesCount,
            isHeterogeneousSet=True,
            excludedViewsDict=testExcludedViews,
            presentTsIds=[TS_03, TS_54])
