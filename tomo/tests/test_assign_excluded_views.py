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
from typing import Optional, Dict, List

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, weakImport
from tomo.objects import SetOfTiltSeries, TiltSeries
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
    excludedViewsAll5 = {
        TS_01: [],
        TS_03: [0, 1, 38, 39],
        TS_43: [0, 1, 2, 39, 40],
        TS_45: [],
        TS_54: [],
    }

    # Tests 2 & 3: subset TS_03 + TS_54; exclude on both (3 and 5 views)
    excludedViewsSubset = {
        TS_03: [0, 38, 39],
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
        cls.sourceAll = cls._runImportTs(objLabel='source: all 5 TS (excl)')
        cls._excludeSetViews(cls.sourceAll, cls.excludedViewsAll5)

        # Target for tests 2 & 3: subset TS_03 + TS_54 (clean)
        cls.tsSubset = cls._runImportTs(
            exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
            objLabel='target: TS_03+TS_54')

        # Source for tests 2 & 3: TS_03 + TS_54 with exclusions
        cls.sourceSubset = cls._runImportTs(
            exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
            objLabel='source: TS_03+TS_54 (excl)')
        cls._excludeSetViews(cls.sourceSubset, cls.excludedViewsSubset)

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
                                 excludedViewsList: Dict[str, List[int]]) -> None:
        tiList = [ti.clone() for ti in obj]
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
                                objLabel: str ='Assign excluded views') -> Optional[SetOfTiltSeries]:
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
    # Test 1: all 5 target TS, exclusions on TS_03 (4 views) and TS_43
    #         (5 views). TS_01, TS_45, TS_54 should remain unmodified.
    # ------------------------------------------------------------------
    def test_excludeViews_all5_exclude_TS03_TS43(self):
        outTsSet = self._runAssignExcludedViews(
            self.sourceAll, self.tsSetAll,
            objLabel='All 5, excl TS_03+TS_43')
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
            excludedViewsDict=self.excludedViewsAll5,
            presentTsIds=[TS_01, TS_03, TS_43, TS_45, TS_54])

    # ------------------------------------------------------------------
    # Test 2: subset target (TS_03 + TS_54), exclusions on both
    #         (TS_03: 3 views, TS_54: 5 views).
    # ------------------------------------------------------------------
    def test_excludeViews_subset_TS03_TS54(self):
        outTsSet = self._runAssignExcludedViews(
            self.sourceSubset, self.tsSubset,
            objLabel='Subset TS_03+TS_54, excl both')
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
            excludedViewsDict=self.excludedViewsSubset,
            presentTsIds=[TS_03, TS_54])

    # ------------------------------------------------------------------
    # Test 3: same exclusions as Test 2, then physically restack with
    #         IMOD's ProtImodExcludeViews to remove disabled images.
    # ------------------------------------------------------------------
    def test_excludeViews_subset_and_restack_TS03_TS54(self):
        # Step 1: assign excluded views
        assignedTsSet = self._runAssignExcludedViews(
            self.sourceSubset, self.tsSubset,
            objLabel='Assign excl before restack')
        self.assertIsNotNone(assignedTsSet,
                             "No output from assign-excluded-views step")

        # Step 2: restack using IMOD excludeviews
        print(magentaStr("\n==> Restacking with IMOD ProtImodExcludeViews:"))
        protRestack = self.newProtocol(
            ProtImodExcludeViews,
            inputSetOfTiltSeries=assignedTsSet)
        protRestack.setObjLabel('Restack TS_03+TS_54')
        self.launchProtocol(protRestack)
        restackedTsSet = getattr(protRestack, OUTPUT_TILTSERIES_NAME, None)
        self.assertIsNotNone(restackedTsSet, "No restacked output produced")

        # Expected image counts after restack (original - excluded)
        restackedAnglesCount = {
            TS_03: self.anglesCountDict2[TS_03] - len(self.excludedViewsSubset[TS_03]),
            TS_54: self.anglesCountDict2[TS_54] - len(self.excludedViewsSubset[TS_54]),
        }

        # Verify set size
        self.assertSetSize(restackedTsSet, self.expectedSetSize2)

        for ts in restackedTsSet:
            tsId = ts.getTsId()
            expectedCount = restackedAnglesCount[tsId]

            # Restacked TS should have fewer images
            self.assertEqual(
                ts.getSize(), expectedCount,
                f"{tsId}: expected {expectedCount} images after restack, "
                f"got {ts.getSize()}")

            # All remaining images must be enabled (disabled ones were removed)
            for ti in ts:
                self.assertTrue(
                    ti.isEnabled(),
                    f"{tsId}: restacked TS should only contain enabled images")

            # TODO: if ProtImodExcludeViews generates new binary files,
            # verify that file dimensions reflect the reduced stack size.
