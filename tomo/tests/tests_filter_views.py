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
import unittest
from typing import Dict, Optional, Set, Tuple

from pwem import ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, weakImport
from tomo.objects import TiltImage, SetOfTiltSeries
from tomo.protocols import ProtImportTs, ProtExclViewFilter
from tomo.protocols.protocol_ts_exclude_views_filter import (
    QualityFilterModes, IN_TS_SET, MIN_TILT, MAX_TILT, MAX_SX, MAX_SY,
    MIN_DOSE, MAX_DOSE, QUALITY_FILTER, DARK_SENSITIVITY,
    DARK_LOCAL_RATION_TH, EXTREME_CONTRAST_TH, MIN_ENTROPY, MIN_EDGE_ENERGY,
    MIN_VIEWS, DO_RESTACK,
)
from tomo.tests import (
    FILTER_EXCLUDED_TS, DataSet_FilterExcludedTs, TS_POS6, TS_POS8,
)
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from tomo.utils import existsPlugin

with weakImport("imod"):
    from imod.constants import OUTPUT_TILTSERIES_NAME
    from imod.protocols import ProtImodImportTransformationMatrix

# ---------------------------------------------------------------------------
# Permissive custom-quality thresholds that effectively disable quality
# filtering.  QualityFilterModes.disabled leaves threshold attributes
# uninitialised (TypeError in _filterByImgQuality), so custom mode with
# permissive values is the workaround.
# ---------------------------------------------------------------------------
_QUALITY_DISABLED = {
    QUALITY_FILTER: QualityFilterModes.custom.value,
    DARK_SENSITIVITY: 6.0,
    DARK_LOCAL_RATION_TH: 0.30,
    EXTREME_CONTRAST_TH: 3.0,
    MIN_ENTROPY: 1.0,
    MIN_EDGE_ENERGY: 0.0,
}


# ===========================================================================
# Base class
# ===========================================================================
class _TestExclViewFilterBase(TestBaseCentralizedLayer):
    ds = None
    expectedSRate = DataSet_FilterExcludedTs.aPixBin4.value
    nTiltImages = DataSet_FilterExcludedTs.dimsTsBin4.value[2]  # 41

    @classmethod
    def setUpClass(cls) -> None:
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(FILTER_EXCLUDED_TS)

    # ------------------------------------------------------------------
    # Import helpers
    # ------------------------------------------------------------------
    @classmethod
    def _runImportTs(cls, motionCorrected: bool = True) -> Optional[SetOfTiltSeries]:
        if motionCorrected:
            filesPath = cls.ds.getFile(
                DataSet_FilterExcludedTs.tsMotionCorrBin4Dir.value)
        else:
            filesPath = cls.ds.getFile(
                DataSet_FilterExcludedTs.tsAliBin4Dir.value)
        print(magentaStr(
            f"\n==> Importing tilt series (mc={motionCorrected}):"))
        prot = cls.newProtocol(
            ProtImportTs,
            filesPath=filesPath,
            filesPattern='*.mdoc',
            voltage=DataSet_FilterExcludedTs.voltage.value,
            magnification=DataSet_FilterExcludedTs.magnification.value,
            sphericalAberration=DataSet_FilterExcludedTs.sphericalAb.value,
            amplitudeContrast=DataSet_FilterExcludedTs.amplitudeContrast.value,
            samplingRate=DataSet_FilterExcludedTs.aPixBin4.value,
            tiltAxisAngle=DataSet_FilterExcludedTs.tiltAxisAngle.value)
        cls.launchProtocol(prot)
        return getattr(prot, prot.OUTPUT_NAME, None)

    @classmethod
    def _runImportTM(cls, inTsSet: SetOfTiltSeries) -> Optional[SetOfTiltSeries]:
        filesPath = cls.ds.getFile(
            DataSet_FilterExcludedTs.tsAliBin4Dir.value)
        print(magentaStr("\n==> Importing transformation matrices:"))
        prot = cls.newProtocol(
            ProtImodImportTransformationMatrix,
            filesPath=filesPath,
            filesPattern='*.xf',
            inputSetOfTiltSeries=inTsSet,
            binningTM=1,
            binningTS=1)
        cls.launchProtocol(prot)
        return getattr(prot, OUTPUT_TILTSERIES_NAME, None)

    # ------------------------------------------------------------------
    # Filter helper
    # ------------------------------------------------------------------
    @classmethod
    def _runFilter(cls, inTsSet: SetOfTiltSeries, objLabel: str = 'filter', **kwargs) -> Tuple[Optional[SetOfTiltSeries], Optional[SetOfTiltSeries]]:
        defaults = {
            MIN_TILT: -100.0, MAX_TILT: 100.0,
            MAX_SX: 0.0, MAX_SY: 0.0,
            MIN_DOSE: 0.0, MAX_DOSE: 300.0,
            MIN_VIEWS: 1,
            DO_RESTACK: False,
        }
        defaults.update(_QUALITY_DISABLED)
        defaults.update(kwargs)
        print(magentaStr(f"\n==> Filter ({objLabel}):"))
        prot = cls.newProtocol(
            ProtExclViewFilter,
            **{IN_TS_SET: inTsSet},
            **defaults)
        prot.setObjLabel(objLabel)
        cls.launchProtocol(prot)
        outTs = getattr(prot, prot._possibleOutputs.tiltSeries.name, None)
        failedTs = getattr(prot, prot._possibleOutputs.failedTiltSeries.name, None)
        return outTs, failedTs

    # ------------------------------------------------------------------
    # Dynamic excluded-views extraction
    # ------------------------------------------------------------------
    @staticmethod
    def _getExcludedViewsDict(outTsSet: SetOfTiltSeries) -> Dict[str, Set[int]]:
        result = {}
        for ts in outTsSet:
            excluded = set()
            for ind, ti in enumerate(ts):
                if not ti.isEnabled():
                    excluded.add(ind)
            result[ts.getTsId()] = excluded
        return result

    # ------------------------------------------------------------------
    # Acquisition reference builders
    #
    # The expected acquisition objects are derived from the dataset's
    # characterised full-set acquisition (DataSet_FilterExcludedTs.tsAcqDict,
    # now holding the corrected tsAcqPos8 values) and from the filter's
    # documented behaviour -- never by reading back the values produced by the
    # protocol under test.
    # ------------------------------------------------------------------
    @classmethod
    def _genFullSetAcqDict(cls) -> dict:
        """Expected per-tilt-series acquisition for a non-re-stacked output.

        A non-re-stacked filter only disables views (the protocol does
        copyInfo, keeping the original acquisition), so the expected acquisition
        equals the full-set acquisition characterised in the dataset
        (DataSet_FilterExcludedTs.tsAcqDict), including the updated per-frame
        dose values now stored for tsAcqPos6/tsAcqPos8.  Clones are returned so
        callers may override fields (e.g. for the re-stacked case) without
        mutating the shared dataset objects.

        :return: dict {tsId: TomoAcquisition}.
        """
        return {tsId: acq.clone()
                for tsId, acq in DataSet_FilterExcludedTs.tsAcqDict.value.items()}

    @classmethod
    def _genReStackedAcqDict(cls, inTsSet: SetOfTiltSeries, excludedViewsDict: dict) -> dict:
        """Expected per-tilt-series acquisition after re-stacking.

        Re-stacking drops the excluded views and updates the acquisition
        accordingly (see ProtExclViewFilter._populateRestackedTs): angleMin and
        angleMax become the min/max tilt angle of the surviving views, accumDose
        the maximum accumulated dose among them and doseInitial the minimum
        initial dose.

        The expectation is computed **independently from the protocol output**:
        the surviving views are those NOT listed in ``excludedViewsDict`` -- the
        deterministic exclusion sets that the non-re-stacked sibling tests assert
        (i.e. the filter's documented behaviour) -- and their tilt angles and
        doses are read from the *input* ground-truth set.  View indices follow
        the tilt-angle order used by both the filter and ``excludedViewsDict``.

        :param inTsSet: input SetOfTiltSeries fed to the filter protocol.
        :param excludedViewsDict: {tsId: set/list of tilt-angle-sorted indices
            expected to be excluded by the filter}.
        :return: dict {tsId: TomoAcquisition}.
        """
        baseAcqDict = cls._genFullSetAcqDict()
        acqDict = {}
        for ts in inTsSet:
            tsId = ts.getTsId()
            excluded = set(excludedViewsDict[tsId])
            tiltAngles = []
            accumDoses = []
            initialDoses = []
            for ind, ti in enumerate(
                    ts.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD)):
                if ind not in excluded:
                    tiAcq = ti.getAcquisition()
                    tiltAngles.append(ti.getTiltAngle())
                    accumDoses.append(tiAcq.getAccumDose())
                    initialDoses.append(tiAcq.getDoseInitial())
            acq = baseAcqDict[tsId]
            acq.setAngleMin(min(tiltAngles))
            acq.setAngleMax(max(tiltAngles))
            acq.setAccumDose(max(accumDoses))
            acq.setDoseInitial(min(initialDoses))
            acqDict[tsId] = acq
        return acqDict


# ===========================================================================
# SCENARIO 1 — Motion-corrected tilt-series (no alignment data)
# ===========================================================================
class TestExclViewFilterMC(_TestExclViewFilterBase):
    """Scenario 1: motion-corrected TS, no alignment.
    Max-shift filter is NOT applicable (no Transform on TiltImages).
    All assertions use checkTiltSeries with excludedViewsDict."""
    importedTs = None

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.importedTs = cls._runImportTs(motionCorrected=True)
        cls.fullSetAcqDict = cls._genFullSetAcqDict()

    # ==================================================================
    # Single-criterion deterministic tests
    # ==================================================================
    def test_01_tiltAngle(self) -> None:
        """Tilt angle filter [-50, 50].  Views at tilt-angle-sorted
        indices 0-6 (angles -70.01° to -52.01°) are excluded."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='tilt [-50,50]',
            **{MIN_TILT: -50.0, MAX_TILT: 50.0})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict={
                TS_POS6: {0, 1, 2, 3, 4, 5, 6},
                TS_POS8: {0, 1, 2, 3, 4, 5, 6},
            }
        )

    def test_02_dose(self) -> None:
        """Dose filter maxDose=65.  Late-acquired views with accumulated
        dose > 65 e/A^2 are excluded.  POS6 and POS8 differ because of
        different per-view doses in the mdoc."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='dose max=65',
            **{MAX_DOSE: 65.0})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict={
                TS_POS6: {0, 1, 2, 3, 4, 36, 37, 38, 39, 40},
                TS_POS8: {0, 1, 2, 38, 39, 40},
            },
        )

    # ==================================================================
    # Quality-filter mode tests (all four modes)
    # ==================================================================
    def test_03_qualityConservative(self) -> None:
        """Quality filter — conservative mode.  Dynamic extraction
        because quality results depend on image content analysis."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='quality conservative',
            **{QUALITY_FILTER: QualityFilterModes.conservative.value})
        evd = self._getExcludedViewsDict(outTs)
        totalExcluded = sum(len(v) for v in evd.values())
        self.assertGreater(totalExcluded, 0,
                           "Conservative quality should exclude at least "
                           "some views (dataset has dark frames)")
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict=evd,
        )

    def test_04_qualityBalanced(self) -> None:
        """Quality filter — balanced mode (default)."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='quality balanced',
            **{QUALITY_FILTER: QualityFilterModes.balanced.value})
        evd = self._getExcludedViewsDict(outTs)
        totalExcluded = sum(len(v) for v in evd.values())
        self.assertGreater(totalExcluded, 0)
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict=evd,
        )

    def test_05_qualityAggressive(self) -> None:
        """Quality filter — aggressive mode."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='quality aggressive',
            **{QUALITY_FILTER: QualityFilterModes.aggressive.value})
        evd = self._getExcludedViewsDict(outTs)
        totalExcluded = sum(len(v) for v in evd.values())
        self.assertGreater(totalExcluded, 0)
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict=evd,
        )

    def test_06_qualityCustom(self) -> None:
        """Quality filter — custom mode with restrictive thresholds
        on all four quality components."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='quality custom',
            **{QUALITY_FILTER: QualityFilterModes.custom.value,
               DARK_SENSITIVITY: 2.0,
               DARK_LOCAL_RATION_TH: 0.65,
               EXTREME_CONTRAST_TH: 1.0,
               MIN_ENTROPY: 3.0,
               MIN_EDGE_ENERGY: 0.005})
        evd = self._getExcludedViewsDict(outTs)
        totalExcluded = sum(len(v) for v in evd.values())
        self.assertGreater(totalExcluded, 0)
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict=evd,
        )

    # ==================================================================
    # Combined criteria
    # ==================================================================
    def test_07_combinedTiltDose(self) -> None:
        """Combined tilt [-50, 50] + dose max=65.  The excluded set is
        the union of each criterion's exclusion set."""
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='tilt+dose',
            **{MIN_TILT: -50.0, MAX_TILT: 50.0, MAX_DOSE: 65.0})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict={
                TS_POS6: {0, 1, 2, 3, 4, 5, 6, 36, 37, 38, 39, 40},
                TS_POS8: {0, 1, 2, 3, 4, 5, 6, 38, 39, 40},
            },
        )

    # ==================================================================
    # Re-stack tests
    # ==================================================================
    def test_08_reStackTiltAngle(self) -> None:
        """Re-stack with tilt [-50, 50].  Both TS keep 34 views
        (41 - 7 excluded).  All output views are enabled."""
        # Same views the tilt [-50, 50] filter excludes in the non-re-stacked
        # sibling test_01 (tilt-angle-sorted indices 0-6, angles < -50).
        restackExcludedViews = {
            TS_POS6: {0, 1, 2, 3, 4, 5, 6},
            TS_POS8: {0, 1, 2, 3, 4, 5, 6},
        }
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='restack tilt [-50,50]',
            **{MIN_TILT: -50.0, MAX_TILT: 50.0, DO_RESTACK: True})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=34,
            testAcqObj=self._genReStackedAcqDict(
                self.importedTs, restackExcludedViews),
            excludedViewsDict={TS_POS6: set(), TS_POS8: set()},
        )

    def test_09_reStackCombined(self) -> None:
        """Re-stack with tilt [-50, 50] + dose max=65.  Heterogeneous
        output: POS6 keeps 29, POS8 keeps 31 views."""
        # Same views the tilt [-50, 50] + dose max=65 filter excludes in the
        # non-re-stacked sibling test_07 (union of the tilt and dose criteria).
        restackExcludedViews = {
            TS_POS6: {0, 1, 2, 3, 4, 5, 6, 36, 37, 38, 39, 40},
            TS_POS8: {0, 1, 2, 3, 4, 5, 6, 38, 39, 40},
        }
        outTs, _ = self._runFilter(
            self.importedTs, objLabel='restack tilt+dose',
            **{MIN_TILT: -50.0, MAX_TILT: 50.0, MAX_DOSE: 65.0,
               DO_RESTACK: True})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount={TS_POS6: 29, TS_POS8: 31},
            testAcqObj=self._genReStackedAcqDict(
                self.importedTs, restackExcludedViews),
            excludedViewsDict={TS_POS6: set(), TS_POS8: set()},
            isHeterogeneousSet=True,
        )

    def test_10_reStackedAsInput(self) -> None:
        """Re-stacked TS as input to a second filter run.
        Step 1: re-stack with tilt [-50, 50] -> 34-view TS.
        Step 2: filter the 34-view TS with balanced quality."""
        restacked, _ = self._runFilter(
            self.importedTs, objLabel='restack step1',
            **{MIN_TILT: -50.0, MAX_TILT: 50.0, DO_RESTACK: True})
        self.assertIsNotNone(restacked)
        # Step 1 re-stacks tilt [-50, 50] (excludes indices 0-6, as in test_08),
        # so its acquisition is derived from the imported set + that exclusion.
        # The second pass does not re-stack, so it carries that acquisition over.
        step1ExcludedViews = {
            TS_POS6: {0, 1, 2, 3, 4, 5, 6},
            TS_POS8: {0, 1, 2, 3, 4, 5, 6},
        }
        reStackedAcqDict = self._genReStackedAcqDict(
            self.importedTs, step1ExcludedViews)
        outTs, _ = self._runFilter(
            restacked, objLabel='restack->quality balanced',
            **{QUALITY_FILTER: QualityFilterModes.balanced.value})
        evd = self._getExcludedViewsDict(outTs)
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            imported=True,
            anglesCount=34,
            testAcqObj=reStackedAcqDict,
            excludedViewsDict=evd,
        )


# ===========================================================================
# SCENARIO 2 — Aligned tilt-series (with alignment via .xf files)
# ===========================================================================
@unittest.skipIf(not existsPlugin('imod'), 'IMOD plugin not available')
class TestExclViewFilterAligned(_TestExclViewFilterBase):
    """Scenario 2: aligned TS with Transform matrices from .xf files.
    All Scenario 1 filters apply, plus max-shift.
    All assertions use checkTiltSeries with excludedViewsDict."""
    importedTs = None
    alignedTs = None

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.importedTs = cls._runImportTs(motionCorrected=False)
        cls.alignedTs = cls._runImportTM(cls.importedTs)
        cls.fullSetAcqDict = cls._genFullSetAcqDict()

    def test_01_maxShift(self) -> None:
        """Max-shift filter at 20% of image dimension.
        POS6: views 0-6 have large alignment shifts (up to 2830 px).
        POS8: views 0-1 have large shifts (3160, 1850 px)."""
        outTs, _ = self._runFilter(
            self.alignedTs, objLabel='shift 20%',
            **{MAX_SX: 0.2, MAX_SY: 0.2})
        evd = self._getExcludedViewsDict(outTs)
        self.assertTrue(
            {0, 1, 2, 3, 4, 5, 6}.issubset(evd[TS_POS6]),
            "POS6 indices 0-6 must be excluded (shifts > 204 px)")
        self.assertTrue(
            {0, 1}.issubset(evd[TS_POS8]),
            "POS8 indices 0-1 must be excluded (shifts > 204 px)")
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            hasAlignment=True,
            alignment=ALIGN_2D,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict=evd,
        )

    def test_02_combinedShiftTilt(self) -> None:
        """Combined shift 5% + tilt [-60, 60].  Deterministic indices.
        POS6: shift excludes 0-11 (dy up to 726 px at threshold 51.2).
        POS8: shift 0-1 + tilt 0-3 + shift idx 22 (-4 deg, dx=52.5)
        and idx 40 (50 deg, dx=-99.4)."""
        outTs, _ = self._runFilter(
            self.alignedTs, objLabel='shift 5% + tilt [-60,60]',
            **{MAX_SX: 0.05, MAX_SY: 0.05,
               MIN_TILT: -60.0, MAX_TILT: 60.0})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            hasAlignment=True,
            alignment=ALIGN_2D,
            anglesCount=self.nTiltImages,
            testAcqObj=self.fullSetAcqDict,
            excludedViewsDict={
                TS_POS6: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                TS_POS8: {0, 1, 2, 3, 22, 40},
            },
        )

    def test_03_reStackFromAligned(self) -> None:
        """Re-stack with shift 5% + tilt [-60, 60].  Heterogeneous
        output: POS6 keeps 29 views, POS8 keeps 35."""
        # Same views the shift 5% + tilt [-60, 60] filter excludes in the
        # non-re-stacked sibling test_02 (union of the shift and tilt criteria).
        restackExcludedViews = {
            TS_POS6: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
            TS_POS8: {0, 1, 2, 3, 22, 40},
        }
        outTs, _ = self._runFilter(
            self.alignedTs, objLabel='restack shift+tilt',
            **{MAX_SX: 0.05, MAX_SY: 0.05,
               MIN_TILT: -60.0, MAX_TILT: 60.0,
               DO_RESTACK: True})
        self.checkTiltSeries(
            outTs,
            expectedSetSize=2,
            expectedSRate=self.expectedSRate,
            hasAlignment=True,
            alignment=ALIGN_2D,
            anglesCount={TS_POS6: 29, TS_POS8: 35},
            testAcqObj=self._genReStackedAcqDict(
                self.alignedTs, restackExcludedViews),
            excludedViewsDict={TS_POS6: set(), TS_POS8: set()},
            isHeterogeneousSet=True,
        )
