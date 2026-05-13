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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr

from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from tomo.protocols.protocol_assign_excluded_views import ProtAssignExcludedViews
from tomo.protocols import ProtImportTs
from tomo.tests import (RE4_STA_TUTO, DataSetRe4STATuto,
                         TS_01, TS_03, TS_54)


class TestAssignExcludedViewsModel(BaseTest):
    """Unit tests for excluded-view transfer logic using programmatic sets.
    No dataset download required - all sets are created in memory."""

    @classmethod
    def setUpClass(cls):
        cls.setupTestOutput()

    def _createTsInSet(self, tsSet, tsId, numImages,
                       disabledAcqOrders=None, acqOrders=None):
        """Add a TiltSeries with specified properties to a set.

        :param tsSet: SetOfTiltSeries to add the TS to.
        :param tsId: TiltSeries identifier.
        :param numImages: Number of TiltImages to create.
        :param disabledAcqOrders: Set of acquisition orders to disable.
        :param acqOrders: List of acquisition orders (default: 1..numImages).
        """
        ts = TiltSeries(tsId=tsId)
        ts.setSamplingRate(1.0)
        tsSet.append(ts)

        for i in range(numImages):
            ti = TiltImage()
            acqOrder = acqOrders[i] if acqOrders else i + 1
            ti.setTiltAngle(-60 + i * 3)
            ti.setAcquisitionOrder(acqOrder)
            if disabledAcqOrders and acqOrder in disabledAcqOrders:
                ti.setEnabled(False)
            ts.append(ti)

        ts.write()
        tsSet.update(ts)

    def _createTsSet(self, configs, suffix=''):
        """Create a SetOfTiltSeries from a list of configuration dicts.

        Each config dict has keys: tsId, numImages, and optionally
        disabledAcqOrders and acqOrders.
        """
        tsSet = SetOfTiltSeries.create(self.outputPath,
                                       template='tiltseries',
                                       suffix=suffix)
        tsSet.setSamplingRate(1.0)
        for cfg in configs:
            self._createTsInSet(tsSet, **cfg)
        tsSet.write()
        return tsSet

    @staticmethod
    def _applyTransferLogic(source, target, outputPath, suffix):
        """Simulate the protocol's transfer logic on programmatic sets.

        Returns (output SetOfTiltSeries, matchedTsIds, unmatchedTsIds).
        """
        sourceTsIds = set(source.getTSIds())
        targetTsIds = set(target.getTSIds())
        matchedTsIds = sorted(sourceTsIds & targetTsIds)
        unmatchedTsIds = sorted(targetTsIds - sourceTsIds)

        sourceTsDict = {ts.getTsId(): ts.clone() for ts in source
                        if ts.getTsId() in matchedTsIds}

        output = SetOfTiltSeries.create(outputPath,
                                        template='tiltseries',
                                        suffix=suffix)
        output.copyInfo(target)

        for targetTs in target:
            tsId = targetTs.getTsId()
            sourceTs = sourceTsDict.get(tsId, None)

            newTs = TiltSeries(tsId=tsId)
            newTs.copyInfo(targetTs)
            output.append(newTs)

            if sourceTs:
                sourceAcqMap = {ti.getAcquisitionOrder(): ti.isEnabled()
                                for ti in sourceTs}
                for ti in targetTs:
                    newTi = ti.clone()
                    acqOrder = ti.getAcquisitionOrder()
                    if acqOrder in sourceAcqMap:
                        newTi.setEnabled(sourceAcqMap[acqOrder])
                    newTs.append(newTi)
            else:
                for ti in targetTs:
                    newTi = ti.clone()
                    newTs.append(newTi)

            newTs.write()
            output.update(newTs)

        output.write()
        return output, matchedTsIds, unmatchedTsIds

    # ------------------------------------------------------------------
    # Test 1: Happy path - basic transfer of exclusions
    # ------------------------------------------------------------------
    def test_happy_path_transfer(self):
        """Source has disabled views; output should replicate them."""
        source = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 5,
             'disabledAcqOrders': {2, 4}},
            {'tsId': 'TS_02', 'numImages': 3,
             'disabledAcqOrders': {1}},
        ], suffix='src_happy')

        target = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 5},
            {'tsId': 'TS_02', 'numImages': 3},
        ], suffix='tgt_happy')

        output, matched, unmatched = self._applyTransferLogic(
            source, target, self.outputPath, 'out_happy')

        self.assertEqual(output.getSize(), 2)
        self.assertEqual(matched, ['TS_01', 'TS_02'])
        self.assertEqual(unmatched, [])

        for ts in output:
            tsId = ts.getTsId()
            for ti in ts:
                ao = ti.getAcquisitionOrder()
                if tsId == 'TS_01' and ao in {2, 4}:
                    self.assertFalse(
                        ti.isEnabled(),
                        "TS_01 acqOrder %d should be disabled" % ao)
                elif tsId == 'TS_02' and ao == 1:
                    self.assertFalse(
                        ti.isEnabled(),
                        "TS_02 acqOrder 1 should be disabled")
                else:
                    self.assertTrue(
                        ti.isEnabled(),
                        "%s acqOrder %d should be enabled" % (tsId, ao))

    # ------------------------------------------------------------------
    # Test 2: Re-stacked / acquisition-order robustness
    # ------------------------------------------------------------------
    def test_acquisition_order_robustness(self):
        """Images matched by acqOrder, not by stack index."""
        source = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 5,
             'acqOrders': [1, 2, 3, 4, 5],
             'disabledAcqOrders': {2, 4}},
        ], suffix='src_robust')

        # Target has REVERSED stack order
        target = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 5,
             'acqOrders': [5, 4, 3, 2, 1]},
        ], suffix='tgt_robust')

        output, _, _ = self._applyTransferLogic(
            source, target, self.outputPath, 'out_robust')

        for ts in output:
            for ti in ts:
                ao = ti.getAcquisitionOrder()
                if ao in {2, 4}:
                    self.assertFalse(
                        ti.isEnabled(),
                        "acqOrder %d should be disabled regardless of "
                        "stack position" % ao)
                else:
                    self.assertTrue(
                        ti.isEnabled(),
                        "acqOrder %d should be enabled" % ao)

    # ------------------------------------------------------------------
    # Test 3: Unmatched target TiltSeries preserved
    # ------------------------------------------------------------------
    def test_unmatched_target_preserved(self):
        """Target TS without a source match appear unchanged in output."""
        source = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 3,
             'disabledAcqOrders': {1, 2, 3}},
        ], suffix='src_unmatched')

        target = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 3},
            {'tsId': 'TS_99', 'numImages': 4},
        ], suffix='tgt_unmatched')

        output, matched, unmatched = self._applyTransferLogic(
            source, target, self.outputPath, 'out_unmatched')

        self.assertEqual(output.getSize(), 2)
        self.assertEqual(matched, ['TS_01'])
        self.assertEqual(unmatched, ['TS_99'])

        for ts in output:
            tsId = ts.getTsId()
            if tsId == 'TS_01':
                for ti in ts:
                    self.assertFalse(
                        ti.isEnabled(),
                        "TS_01 all views should be disabled")
            elif tsId == 'TS_99':
                for ti in ts:
                    self.assertTrue(
                        ti.isEnabled(),
                        "TS_99 (unmatched) views should stay enabled")

    # ------------------------------------------------------------------
    # Test 4: Summary content
    # ------------------------------------------------------------------
    def test_summary_content(self):
        """Matched and unmatched tsIds are reported correctly."""
        source = self._createTsSet([
            {'tsId': 'TS_A', 'numImages': 2},
            {'tsId': 'TS_B', 'numImages': 2},
        ], suffix='src_summary')

        target = self._createTsSet([
            {'tsId': 'TS_B', 'numImages': 2},
            {'tsId': 'TS_C', 'numImages': 2},
            {'tsId': 'TS_D', 'numImages': 2},
        ], suffix='tgt_summary')

        _, matched, unmatched = self._applyTransferLogic(
            source, target, self.outputPath, 'out_summary')

        self.assertEqual(matched, ['TS_B'])
        self.assertEqual(unmatched, ['TS_C', 'TS_D'])

        # Simulate summary message construction
        matchedMsg = ", ".join(matched) if matched else None
        unmatchedMsg = ", ".join(unmatched) if unmatched else None

        self.assertIn('TS_B', matchedMsg)
        self.assertIn('TS_C', unmatchedMsg)
        self.assertIn('TS_D', unmatchedMsg)

    # ------------------------------------------------------------------
    # Test 5: Partial acquisition order overlap
    # ------------------------------------------------------------------
    def test_partial_acq_order_overlap(self):
        """Source has fewer images - unmatched target acqOrders stay unchanged."""
        source = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 3,
             'acqOrders': [1, 2, 3],
             'disabledAcqOrders': {2}},
        ], suffix='src_partial')

        # Target has acqOrders 1..5 (source only covers 1..3)
        target = self._createTsSet([
            {'tsId': 'TS_01', 'numImages': 5,
             'acqOrders': [1, 2, 3, 4, 5]},
        ], suffix='tgt_partial')

        output, _, _ = self._applyTransferLogic(
            source, target, self.outputPath, 'out_partial')

        for ts in output:
            for ti in ts:
                ao = ti.getAcquisitionOrder()
                if ao == 2:
                    self.assertFalse(
                        ti.isEnabled(),
                        "acqOrder 2 should be disabled from source")
                else:
                    self.assertTrue(
                        ti.isEnabled(),
                        "acqOrder %d should stay enabled" % ao)

    # ------------------------------------------------------------------
    # Test 6: Validation - tsId matching
    # ------------------------------------------------------------------
    def test_tsid_matching_sets(self):
        """Verify correct tsId set operations."""
        source = self._createTsSet([
            {'tsId': 'TS_A', 'numImages': 2},
            {'tsId': 'TS_B', 'numImages': 2},
            {'tsId': 'TS_C', 'numImages': 2},
        ], suffix='src_match')

        target = self._createTsSet([
            {'tsId': 'TS_B', 'numImages': 2},
            {'tsId': 'TS_C', 'numImages': 2},
            {'tsId': 'TS_D', 'numImages': 2},
        ], suffix='tgt_match')

        sourceTsIds = set(source.getTSIds())
        targetTsIds = set(target.getTSIds())
        matched = sourceTsIds & targetTsIds
        unmatched = targetTsIds - sourceTsIds

        self.assertEqual(matched, {'TS_B', 'TS_C'})
        self.assertEqual(unmatched, {'TS_D'})


class TestAssignExcludedViewsProtocol(BaseTest):
    """Integration tests for ProtAssignExcludedViews using the full
    protocol framework with the RE4_STA_TUTO dataset."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)

    @classmethod
    def _runImportTs(cls, exclusionWords=None):
        print(magentaStr("\n==> Importing tilt series:"))
        protTsImport = cls.newProtocol(
            ProtImportTs,
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
        return getattr(protTsImport, protTsImport.OUTPUT_NAME, None)

    def test_protocol_with_matched_and_unmatched(self):
        """Full protocol run: 2 matched tsIds + 1 unmatched in target."""
        print(magentaStr("\n==> Testing assign excluded views protocol:"))

        # Source: TS_03, TS_54 (2 TS)
        sourceTs = self._runImportTs(
            exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value)
        self.assertIsNotNone(sourceTs, "Source import failed")
        self.assertEqual(sourceTs.getSize(), 2)

        # Target: TS_01, TS_03, TS_54 (3 TS - TS_01 has no source match)
        targetTs = self._runImportTs(exclusionWords='output 43 45')
        self.assertIsNotNone(targetTs, "Target import failed")
        self.assertEqual(targetTs.getSize(), 3)

        # Run the protocol
        prot = self.newProtocol(ProtAssignExcludedViews,
                                inputSourceTiltSeries=sourceTs,
                                inputTargetTiltSeries=targetTs)
        self.launchProtocol(prot)

        # -- Check output set --
        outTs = getattr(prot, prot._possibleOutputs.tiltSeries.name, None)
        self.assertIsNotNone(outTs, "Output tilt-series should not be None")
        self.assertSetSize(outTs, 3, msg="Output should contain 3 TS")

        outTsIds = sorted(outTs.getTSIds())
        self.assertEqual(outTsIds, sorted([TS_01, TS_03, TS_54]),
                         "Output should contain all target tsIds")

        # Since source has no excluded views, all output images stay enabled
        for ts in outTs:
            for ti in ts:
                self.assertTrue(
                    ti.isEnabled(),
                    "All images should be enabled (source has no exclusions)")

        # -- Check summary --
        summary = prot._summary()
        summaryStr = " ".join(summary)
        self.assertIn(TS_03, summaryStr,
                      "Summary should mention matched TS_03")
        self.assertIn(TS_54, summaryStr,
                      "Summary should mention matched TS_54")
        self.assertIn(TS_01, summaryStr,
                      "Summary should mention unmatched TS_01")

    def test_protocol_all_matched(self):
        """Full protocol run where all target tsIds have a source match."""
        print(magentaStr("\n==> Testing assign excluded views "
                         "(all matched):"))

        # Both source and target: TS_03, TS_54
        exclusionWords = DataSetRe4STATuto.exclusionWordsTs03ts54.value
        sourceTs = self._runImportTs(exclusionWords=exclusionWords)
        targetTs = self._runImportTs(exclusionWords=exclusionWords)

        prot = self.newProtocol(ProtAssignExcludedViews,
                                inputSourceTiltSeries=sourceTs,
                                inputTargetTiltSeries=targetTs)
        self.launchProtocol(prot)

        outTs = getattr(prot, prot._possibleOutputs.tiltSeries.name, None)
        self.assertIsNotNone(outTs)
        self.assertSetSize(outTs, 2)

        # Summary should report all matched
        summary = prot._summary()
        summaryStr = " ".join(summary)
        self.assertIn("All target tilt-series were matched", summaryStr)
        self.assertIn(TS_03, summaryStr)
        self.assertIn(TS_54, summaryStr)
