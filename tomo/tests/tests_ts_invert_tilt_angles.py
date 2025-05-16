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
import glob
from os.path import join
from typing import Union

import numpy as np

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr
from tomo.objects import SetOfTiltSeries
from tomo.protocols import ProtImportTs, ProtInverTiltAngles
from tomo.protocols.protocol_invert_tilt_angles import IN_TS_SET
from tomo.tests import DataSetRe4STATuto, RE4_STA_TUTO, TS_03, TS_54
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestTsInvertTiltAngles(TestBaseCentralizedLayer):
    ds = None
    importedTs = None
    expectedSRate = None
    expectedSetSize = None
    expectedAcqDict = None
    expectedTsDims = None
    expectedAnglesCount = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        tsIds = (TS_03, TS_54)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.expectedSetSize = len(tsIds)
        cls.expectedSRate = DataSetRe4STATuto.unbinnedPixSize.value
        cls.expectedAcqDict, cls.expectedTsDims, cls.expectedAnglesCount = DataSetRe4STATuto.genTestTsDicts(tsIds)
        cls.expectedTsDims = DataSetRe4STATuto.dimsTsBin1Dict.value
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        cls._runImportTs()
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

    @classmethod
    def _runImportTs(cls) -> None:
        print(magentaStr("\n==> Importing the tilt-series:"))
        protTsImport = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
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
        cls.importedTs = getattr(protTsImport, protTsImport.OUTPUT_NAME, None)

    @classmethod
    def _runInvertTiltAngles(cls) -> Union[SetOfTiltSeries, None]:
        print(magentaStr("\n==> Inverting the tilt-angles:"))
        argsDict = {IN_TS_SET: cls.importedTs}
        protInvertTa = cls.newProtocol(ProtInverTiltAngles, **argsDict)
        cls.launchProtocol(protInvertTa)
        return getattr(protInvertTa, protInvertTa._possibleOutputs.tiltSeries.name, None)

    @classmethod
    def _readAnglesFromTlt(cls, tltFile: str) -> np.array:
        tiltAngles = []
        with open(tltFile, 'r') as file:
            for line in file:
                values = line.split()  # The first column is the tilt angle, while the second is the dose: -57.00 117
                tiltAngles.append(float(values[0]))

        return np.array(tiltAngles)

    def testInvertTiltAngles(self) -> None:
        invTsSet = self._runInvertTiltAngles()
        # Check the tilt-series
        self.checkTiltSeries(invTsSet,
                             expectedSetSize=self.expectedSetSize,
                             expectedSRate=self.expectedSRate,
                             imported=True,
                             expectedDimensions=self.expectedTsDims,
                             testAcqObj=self.expectedAcqDict,
                             anglesCount=self.expectedAnglesCount,
                             isHeterogeneousSet=True)
        # Check the tilt angles
        tsPath = self.ds.getFile(DataSetRe4STATuto.tsPath.value)
        ts03TltFile = glob.glob(join(tsPath, TS_03, '03.tlt'))[0]
        ts54TltFile = glob.glob(join(tsPath, TS_54, '54.tlt'))[0]
        expectedTiltsDict = {
            TS_03: -1 * self._readAnglesFromTlt(ts03TltFile),
            TS_54: -1 * self._readAnglesFromTlt(ts54TltFile)
        }
        resultingTiltsDict = {}
        for ts in invTsSet:
            resultingTiltsDict[ts.getTsId()] = np.array([ti.getTiltAngle() for ti in ts])

        self.assertEqual(expectedTiltsDict.keys(), resultingTiltsDict.keys())
        for tsId in expectedTiltsDict.keys():
            self.assertTrue(np.all(abs(expectedTiltsDict[tsId] - resultingTiltsDict[tsId]) <= 1e-4))





