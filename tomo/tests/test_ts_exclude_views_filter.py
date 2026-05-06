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
from typing import Optional
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.objects import SetOfTiltSeries
from tomo.protocols import ProtImportTs, ProtExclViewFilter
from tomo.tests import FILTER_EXCLUDED_TS, DataSet_FilterExcludedTs
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestTsExcludeViewsFilterBase(TestBaseCentralizedLayer):
    ds = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(FILTER_EXCLUDED_TS)

    @classmethod
    def _runImportTs(cls,
                     motionCorrected: bool = True,
                     exclusionWords: Optional[str] = None) -> Optional[SetOfTiltSeries]:
        if motionCorrected:
            strMsg = 'motion-corrected'
            filesPath = cls.ds.getFile(DataSet_FilterExcludedTs.tsMotionCorrBin4Dir.value)
        else:
            strMsg = 'with alignment metadata'
            filesPath = cls.ds.getFile(DataSet_FilterExcludedTs.tsAliBin4Dir.value)

        print(magentaStr(f"\n==> Importing the tilt series ({strMsg}):"))
        protTsImport = cls.newProtocol(ProtImportTs,
                                       objLabel=f'import ts {strMsg}',
                                       filesPath=filesPath,
                                       filesPattern='*.mdoc',
                                       exclusionWords=exclusionWords,
                                       voltage=DataSet_FilterExcludedTs.voltage.value,
                                       magnification=DataSet_FilterExcludedTs.magnification.value,
                                       sphericalAberration=DataSet_FilterExcludedTs.sphericalAb.value,
                                       amplitudeContrast=DataSet_FilterExcludedTs.amplitudeContrast.value,
                                       samplingRate=DataSet_FilterExcludedTs.aPixBin4.value,
                                       tiltAxisAngle=DataSet_FilterExcludedTs.tiltAxisAngle.value)

        cls.launchProtocol(protTsImport)
        tsImported = getattr(protTsImport, protTsImport.OUTPUT_NAME, None)
        return tsImported

    @classmethod
    def _runTsExcludeViewsFilter(cls,
                                 inTsSet: Optional[SetOfTiltSeries],
                                 maxSx: float = 0.,
                                 maxSy: float = 0.,
                                 minTilt: float = -100.,
                                 maxTilt: float = 100.,
                                 minDose: float = 0.,
                                 maxDose: float = 300.,
                                 darkFactor: float = 2.,
                                 minNoTilts: int = 50,
                                 doReStack: bool = False) -> Optional[SetOfTiltSeries]:
        print(magentaStr(f"\n==> Excluding the views:"
                         f"\n\t- MaxSx [%] = {maxSx}"
                         f"\n\t- MaxSy [%] = {maxSy}"
                         f"\n\t- MinTilt [deg.] = {minTilt}"
                         f"\n\t- MaxTilt [deg.] = {maxTilt}"
                         f"\n\t- MinDose [e/A²] = {minDose}"
                         f"\n\t- MaxDose [e/A²] = {maxDose}"
                         f"\n\t- Dark sensitivity = {darkFactor}"
                         f"\n\t- Min no. tilts = {minNoTilts}"
                         f"\n\t- Re-stack = {doReStack}",
                         ))
        protExcViewsFilter = cls.newProtocol(ProtExclViewFilter,
                                             inTsSet=inTsSet,
                                             minTilt=minTilt,
                                             maxTilt=maxTilt,
                                             maxShiftX=maxSx,
                                             maxShiftY=maxSy,
                                             minDose=minDose,
                                             maxDose=maxDose,
                                             darkSensitivity=darkFactor,
                                             minViews=minNoTilts,
                                             doReStack=doReStack)
        cls.launchProtocol(protExcViewsFilter)
        outTsSet = getattr(protExcViewsFilter, protExcViewsFilter._possibleOutputs.tiltSeries.name, None)
        return outTsSet


class TestTsExcludeViewsFilterMC(TestTsExcludeViewsFilterBase):

    def testTsExcludeViewsFilterMC_01(self):
        importedTsSet = self._runImportTs()
        self._runTsExcludeViewsFilter(importedTsSet)


