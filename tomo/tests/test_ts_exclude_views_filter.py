# **************************************************************************
# *
# * Authors:     Alberto García Mena (alberto.garcia@cnb.csic.es)
# *# *
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

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSetRe4STATuto, RE4_STA_TUTO, TS_03, TS_54
from tomo.objects import SetOfTiltSeries
from . import DataSet

from tomo.protocols.protocol_ts_exclude_views_filter import ProtExclViewFilter
from tomo.protocols.protocol_ts_import import ProtImportTs


class TestExclViewsFilter(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        tsIds = (TS_03, TS_54)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.expectedSetSize = len(tsIds)
        cls.expectedSRate = DataSetRe4STATuto.unbinnedPixSize.value
        cls.expectedAcqDict, cls.expectedTsDims, cls.expectedAnglesCount = DataSetRe4STATuto.genTestTsDicts(tsIds)
        cls.expectedTsDims = DataSetRe4STATuto.dimsTsBin1Dict.value
        cls.protImportTs = cls._runImportTs()

    @classmethod
    def _runImportTs(cls,
                     filesPattern: str = DataSetRe4STATuto.tsPattern.value,
                     exclusionWords: str = DataSetRe4STATuto.exclusionWordsTs03ts54.value) -> SetOfTiltSeries:
        print(magentaStr("\n==> Importing the tilt-series:"))
        protImport = cls.newProtocol(ProtImportTs,
                                     filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                     filesPattern=filesPattern,
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
        cls.launchProtocol(protImport)
        tsImported = getattr(protImport, 'outputTiltSeries', None)
        print(tsImported)
        return tsImported
    
    @classmethod
    def _runExcludeViewsFilter(cls, inTS, **kargs) -> None:
        print(magentaStr("\n==> Excluding views in the tilt-series:"))
        protExclViewsFilter = cls.newProtocol(ProtExclViewFilter,
                                              inTsSet=inTS, **kargs)
        cls.launchProtocol(protExclViewsFilter)
        return protExclViewsFilter

    def test_exclude_by_tilt(self) -> None:
        print(magentaStr("\n==> Excluding views by tilting range:"))
        maxShiftX = 200
        maxShiftY = 200
        mintilt = 200
        maxtilt = 200
        maxDose = 200
        minViews = 200
        inTs = self.protImportTs
        print(inTs)
        protExclViewsFilter = self._runExcludeViewsFilter(inTs,
                                                          maxShiftX=maxShiftX,
                                                          maxShiftY=maxShiftY,
                                                          mintilt=mintilt,
                                                          maxtilt=maxtilt,
                                                          minDose=maxDose,
                                                          maxDose=minViews,
                                                          minViews=minViews)
        #self._checkTs(protExclViewsFilter)
        #self.assertIsNone(protExclViewsFilter)
    
    # def test_exclude_by_shift(self) -> None:
    #     print(magentaStr("\n==> Excluding views by tilting range:"))
    #     invTsSet, ctfTomoSet = self._runExcludeViewsFilter()
    #     self._checkTs(invTsSet)
    #     self.assertIsNone(ctfTomoSet)
    
    # def test_exclude_by_dose(self) -> None:
    #     print(magentaStr("\n==> Excluding views by tilting range:"))
    #     invTsSet, ctfTomoSet = self._runExcludeViewsFilter()
    #     self._checkTs(invTsSet)
    #     self.assertIsNone(ctfTomoSet)
    
    # def test_exclude_by_multipleChoice(self) -> None:
    #     print(magentaStr("\n==> Excluding views by tilting range:"))
    #     invTsSet, ctfTomoSet = self._runExcludeViewsFilter()
    #     self._checkTs(invTsSet)
    #     self.assertIsNone(ctfTomoSet)

        
