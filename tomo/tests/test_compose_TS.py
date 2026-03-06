# **************************************************************************
# *
# * Authors:     Alberto García Mena (alberto.garcia@cnb.csic.es)
# *              Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
# *
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
from typing import Optional
from pwem.objects import SetOfMovies, SetOfMicrographs
from pyworkflow.tests import setupTestProject
from pyworkflow.utils import magentaStr, cyanStr
from pwem.protocols import ProtImportMovies
from . import DataSet, RE_STA_TUTO_MOVIES, DataSetRe4STATuto, DataSet_RE_STA_TUTO_MOVIES, TS_03, TS_54
from .test_base_centralized_layer import TestBaseCentralizedLayer
from tomo.protocols.protocol_compose_TS import ProtComposeTS, OUT_TS_SET
from motioncorr.protocols import ProtMotionCorrNewStreaming
from ..objects import SetOfTiltSeries


class TestTestTomoComposeTS(TestBaseCentralizedLayer):
    binFactor = 2

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE_STA_TUTO_MOVIES)
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        cls._runPreviousProtocols()
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

    @classmethod
    def _runPreviousProtocols(cls):
        importedMovies = cls._runImportMovies()
        cls.mcMovies = cls._runAlignMovies(importedMovies)

    @classmethod
    def _runImportMovies(cls, blackList=None) -> Optional[SetOfMovies]:
        print(magentaStr(f"\n==> Importing the movies: \n"))
        protMovieImport = cls.newProtocol(ProtImportMovies,
                                          importFrom=ProtImportMovies.IMPORT_FROM_FILES,
                                          filesPath=cls.ds.getFile(DataSet_RE_STA_TUTO_MOVIES.framesDir.name),
                                          filesPattern='*.mrc',
                                          blacklistSet=blackList,
                                          voltage=DataSetRe4STATuto.voltage.value,
                                          magnification=DataSetRe4STATuto.magnification.value,
                                          sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                          amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                          samplingRate=DataSetRe4STATuto.unbinnedPixSize.value,
                                          doseInitial=DataSetRe4STATuto.initialDose.value,
                                          dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value)

        cls.launchProtocol(protMovieImport)
        return getattr(protMovieImport, 'outputMovies', None)

    @classmethod
    def _runAlignMovies(cls, movies: SetOfMovies) -> Optional[SetOfMicrographs]:
        print(magentaStr(f"\n==> Running the motion correction with Motioncorr: \n"))
        protMc = cls.newProtocol(ProtMotionCorrNewStreaming,
                                 inputMovies=movies,
                                 binFactor=2,
                                 doApplyDoseFilter=True,
                                 splitEvenOdd=True)
        cls.launchProtocol(protMc)
        return getattr(protMc, protMc._possibleOutputs.micrographsDW.name, None)

    @classmethod
    def _runComposeTS(cls,
                      doEvenOdd: bool = False,
                      percentTiltsRequired: int = 80,
                      time4NextTilt: int = 20) -> Optional[SetOfTiltSeries]:
        print(magentaStr(f"\n==> Running the composeTS:"))
        print(magentaStr(f"\n\t- Odd/Even: {doEvenOdd}\n"))
        protComposeTS = cls.newProtocol(ProtComposeTS,
                                        objLabel=f'Compose ts, oe = {doEvenOdd}',
                                        inputMicrographs=cls.mcMovies,
                                        filesPath=cls.ds.getFile(DataSet_RE_STA_TUTO_MOVIES.framesDir.name),
                                        mdocPattern='*mrc.mdoc',
                                        doEvenOdd=doEvenOdd,
                                        percentTiltsRequired=percentTiltsRequired,
                                        time4NextTilt=time4NextTilt)

        cls.launchProtocol(protComposeTS)
        return getattr(protComposeTS, OUT_TS_SET, None)

    def testComposeTs01(self):
        tsSet = self._runComposeTS()
        self._checkTs(tsSet)

    def testComposeTs02(self):
        tsSet = self._runComposeTS(doEvenOdd=True)
        self._checkTs(tsSet, hasOddEven=True)

    def _checkTs(self,
                 tsSet: SetOfTiltSeries,
                 hasOddEven: bool = False):
        self.checkTiltSeries(tsSet,
                             expectedSetSize=2,
                             expectedSRate=DataSet_RE_STA_TUTO_MOVIES.unbinnedPixSize.value * self.binFactor,
                             hasAlignment=False,
                             isHeterogeneousSet=True,
                             hasOddEven=hasOddEven,
                             imported=True,
                             expectedDimensions=DataSet_RE_STA_TUTO_MOVIES.dimsTsBin2Dict.value,
                             testAcqObj=DataSet_RE_STA_TUTO_MOVIES.tsAcqDict.value,
                             anglesCount={TS_03: 5, TS_54: 6})

    # def test_composeTSBasic(self):
    #     print(magentaStr(f"\n==> Running the basic Test: \n"))
    #     outputMovies = self._runImportMovies()
    #     # outputMicrographs = self._runAlignMoviesFlexAlign(outputMovies)
    #     outputMicrographs = self._runAlignMovies(outputMovies)
    #
    #     mdocPattern = '*mrc.mdoc'
    #     filesPath = self.ds.getFile(DataSet_RE_STA_TUTO_MOVIES.framesDir.name)
    #     TiltSeries = self._runComposeTS(outputMicrographs, filesPath, mdocPattern, percentTiltsRequired='100')
    #
    #     # TEST VALUES
    #     expectedSetSize = 2
    #     anglesCount = {TS_03: 5, TS_54: 6}
    #
    #     print(magentaStr(f"\n==> Checking Tilt Series: \n"))
    #     self.checkTiltSeries(TiltSeries,
    #                          expectedSetSize=expectedSetSize,
    #                          expectedSRate=DataSet_RE_STA_TUTO_MOVIES.unbinnedPixSize.value,
    #                          hasAlignment=False,
    #                          isHeterogeneousSet=False,
    #                          imported=True,
    #                          expectedDimensions=DataSet_RE_STA_TUTO_MOVIES.dimsTsBin1Dict.value,
    #                          testAcqObj=DataSet_RE_STA_TUTO_MOVIES.tsAcqDict.value,
    #                          anglesCount=anglesCount)
    #
    #     print(magentaStr(f"\n==> Running the rejected mics Test: \n"))
    #     mdocPattern = '*rejecting.mdoc'
    #     TiltSeries = self._runComposeTS(outputMicrographs, filesPath, mdocPattern, percentTiltsRequired='80')
    #     expectedSetSize = 1
    #     anglesCount = {TS_54: 5}
    #     self.assertSetSize(TiltSeries, expectedSetSize)
    #     print(magentaStr(f"\n==> Checking Tilt Series: \n"))
    #     self.checkTiltSeries(TiltSeries,
    #                          expectedSetSize=expectedSetSize,
    #                          expectedSRate=DataSet_RE_STA_TUTO_MOVIES.unbinnedPixSize.value,
    #                          hasAlignment=False,
    #                          isHeterogeneousSet=False,
    #                          imported=True,
    #                          expectedDimensions=DataSet_RE_STA_TUTO_MOVIES.dimsTs54Bin1Dict.value,
    #                          testAcqObj=DataSet_RE_STA_TUTO_MOVIES.testAcq54_rejectDict.value,
    #                          anglesCount=anglesCount)
