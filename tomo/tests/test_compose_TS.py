# **************************************************************************
# *
# * Authors:    alberto garcia (alberto.garcia@cnb.csic.es)(scipion@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.tests import setupTestProject
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportMovies
from . import DataSet, RE_STA_TUTO_MOVIES, DataSetRe4STATuto, TS_03, TS_54
from .test_base_centralized_layer import TestBaseCentralizedLayer
from pyworkflow.plugin import Domain
from tomo.protocols.protocol_compose_TS import ProtComposeTS
from motioncorr.protocols import ProtMotionCorr

import numpy as np

class TestTestTomoComposeTS(TestBaseCentralizedLayer):
	""" This class check if the protocol to compose TiltSeries works properly."""
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.ds = DataSet.getDataSet(RE_STA_TUTO_MOVIES)


	@classmethod
	def _runImportMovies(cls, blackList=None):
		#testDataTS = testData['stack10'] or testData['stack31']
		print(magentaStr(f"\n==> Running the import movies preprocessing: \n" ))
		protMovieImport = cls.newProtocol(ProtImportMovies,
		                                   objLabel='Import movies (SPA)',
		                                   importFrom=ProtImportMovies.IMPORT_FROM_FILES,
		                                   filesPath=cls.ds.getFile(RE_STA_TUTO_MOVIES.framesDir.name),
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
		cls.assertIsNotNone(protMovieImport.outputMovies, 'Movies not imported')
		return protMovieImport.outputMovies

	def _runAlignMovies(cls, movies):
		print(magentaStr(f"\n==> Running the align movies preprocessing: \n" ))
		protMc = cls.newProtocol(ProtMotionCorr,
									objLabel='Movie Alignment (SPA)')
		protMc.inputMovies.set(movies)
		cls.launchProtocol(protMc)
		cls.assertIsNotNone(protMc.outputMovies, 'Micrograph not generated')
		return getattr(protMc, 'outputMicrographs', None)

	def _runComposeTS(cls, outputMicrographs, filesPath, mdocPattern, isTomo5=False, mdoc_bug_Correction=False, percentTiltsRequired='80', time4NextTilt='20'):
		print(magentaStr(f"\n==> Running the composeTS: \n"))

		protComposeTS = cls.newProtocol(ProtComposeTS,
	                                   objLabel='Compose TiltSeries',
	                                   inputMicrographs=outputMicrographs,
	                                   filesPath=filesPath,
	                                   mdocPattern=mdocPattern,
		                               isTomo5=isTomo5,
		                               mdoc_bug_Correction=mdoc_bug_Correction,
		                               percentTiltsRequired=percentTiltsRequired,
		                               time4NextTilt=time4NextTilt)

		cls.launchProtocol(protComposeTS)
		return getattr(protComposeTS, 'TiltSeries', None)


	def test_composeTSBasic(self):
		print(magentaStr(f"\n==> Running the basic Test: \n"))
		outputMovies = self._runImportMovies()
		outputMicrographs = self._runAlignMovies(outputMovies)
		mdocPattern = '*mrc.mdoc'
		filesPath = self.ds.getFile(RE_STA_TUTO_MOVIES.framesDir.name)
		TiltSeries = self._runComposeTS(outputMicrographs, filesPath, mdocPattern, percentTiltsRequired='100')

		#TEST VALUES
		expectedSetSize = 2
		anglesCount = {TS_03: 5, TS_54: 6}

		print(f'testSetAcqObj=RE_STA_TUTO_MOVIES.tsAcqDict.value: {RE_STA_TUTO_MOVIES.tsAcqDict.value}')
		self.checkTiltSeries(TiltSeries,
		                     expectedSetSize=expectedSetSize,
		                     expectedSRate=RE_STA_TUTO_MOVIES.unbinnedPixSize.value,
		                     hasAlignment=False,
		                     isHeterogeneousSet=False,
		                     imported=True,
		                     expectedDimensions=RE_STA_TUTO_MOVIES.dimsTsBin1Dict.value,
		                     testAcqObj=RE_STA_TUTO_MOVIES.tsAcqDict.value,
		                     anglesCount=anglesCount)


		print(magentaStr(f"\n==> Running the rejected mics Test: \n"))
		mdocPattern = '*rejecting.mdoc'
		TiltSeries = self._runComposeTS(outputMicrographs, filesPath, mdocPattern, percentTiltsRequired='80')
		expectedSetSize = 1
		anglesCount = {TS_54: 5}

		self.assertSetSize(TiltSeries, expectedSetSize)
		self.checkTiltSeries(TiltSeries,
		                     expectedSetSize=expectedSetSize,
		                     expectedSRate=RE_STA_TUTO_MOVIES.unbinnedPixSize.value,
		                     hasAlignment=False,
		                     isHeterogeneousSet=False,
		                     imported=True,
		                     expectedDimensions=RE_STA_TUTO_MOVIES.dimsTs54Bin1Dict.value,
		                     testAcqObj=RE_STA_TUTO_MOVIES.tsAcq54Dict.value,
		                     anglesCount=anglesCount)



