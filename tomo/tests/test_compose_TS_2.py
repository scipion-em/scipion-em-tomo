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
from . import DataSet, RE_STA_TUTO_MOVIES, DataSetRE_STA_TUTO_MOVIES, DataSetRe4STATuto, tsAcqDict, testAcq03, testAcq54, TS_03, TS_54
from .test_base_centralized_layer import TestBaseCentralizedLayer
from pyworkflow.plugin import Domain
from tomo.protocols.protocol_compose_TS import ProtComposeTS


import numpy as np

class TestTestTomoComposeTS2(TestBaseCentralizedLayer):
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
		                                   filesPath=cls.ds.getFile(DataSetRE_STA_TUTO_MOVIES.framesDir.name),
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

		xmipp3 = Domain.importFromPlugin('xmipp3.protocols', doRaise=True)
		protAlign = cls.newProtocol(xmipp3.XmippProtFlexAlign,
		                            objLabel='Movie Alignment (SPA)',
		                            alignFrame0=1,
		                            alignFrameN=0,
		                            useAlignToSum=True,
		                            doLocalAlignment=False)
		protAlign.inputMovies.set(movies)
		cls.launchProtocol(protAlign)
		cls.assertIsNotNone(protAlign.outputMovies, 'Micrograph not generated')
		return getattr(protAlign, 'outputMicrographs', None)

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
		#DYNAMIC TEMPLATE STARTS
		import os
		fname = "/home/agarcia/Documents/test_DEBUGALBERTO.txt"
		if os.path.exists(fname):
		    os.remove(fname)
		fjj = open(fname, "a+")
		fjj.write('ALBERTO--------->onDebugMode PID {}'.format(os.getpid()))
		fjj.close()
		print('ALBERTO--------->onDebugMode PID {}'.format(os.getpid()))
		import time
		time.sleep(10)

		#DYNAMIC TEMPLATE ENDS
		outputMovies = self._runImportMovies()
		outputMicrographs = self._runAlignMovies(outputMovies)
		mdocPattern = '*.mdoc'
		filesPath = self.ds.getFile(DataSetRE_STA_TUTO_MOVIES.framesDir.name)



		TiltSeries = self._runComposeTS(outputMicrographs, filesPath, mdocPattern)

		#TEST VALUES
		testAcq03.setAngleMin(-6.0)
		testAcq54.setAngleMin(-6.0)
		expectedSetSize = 2
		dimsTsBin1Dict = {TS_03: [3710, 3838, 5], TS_54: [3710, 3838, 6]}
		sampligRate = DataSetRE_STA_TUTO_MOVIES.unbinnedPixSize.value
		anglesCount = {TS_03: 5, TS_54: 6}

		print(f'testSetAcqObj=DataSetRE_STA_TUTO_MOVIES.tsAcqDict.value: {DataSetRE_STA_TUTO_MOVIES.tsAcqDict.value}')
		self.checkTiltSeries(TiltSeries,
		                     expectedSetSize=expectedSetSize,
		                     expectedSRate=sampligRate,
		                     imported=True,
		                     expectedDimensions=dimsTsBin1Dict,
		                     testAcqObj=DataSetRE_STA_TUTO_MOVIES.tsAcqDict.value,
		                     anglesCount=anglesCount)


	def test_composeTS_TiltsRejected(self):
		print(magentaStr(f"\n==> Running the rejected mics Test: \n"))

		reg_TS03_03= 'TS_03_003_-6.0.mrc'
		reg_TS03_04 = 'TS_03_004_6.0.mrc'
		reg_TS54_05 = 'TS_54_005_9.0.mrc'
		reg_TS54_04 = 'TS_54_004_6.0.mrc'

		with open(self.proj.getTmpPath('blacklist_regex.txt'), 'w') as f:
			f.write(reg_TS03_03)
			f.write(reg_TS03_04)
			f.write(reg_TS54_05)
			f.write(reg_TS54_04)

		#outputMovies = self._runImportMovies(blackList='blacklist_regex.txt')


