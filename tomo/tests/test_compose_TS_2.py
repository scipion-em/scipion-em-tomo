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


from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr, createLink, makePath, copyFile
from pwem.protocols import ProtImportMovies
from . import DataSet, RE4_STA_TUTO, DataSetRe4STATuto, TS_01, TS_43, TS_45, TS_03, TS_54
from .test_base_centralized_layer import TestBaseCentralizedLayer
from pyworkflow.plugin import Domain
from tomo.protocols.protocol_compose_TS import ProtComposeTS

import numpy as np

class TestTestTomoComposeTS2(TestBaseCentralizedLayer, BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.dataset = DataSet.getDataSet('tomo-em')
		cls.parentDir = cls.dataset.getFile('tsMParentFolder')
		cls.ts10Dir = 'tsM10Dir'
		cls.ts31Dir = 'tsM31Dir'
		cls.ts10_31 = ['tsM10Dir', 'tsM31Dir']
		cls.path = cls.dataset.getPath()
		cls.pattern = '*/stack*.mdoc'
		cls.sRate = 1.716

	@staticmethod
	def _getListOfFileNames(tomoNum, isTsMovie):
		angles = [0, 3, -3, -6, 6]
		if isTsMovie:
			# Each slice will be a stack of frames, which means a different file
			pattern = 'SKvesicles_Pertuzumab_tomo%2.0f_%03d_%1.1f.mrc'

			return [pattern % (tomoNum, i, angle)
			        for i, angle in enumerate(angles)]
		else:
			# Each slice will be referred to the same stack file, the one which represents the tilt series
			return ['stack%2.0f.mrcs' % tomoNum] * len(angles)

	@classmethod
	def _genTestDict(cls, **kwargs):
		return {
			cls.VOLTAGE: kwargs.get(cls.VOLTAGE, None),
			cls.MAGNIFICATION: kwargs.get(cls.MAGNIFICATION, None),
			cls.SPH_ABERRATION: kwargs.get(cls.SPH_ABERRATION, None),
			cls.AMP_CONTRAST: kwargs.get(cls.AMP_CONTRAST, None),
			cls.DOSE_PER_FRAME: kwargs.get(cls.DOSE_PER_FRAME, None),
			cls.ACCUM_DOSE: kwargs.get(cls.ACCUM_DOSE, None),
			cls.PIXEL_SIZE: kwargs.get(cls.PIXEL_SIZE, None),
			cls.SIZE: kwargs.get(cls.SIZE, None),
			cls.FILENAME_LIST: kwargs.get(cls.FILENAME_LIST, None),
			cls.INCOMING_DOSE_LIST: kwargs.get(cls.INCOMING_DOSE_LIST, None),
			cls.ACCUM_DOSE_LIST: kwargs.get(cls.ACCUM_DOSE_LIST, None),
			cls.ANGLE_LIST: kwargs.get(cls.ANGLE_LIST, None),
			cls.ACQ_ORDER_LIST: kwargs.get(cls.ACQ_ORDER_LIST, None),
			cls.ANGLE_MAX: kwargs.get(cls.ANGLE_MAX, None),
			cls.ANGLE_MIN: kwargs.get(cls.ANGLE_MIN, None),
			cls.ANGLE_STEP: kwargs.get(cls.ANGLE_STEP, None),
		}

	def sortTestDataByAngle(self, testDataDict):

		for key in testDataDict:
			testData = testDataDict[key]
			ind = np.argsort(testData[self.ANGLE_LIST])

			incDoses = []
			accumDoses = []
			newAcqOrders = []
			newAngles = []
			newFiles = []
			for index in ind:
				newAcqOrders.append(testData[self.ACQ_ORDER_LIST][index])
				incDoses.append(testData[self.INCOMING_DOSE_LIST][index])
				accumDoses.append(testData[self.ACCUM_DOSE_LIST][index])
				newAngles.append(testData[self.ANGLE_LIST][index])
				newFiles.append(testData[self.FILENAME_LIST][index])

			testData[self.ANGLE_LIST] = newAngles
			testData[self.INCOMING_DOSE_LIST] = incDoses
			testData[self.ACCUM_DOSE_LIST] = accumDoses
			testData[self.ACQ_ORDER_LIST] = newAcqOrders

	def _genTestData(self, isTsMovie):
		anglesStack31 = [0.0036, 2.9683, -3.0250, -6.0251, 5.9684]
		testData = {
		    'stack': self._genTestDict(
		        voltage=300,
		        magnification=53000,
		        sphAberration=2.7,
		        ampContrast=0.1,
		        dosePerFrame=2.2096,
		        accumDose=11.0479,
		        pixelSize=self.sRate,
		        size=5,
		        filenames=self._getListOfFileNames(tomoNum=31, isTsMovie=isTsMovie),
		        incDoses=[2.2088, 2.1975, 2.2245, 2.2275, 2.1895],
		        accumDoses=[2.2088, 4.4063, 6.6308, 8.8583, 11.0479],
		        angles=anglesStack31,
		        acqOrder=[1, 2, 3, 4, 5],
		        angleMax=max(anglesStack31),
		        angleMin=min(anglesStack31),
		        step=3
		    ),

		}

		if not isTsMovie:
			self.sortTestDataByAngle(testData)

		return testData


	@classmethod
	def _runImportMovies(cls, pathFolder, testDataTS):
		#testDataTS = testData['stack10'] or testData['stack31']
		print(magentaStr(f"\n==> Running the import movies preprocessing: \n" ))
		protMovieImport = cls.newProtocol(ProtImportMovies,
		                                   objLabel='Import movies (SPA)',
		                                   importFrom=ProtImportMovies.IMPORT_FROM_FILES,
		                                   filesPath=pathFolder,
		                                   filesPattern=cls.pattern,
		                                   voltage=300,
		                                   magnification=testDataTS['magnification'],
		                                   sphAberration=testDataTS['sphAberration'],
		                                   ampContrast=testDataTS['ampContrast'],
		                                   dosePerFrame=testDataTS['dosePerFrame'],
		                                   accumDose=testDataTS['accumDose'],
		                                   pixelSize=testDataTS['pixelSize'])
		cls.launchProtocol(protMovieImport)
		cls.assertIsNotNone(protMovieImport.outputMovies, 'Movies not imported')
		return protMovieImport.outputMovies


	def _runAlignMovies(cls, movies):
		print(magentaStr(f"\n==> Running the import movies preprocessing: \n" ))

		xmipp3 = Domain.importFromPlugin('xmipp3.protocols', doRaise=True)
		protAlign = cls.newProtocol(xmipp3.XmippProtFlexAlign,
		                            objLabel='Movie Alignment (SPA)',
		                            alignFrame0=1,
		                            alignFrameN=0,
		                            useAlignToSum=True,
		                            doLocalAlignment=False)
		protAlign.inputMovies.set(movies)
		cls.launchProtocol(protAlign)
		cls.assertIsNotNone(protAlign.outputMicrographs, 'Micrograph not generated')
		return getattr(protAlign, 'outputMicrographs', None)

	def _runComposeTS(self, outputMicrographs, filesPath, mdocPattern):
		print(magentaStr(f"\n==> Running the composeTS: \n"))

		protComposeTS = self.newProtocol(ProtComposeTS,
	                                   objLabel='Import movies (SPA)',
	                                   importFrom=outputMicrographs,
	                                   filesPath=self.partFolderPath,
	                                   mdocPattern=self.pattern,
		                               isTomo5=False,
		                               mdoc_bug_Correction=False,
		                               percentTiltsRequired='80',
		                               time4NextTilt='20s')
		protComposeTS.Micrographs.set(outputMicrographs)
		return protComposeTS.TiltSeries

	def test_composeTSBasic(self):
		pathFolder= ''
		testDataTS = ''
		outputMovies = self._runImportMovies(pathFolder, testDataTS)

		outputMicrographs = self._runAlignMovies(outputMovies)

		filesPath = ''
		mdocPattern = ''
		TiltSeries = self._runComposeTS(outputMicrographs, filesPath, mdocPattern)
