# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
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
import glob
import os
import random
import string
import tempfile
from os.path import join, exists, abspath, basename
from typing import List, Set, Tuple

import numpy as np

from tomo.convert import getOrderFromList
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr, createLink, makePath, copyFile
from pyworkflow.object import Pointer
from pwem.protocols import ProtSplitSet, ProtSetFilter, ProtSetEditor

from tomo.protocols.protocol_ts_import import MDoc, ProtImportTs
from . import DataSet, RE4_STA_TUTO, DataSetRe4STATuto, TS_01, TS_43, TS_45, TS_03, TS_54
from .test_base_centralized_layer import TestBaseCentralizedLayer
from ..constants import BOTTOM_LEFT_CORNER, TOP_LEFT_CORNER, ERR_COORDS_FROM_SQLITE_NO_MATCH, ERR_NO_TOMOMASKS_GEN, \
    ERR_NON_MATCHING_TOMOS, SCIPION
import tomo.protocols
from ..objects import TomoAcquisition
from ..protocols import ProtImportTomograms, ProtImportTomomasks
from ..protocols.protocol_base import ProtTomoImportAcquisition
from ..protocols.protocol_import_coordinates import IMPORT_FROM_AUTO, ProtImportCoordinates3D
from ..protocols.protocol_import_coordinates_from_scipion import ProtImportCoordinates3DFromScipion, outputObjs
from ..protocols.protocol_import_ctf import ImportChoice, ProtImportTsCTF
from ..protocols.protocol_import_tomograms import OUTPUT_NAME
from ..utils import existsPlugin

from imod.protocols import ProtImodTomoNormalization


class TestTomoImportTomogramsFromPattern(TestBaseCentralizedLayer):

    ds = None
    bin4SRate = DataSetRe4STATuto.unbinnedPixSize.value * 4

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)

    @classmethod
    def _runImportTomograms(cls, filesPattern=None, exclusionWords=None, objLabel=None):
        print(magentaStr(f"\n==> Importing the tomograms:"
                         f"\n\t- Files pattern = {filesPattern}"
                         f"\n\t- Excluded words = {exclusionWords}"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=filesPattern,
                                          exclusionWords=exclusionWords,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value,  # Bin 4
                                          importAcquisitionFrom=ProtTomoImportAcquisition.MANUAL_IMPORT,
                                          oltage=DataSetRe4STATuto.voltage.value,
                                          sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                          amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value)
        if objLabel:
            protImportTomos.setObjLabel(objLabel)
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    def _checkTomos(self, inTomos, testAcqObjDict=None, expectedDimensionsDict=None, isHeterogeneousSet=True):
        self.checkTomograms(inTomos,
                            expectedSetSize=len(expectedDimensionsDict),
                            expectedSRate=self.bin4SRate,
                            testAcqObj=testAcqObjDict,
                            expectedDimensions=expectedDimensionsDict,
                            isHeterogeneousSet=isHeterogeneousSet)

    def testImportTomos01(self):
        filesPattern = DataSetRe4STATuto.tomosPattern.value
        importedTomos = self._runImportTomograms(filesPattern=filesPattern, objLabel='testImportTomos01')
        # Check the results
        testAcqObjDict, expectedDimensionsDict = DataSetRe4STATuto.genTestTomoDicts()
        self._checkTomos(importedTomos, testAcqObjDict=testAcqObjDict, expectedDimensionsDict=expectedDimensionsDict)

    def testImportTomos02(self):
        filesPattern = '*{TS}.mrc'
        importedTomos = self._runImportTomograms(filesPattern=filesPattern, objLabel='testImportTomos02')
        # Check the results
        testAcqObjDict, expectedDimensionsDict = DataSetRe4STATuto.genTestTomoDicts()
        self._checkTomos(importedTomos, testAcqObjDict=testAcqObjDict, expectedDimensionsDict=expectedDimensionsDict)

    def testImportTomos03(self):
        filesPattern = '*TS_{TS}.mrc'  # With this pattern, the tsId become numeric, but the protocol should fix them
        exclusionWords = 'TS_03 TS_54'
        importedTomos = self._runImportTomograms(filesPattern=filesPattern,
                                                 exclusionWords=exclusionWords,
                                                 objLabel='testImportTomos03')
        # Check the results
        tsIdList = (TS_01, TS_43, TS_45)
        testAcqObjDict, expectedDimensionsDict = DataSetRe4STATuto.genTestTomoDicts(tsIdList=tsIdList)
        self._checkTomos(importedTomos, testAcqObjDict=testAcqObjDict, expectedDimensionsDict=expectedDimensionsDict)

    def testImportTomos04(self):
        filesPattern = '*{TS}.mrc'
        exclusionWords = DataSetRe4STATuto.exclusionWordsTs03ts54.value
        importedTomos = self._runImportTomograms(filesPattern=filesPattern,
                                                 exclusionWords=exclusionWords,
                                                 objLabel='testImportTomos04')
        # Check the results
        tsIdList = (TS_03,  TS_54)
        testAcqObjDict, expectedDimensionsDict = DataSetRe4STATuto.genTestTomoDicts(tsIdList=tsIdList)
        self._checkTomos(importedTomos, testAcqObjDict=testAcqObjDict,
                         expectedDimensionsDict=expectedDimensionsDict,
                         isHeterogeneousSet=False)  # Now the set is homogeneous


class TestTomoImportTsFromPattern(TestBaseCentralizedLayer):
    testAcq = None
    samplingRate = 1.35

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.getFile = cls.dataset.getFile('etomo')
        cls.getFileM = cls.dataset.getFile('empiar')
        # Test acquisition: set the common attributes
        testAcq = TomoAcquisition()
        testAcq.setVoltage(300)
        testAcq.setMagnification(105000)
        testAcq.setSphericalAberration(2.7)
        testAcq.setAmplitudeContrast(0.1)
        testAcq.setDoseInitial(0)
        testAcq.setDosePerFrame(0.3)
        testAcq.setTiltAxisAngle(84.1)
        # EmpiarTestData specific parameters
        testAcqEmpiarTestData = testAcq.clone()
        testAcqEmpiarTestData.setAngleMin(-6)
        testAcqEmpiarTestData.setAngleMax(6)
        testAcqEmpiarTestData.setStep(3)
        cls.testAcqEmpiarTestData = testAcqEmpiarTestData
        # EmpiarTestData2 specific parameters
        testAcqEmpiarTestData2 = testAcq.clone()
        testAcqEmpiarTestData2.setAngleMin(-60)
        testAcqEmpiarTestData2.setAngleMax(60)
        testAcqEmpiarTestData2.setStep(60)
        cls.testAcqEmpiarTestData2 = testAcqEmpiarTestData2
        # ImodTestDat specific parameters
        testAcqImodTestData = testAcq.clone()
        testAcqImodTestData.setAngleMin(-55)
        testAcqImodTestData.setAngleMax(65)
        testAcqImodTestData.setStep(2)
        cls.testAcqImodTestData = testAcqImodTestData

    def _runImportTiltSeriesM(self, acq, filesPath, filesPattern='{TS}_{TO}_{TA}.mrc'):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTsMovies,
            filesPath=filesPath,
            filesPattern=filesPattern,
            voltage=acq.getVoltage(),
            magnification=acq.getMagnification(),
            sphericalAberration=acq.getSphericalAberration(),
            amplitudeContrast=acq.getAmplitudeContrast(),
            samplingRate=self.samplingRate,
            doseInitial=acq.getDoseInitial(),
            dosePerFrame=acq.getDosePerFrame(),
            tiltAxisAngle=acq.getTiltAxisAngle()
        )
        self.launchProtocol(protImport)
        return protImport

    def _runImportTiltSeries(self, acq):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTs,
            filesPath=self.getFile,
            filesPattern='BB{TS}.st',
            minAngle=acq.getAngleMin(),
            maxAngle=acq.getAngleMax(),
            stepAngle=acq.getStep(),
            voltage=acq.getVoltage(),
            magnification=acq.getMagnification(),
            sphericalAberration=acq.getSphericalAberration(),
            amplitudeContrast=acq.getAmplitudeContrast(),
            samplingRate=self.samplingRate,
            doseInitial=acq.getDoseInitial(),
            dosePerFrame=acq.getDosePerFrame(),
            tiltAxisAngle=acq.getTiltAxisAngle()
        )
        self.launchProtocol(protImport)
        return protImport

    def test_importTiltSeriesM(self):
        # Expected values
        testAcq = self.testAcqEmpiarTestData2
        expectedSetSize = 2
        expectedDimensions = [7420, 7676, 8]
        expectedAnglesCount = 3
        presentAcqOrders = [0, 39, 40]
        testAcq.setAccumDose(testAcq.getDosePerFrame() * max(presentAcqOrders))

        # Run teh protocol
        protImport = self._runImportTiltSeriesM(testAcq, filesPath=self.getFileM)

        # Check the results
        # self.checkTSSet(protImport.outputTiltSeriesM, 2, 3)
        self.checkTiltSeriesM(protImport.outputTiltSeriesM,
                              expectedSetSize=expectedSetSize,
                              expectedSRate=self.samplingRate,
                              expectedDimensions=expectedDimensions,
                              testAcqObj=testAcq,
                              anglesCount=expectedAnglesCount,
                              checkIds=True)
        return protImport

    def test_importTiltSeriesM_withBrackets(self):
        # Make a tmp dir with the angular stacks renamed, so they contain bracket characters in their base names
        tmpBaseName = ''.join(random.choices(string.ascii_letters, k=5))
        tmpDir = join(tempfile.gettempdir(), tmpBaseName)
        makePath(tmpDir)
        tsMDirKeys = ['tsM10Dir', 'tsM31Dir']
        str2replace = 'Pertuzumab_'
        newStr = str2replace + '015[16]_'
        for tsMDirKey in tsMDirKeys:
            dirName = self.dataset.getFile(tsMDirKey)
            newPath = join(tmpDir, basename(dirName))
            makePath(dirName, newPath)
            origFiles = glob.glob(join(dirName, '*.mrc'))
            [copyFile(origFile, join(newPath, basename(origFile).replace(str2replace, newStr))) for origFile in
             origFiles]

        # Expected values
        testAcq = self.testAcqEmpiarTestData
        expectedSetSize = 2
        expectedDimensions = [1152, 1152, 6]
        expectedAnglesCount = 5
        presentAcqOrders = [0, 1, 2, 3, 4]
        testAcq.setAccumDose(testAcq.getDosePerFrame() * max(presentAcqOrders))

        # 1: Brackets in out of the {} labels
        filesPattern = '*/SKvesicles_Pertuzumab_015[16]_{TS}_{TO}_{TA}.mrc'
        protImport = self._runImportTiltSeriesM(testAcq, filesPath=tmpDir, filesPattern=filesPattern)
        self.checkTiltSeriesM(protImport.outputTiltSeriesM,
                              expectedSetSize=expectedSetSize,
                              expectedSRate=self.samplingRate,
                              expectedDimensions=expectedDimensions,
                              testAcqObj=testAcq,
                              anglesCount=expectedAnglesCount,
                              checkIds=True)

        # 2: Brackets in the {TS} label
        filesPattern = '*/SKvesicles_Pertuzumab_{TS}_{TO}_{TA}.mrc'
        protImport = self._runImportTiltSeriesM(testAcq, filesPath=tmpDir, filesPattern=filesPattern)
        self.checkTiltSeriesM(protImport.outputTiltSeriesM,
                              expectedSetSize=expectedSetSize,
                              expectedSRate=self.samplingRate,
                              expectedDimensions=expectedDimensions,
                              testAcqObj=testAcq,
                              anglesCount=expectedAnglesCount,
                              checkIds=True)

    def test_importTiltSeries(self):
        testAcq = self.testAcqImodTestData
        expectedSetSize = 2
        expectedDimensions = [512, 512, 61]
        expectedAnglesCount = 61
        testAcq.setAccumDose(testAcq.getDosePerFrame() * expectedAnglesCount)

        # Run the protocol
        protImport = self._runImportTiltSeries(testAcq)

        # Check teh results
        self.checkTiltSeries(protImport.outputTiltSeries,
                             expectedSetSize=expectedSetSize,
                             expectedSRate=self.samplingRate,
                             imported=True,
                             expectedDimensions=expectedDimensions,
                             testSetAcqObj=testAcq,
                             testAcqObj=testAcq,
                             anglesCount=expectedAnglesCount)

    def checkTSSet(self, tsSet, expectedSetSize=-1, expectedSRate=-1, expectedDimensions=None,
                   testAcqObj=None, anglesCount=None, checkIds=False):

        self.checkTiltSeries(tsSet,
                             expectedSetSize=expectedSetSize,
                             expectedSRate=expectedSRate,
                             imported=True,
                             expectedDimensions=expectedDimensions,
                             testSetAcqObj=testAcqObj,
                             testAcqObj=testAcqObj,
                             anglesCount=anglesCount)
        # Other checks
        for ts in tsSet:
            for i, ti in enumerate(ts):
                self.assertFalse(os.path.isabs(ti.getFileName()),
                                 "Tilt image file %s is not absolute!. Should be relative.")
                self.assertTrue(os.path.islink(ti.getFileName()),
                                "Tilt series file %s is not a link." % ti.getFileName())

                # objId must start with 1 --> i +1
                # When using {TO} the objId ends up with that value. For Movies it is 0, 39, 40
                if checkIds:
                    self.assertEqual(i + 1, ti.getObjId(), "Tilt image Movie objId is incorrect")


class TestTomoImportTsFromMdoc(BaseTest):
    VOLTAGE = 'voltage'
    MAGNIFICATION = 'magnification'
    SPH_ABERRATION = 'sphAberration'
    AMP_CONTRAST = 'ampContrast'
    DOSE_PER_FRAME = 'dosePerFrame'
    ACCUM_DOSE = 'accumDose'
    PIXEL_SIZE = 'pixelSize'
    SIZE = 'size'
    FILENAME_LIST = 'filenames'
    INCOMING_DOSE_LIST = 'incDoses'
    ACCUM_DOSE_LIST = 'accumDoses'
    ANGLE_LIST = 'angles'
    ACQ_ORDER_LIST = 'acqOrder'
    ANGLE_MAX = 'angleMax'
    ANGLE_MIN = 'angleMin'
    ANGLE_STEP = 'step'

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.parentDir = cls.dataset.getFile('tsMParentFolder')
        cls.ts10Dir = cls.dataset.getFile('tsM10Dir')
        cls.ts31Dir = cls.dataset.getFile('tsM31Dir')
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
        anglesStack10 = [0.0026, 2.9708, -3.0255, -6.0241, 5.9684]
        testData = {
            'stack31': self._genTestDict(
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
            'stack10': self._genTestDict(
                voltage=300,
                magnification=53000,
                sphAberration=2.7,
                ampContrast=0.1,
                dosePerFrame=2.3406,
                accumDose=11.7028,
                pixelSize=self.sRate,
                size=5,
                filenames=self._getListOfFileNames(tomoNum=10, isTsMovie=isTsMovie),
                incDoses=[2.3443, 2.3421, 2.3425, 2.3373, 2.3366],
                accumDoses=[2.3443, 4.6864, 7.0289, 9.3662, 11.7028],
                angles=anglesStack10,
                acqOrder=[1, 2, 3, 4, 5],
                angleMax=max(anglesStack10),
                angleMin=min(anglesStack10),
                step=3
            )
        }

        if not isTsMovie:
            self.sortTestDataByAngle(testData)

        return testData

    def _runImportTiltSeries(self, filesPath, pattern, isTsMovie=False, exclusionWords=None,
                             label="import tilt series"):
        prot = tomo.protocols.ProtImportTsMovies if isTsMovie else tomo.protocols.ProtImportTs
        attribDict = {
            'filesPath': filesPath,
            'filesPattern': pattern,
            'objLabel': label
        }
        if exclusionWords:
            attribDict['exclusionWords'] = exclusionWords
        protImport = self.newProtocol(
            prot,
            **attribDict)
        self.launchProtocol(protImport)
        return protImport

    def test_importTiltSeries_numeric_basename(self):
        # Make a tmp dir with the TS, renamed to have a numeric base name, and also edit the corresponding mdoc files
        # for consistency
        tmpBaseName = ''.join(random.choices(string.ascii_letters, k=5))
        tmpDir = join(tempfile.gettempdir(), tmpBaseName)
        makePath(tmpDir)
        tsMDirKeys = ['tsM10Dir', 'tsM31Dir']
        fileExt2Copy = ['*.mrcs', '*.mdoc']
        str2replace = 'stack'
        origFiles = []
        labelImageFileLineNum = 2
        for tsMDirKey in tsMDirKeys:
            dirName = self.dataset.getFile(tsMDirKey)
            newPath = join(tmpDir, basename(dirName))
            makePath(dirName, newPath)
            # Copy the mrcs and mdoc files to a temp location, removing the 'stack' word from the name, so they become
            # numerical names: stack10.mdoc --> 10.mdoc
            for ext in fileExt2Copy:
                origFiles += glob.glob(join(dirName, ext))
            [copyFile(origFile, join(newPath, basename(origFile).replace(str2replace, ''))) for origFile in origFiles]
            copiedMdoc = glob.glob(join(newPath, '*.mdoc'))[0]
            # Update the Image File label in the mdoc to point to the numerical renamed TS
            with open(copiedMdoc, "r") as f:
                lines = f.readlines()
            lines[labelImageFileLineNum] = lines[labelImageFileLineNum].replace(str2replace, '')
            with open(copiedMdoc, "w") as f:
                f.writelines(lines)

        isTsMovie = False
        label = 'Import TS numeric base name'
        # Run protocol
        protImport = self._runImportTiltSeries(tmpDir, '*/*.mdoc', isTsMovie=isTsMovie)
        protImport.setObjLabel(label)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeries', None)
        self._checkResults(outputSet, isTsMovie, (1440, 1024, 1), protImport, size=2, numericBaseName=True)
        return protImport

    def _checkResults(self, outputSet, isTsMovie, dimensions, prot, size=None, numericBaseName=False):
        # Generate the test data
        testDict = self._genTestData(isTsMovie)
        # Check set properties
        self.assertSetSize(outputSet, size=size)
        self.assertAlmostEqual(outputSet.getSamplingRate(), self.sRate, delta=0.001)
        self.assertEqual(outputSet.getDimensions(), dimensions)

        # Check tilt series movies
        for tsM in outputSet:
            # If numeric base name, the expected results need to be adapted as new files are generated in execution
            # time to simulate the scenario for that test without having to add more files to the test dataset
            tsId = tsM.getTsId().replace('TS_', 'stack') if numericBaseName else tsM.getTsId()
            testDataDict = testDict[tsId]
            # Check set acquisition
            acq = tsM.getAcquisition()
            self.assertAlmostEqual(tsM.getSamplingRate(), testDataDict[self.PIXEL_SIZE], delta=0.001)
            self.assertAlmostEqual(acq.getVoltage(), testDataDict[self.VOLTAGE], delta=0.1)
            self.assertAlmostEqual(acq.getSphericalAberration(), testDataDict[self.SPH_ABERRATION], delta=0.01)
            self.assertAlmostEqual(acq.getAmplitudeContrast(), testDataDict[self.AMP_CONTRAST], delta=0.001)
            self.assertAlmostEqual(acq.getDosePerFrame(), testDataDict[self.DOSE_PER_FRAME], delta=0.0001)
            self.assertAlmostEqual(acq.getAccumDose(), testDataDict[self.ACCUM_DOSE], delta=0.0001)
            self.assertAlmostEqual(acq.getAngleMax(), testDataDict[self.ANGLE_MAX], delta=0.001)
            self.assertAlmostEqual(acq.getAngleMin(), testDataDict[self.ANGLE_MIN], delta=0.001)
            self.assertAlmostEqual(acq.getStep(), testDataDict[self.ANGLE_STEP], delta=0.1)
            # Check angles and accumulated dose per angular acquisition (tilt series image)
            filesList = testDataDict[self.FILENAME_LIST]
            if numericBaseName:
                # If numeric base name, the expected results need to be adapted as new files are generated in execution
                # time to simulate the scenario for that test without having to add more files to the test dataset
                filesList = [file.replace('stack', '') for file in filesList]
            incDoseList = testDataDict[self.INCOMING_DOSE_LIST]
            accumDoseList = testDataDict[self.ACCUM_DOSE_LIST]
            angleList = testDataDict[self.ANGLE_LIST]
            acqOrderList = testDataDict[self.ACQ_ORDER_LIST]

            previousAngle = None

            for i, tiM in enumerate(tsM):
                self.assertEqual(tiM.getFileName(), prot._getExtraPath(filesList[i]))
                self.assertTrue(os.path.islink(tiM.getFileName()),
                                "Tilt series file %s is not a link." % tiM.getFileName())
                self.assertAlmostEqual(tiM.getAcquisition().getDosePerFrame(), incDoseList[i], delta=0.0001)
                self.assertAlmostEqual(tiM.getAcquisition().getAccumDose(), accumDoseList[i], delta=0.0001)
                self.assertAlmostEqual(tiM.getTiltAngle(), angleList[i], delta=0.0001)
                # objId must start with 1 --> i +1
                self.assertEqual(i + 1, tiM.getObjId(), "Tilt image Movie objId is incorrect")
                # Acquisition order
                self.assertEqual(tiM.getAcquisitionOrder(), acqOrderList[i])

                if not isTsMovie:
                    self.assertEqual(i + 1, tiM.getIndex())
                    if previousAngle is not None:
                        self.assertTrue(previousAngle < tiM.getTiltAngle(), "Tilt images are not sorted by angle.")
                    previousAngle = tiM.getTiltAngle()

    def _runTestImportTsM(self, exclusionWords=None, outputSize=None):
        isTsMovie = True
        label = 'Import TsM'
        # Run protocol
        if exclusionWords:
            label += ' excluding words'

        protImport = self._runImportTiltSeries(self.parentDir, self.pattern, label=label, isTsMovie=isTsMovie,
                                               exclusionWords=exclusionWords)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeriesM', None)
        self._checkResults(outputSet, isTsMovie, (1152, 1152, 6), protImport, size=outputSize)  # ts and tsM have
        # different dimensions because they have been downsampled separately in order to get a lighter test dataset
        return protImport

    def test_importTiltSeriesM(self):
        self._runTestImportTsM(outputSize=2)

    def test_importTiltSeriesM_excludingWords(self):
        """There are 2 mdocs (which means 2 TsM), with tsIds = [stack10, stack31]. Let's exclude the first one
        using the functionality exclusionWords offered by the protocol"""
        protImport = self._runTestImportTsM(exclusionWords='stack10', outputSize=1)
        outputSet = getattr(protImport, 'outputTiltSeriesM', None)
        self.assertEqual(outputSet.getFirstItem().getTsId(), 'stack31')

    def test_importTiltSeries(self):
        isTsMovie = False
        # Run protocol
        protImport = self._runImportTiltSeries(self.parentDir, self.pattern, isTsMovie=isTsMovie)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeries', None)
        self._checkResults(outputSet, isTsMovie, (1440, 1024, 1), protImport, size=2)  # ts and tsM have different dims
        # because they have been downsampled separately in order to get a lighter test dataset

    def test_mdocFileCheckerRealOk(self):
        # Real data from EMPIAR --> OK files
        mdocList = glob.glob(join(self.dataset.getFile('empiarMdocDirOk'), '*.mdoc'))
        expectedErrorKeyWordList = None
        self._checkMDocParsingErrorMsg(mdocList, expectedErrorKeyWordList)

    def test_mdocFileCheckerRealNoOk(self):
        # Real data from EMPIAR --> No OK files
        dataSet = self.dataset.filesDict
        noOkMdocDir = self.dataset.getFile('empiarMdocDirNoOk')
        mdocList = [
            join(noOkMdocDir, dataSet['realFileNoVoltage1']),
            join(noOkMdocDir, dataSet['realFileNoVoltage2'])
        ]
        VOLTAGE = 'Voltage'
        expectedErrorKeyWordList = [
            VOLTAGE,  # Missing voltage
            VOLTAGE  # Missing voltage
        ]
        self._checkMDocParsingErrorMsg(mdocList, expectedErrorKeyWordList)

    def test_mdocFileCheckerSimpleError(self):
        # Edited data to simulate simple errors
        dataSet = self.dataset.filesDict
        simErrorMdocDir = self.dataset.getFile('simErrorMdocDir')
        mdocList = [
            join(simErrorMdocDir, dataSet['noMaginficationMdoc']),
            join(simErrorMdocDir, dataSet['noSamplingRateMdoc']),
            join(simErrorMdocDir, dataSet['noDoseMdoc'])
        ]
        expectedErrorKeyWordList = [
            'Magnification',  # Missing Magnification
            'PixelSpacing',  # Missing Sampling Rate
            'Dose'  # Not able to get the dose
        ]
        self._checkMDocParsingErrorMsg(mdocList, expectedErrorKeyWordList)

    def test_mdocFileCheckerMultipleErrors(self):
        # Edited data to simulate multiple errors, which are more than one acquisition magnitude or an angle specific
        # error present in more than one angular slice
        dataSet = self.dataset.filesDict
        simErrorMdocDir = self.dataset.getFile('simErrorMdocDir')
        mdocList = [
            join(simErrorMdocDir, dataSet['noVoltagenoSRateMdoc']),
            join(simErrorMdocDir, dataSet['someMissingAnglesMdoc'])
        ]
        expectedErrorKeyWordList = [
            ['Voltage', 'PixelSpacing'],  # Missing voltage and sampling rate
            'TiltAngle for Z values: 1, 7, 48'  # Missing tilt angles in slices 1, 7 and 48
        ]
        self._checkMDocParsingErrorMsg(mdocList, expectedErrorKeyWordList)

    def _checkMDocParsingErrorMsg(self, mdocList, expectedErrorKeyWordList):
        for i, mdoc in enumerate(mdocList):
            mdocObj = MDoc(mdoc)
            errorMsg = mdocObj.read(
                isImportingTsMovies=True,
                ignoreFilesValidation=True)  # Ignore files validation in order to make the dataset lighter

            if expectedErrorKeyWordList:  # There can be more than one error keyword to be checked per file
                keywords = expectedErrorKeyWordList[i]
                if type(keywords) is str:
                    keywords = [keywords]
                for errorKeyword in keywords:
                    self.assertTrue(errorKeyword in errorMsg,
                                    msg="%s not found in error message %s after validating %s."
                                        % (errorKeyword, errorMsg, mdoc))
            else:
                self.assertTrue(not errorMsg, "There are errors unexpected when validating the Mdocs: %s" % errorMsg)


class TestImportTomoCtfPattern(TestBaseCentralizedLayer):
    ds = None
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()  # 5 tilt-series are imported with the test default pattern

    @classmethod
    def _runImportTs(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords='output',
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runImportCtf(cls, filesPattern=None, exclusionWords=None, objLabel=None):
        print(magentaStr(f"\n==> Importing the CTFs:"
                         f"\n\t- Files pattern = {filesPattern}"
                         f"\n\t- Exclusion words = {exclusionWords}"))
        protImportCtf = cls.newProtocol(ProtImportTsCTF,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                        filesPattern=filesPattern,
                                        exclusionWords=exclusionWords,
                                        importFrom=ImportChoice.CTFFIND.value,
                                        inputSetOfTiltSeries=cls.importedTs)
        if objLabel:
            protImportCtf.setObjLabel(objLabel)
        cls.launchProtocol(protImportCtf)
        outputMask = getattr(protImportCtf, protImportCtf._possibleOutputs.CTFs.name, None)
        return outputMask


    def testImportCtfFromPattern01(self):
        exclusionWords = '01 43 45'  # Imported CTF should be TS_03 and TS_54
        ctfSet = self._runImportCtf(filesPattern='*/{TS}.defocus',
                                    exclusionWords=exclusionWords,
                                    objLabel='testImportCtfFromPattern01')
        self.checkCTFs(ctfSet, expectedSetSize=2)

    def testImportCtfFromPattern02(self):
        exclusionWords = '03 54'  # Imported CTF should be TS_01, TS_43, TS_45
        ctfSet = self._runImportCtf(filesPattern='*/{TS}.defocus',
                                    exclusionWords=exclusionWords,
                                    objLabel='testImportCtfFromPattern02')
        self.checkCTFs(ctfSet, expectedSetSize=3)
