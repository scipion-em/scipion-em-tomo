# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
from os.path import join

import numpy
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.protocols.protocol_ts_import import MDoc

from . import DataSet
from ..constants import BOTTOM_LEFT_CORNER, ERR_COORDS_FROM_SQLITE_NO_MATCH
from ..protocols.protocol_import_coordinates import IMPORT_FROM_AUTO, ProtImportCoordinates3D
from ..protocols.protocol_import_coordinates_from_scipion import ProtImportCoordinates3DFromScipion
from ..utils import existsPlugin
import tomo.protocols


class TestTomoImportSubTomograms(BaseTest):
    """ This class check if the protocol to import sub tomograms works
     properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        # cls.tomogram = cls.dataset.getFile('tomo1')
        # cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
        cls.table = cls.dataset.getFile('initial.tbl')
        cls.path = cls.dataset.getPath()
        cls.subtomos = cls.dataset.getFile('basename.hdf')

    def _runImportSubTomograms(self):
        # protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
        #                                       filesPath=self.tomogram,
        #                                       samplingRate=5)
        # self.launchProtocol(protImportTomogram)
        #
        # protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
        #                          auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
        #                          filesPath=self.coords3D,
        #                          importTomograms=protImportTomogram.outputTomograms,
        #                          filesPattern='', boxSize=32,
        #                          samplingRate=5)
        # self.launchProtocol(protImportCoordinates3d)

        protImport = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.subtomos,
                                      samplingRate=1.35)
        # importCoordinates=protImportCoordinates3d.outputCoordinates)
        self.launchProtocol(protImport)
        return protImport

    def _runImportSubTomograms2(self):
        # protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
        #                                       filesPath=self.tomogram,
        #                                       filesPattern='',
        #                                       samplingRate=5)
        # self.launchProtocol(protImportTomogram)
        #
        # protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
        #                          auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
        #                          filesPath=self.coords3D,
        #                          importTomograms=protImportTomogram.outputTomograms,
        #                          filesPattern='', boxSize=32,
        #                          samplingRate=5)
        # self.launchProtocol(protImportCoordinates3d)

        protImport = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.subtomos,
                                      samplingRate=1.35,
                                      acquisitionAngleMax=40,
                                      acquisitionAngleMin=-40)
        # importCoordinates=protImportCoordinates3d.outputCoordinates)
        self.launchProtocol(protImport)
        return protImport

    def test_import_sub_tomograms(self):
        protImport = self._runImportSubTomograms()
        output = getattr(protImport, 'outputSubTomograms', None)
        self.assertTrue(output.getSamplingRate() == 1.35)
        self.assertTrue(output.getFirstItem().getSamplingRate() == 1.35)
        self.assertTrue(output.getDim()[0] == 32)
        self.assertTrue(output.getDim()[1] == 32)
        self.assertTrue(output.getDim()[2] == 32)
        self.assertEqual(output.getFirstItem().getAcquisition().getAngleMax(), 60)
        self.assertEqual(output.getFirstItem().getAcquisition().getAngleMin(), -60)
        # self.assertTrue(output.getFirstItem().getCoordinate3D().getX() == 314)
        # self.assertTrue(output.getFirstItem().getCoordinate3D().getY() == 350)
        # self.assertTrue(output.getFirstItem().getCoordinate3D().getZ() == 256)
        self.assertIsNotNone(output,
                             "There was a problem with Import SubTomograms protocol")

        protImport2 = self._runImportSubTomograms2()
        output2 = getattr(protImport2, 'outputSubTomograms', None)
        self.assertIsNotNone(output2,
                             "There was a problem with Import SubTomograms protocol")
        self.assertTrue(output2.getSamplingRate() == 1.35)
        self.assertTrue(output2.getDim()[0] == 32)
        self.assertTrue(output2.getDim()[1] == 32)
        self.assertTrue(output2.getDim()[2] == 32)
        self.assertEqual(output2.getFirstItem().getAcquisition().getAngleMax(), 40)
        self.assertEqual(output2.getFirstItem().getAcquisition().getAngleMin(), -40)
        # for i, subtomo in enumerate(output2.iterItems()):
        #     if i == 1:
        #         self.assertTrue(subtomo.getCoordinate3D().getX() == 174)
        #         self.assertTrue(subtomo.getCoordinate3D().getY() == 172)
        #         self.assertTrue(subtomo.getCoordinate3D().getZ() == 256)
        #     if i == 0:
        #         self.assertTrue(subtomo.getCoordinate3D().getX() == 314)
        #         self.assertTrue(subtomo.getCoordinate3D().getY() == 350)
        #         self.assertTrue(subtomo.getCoordinate3D().getZ() == 256)

        return output2


class TestTomoImportSetOfCoordinates3D(BaseTest):
    """This class check if the protocol to import set of coordinates 3d works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.tomoDs = DataSet.getDataSet('tomo-em')
        cls.emdDs = DataSet.getDataSet('emd_10439')

    def _importTomograms(self, filesPath, samplingRate):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=filesPath,
                                              samplingRate=samplingRate)

        self.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, 'outputTomograms', None)
        self.assertIsNotNone(outputTomos, 'No tomograms were genetated.')

        return outputTomos

    def _runTomoImportSetOfCoordinates(self, pattern, program, ext):
        sRate = 5.5
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   objLabel='Import from %s - %s' % (program, ext),
                                                   auto=IMPORT_FROM_AUTO,
                                                   filesPath=self.tomoDs.getPath(),
                                                   importTomograms=self._importTomograms(self.tomoDs.getFile('tomo1'), sRate),
                                                   filesPattern=pattern,
                                                   boxSize=32,
                                                   samplingRate=sRate)
        self.launchProtocol(protImportCoordinates3d)
        return protImportCoordinates3d

    def _runImportSetoOfCoordsFromScipionSqlite(self, inTomos, boxSize):
        protImportCoordsFromSqlite = self.newProtocol(ProtImportCoordinates3DFromScipion,
                                                      objLabel='Import from Scipion sqlite',
                                                      sqliteFile=self.emdDs.getFile('scipionSqlite3dCoords'),
                                                      importTomograms=inTomos,
                                                      boxSize=boxSize)

        self.launchProtocol(protImportCoordsFromSqlite)
        return getattr(protImportCoordsFromSqlite, 'outputSetOfCoordinates', None)

    def testImport3dCoordsFromSqlite_AllMatch(self):
        boxSize = 20
        samplingRate = 13.68
        setSize = 2339
        inTomos = self._importTomograms(self.emdDs.getFile('tomoEmd10439'), samplingRate)
        outputCoordsSet = self._runImportSetoOfCoordsFromScipionSqlite(inTomos, boxSize)
        self.assertCoordinates(outputCoordsSet, setSize, boxSize, samplingRate)

    def testImport3dCoordsFromSqlite_NoneMatch(self):
        inTomos = self._importTomograms(self.tomoDs.getFile('tomo1'), 13.68)
        with self.assertRaises(Exception) as eType:
            self._runImportSetoOfCoordsFromScipionSqlite(inTomos, 20)
            self.assertEqual(str(eType.exception), ERR_COORDS_FROM_SQLITE_NO_MATCH)

    def test_import_set_of_coordinates_3D(self):
        boxSize = 32
        samplingRate = 5.5

        # From txt
        protCoordinates = self._runTomoImportSetOfCoordinates('*.txt', 'TOMO', 'TXT')
        output = getattr(protCoordinates, 'outputCoordinates', None)
        self.assertCoordinates(output, 5, boxSize, samplingRate)

        # From emantomo file
        if existsPlugin('emantomo'):
            protCoordinates = self._runTomoImportSetOfCoordinates('*.json', 'EMAN', 'JSON')
            output = getattr(protCoordinates, 'outputCoordinates', None)
            self.assertCoordinates(output, 19, boxSize, samplingRate)

            # Check se are not loosing precision
            firstCoord = output.getFirstItem()
            # First row --> 224, 316, 260
            emanCoords = [224, 316, 260]
            self.assertEqual(firstCoord.getX(BOTTOM_LEFT_CORNER), emanCoords[0], "eman coordinate x has a wrong value")
            self.assertEqual(firstCoord.getY(BOTTOM_LEFT_CORNER), emanCoords[1], "eman coordinate y has a wrong value")
            self.assertEqual(firstCoord.getZ(BOTTOM_LEFT_CORNER), emanCoords[2], "eman coordinate z has a wrong value")

        # From dynamo file
        if existsPlugin('dynamo'):
            protCoordinates = self._runTomoImportSetOfCoordinates('*.tbl', 'DYNAMO', 'TBL')
            output = getattr(protCoordinates, 'outputCoordinates', None)
            self.assertCoordinates(output, 5, boxSize, samplingRate)

    def assertCoordinates(self, coordSet, size, boxSize, samplingRate):
        self.assertTrue(coordSet, 'No 3d coordinates were generated.')
        self.assertSetSize(coordSet, size=size)
        self.assertTrue(coordSet.getBoxSize() == boxSize)
        self.assertTrue(coordSet.getSamplingRate() == samplingRate)


class TestTomoImportTomograms(BaseTest):
    """This class check if the protocol to import tomograms works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('*.em')

    def _runImportTomograms(self):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTomograms,
            filesPath=self.tomogram,
            filesPattern='',
            acquisitionAngleMax=40,
            acquisitionAngleMin=-40,
            samplingRate=1.35)
        self.launchProtocol(protImport)
        return protImport

    def _runImportTomograms2(self):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTomograms,
            filesPath=self.tomogram,
            filesPattern='',
            samplingRate=1.35)
        self.launchProtocol(protImport)
        return protImport

    def test_importTomograms(self):
        protImport = self._runImportTomograms()
        output = getattr(protImport, 'outputTomograms', None)
        self.assertIsNotNone(output,
                             "There was a problem with Import Tomograms protocol")

        for tomogram in protImport.outputTomograms.iterItems():
            self.assertTrue(tomogram.getXDim() == 1024,
                            "There was a problem with Import Tomograms protocol")
            self.assertIsNotNone(tomogram.getYDim() == 1024,
                                 "There was a problem with Import Tomograms protocol")

            self.assertTrue(tomogram.getAcquisition().getAngleMax() == 40,
                            "There was a problem with the acquisition angle max")
            self.assertTrue(tomogram.getAcquisition().getAngleMin() == -40,
                            "There was a problem with the acquisition angle min")

            break

        protImport2 = self._runImportTomograms2()
        output2 = getattr(protImport2, 'outputTomograms', None)
        self.assertIsNotNone(output2,
                             "There was a problem with Import Tomograms protocol")

        for tomogram in protImport2.outputTomograms.iterItems():
            self.assertTrue(tomogram.getXDim() == 1024,
                            "There was a problem with Import Tomograms protocol")
            self.assertIsNotNone(tomogram.getYDim() == 1024,
                                 "There was a problem with Import Tomograms protocol")

            self.assertTrue(tomogram.getAcquisition().getAngleMax() == 60,
                            "There was a problem with the acquisition angle max")
            self.assertTrue(tomogram.getAcquisition().getAngleMin() == -60,
                            "There was a problem with the acquisition angle min")

            break

class TestTomoBaseProtocols(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.getFile = cls.dataset.getFile('etomo')
        cls.getFileM = cls.dataset.getFile('empiar')

    def _runImportTiltSeriesM(self, filesPattern='{TS}_{TO}_{TA}.mrc'):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTsMovies,
            filesPath=self.getFileM,
            filesPattern=filesPattern,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=0.675,
            doseInitial=0,
            dosePerFrame=0.375,
            tiltAxisAngle=84.1
        )
        self.launchProtocol(protImport)
        return protImport

    def test_importTiltSeriesM(self):
        protImport = self._runImportTiltSeriesM()
        self.assertSetSize(protImport.outputTiltSeriesM, 2)

        return protImport

    def test_motioncorTiltSeriesM(self):
        protImport = self.test_importTiltSeriesM()

        # --------- Motion correction with motioncor2 for Tilt-series ------
        protMc = self.newProtocol(
            tomo.protocols.ProtTsAverage,
            binFactor=2.0
        )

        protMc.inputTiltSeriesM.set(protImport.outputTiltSeriesM)
        self.launchProtocol(protMc)


class TestTomoImportTsFromPattern(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.getFile = cls.dataset.getFile('etomo')
        cls.getFileM = cls.dataset.getFile('empiar')
        cls.tiltAxisAngle = 84.1

    def _runImportTiltSeriesM(self, filesPattern='{TS}_{TO}_{TA}.mrc'):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTsMovies,
            filesPath=self.getFileM,
            filesPattern=filesPattern,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=0.675,
            doseInitial=0,
            dosePerFrame=0.375,
            tiltAxisAngle=self.tiltAxisAngle
        )
        self.launchProtocol(protImport)
        return protImport

    def _runImportTiltSeries(self):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTs,
            filesPath=self.getFile,
            filesPattern='BB{TS}.st',
            minAngle=-55,
            maxAngle=65,
            stepAngle=2,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=1.35,
            doseInitial=0,
            dosePerFrame=0.3,
            tiltAxisAngle=self.tiltAxisAngle)
        self.launchProtocol(protImport)
        return protImport

    def test_importTiltSeriesM(self):
        protImport = self._runImportTiltSeriesM()
        self.checkTSSet(protImport.outputTiltSeriesM, 2, 3)

        return protImport

    def test_importTiltSeries(self):
        protImport = self._runImportTiltSeries()
        self.checkTSSet(protImport.outputTiltSeries, 2, 61, checkIds=True)

    def checkTSSet(self, set, size, anglesCount, checkIds=False):
        """
            Check basic attributes of a TS set

            :param set: TiltSeries set (Movies or Images)
            :param size: Expected size
            :param anglesCount: Expected number of tilts
            :param checkIds: check if ids start with 1 and increments by one
            :return: None
            """

        self.assertSetSize(set, size)
        for ts in set:
            self.assertEquals(ts.getSize(), anglesCount, "Size of tilt images is wrong.")
            self.assertEqual(ts.getAcquisition().getTiltAxisAngle(), self.tiltAxisAngle)
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
            cls.ACQ_ORDER_LIST: kwargs.get(cls.ACQ_ORDER_LIST, None)
        }

    def sortTestDataByAngle(self, testDataDict):

        for key in testDataDict:
            testData = testDataDict[key]
            ind = numpy.argsort(testData[self.ANGLE_LIST])

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
                angles=[0.0036, 2.9683, -3.0250, -6.0251, 5.9684],
                acqOrder=[1, 2, 3, 4, 5]
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
                angles=[0.0026, 2.9708, -3.0255, -6.0241, 5.9684],
                acqOrder=[1, 2, 3, 4, 5]
            )
        }

        if not isTsMovie:
            self.sortTestDataByAngle(testData)

        return testData

    def _runImportTiltSeries(self, isTsMovie=False):
        prot = tomo.protocols.ProtImportTsMovies if isTsMovie else tomo.protocols.ProtImportTs
        protImport = self.newProtocol(
            prot,
            filesPath=self.parentDir,
            filesPattern=self.pattern)
        self.launchProtocol(protImport)
        return protImport

    def _checkResults(self, outputSet, isTsMovie, dimensions, prot):
        # Generate the test data
        testDict = self._genTestData(isTsMovie)
        # Check set properties
        self.assertSetSize(outputSet, size=2)
        self.assertAlmostEqual(outputSet.getSamplingRate(), self.sRate, delta=0.001)
        self.assertEqual(outputSet.getDimensions(), dimensions)

        # Check tilt series movies
        for tsM in outputSet:
            testDataDict = testDict[tsM.getTsId()]
            # Check set acquisition
            acq = tsM.getAcquisition()
            self.assertAlmostEqual(tsM.getSamplingRate(), testDataDict[self.PIXEL_SIZE], delta=0.001)
            self.assertAlmostEqual(acq.getVoltage(), testDataDict[self.VOLTAGE], delta=0.1)
            self.assertAlmostEqual(acq.getMagnification(), testDataDict[self.MAGNIFICATION], delta=0.1)
            self.assertAlmostEqual(acq.getSphericalAberration(), testDataDict[self.SPH_ABERRATION], delta=0.01)
            self.assertAlmostEqual(acq.getAmplitudeContrast(), testDataDict[self.AMP_CONTRAST], delta=0.001)
            self.assertAlmostEqual(acq.getDosePerFrame(), testDataDict[self.DOSE_PER_FRAME], delta=0.0001)
            self.assertAlmostEqual(acq.getAccumDose(), testDataDict[self.ACCUM_DOSE], delta=0.0001)
            # Check angles and accumulated dose per angular acquisition (tilt series image)
            filesList = testDataDict[self.FILENAME_LIST]
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
                    self.assertEqual(i+1, tiM.getIndex())
                    if previousAngle is not None:
                        self.assertTrue(previousAngle < tiM.getTiltAngle(), "Tilt images are not sorted by angle.")
                    previousAngle = tiM.getTiltAngle()

    def test_importTiltSeriesM(self):
        isTsMovie = True
        # Run protocol
        protImport = self._runImportTiltSeries(isTsMovie=isTsMovie)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeriesM', None)
        self._checkResults(outputSet, isTsMovie, (1152, 1152, 6), protImport)  # ts and tsM has different dimensions
        # because they have been downsampled separately in order to get a lighter test dataset

    def test_importTiltSeries(self):
        isTsMovie = False
        # Run protocol
        protImport = self._runImportTiltSeries(isTsMovie=isTsMovie)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeries', None)
        self._checkResults(outputSet, isTsMovie, (1440, 1024, 1), protImport)  # ts and tsM have different dimensions
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
        VOLTAGE = '*Voltage*'
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
            '*Magnification*',  # Missing Magnification
            '*PixelSpacing*',  # Missing Sampling Rate
            '*dose*'  # Not able to get the dose
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
            ['*Voltage*', '*PixelSpacing*'],  # Missing voltage and sampling rate
            '*TiltAngle*: 1 7 48'  # Missing tilt angles in slices 1, 7 and 48
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
                    self.assertTrue(errorKeyword in errorMsg)
            else:
                self.assertTrue(not errorMsg)
