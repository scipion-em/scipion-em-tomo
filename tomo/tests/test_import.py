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
from os.path import join, exists, abspath

import numpy

from tomo.convert import getOrderFromList
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr, createLink
from pyworkflow.object import Pointer
from pwem.protocols import ProtSplitSet, ProtSetFilter, ProtSetEditor

from tomo.protocols.protocol_ts_import import MDoc
from . import DataSet
from ..constants import BOTTOM_LEFT_CORNER, TOP_LEFT_CORNER, ERR_COORDS_FROM_SQLITE_NO_MATCH, ERR_NO_TOMOMASKS_GEN, \
    ERR_NON_MATCHING_TOMOS, SCIPION
import tomo.protocols
from ..protocols import ProtImportTomograms, ProtImportTomomasks
from ..protocols.protocol_import_coordinates import IMPORT_FROM_AUTO, ProtImportCoordinates3D
from ..protocols.protocol_import_coordinates_from_scipion import ProtImportCoordinates3DFromScipion, outputObjs
from ..utils import existsPlugin

from imod.protocols import ProtImodTomoNormalization


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
        cls.sqBoxSize = 20
        cls.sqSamplingRate = 13.68

    def _importTomograms(self, filesPath, samplingRate, pattern=None):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=filesPath,
                                              filesPattern=pattern,
                                              samplingRate=samplingRate)

        self.launchProtocol(protImportTomogram)

        self.assertSetSize(protImportTomogram.Tomograms, size=1, msg='No tomograms were generated.')

        return protImportTomogram.Tomograms

    def _runTomoImportSetOfCoordinates(self, pattern, program, ext):
        sRate = 5.5
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   objLabel='Import from %s - %s' % (program, ext),
                                                   auto=IMPORT_FROM_AUTO,
                                                   filesPath=self.tomoDs.getPath(),
                                                   importTomograms=self._importTomograms(self.tomoDs.getFile('tomo1'),
                                                                                         sRate),
                                                   filesPattern=pattern,
                                                   boxSize=32,
                                                   samplingRate=sRate)
        self.launchProtocol(protImportCoordinates3d)
        return protImportCoordinates3d

    def _runImportSetoOfCoordsFromScipionSqlite(self, sqliteFile, inTomos, boxSize, objLabel):
        protImportCoordsFromSqlite = self.newProtocol(ProtImportCoordinates3DFromScipion,
                                                      objLabel=objLabel,
                                                      sqliteFile=sqliteFile,
                                                      importTomograms=inTomos,
                                                      boxSize=boxSize)

        self.launchProtocol(protImportCoordsFromSqlite)
        return getattr(protImportCoordsFromSqlite, outputObjs.coordinates.name, None), \
               getattr(protImportCoordsFromSqlite, 'outputTomograms', None)

    def testImport3dCoordsFromSqlite_FullMatch(self):
        setSize = 2339
        inTomos = self._importTomograms(self.emdDs.getFile('tomoEmd10439'), self.sqSamplingRate)
        outputCoordsSet, outputTomoSet = self._runImportSetoOfCoordsFromScipionSqlite(
            self.emdDs.getFile('scipionSqlite3dCoords'), inTomos, self.sqBoxSize, 'Scipion - Full match')
        self.assertCoordinates(outputCoordsSet, setSize, self.sqBoxSize, self.sqSamplingRate)
        self.assertFalse(outputTomoSet)

        # Subset of coordinates 3d
        splitSet = self.newProtocol(ProtSplitSet,
                                    inputSet=outputCoordsSet)

        # Launch the split set protocol
        self.launchProtocol(splitSet)
        # test we have teh tomograms associatted
        self.assertCoordinates(splitSet.outputCoordinates3D01, 1170, self.sqBoxSize, self.sqSamplingRate)

        # Launch the filter set
        filterSet = self.newProtocol(ProtSetFilter, formula="True")
        filterSet.inputSet = Pointer(splitSet, extended="outputCoordinates3D01")

        self.launchProtocol(filterSet)
        # test the set is correct
        self.assertCoordinates(filterSet.outputCoordinates3D01, 1170, self.sqBoxSize, self.sqSamplingRate)

        # Launch the edit set
        editSet = self.newProtocol(ProtSetEditor, formula="True")
        editSet.inputSet = Pointer(filterSet, extended="outputCoordinates3D01")

        self.launchProtocol(editSet)
        # test the set is correct
        self.assertCoordinates(editSet.outputCoordinates3D01, 1170, self.sqBoxSize, self.sqSamplingRate)

    def testImport3dCoordsFromSqlite_SomeCoordsExcluded(self):
        setSize = 2335
        inTomos = self._importTomograms(self.emdDs.getFile('tomoEmd10439'), self.sqSamplingRate)
        outputCoordsSet, outputTomoSet = self._runImportSetoOfCoordsFromScipionSqlite(
            self.emdDs.getFile('scipionSqlite3dCoordsSomeBad'), inTomos, self.sqBoxSize,
            'Scipion - Some coords excluded')
        self.assertCoordinates(outputCoordsSet, setSize, self.sqBoxSize, self.sqSamplingRate)
        self.assertFalse(outputTomoSet)

    def testImport3dCoordsFromSqlite_SomeTomosAndCoordsExcluded(self):
        # Generate a symbolic link to another tomogram to have a set of two, and the coordinates referred only to one
        # of them
        setSize = 2335
        workingPath = self.getOutputPath()
        emTomogram = self.tomoDs.getFile('tomo1')
        mrcTomogram = self.emdDs.getFile('tomoEmd10439')

        for tomgramPath in [emTomogram, mrcTomogram]:

            fileName = os.path.basename(tomgramPath).replace(".ed", ".mrc")
            linkedTomo = join(workingPath, fileName)
            if not exists(linkedTomo):
                os.symlink(abspath(tomgramPath), abspath(linkedTomo))

        inTomos = self._importTomograms(workingPath, self.sqSamplingRate, pattern='*.mrc')
        outputCoordsSet, outputTomoSet = self._runImportSetoOfCoordsFromScipionSqlite(
            self.emdDs.getFile('scipionSqlite3dCoordsSomeBad'), inTomos, self.sqBoxSize,
            'Scipion - Some tomos and coords excluded')
        self.assertCoordinates(outputCoordsSet, setSize, self.sqBoxSize, self.sqSamplingRate)

    def testImport3dCoordsFromSqlite_NoneMatch(self):
        inTomos = self._importTomograms(self.tomoDs.getFile('tomo1'), self.sqSamplingRate)
        with self.assertRaises(Exception) as eType:
            self._runImportSetoOfCoordsFromScipionSqlite(
                self.emdDs.getFile('scipionSqlite3dCoords'), inTomos, self.sqBoxSize, 'Scipion - No match')
            self.assertEqual(str(eType.exception), ERR_COORDS_FROM_SQLITE_NO_MATCH)

    def test_import_set_of_coordinates_3D(self):
        boxSize = 32
        samplingRate = 5.5

        # From txt
        protCoordinates = self._runTomoImportSetOfCoordinates('*.txt', 'TOMO', 'TXT')
        output = getattr(protCoordinates, 'outputCoordinates', None)
        self.assertCoordinates(output, 5, boxSize, samplingRate)

        # From txt
        protCoordinates = self._runTomoImportSetOfCoordinates('*.cbox', 'CRYOLO', 'CBOX')
        output = getattr(protCoordinates, 'outputCoordinates', None)
        self.assertCoordinates(output, 3, boxSize, samplingRate)

        def checkCoordinates(expectedValues, coordSet):

            for index, row in enumerate(expectedValues):
                coord = coordSet[index + 1]
                self.assertEqual(coord.getX(SCIPION), row[0], "X coordinate not converted properly")
                self.assertEqual(coord.getY(SCIPION), row[1], "Y coordinate not converted properly")
                self.assertEqual(coord.getZ(SCIPION), row[2], "Z coordinate not converted properly")

        expectedCoords = [
            [-512, -512, -256],
            [0, 0, 0],
            [512, 512, 256],
        ]

        checkCoordinates(expectedCoords, output)

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

        self.assertIsNotNone(coordSet.getPrecedents(), "Tomograms not associated in the output set")

        for coord in coordSet.iterCoordinates():
            # Access the coordinate, this should work if tomograms are associated
            X = coord.getY(BOTTOM_LEFT_CORNER)
            Y = coord.getY(TOP_LEFT_CORNER)

            break


class TestTomoImportTomograms(BaseTest):
    """This class check if the protocol to import tomograms works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('*.em')

    def _runImportTomograms(self)-> tomo.protocols.ProtImportTomograms:
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTomograms,
            filesPath=self.tomogram,
            filesPattern='',
            acquisitionAngleMax=40,
            acquisitionAngleMin=-40,
            samplingRate=1.35)
        self.launchProtocol(protImport)
        return protImport

    def _runImportTomograms2(self)->tomo.protocols.ProtImportTomograms:
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTomograms,
            filesPath=self.tomogram,
            filesPattern='',
            samplingRate=1.35)
        self.launchProtocol(protImport)
        return protImport

    def test_importTomograms(self):
        protImport = self._runImportTomograms()
        self.assertSetSize(protImport.Tomograms, size=2, msg="There was a problem with Import Tomograms protocol")

        for tomogram in protImport.Tomograms.iterItems():
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
        self.assertSetSize(protImport.Tomograms, size=2, msg=
                             "There was a problem with Import Tomograms protocol")

        for tomogram in protImport2.Tomograms.iterItems():
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
        self.assertSetSize(protImport.outputTiltSeriesM, size=2, msg="Importing tilt series movies is failing")

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

    def _runImportTiltSeries(self, isTsMovie=False, exclusionWords=None):
        prot = tomo.protocols.ProtImportTsMovies if isTsMovie else tomo.protocols.ProtImportTs
        attribDict = {
            'filesPath': self.parentDir,
            'filesPattern': self.pattern,
        }
        if exclusionWords:
            attribDict['exclusionWords'] = exclusionWords
        protImport = self.newProtocol(
            prot,
            **attribDict)
        self.launchProtocol(protImport)
        return protImport

    def _checkResults(self, outputSet, isTsMovie, dimensions, prot, size=None):
        # Generate the test data
        testDict = self._genTestData(isTsMovie)
        # Check set properties
        self.assertSetSize(outputSet, size=size)
        self.assertAlmostEqual(outputSet.getSamplingRate(), self.sRate, delta=0.001)
        self.assertEqual(outputSet.getDimensions(), dimensions)

        # Check tilt series movies
        for tsM in outputSet:
            testDataDict = testDict[tsM.getTsId()]
            # Check set acquisition
            acq = tsM.getAcquisition()
            self.assertAlmostEqual(tsM.getSamplingRate(), testDataDict[self.PIXEL_SIZE], delta=0.001)
            self.assertAlmostEqual(acq.getVoltage(), testDataDict[self.VOLTAGE], delta=0.1)
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
                    self.assertEqual(i + 1, tiM.getIndex())
                    if previousAngle is not None:
                        self.assertTrue(previousAngle < tiM.getTiltAngle(), "Tilt images are not sorted by angle.")
                    previousAngle = tiM.getTiltAngle()

    def _runTestImportTsM(self, exclusionWords=None, outputSize=None):
        isTsMovie = True
        label = 'Import TsM'
        # Run protocol
        protImport = self._runImportTiltSeries(isTsMovie=isTsMovie, exclusionWords=exclusionWords)
        if exclusionWords:
            label += ' excluding words'
        protImport.setObjLabel(label)
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
        protImport = self._runImportTiltSeries(isTsMovie=isTsMovie)
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


class TestImportTomoMasks(BaseTest):
    outputPath = None
    ds = None
    samplingRate = 13.68
    inTomoSet = None
    inTomoSetBinned = None
    inNotMatchingTomoSet = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        ds = DataSet.getDataSet('emd_10439')
        cls.ds = ds
        cls.inTomoSet = cls._importTomograms(ds.getFile('tomoEmd10439'))
        cls.inTomoSetBinned = cls._normalizeTomo()
        # Create a link pointing to the tomogram but renaming it in a way that it won't be possible to match it with
        # the tomomask
        nonMatchingTomoName = join(ds.getPath(), 'nonMatchingTomo.mrc')
        createLink(ds.getFile('tomoEmd10439'), nonMatchingTomoName)
        cls.inNotMatchingTomoSet = cls._importTomograms(nonMatchingTomoName)

    @classmethod
    def _importTomograms(cls, filesPath):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=filesPath,
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        cls.assertIsNotNone(protImportTomogram.Tomograms, 'No tomograms were generated importing %s.' % filesPath)

        return protImportTomogram.Tomograms

    @classmethod
    def _normalizeTomo(cls):
        print(magentaStr("\n==> Tomogram normalization:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=cls.inTomoSet,
                                                binning=2)

        cls.launchProtocol(protTomoNormalization)
        outputTomos = getattr(protTomoNormalization, 'Tomograms', None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were generated in tomo normalization.')

        return outputTomos

    def testImportTomoMasksAllGood(self):
        print(magentaStr("\n==> Importing data - tomoMasks:"))
        protImportTomomasks = self.newProtocol(ProtImportTomomasks,
                                               filesPath=self.ds.getFile('tomoMaskAnnotated'),
                                               inputTomos=self.inTomoSetBinned)

        self.launchProtocol(protImportTomomasks)
        tomoMaskSet = getattr(protImportTomomasks, 'outputTomoMasks', None)

        self.assertSetSize(tomoMaskSet, size=1, msg="Importing tomoMasks failed.")
        self.assertEqual(tomoMaskSet.getDimensions(), (464, 464, 250))
        self.assertEqual(tomoMaskSet.getSamplingRate(), 2 * self.samplingRate)
        self.assertFalse(protImportTomomasks.warnMsg)

    # The tomogram and the tomomask have different sizes
    def testImportTomoMasksDiffSize(self):
        print(magentaStr("\n==> Importing data - tomoMasks:"))
        protImportTomomasks = self.newProtocol(ProtImportTomomasks,
                                               filesPath=self.ds.getFile('tomoMaskAnnotated'),
                                               inputTomos=self.inTomoSet)

        with self.assertRaises(Exception) as eType:
            self.launchProtocol(protImportTomomasks)
            self.assertEqual(str(eType.exception), ERR_NO_TOMOMASKS_GEN)

    # Tomomask and tomogram doesn't correspond
    def testImportTomoMasksNoneMatch(self):
        print(magentaStr("\n==> Importing data - tomoMasks:"))
        protImportTomomasks = self.newProtocol(ProtImportTomomasks,
                                               filesPath=self.ds.getFile('tomoMaskAnnotated'),
                                               inputTomos=self.inNotMatchingTomoSet)

        with self.assertRaises(Exception) as eType:
            self.launchProtocol(protImportTomomasks)
            self.assertEqual(str(eType.exception), ERR_NON_MATCHING_TOMOS)

class TestConvert(BaseTest):

    def test_list_order(self):

        myList = [20,40,50,10,30]
        orders = getOrderFromList(myList)
        self.assertEqual(orders, [2,4,5,1,3])

        myList = [20, 40, 50, 10, 30]
        orders = getOrderFromList(myList)
        self.assertEqual(orders, [2, 4, 5, 1, 3])

        # Simulate dose list in symmetric dose scheme
        myList = [24,21,12,9,3,6,15,18,27]
        orders = getOrderFromList(myList)
        self.assertEqual(orders, [8,7,4,3,1,2,5,6,9])



