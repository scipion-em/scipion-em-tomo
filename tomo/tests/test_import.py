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
import os
from os.path import join, exists, abspath
import numpy as np
from tomo.convert import getOrderFromList
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr, createLink
from pyworkflow.object import Pointer
from pwem.protocols import ProtSplitSet, ProtSetFilter, ProtSetEditor
from tomo.protocols.protocol_ts_import import MDoc, ProtImportTs
from . import DataSet, RE4_STA_TUTO, DataSetRe4STATuto, EMD_10439, DataSetEmd10439
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
        cls.sRate = 5.5
        cls.boxSize = 32

    def _importTomograms(self, filesPath, samplingRate, pattern=None):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=filesPath,
                                              filesPattern=pattern,
                                              samplingRate=samplingRate)

        self.launchProtocol(protImportTomogram)

        self.assertSetSize(protImportTomogram.Tomograms, size=1, msg='No tomograms were generated.')

        return protImportTomogram.Tomograms

    def _runImportSetOfCoordinates(self, pattern, program, ext, importTomograms=None, coordsSRate=None,
                                   boxSize=None, extraText=''):
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_AUTO,
                                                   filesPath=self.tomoDs.getPath(),
                                                   importTomograms=importTomograms,
                                                   filesPattern=pattern,
                                                   boxSize=boxSize,
                                                   samplingRate=coordsSRate)
        protImportCoordinates3d.setObjLabel('Import from %s - %s ' % (program, ext) + extraText)
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
        OUT_COORDS = 'outputCoordinates'
        importedTomos = self._importTomograms(self.tomoDs.getFile('tomo1'), self.sRate)

        def checkCoordinates(expectedValues, coordSet):

            for index, row in enumerate(expectedValues):
                coord = coordSet[index + 1]
                self.assertEqual(coord.getX(SCIPION), row[0], "X coordinate not converted properly")
                self.assertEqual(coord.getY(SCIPION), row[1], "Y coordinate not converted properly")
                self.assertEqual(coord.getZ(SCIPION), row[2], "Z coordinate not converted properly")

        # From txt -----------------------------------------------------------------------------------------------------
        nCoords = 5
        tomoDim = [1024, 1024, 512]

        def getExpectedCoords():
            expCoords = np.array([
                [314, 350, 256],
                [174, 172, 256],
                [413, 346, 255],
                [46, 189, 255],
                [167, 224, 255]
            ], dtype=float)
            # The TXT case is coded with an origin of type BOTTOM_LEFT_CORNER
            for i in range(len(tomoDim)):
                expCoords[:, i] = expCoords[:, i] * scaleFactor - tomoDim[i] / 2
            return expCoords

        def runImportTxtCoords(sRate, boxSize, extraText):
            protImportCoordinates = self._runImportSetOfCoordinates('*.txt', 'TOMO', 'TXT',
                                                                    importTomograms=importedTomos,
                                                                    coordsSRate=sRate,
                                                                    boxSize=boxSize,
                                                                    extraText=extraText)
            return getattr(protImportCoordinates, OUT_COORDS, None)

        # 1) Coordinates with a sampling rate equal to the tomograms' sampling rate
        scaleFactor = 1
        output = runImportTxtCoords(self.sRate, self.boxSize, 'sRate equal')
        self.assertCoordinates(output, nCoords, self.boxSize, self.sRate)
        expectedValues = getExpectedCoords()
        checkCoordinates(expectedValues, output)

        # 2) Coordinates' sampling rate not provided --> assume it's the same as the tomograms' sampling rate
        scaleFactor = 1
        output = runImportTxtCoords(self.sRate, self.boxSize, 'sRate not provided')
        self.assertCoordinates(output, nCoords, self.boxSize, self.sRate)
        expectedValues = getExpectedCoords()
        checkCoordinates(expectedValues, output)

        # 3) Coordinates' sampling rate is half of the tomograms' sampling rate
        scaleFactor = 0.5
        output = runImportTxtCoords(scaleFactor * self.sRate, self.boxSize, 'sRate half')
        self.assertCoordinates(output, nCoords, scaleFactor * self.boxSize, self.sRate)
        expectedValues = getExpectedCoords()
        checkCoordinates(expectedValues, output)

        # 4) Coordinates' sampling rate is double of the tomograms' sampling rate
        scaleFactor = 2
        output = runImportTxtCoords(scaleFactor * self.sRate, self.boxSize, 'sRate double')
        self.assertCoordinates(output, nCoords, scaleFactor * self.boxSize, self.sRate)
        expectedValues = getExpectedCoords()
        checkCoordinates(expectedValues, output)

        # From cbox ----------------------------------------------------------------------------------------------------
        # Cryolo's test data cbox file provides a box of [20, 20, 1] --> expected coords = coords + boxSize/2
        expectedCoords = [
            [-502, -502, -256],
            [10, 10, 0],
            [522, 522, 256],
        ]
        protCoordinates = self._runImportSetOfCoordinates('*.cbox', 'CRYOLO', 'CBOX',
                                                          importTomograms=importedTomos,
                                                          coordsSRate=self.sRate,
                                                          boxSize=self.boxSize)
        output = getattr(protCoordinates, OUT_COORDS, None)
        self.assertCoordinates(output, 3, self.boxSize, self.sRate)
        checkCoordinates(expectedCoords, output)

        # From emantomo file -------------------------------------------------------------------------------------------
        if existsPlugin('emantomo'):
            protCoordinates = self._runImportSetOfCoordinates('*.json', 'EMAN', 'JSON',
                                                              importTomograms=importedTomos,
                                                              coordsSRate=self.sRate,
                                                              boxSize=self.boxSize)
            output = getattr(protCoordinates, OUT_COORDS, None)
            self.assertCoordinates(output, 19, self.boxSize, self.sRate)

            # Check se are not loosing precision
            firstCoord = output.getFirstItem()
            # First row --> 224, 316, 260
            emanCoords = [224, 316, 260]
            self.assertEqual(firstCoord.getX(BOTTOM_LEFT_CORNER), emanCoords[0], "eman coordinate x has a wrong value")
            self.assertEqual(firstCoord.getY(BOTTOM_LEFT_CORNER), emanCoords[1], "eman coordinate y has a wrong value")
            self.assertEqual(firstCoord.getZ(BOTTOM_LEFT_CORNER), emanCoords[2], "eman coordinate z has a wrong value")

        # From dynamo file ---------------------------------------------------------------------------------------------
        if existsPlugin('dynamo'):
            protCoordinates = self._runImportSetOfCoordinates('*.tbl', 'DYNAMO', 'TBL',
                                                              importTomograms=importedTomos,
                                                              coordsSRate=self.sRate,
                                                              boxSize=self.boxSize)
            output = getattr(protCoordinates, OUT_COORDS, None)
            self.assertCoordinates(output, 5, self.boxSize, self.sRate)

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

    def _runImportTomograms(self) -> tomo.protocols.ProtImportTomograms:
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTomograms,
            filesPath=self.tomogram,
            filesPattern='',
            acquisitionAngleMax=40,
            acquisitionAngleMin=-40,
            samplingRate=1.35)
        self.launchProtocol(protImport)
        return protImport

    def _runImportTomograms2(self) -> tomo.protocols.ProtImportTomograms:
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
        ds = DataSet.getDataSet(EMD_10439)
        tomoFile = ds.getFile(DataSetEmd10439.tomoEmd10439.name)
        cls.ds = ds
        cls.inTomoSet = cls._importTomograms(tomoFile)
        cls.inTomoSetBinned = cls._normalizeTomo()
        # Create a link pointing to the tomogram but renaming it in a way that it won't be possible to match it with
        # the tomomask
        nonMatchingTomoName = join(ds.getPath(), 'nonMatchingTomo.mrc')
        createLink(ds.getFile(tomoFile), nonMatchingTomoName)
        cls.inNotMatchingTomoSet = cls._importTomograms(nonMatchingTomoName)

    @classmethod
    def _importTomograms(cls, filesPath):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=filesPath,
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        tomograms = getattr(protImportTomogram, OUTPUT_NAME, None)
        cls.assertIsNotNone(tomograms, 'No tomograms were generated importing %s.' % filesPath)
        return tomograms

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

    @classmethod
    def _importTomoMasks(cls, filesPath, inputTomos, filesPattern=None,
                         descriptionMsg='', protLabel=''):
        print(magentaStr(f"\n==> Importing data - tomoMasks:\n{descriptionMsg}"))
        protImportTomomasks = cls.newProtocol(ProtImportTomomasks,
                                              filesPath=filesPath,
                                              filesPattern=filesPattern,
                                              inputTomos=inputTomos)
        if protLabel:
            protImportTomomasks.setObjLabel(protLabel)
        cls.launchProtocol(protImportTomomasks)
        return protImportTomomasks

    def _checkResults(self, protImportTomomasks):
        tomoMaskSet = getattr(protImportTomomasks, protImportTomomasks._possibleOutputs.tomomasks.name, None)
        self.assertSetSize(tomoMaskSet, size=1, msg="Importing tomoMasks failed.")
        self.assertEqual(tomoMaskSet.getDimensions(), (464, 464, 250))
        self.assertEqual(tomoMaskSet.getSamplingRate(), 2 * self.samplingRate)
        self.assertFalse(protImportTomomasks.warnMsg.get())

    def testImportTomoMasksAllGood(self):
        protImportTomomasks = self._importTomoMasks(filesPath=self.ds.getFile(DataSetEmd10439.tomoMaskAnnotated.value),
                                                    inputTomos=self.inTomoSetBinned)
        self._checkResults(protImportTomomasks)

    def testImportTomoMasksAllGood_pattern(self):
        protImportTomomasks = self._importTomoMasks(filesPath=self.ds.getFile(DataSetEmd10439.tomomasksAnnotatedDir.value),
                                                    filesPattern='{TS}_materials.mrc',
                                                    inputTomos=self.inTomoSetBinned)
        self._checkResults(protImportTomomasks)


    # The tomogram and the tomomask have different sizes
    def testImportTomoMasksDiffSize(self):
        with self.assertRaises(Exception) as eType:
            self._importTomoMasks(filesPath=self.ds.getFile(DataSetEmd10439.tomoMaskAnnotated.value),
                                  inputTomos=self.inTomoSet)
            self.assertEqual(str(eType.exception), ERR_NO_TOMOMASKS_GEN)

    # Tomomask and tomogram doesn't correspond
    def testImportTomoMasksNoneMatch(self):
        with self.assertRaises(Exception) as eType:
            self._importTomoMasks(filesPath=self.ds.getFile(DataSetEmd10439.tomoMaskAnnotated.value),
                                  inputTomos=self.inNotMatchingTomoSet)
            self.assertEqual(str(eType.exception), ERR_NON_MATCHING_TOMOS)


class TestConvert(BaseTest):

    def test_list_order(self):
        myList = [20, 40, 50, 10, 30]
        orders = getOrderFromList(myList)
        self.assertEqual(orders, [2, 4, 5, 1, 3])

        myList = [20, 40, 50, 10, 30]
        orders = getOrderFromList(myList)
        self.assertEqual(orders, [2, 4, 5, 1, 3])

        # Simulate dose list in symmetric dose scheme
        myList = [24, 21, 12, 9, 3, 6, 15, 18, 27]
        orders = getOrderFromList(myList)
        self.assertEqual(orders, [8, 7, 4, 3, 1, 2, 5, 6, 9])


class TestImportTomoCtf(TestBaseCentralizedLayer):
    ds = None
    nTiltSeries = 2

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)

    @classmethod
    def _runImportTs(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
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
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protTsImport)
        tsImported = getattr(protTsImport, protTsImport.OUTPUT_NAME, None)
        return tsImported

    def test_import_ctf_aretomo_01(self):
        inTsSet = self._runImportTs()

        print(magentaStr("\n==> Testing tomo CTF import for AreTomo:"
                         "\n\t- TS without excluded views"))
        prot = self.newProtocol(ProtImportTsCTF,
                                importFrom=ImportChoice.ARETOMO.value,
                                filesPath=self.ds.getFile(DataSetRe4STATuto.aretomoCtfFilesPath.value),
                                filesPattern='*.txt',
                                inputSetOfTiltSeries=inTsSet)
        self.launchProtocol(prot)
        self._checkCTFs(getattr(prot, prot._possibleOutputs.CTFs.name, None))

    def test_import_ctf_aretomo_02(self):
        inTsSet = self._runImportTs()

        print(magentaStr("\n==> Testing tomo CTF import for AreTomo:"
                         "\n\t- TS with some excluded views"))

        # Exclude some views in each TS
        exludedViews = {'TS_03': [0, 1, 2, 34, 35, 36, 37, 38, 39],
                        'TS_54': [0, 1, 2, 39, 40]}
        nTiltImages = {'TS_03': 40,
                       'TS_54': 41}
        # To do this properly, the sets must be iterated by direct indexing of the elements, as the iterators
        # stop iterating after having executed a mapper update operation
        for n in range(self.nTiltSeries):
            ts = inTsSet[n + 1]
            tsId = ts.getTsId()
            exclusionList = exludedViews[tsId]
            for i in range(nTiltImages[tsId]):
                ti = ts[i + 1]
                if i in exclusionList:
                    ti.setEnabled(False)
                    ts.update(ti)
            ts.write()
            inTsSet.update(ts)

        # Launch the protocol and check the results
        prot = self.newProtocol(ProtImportTsCTF,
                                importFrom=ImportChoice.ARETOMO.value,
                                filesPath=self.ds.getFile(DataSetRe4STATuto.aretomoCtfFilesPath.value),
                                filesPattern='*.txt',
                                inputSetOfTiltSeries=inTsSet)
        prot.setObjLabel('Import CTFs, some views excluded')
        self.launchProtocol(prot)
        self._checkCTFs(getattr(prot, prot._possibleOutputs.CTFs.name), excludedViewsDict=exludedViews)

    def _checkCTFs(self, ctfSet, excludedViewsDict=None):
        self.checkCTFs(ctfSet,
                       expectedSetSize=self.nTiltSeries,
                       excludedViewsDict=excludedViewsDict,
                       expectedPsdFile=True)
