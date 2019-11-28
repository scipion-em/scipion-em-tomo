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
from shutil import copyfile

import pyworkflow as pw
from pyworkflow.tests import BaseTest, setupTestOutput, setupTestProject
from pyworkflow.em import Domain, CTFModel

from tomo.tests import DataSet
from tomo.objects import SetOfTiltSeriesM, SetOfTiltSeries
import tomo.protocols


class TestTomoBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('tomo')
        cls.getFile = cls.dataset.getFile

    def test_plugin(self):
        # Really stupid test to check that tomo plugin is defined
        tomo = Domain.getPlugin('tomo')

        self.assertFalse(tomo is None)
        self.assertTrue(hasattr(tomo, 'Plugin'))

        # Check that defined objects here are found
        objects = Domain.getObjects()

        expected = ['TiltImage', 'TiltSeries', 'SetOfTiltSeries',
                    'TiltImageM', 'TiltSeriesM', 'SetOfTiltSeriesM']
        for e in expected:
            self.assertTrue(e in objects,
                            "%s should be in Domain.getObjects" % e)

    def _create_tiltseries(self, tiltSeriesClass):
        setFn = self.getOutputPath('%s.sqlite' % tiltSeriesClass.__name__)
        pw.utils.cleanPath(setFn)

        testSet = tiltSeriesClass(filename=setFn)

        for tsId in ['TS01', 'TS02', 'TS03']:
            tsObj = testSet.ITEM_TYPE(tsId=tsId)
            testSet.append(tsObj)
            for i, a in enumerate(range(-60, 60, 5)):
                tsFn = '%s_%02d_%s.mrc' % (tsId, i, a)
                tim = tsObj.ITEM_TYPE(location=tsFn,
                                      acquisitionOrder=i,
                                      tiltAngle=a)
                tsObj.append(tim)

        testSet.write()
        testSet.close()

        testSet2 = SetOfTiltSeriesM(filename=setFn)
        self.assertEqual(testSet2.getSize(), 3)
        for tsObj in testSet2:
            self.assertEqual(tsObj.getSize(), 24)

        testSet2.close()

        return testSet

    def test_create_tiltseries(self):
        self._create_tiltseries(SetOfTiltSeries)

    def test_create_tiltseriesM(self):
        self._create_tiltseries(SetOfTiltSeriesM)

    def test_copyItems(self):
        testSet = self._create_tiltseries(SetOfTiltSeries)

        fn = testSet.getFileName()
        print("Input sqlite: %s" % fn)

        fn2 = fn.replace('.sqlite', '-withctf.sqlite')
        print("Writing new set to: %s" % fn2)

        def _updateCtf(j, ts, ti, tsOut, tiOut):
            tiOut.setCTF(CTFModel(defocusU=j*1000, defocusAngle=j/2.0))

        testSet2 = SetOfTiltSeries(filename=fn2)
        testSet2.copyItems(testSet, updateTiCallback=_updateCtf)

        testSet2.write()
        testSet2.close()


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
            filesPath=os.path.join(self.dataPath, 'data', 'frames'),
            filesPattern=filesPattern,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=1.35,
            doseInitial=0,
            dosePerFrame=0.3)
        self.launchProtocol(protImport)
        return protImport

    def test_importTiltSeriesM(self):
        protImport = self._runImportTiltSeriesM()
        output = getattr(protImport, 'outputTiltSeriesM', None)
        self.assertFalse(output is None)

        return protImport

    def test_motioncorTiltSeriesM(self):
        protImport = self.test_importTiltSeriesM()

        # --------- Motion correction with motioncor2 for Tilt-series ------
        import motioncorr.protocols
        protMc = self.newProtocol(
            motioncorr.protocols.ProtTsMotionCorr,
            binFactor=2.0
        )

        protMc.inputTiltSeriesM.set(protImport.outputTiltSeriesM)
        self.launchProtocol(protMc)


class TestTomoImportTs(BaseTest):
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
                samplingRate=1.35,
                doseInitial=0,
                dosePerFrame=0.3)
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
                dosePerFrame=0.3)
            self.launchProtocol(protImport)
            return protImport

        def test_importTiltSeriesM(self):
            protImport = self._runImportTiltSeriesM()
            output = getattr(protImport, 'outputTiltSeriesM', None)
            self.assertFalse(output is None)
            self.assertEqual(output.getSize(), 2)

            return protImport

        def test_importTiltSeries(self):
            protImport = self._runImportTiltSeries()
            output = getattr(protImport, 'outputTiltSeries', None)
            self.assertFalse(output is None)
            self.assertEqual(output.getSize(), 2)


class TestTomoImportTomograms(BaseTest):
    """This class check if the protocol to import tomograms works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo2')

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

    def test_importTomograms(self):
        protImport = self._runImportTomograms()
        output = getattr(protImport, 'outputTomograms', None)
        self.assertIsNotNone(output,
                             "There was a problem with Import Tomograms protocol")

        for tomo in protImport.outputTomograms.iterItems():

            self.assertTrue(tomo.getXDim() == 1024,
                                "There was a problem with Import Tomograms protocol")
            self.assertIsNotNone(tomo.getYDim() == 1024,
                                "There was a problem with Import Tomograms protocol")

            self.assertTrue(tomo.getAcquisition().getAngleMax() == 40, "There was a problem with the aquisition angle max")
            self.assertTrue(tomo.getAcquisition().getAngleMin() == -40, "There was a problem with the aquisition angle min")

            break



class TestTomoImportSubTomograms(BaseTest):
     """ This class check if the protocol to import sub tomograms works
     properly."""
     @classmethod
     def setUpClass(cls):
         setupTestProject(cls)
         cls.dataset = DataSet.getDataSet('tomo-em')
         cls.tomogram = cls.dataset.getFile('tomo1')
         cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
         cls.path = cls.dataset.getPath()

     def _runImportSubTomograms(self):

        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)
        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                 auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                 filesPath= self.coords3D,
                                 importTomograms=protImportTomogram.outputTomograms,
                                 filesPattern='', boxSize=32,
                                 samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)

        protImport = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.tomogram,
                                      filesPattern='',
                                      samplingRate=1.35,
                                      importCoordinates=protImportCoordinates3d.outputCoordinates)
        self.launchProtocol(protImport)
        return protImport

     def _runImportSubTomograms2(self):

        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              filesPattern='',
                                              samplingRate=5)
        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                 auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                 filesPath= self.coords3D,
                                 importTomograms=protImportTomogram.outputTomograms,
                                 filesPattern='', boxSize=32,
                                 samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)

        protImport = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.path,
                                      filesPattern='*.em',
                                      samplingRate=1.35,
                                      importCoordinates=protImportCoordinates3d.outputCoordinates)
        self.launchProtocol(protImport)
        return protImport

     def test_import_sub_tomograms(self):
         protImport = self._runImportSubTomograms()
         output = getattr(protImport, 'outputSubTomograms', None)
         self.assertTrue(output.getSamplingRate() == 1.35)
         self.assertTrue(output.getFirstItem().getSamplingRate() == 1.35)
         self.assertTrue(output.getDim()[0] == 1024)
         self.assertTrue(output.getDim()[1] == 1024)
         self.assertTrue(output.getDim()[2] == 512)
         self.assertTrue(output.getFirstItem().getCoordinate3D().getX() == 314)
         self.assertTrue(output.getFirstItem().getCoordinate3D().getY() == 350)
         self.assertTrue(output.getFirstItem().getCoordinate3D().getZ() == 256)
         self.assertIsNotNone(output,
                             "There was a problem with Import SubTomograms protocol")

         protImport2 = self._runImportSubTomograms2()
         output2 = getattr(protImport2, 'outputSubTomograms', None)
         self.assertIsNotNone(output2,
                              "There was a problem with Import SubTomograms protocol")
         self.assertTrue(output2.getSamplingRate() == 1.35)
         self.assertTrue(output2.getDim()[0] == 1024)
         self.assertTrue(output2.getDim()[1] == 1024)
         self.assertTrue(output2.getDim()[2] == 512)
         for i, subtomo in enumerate(output2.iterItems()):
             if i == 1:
                 self.assertTrue(subtomo.getCoordinate3D().getX() == 174)
                 self.assertTrue(subtomo.getCoordinate3D().getY() == 172)
                 self.assertTrue(subtomo.getCoordinate3D().getZ() == 256)
             if i == 0:
                 self.assertTrue(subtomo.getCoordinate3D().getX() == 314)
                 self.assertTrue(subtomo.getCoordinate3D().getY() == 350)
                 self.assertTrue(subtomo.getCoordinate3D().getZ() == 256)

         return output2


class TestTomoImportSetOfCoordinates3D(BaseTest):
    """This class check if the protocol to import set of coordinates 3d works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')

    def _runTomoImportSetOfCoordinates(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                 filesPath=self.tomogram,
                                 samplingRate=5)

        self.launchProtocol(protImportTomogram)

        output = getattr(protImportTomogram, 'outputTomograms', None)

        self.assertIsNotNone(output,
                             "There was a problem with tomogram output")

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                 auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                 filesPath= self.coords3D,
                                 importTomograms=protImportTomogram.outputTomograms,
                                 filesPattern='', boxSize=32,
                                 samplingRate=5.5)
        self.launchProtocol(protImportCoordinates3d)

        return protImportCoordinates3d

    def test_import_set_of_coordinates_3D(self):
        protCoordinates = self._runTomoImportSetOfCoordinates()
        output = getattr(protCoordinates, 'outputCoordinates', None)
        self.assertTrue(output,
                             "There was a problem with coordinates 3d output")
        self.assertTrue(output.getSize() == 5)
        self.assertTrue(output.getBoxSize() == 32)
        self.assertTrue(output.getSamplingRate() == 5.5)
        return output


class TestTomoPreprocessing(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataPath = os.environ.get('SCIPION_TOMO_EMPIAR10164', '')

        if not os.path.exists(cls.dataPath):
            raise Exception("Can not run tomo tests, "
                            "SCIPION_TOMO_EMPIAR10164 variable not defined. ")

    def _runImportTiltSeriesM(self, filesPattern='{TS}_{TO}_{TA}.mrc'):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTsMovies,
            filesPath=os.path.join(self.dataPath, 'data', 'frames'),
            filesPattern=filesPattern,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=1.35,
            doseInitial=0,
            dosePerFrame=0.3)
        self.launchProtocol(protImport)
        return protImport

    def test_preprocess1(self):
        """ Run the basic preprocessing pipeline for just one TS. """
        protImport = self._runImportTiltSeriesM(
            filesPattern='{TS}1_{TO}_{TA}.mrc')
        output = getattr(protImport, 'outputTiltSeriesM', None)
        self.assertFalse(output is None)
        self.assertEqual(protImport.outputTiltSeriesM.getSize(), 1)

        gpuList = os.environ.get('SCIPION_TEST_GPULIST', '0')
        threads = len(gpuList.split()) + 1

        # --------- Motion correction with motioncor2 for Tilt-series ------
        import motioncorr.protocols
        protMc = self.newProtocol(
            motioncorr.protocols.ProtTsMotionCorr,
            inputTiltSeriesM=protImport.outputTiltSeriesM,
            binFactor=2.0,
            gpuList=gpuList,
            numberOfThreads=threads
        )
        self.launchProtocol(protMc)

        # --------- CTF estimation with Gctf ------
        import gctf.protocols
        protGctf = self.newProtocol(
            gctf.protocols.ProtTsGctf,
            inputTiltSeries=protMc.outputTiltSeries,
            gpuList=gpuList,
            numberOfThreads=threads,
        )
        self.launchProtocol(protGctf)

        # -------- Basic alignment and reconstruction with IMOD ------
        import imod.protocols
        protImodAuto = self.newProtocol(
            imod.protocols.ProtImodAuto3D,
            inputTiltSeries=protGctf.outputTiltSeries,
            excludeList=1,
            rotationAngle=90,
            zWidth=400,
            useRaptor=True,
            markersDiameter=20,
            markersNumber=20
        )
        self.launchProtocol(protImodAuto)


if __name__ == 'main':
    pass
