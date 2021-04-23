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
import random
from os.path import join

import pyworkflow.utils as pwutils
import pwem.protocols as emprot
from pyworkflow.tests import BaseTest, setupTestOutput, setupTestProject
from pwem import Domain
from pwem.objects import CTFModel
from tomo.protocols.protocol_ts_import import MDoc

from . import DataSet
from ..objects import SetOfTiltSeriesM, SetOfTiltSeries
from ..utils import existsPlugin
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
        pwutils.cleanPath(setFn)

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

    def test_importTiltSeriesM(self):
        protImport = self._runImportTiltSeriesM()
        output = getattr(protImport, 'outputTiltSeriesM', None)
        self.assertFalse(output is None)

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


class TestTomoImportTsFromMdoc(BaseTest):

    VOLTAGE = 'voltage'
    MAGNIFICATION = 'magnification'
    SPH_ABERRATION = 'sphAberration'
    AMP_CONTRAST = 'ampContrast'
    DOSE_PER_FRAME = 'dosePerFrame'
    PIXEL_SIZE = 'pixelSize'
    SIZE = 'size'
    FILENAME_LIST = 'filenames'
    ACCUM_DOSE_LIST = 'doses'
    ANGLE_LIST = 'angles'

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

    def _getListOfFileNames(self, tomoNum, isTsMovie):
        fpath = self.ts31Dir if tomoNum == 31 else self.ts10Dir
        angles = [0, 3, -3, -6, 6]
        if isTsMovie:
            # Each slice will be a stack of frames, which means a different file
            pattern = 'SKvesicles_Pertuzumab_tomo%2.0f_%03d_%1.1f.mrc'

            return [join(self.dataset.getPath(), fpath, pattern % (tomoNum, i, angle))
                    for i, angle in enumerate(angles)]
        else:
            # Each slice will be referred to the same stack file, the one which represents the tilt series
            return [join(self.dataset.getPath(), fpath, 'stack%2.0f.mrcs' % tomoNum)] * len(angles)

    @classmethod
    def _genTestDict(cls, **kwargs):
        return {
            cls.VOLTAGE: kwargs.get(cls.VOLTAGE, None),
            cls.MAGNIFICATION: kwargs.get(cls.MAGNIFICATION, None),
            cls.SPH_ABERRATION: kwargs.get(cls.SPH_ABERRATION, None),
            cls.AMP_CONTRAST: kwargs.get(cls.AMP_CONTRAST, None),
            cls.DOSE_PER_FRAME: kwargs.get(cls.DOSE_PER_FRAME, None),
            cls.PIXEL_SIZE: kwargs.get(cls.PIXEL_SIZE, None),
            cls.SIZE: kwargs.get(cls.SIZE, None),
            cls.FILENAME_LIST: kwargs.get(cls.FILENAME_LIST, None),
            cls.ACCUM_DOSE_LIST: kwargs.get(cls.ACCUM_DOSE_LIST, None),
            cls.ANGLE_LIST: kwargs.get(cls.ANGLE_LIST, None)
        }

    def _genTestData(self, isTsMovie):
        return {
            'stack31': self._genTestDict(
                voltage=300,
                magnification=53000,
                sphAberration=2.7,
                ampContrast=0.1,
                dosePerFrame=0.0891,
                pixelSize=self.sRate,
                size=5,
                filenames=self._getListOfFileNames(tomoNum=31, isTsMovie=isTsMovie),
                doses=[0.0424, 0.2121, 0.2546, 0.3394, 0.4455],
                angles=[0.0036, 2.9683, -3.0250, -6.0251, 5.9684]
            ),
            'stack10': self._genTestDict(
                voltage=300,
                magnification=53000,
                sphAberration=2.7,
                ampContrast=0.1,
                dosePerFrame=0.0806,
                pixelSize=self.sRate,
                size=5,
                filenames=self._getListOfFileNames(tomoNum=10, isTsMovie=isTsMovie),
                doses=[0.1909, 0.2333, 0.2546, 0.2970, 0.4030],
                angles=[0.0026, 2.9708, -3.0255, -6.0241, 5.9684]
            )
        }

    def _runImportTiltSeries(self, isTsMovie=False):
        prot = tomo.protocols.ProtImportTsMovies if isTsMovie else tomo.protocols.ProtImportTs
        protImport = self.newProtocol(
            prot,
            filesPath=self.parentDir,
            filesPattern=self.pattern)
        self.launchProtocol(protImport)
        return protImport

    def _checkResults(self, outputSet, isTsMovie, dimensions):
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
            # Check angles and accumulated dose per angular acquisition (tilt series image)
            filesList = testDataDict[self.FILENAME_LIST]
            accumDoseList = testDataDict[self.ACCUM_DOSE_LIST]
            angleList = testDataDict[self.ANGLE_LIST]
            for i, tiM in enumerate(tsM):
                self.assertEqual(tiM.getFileName(), filesList[i])
                self.assertAlmostEqual(tiM.getAcquisition().getDosePerFrame(), accumDoseList[i], delta=0.0001)
                self.assertAlmostEqual(tiM.getTiltAngle(), angleList[i], delta=0.0001)

    def test_importTiltSeriesM(self):
        isTsMovie = True
        # Run protocol
        protImport = self._runImportTiltSeries(isTsMovie=isTsMovie)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeriesM', None)
        self._checkResults(outputSet, isTsMovie, dimensions=(1152, 1152, 6))  # ts and tsM has different dimensions
        # because they have been downsampled separately in order to get a lighter test dataset

    def test_importTiltSeries(self):
        isTsMovie = False
        # Run protocol
        protImport = self._runImportTiltSeries(isTsMovie=isTsMovie)
        # Check results
        outputSet = getattr(protImport, 'outputTiltSeries', None)
        self._checkResults(outputSet, isTsMovie, dimensions=(1440, 1024, 1))  # ts and tsM have different dimensions
        # because they have been downsampled separately in order to get a lighter test dataset

    def test_mdocFileChecker(self):
        # JORGE
        import os
        fname = "/home/jjimenez/Desktop/test_JJ.txt"
        if os.path.exists(fname):
            os.remove(fname)
        fjj = open(fname, "a+")
        fjj.write('JORGE--------->onDebugMode PID {}'.format(os.getpid()))
        fjj.close()
        print('JORGE--------->onDebugMode PID {}'.format(os.getpid()))
        import time
        time.sleep(10)
        # JORGE_END
        mdocList = glob.glob(join(self.dataset.getFile('empiarMdocDir'), '*.mdoc'))
        self.assertEqual(len(mdocList), 10)

        errorMsgList = []
        for mdoc in mdocList:
            mdocObj = MDoc(mdoc)
            errorMsgList.append(mdocObj.read(isImportingTsMovies=True, ignoreFilesValidation=True))

        a = 1


class TestTomoSubSetsTs(BaseTest):
    """This class check if the protocol to import Tilt series, create subsets,
     join subset and create a intersection of subsets of tomogram
     works properly."""

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

    def split(self, tsSet, n, randomize):
        """Return a run split protocol over input set tomoSet."""

        pSplit = self.proj.newProtocol(emprot.ProtSplitSet, numberOfSets=n)
        pSplit.inputSet.set(tsSet)
        pSplit.randomize.set(randomize)
        self.proj.launchProtocol(pSplit, wait=True)
        return pSplit

    @staticmethod
    def outputs(p):
        """Iterator over all the elements in the outputs of protocol p."""
        for key, output in p.iterOutputAttributes():
            yield output

    def test_createSubSetOfTs(self):
        """This class check if the protocol to import tilt series, create subsets,
             join subset and create a intersection of subsets of tomogram
             works properly."""
        print(pwutils.magentaStr("\n==> Importing a Set of TiltSeriesM..."))

        protImport = self._runImportTiltSeriesM()
        self.assertFalse(protImport.outputTiltSeriesM is None)
        self.assertEqual(protImport.outputTiltSeriesM.getSize(), 2)

        print(pwutils.magentaStr("\n==> Importing a Set of TiltSeries..."))
        protImport2 = self._runImportTiltSeries()
        self.assertFalse(protImport2.outputTiltSeries is None)
        self.assertEqual(protImport2.outputTiltSeries.getSize(), 2)

        # Create a subset 1 of TiltSeriesM
        print(pwutils.magentaStr("\n==> Create the subset 1 of TiltSeriesM"))
        tsSubset1 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 1',
                                       chooseAtRandom=True,
                                       nElements=1)

        tsSubset1.inputFullSet.set(protImport.outputTiltSeriesM)
        self.launchProtocol(tsSubset1)

        # Create a subset subset 2 of TiltSeriesM
        print(pwutils.magentaStr("\n==> Create the subset 2 of TiltSeriesM"))
        tsSubset2 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 2',
                                       chooseAtRandom=True,
                                       nElements=1)

        tsSubset2.inputFullSet.set(protImport.outputTiltSeriesM)
        self.launchProtocol(tsSubset2)

        # Create a subset subset 3 of TiltSeriesM
        print(pwutils.magentaStr("\n==> Create the subset 3 of TiltSeriesM"))
        tsSubset3 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 3',
                                       chooseAtRandom=True,
                                       nElements=2)

        tsSubset3.inputFullSet.set(protImport.outputTiltSeriesM)
        self.launchProtocol(tsSubset3)

        # create merge protocol
        print(pwutils.magentaStr("\n==> Join subsets 1 and 2 "))
        joinTsSets = self.newProtocol(emprot.ProtUnionSet,
                                       objLabel='join TiltSeriesM sets',
                                       ignoreExtraAttributes=True)
        joinTsSets.inputSets.append(tsSubset1.outputTiltSeriesM)
        joinTsSets.inputSets.append(tsSubset2.outputTiltSeriesM)
        self.launchProtocol(joinTsSets, wait=True)

        print(pwutils.magentaStr("\n==> Split the subset 3 "))

        tomoSplit1 = self.split(tsSubset3.outputTiltSeriesM, n=2,
                                randomize=True)
        tomoSplit2 = self.split(tsSubset3.outputTiltSeriesM, n=1,
                                randomize=True)

        setFull = random.choice(list(self.outputs(tomoSplit1)))
        setSub = random.choice(list(self.outputs(tomoSplit2)))

        # Launch intersection subset
        print(pwutils.magentaStr(
            "\n==> Intersection subsets of TiltSeriesM"))
        label = '%s - %s,%s ' % (
        tsSubset3.outputTiltSeriesM.getClassName(), 1, 2)
        tsSubset = self.newProtocol(emprot.ProtSubSet)
        tsSubset.setObjLabel(label + 'intersection')
        tsSubset.inputFullSet.set(setFull)
        tsSubset.inputSubSet.set(setSub)
        self.launchProtocol(tsSubset)

        setFullIds = setFull.getIdSet()
        setSubIds = setSub.getIdSet()

        # Check intersection
        outputs = [o for o in self.outputs(tsSubset)]
        if outputs:
            outputTs = outputs[0]

            # Check properties
            self.assertTrue(
                tsSubset2.outputTiltSeriesM.equalAttributes(outputTs,
                                                            ignore=[
                                                                '_mapperPath',
                                                                '_size'],
                                                            verbose=True),
                "Intersection subset attributes are wrong")

            for elem in outputTs:
                self.assertTrue(elem.getObjId() in setFullIds)
                self.assertTrue(elem.getObjId() in setSubIds,
                                'object id %s not in set: %s'
                                % (elem.getObjId(), setSubIds))


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


class TestTomoSubSetsTomograms(BaseTest):
    """This class check if the protocol to import tomograms, create subsets,
     join subset and create a intersection of subsets of tomogram
     works properly."""
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

    def split(self, tomoSet, n, randomize):
        """Return a run split protocol over input set tomoSet."""

        pSplit = self.proj.newProtocol(emprot.ProtSplitSet, numberOfSets=n)
        pSplit.inputSet.set(tomoSet)
        pSplit.randomize.set(randomize)
        self.proj.launchProtocol(pSplit, wait=True)
        return pSplit

    @staticmethod
    def outputs(p):
        """Iterator over all the elements in the outputs of protocol p."""
        for key, output in p.iterOutputAttributes():
            yield output

    def test_createSubSetOfTomogram(self):
        print(pwutils.magentaStr("\n==> Importing a Set of Tomograms..."))
        protImport = self._runImportTomograms()

        # Check importing
        for tomo in protImport.outputTomograms.iterItems():
            self.assertTrue(tomo.getXDim() == 1024,
                            "There was a problem with Import Tomograms protocol")
            self.assertIsNotNone(tomo.getYDim() == 1024,
                                 "There was a problem with Import Tomograms "
                                 "protocol")

            self.assertTrue(tomo.getAcquisition().getAngleMax() == 40,
                            "There was a problem with the aquisition angle max")
            self.assertTrue(tomo.getAcquisition().getAngleMin() == -40,
                            "There was a problem with the aquisition angle min")

        # Create a subset with 1 tomograms
        print(pwutils.magentaStr("\n==> Create the subset 1 of Tomograms"))
        tomoSubset1 = self.newProtocol(emprot.ProtSubSet,
                                      objLabel='subset 1',
                                      chooseAtRandom=True,
                                      nElements=1)

        tomoSubset1.inputFullSet.set(protImport.outputTomograms)
        self.launchProtocol(tomoSubset1)

        # Create a subset with 2 tomograms
        print(pwutils.magentaStr("\n==> Create the subset 2 of Tomograms"))
        tomoSubset2 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 2',
                                       chooseAtRandom=True,
                                       nElements=1)

        tomoSubset2.inputFullSet.set(protImport.outputTomograms)
        self.launchProtocol(tomoSubset2)

        # Create a subset with 3 tomograms
        print(pwutils.magentaStr("\n==> Create the subset 3 of Tomograms"))
        tomoSubset3 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 3',
                                       chooseAtRandom=True,
                                       nElements=2)

        tomoSubset3.inputFullSet.set(protImport.outputTomograms)
        self.launchProtocol(tomoSubset3)

        # create merge protocol
        print(pwutils.magentaStr("\n==> Join subsets 1 and 2 "))
        joinTomoSet = self.newProtocol(emprot.ProtUnionSet,
                                   objLabel='join Tomograms sets',
                                   ignoreExtraAttributes=True)
        joinTomoSet.inputSets.append(tomoSubset1.outputTomograms)
        joinTomoSet.inputSets.append(tomoSubset2.outputTomograms)
        self.launchProtocol(joinTomoSet, wait=True)

        print(pwutils.magentaStr("\n==> Split the subset 3 "))

        tomoSplit1 = self.split(tomoSubset3.outputTomograms, n=2, randomize=True)
        tomoSplit2 = self.split(tomoSubset3.outputTomograms, n=1, randomize=True)

        setFull = random.choice(list(self.outputs(tomoSplit1)))
        setSub = random.choice(list(self.outputs(tomoSplit2)))

        # Launch intersection subset
        print(pwutils.magentaStr(
            "\n==> Intersection subsets of Tomograms"))
        label = '%s - %s,%s ' % (tomoSubset3.outputTomograms.getClassName(), 1, 2)
        tomoSubset = self.newProtocol(emprot.ProtSubSet)
        tomoSubset.setObjLabel(label + 'intersection')
        tomoSubset.inputFullSet.set(setFull)
        tomoSubset.inputSubSet.set(setSub)
        self.launchProtocol(tomoSubset)

        setFullIds = setFull.getIdSet()
        setSubIds = setSub.getIdSet()

        # Check intersection
        outputs = [o for o in self.outputs(tomoSubset)]
        if outputs:
            outputTomo = outputs[0]

            # Check properties
            self.assertTrue(tomoSubset2.outputTomograms.equalAttributes(outputTomo,
                                                 ignore=['_mapperPath',
                                                         '_size'],
                                                 verbose=True),
                            "Intersection subset attributes are wrong")

            for elem in outputTomo:
                self.assertTrue(elem.getObjId() in setFullIds)
                self.assertTrue(elem.getObjId() in setSubIds,
                                'object id %s not in set: %s'
                                % (elem.getObjId(), setSubIds))


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
                                      samplingRate=1.35)
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


class TestTomoSubSetsSubTomograms(BaseTest):
    """This class check if the protocol to import subtomograms, create subsets,
     join subset and create a intersection of subsets of tomogram
     works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
        cls.table = cls.dataset.getFile('initial.tbl')
        cls.path = cls.dataset.getPath()

    def _runImportSubTomograms(self):

        protImportSubTomogram = self.newProtocol(
            tomo.protocols.ProtImportTomograms,
            filesPath=self.tomogram,
            samplingRate=5)
        self.launchProtocol(protImportSubTomogram)

        protImportCoordinates3d = self.newProtocol(
            tomo.protocols.ProtImportCoordinates3D,
            auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_EMAN,
            filesPath=self.coords3D,
            importTomograms=protImportSubTomogram.outputTomograms,
            filesPattern='', boxSize=32,
            samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)

        protImport = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.path,
                                      filesPattern='base*.hdf',
                                      samplingRate=1.35,
                                      importCoordinates=protImportCoordinates3d.outputCoordinates)
        self.launchProtocol(protImport)
        return protImport

    def split(self, subtomoSet, n, randomize):
        """Return a run split protocol over input set tomoSet."""

        pSplit = self.proj.newProtocol(emprot.ProtSplitSet, numberOfSets=n)
        pSplit.inputSet.set(subtomoSet)
        pSplit.randomize.set(randomize)
        self.proj.launchProtocol(pSplit, wait=True)
        return pSplit

    @staticmethod
    def outputs(p):
        """Iterator over all the elements in the outputs of protocol p."""
        for key, output in p.iterOutputAttributes():
            yield output

    def test_createSubSetOfTomogram(self):
        print(pwutils.magentaStr("\n==> Importing a Set of protImportSubTomogram..."))
        protImport = self._runImportSubTomograms()

        # Check importing
        self.assertTrue(protImport.outputSubTomograms.getSamplingRate() == 1.35)
        self.assertTrue(protImport.outputSubTomograms.getFirstItem().getSamplingRate() == 1.35)
        self.assertTrue(protImport.outputSubTomograms.getDim()[0] == 32)
        self.assertTrue(protImport.outputSubTomograms.getDim()[1] == 32)
        self.assertTrue(protImport.outputSubTomograms.getDim()[2] == 32)
        self.assertSetSize(protImport.outputSubTomograms, 4,
                             "There was a problem with Import SubTomograms protocol")

        # Create a subset with 1 tomograms
        print(pwutils.magentaStr("\n==> Create the subset 1 of SubTomograms"))
        tomoSubset1 = self.newProtocol(emprot.ProtSubSet,
                                      objLabel='subset 1',
                                      chooseAtRandom=True,
                                      nElements=2)

        tomoSubset1.inputFullSet.set(protImport.outputSubTomograms)
        self.launchProtocol(tomoSubset1)

        # Create a subset with 2 tomograms
        print(pwutils.magentaStr("\n==> Create the subset 2 of SubTomograms"))
        tomoSubset2 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 2',
                                       chooseAtRandom=True,
                                       nElements=2)

        tomoSubset2.inputFullSet.set(protImport.outputSubTomograms)
        self.launchProtocol(tomoSubset2)

        # Create a subset with 3 tomograms
        print(pwutils.magentaStr("\n==> Create the subset 3 of SubTomograms"))
        tomoSubset3 = self.newProtocol(emprot.ProtSubSet,
                                       objLabel='subset 3',
                                       chooseAtRandom=True,
                                       nElements=2)

        tomoSubset3.inputFullSet.set(protImport.outputSubTomograms)
        self.launchProtocol(tomoSubset3)

        # create merge protocol
        print(pwutils.magentaStr("\n==> Join subsets 1 and 2 "))
        joinTomoSet = self.newProtocol(emprot.ProtUnionSet,
                                   objLabel='join SubTomograms sets',
                                   ignoreExtraAttributes=True)
        joinTomoSet.inputSets.append(tomoSubset1.outputSubTomograms)
        joinTomoSet.inputSets.append(tomoSubset2.outputSubTomograms)
        self.launchProtocol(joinTomoSet, wait=True)

        print(pwutils.magentaStr("\n==> Split the subset 3 "))

        tomoSplit1 = self.split(tomoSubset3.outputSubTomograms, n=2, randomize=True)
        tomoSplit2 = self.split(tomoSubset3.outputSubTomograms, n=1, randomize=True)

        setFull = random.choice(list(self.outputs(tomoSplit1)))
        setSub = random.choice(list(self.outputs(tomoSplit2)))

        # Launch intersection subset
        print(pwutils.magentaStr(
            "\n==> Intersection subsets of SubTomograms"))
        label = '%s - %s,%s ' % (tomoSubset3.outputSubTomograms.getClassName(), 1, 2)
        tomoSubset = self.newProtocol(emprot.ProtSubSet)
        tomoSubset.setObjLabel(label + 'intersection')
        tomoSubset.inputFullSet.set(setFull)
        tomoSubset.inputSubSet.set(setSub)
        self.launchProtocol(tomoSubset)

        setFullIds = setFull.getIdSet()
        setSubIds = setSub.getIdSet()

        # Check intersection
        outputs = [o for o in self.outputs(tomoSubset)]
        if outputs:
            outputTomo = outputs[0]

            # Check properties
            self.assertTrue(tomoSubset2.outputSubTomograms.equalAttributes(outputTomo,
                                                 ignore=['_mapperPath',
                                                         '_size'],
                                                 verbose=True),
                            "Intersection subset attributes are wrong")

            for elem in outputTomo:
                self.assertTrue(elem.getObjId() in setFullIds)
                self.assertTrue(elem.getObjId() in setSubIds,
                                'object id %s not in set: %s'
                                % (elem.getObjId(), setSubIds))

class TestTomoImportSetOfCoordinates3D(BaseTest):
    """This class check if the protocol to import set of coordinates 3d works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')

    def _runTomoImportSetOfCoordinates(self, pattern, program, ext):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        output = getattr(protImportTomogram, 'outputTomograms', None)

        self.assertIsNotNone(output,
                             "There was a problem with tomogram output")

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                 objLabel='Import from %s - %s' %(program, ext),
                                 auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_AUTO,
                                 filesPath=self.dataset.getPath(),
                                 importTomograms=protImportTomogram.outputTomograms,
                                 filesPattern=pattern, boxSize=32,
                                 samplingRate=5.5)
        self.launchProtocol(protImportCoordinates3d)

        return protImportCoordinates3d

    def test_import_set_of_coordinates_3D(self):
        outputs = []
        if existsPlugin('emantomo'):
            protCoordinates = self._runTomoImportSetOfCoordinates('*.json', 'EMAN', 'JSON')
            output = getattr(protCoordinates, 'outputCoordinates', None)
            self.assertTrue(output,
                                 "There was a problem with coordinates 3d output")
            self.assertTrue(output.getSize() == 19)
            self.assertTrue(output.getBoxSize() == 32)
            self.assertTrue(output.getSamplingRate() == 5.5)
            outputs.append(output)

        protCoordinates2 = self._runTomoImportSetOfCoordinates('*.txt', 'TOMO', 'TXT')
        output2 = getattr(protCoordinates2, 'outputCoordinates', None)
        self.assertTrue(output2,
                             "There was a problem with coordinates 3d output")
        self.assertTrue(output2.getSize() == 5)
        self.assertTrue(output2.getBoxSize() == 32)
        self.assertTrue(output2.getSamplingRate() == 5.5)
        outputs.append(output2)

        if existsPlugin('dynamo'):
            protCoordinates2 = self._runTomoImportSetOfCoordinates('*.tbl', 'DYNAMO', 'TBL')
            output3 = getattr(protCoordinates2, 'outputCoordinates', None)
            self.assertTrue(output2,
                                 "There was a problem with coordinates 3d output")
            self.assertTrue(output2.getSize() == 5)
            self.assertTrue(output2.getBoxSize() == 32)
            self.assertTrue(output2.getSamplingRate() == 5.5)
            outputs.append(output3)

        return outputs


class TestTomoPreprocessing(BaseTest):
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

    def test_preprocess1(self):
        try:
            import gctf.protocols
            import imod.protocols
        except ImportError:
            print('To run this test scipion-em-gctf and scipion-em-imod plugins are required. '
                  'This test will not be executed unless this plugins are installed')
            return

        """ Run the basic preprocessing pipeline for just one TS. """
        protImport = self._runImportTiltSeriesM(
            filesPattern='{TS}1_{TO}_{TA}.mrc')
        output = getattr(protImport, 'outputTiltSeriesM', None)
        self.assertFalse(output is None)
        self.assertEqual(protImport.outputTiltSeriesM.getSize(), 1)

        gpuList = os.environ.get('SCIPION_TEST_GPULIST', '0')
        threads = len(gpuList.split()) + 1

        # --------- Motion correction with motioncor2 for Tilt-series ------
        protMc = self.newProtocol(
            tomo.protocols.ProtTsAverage,
            inputTiltSeriesM=protImport.outputTiltSeriesM,
            binFactor=2.0,
            gpuList=gpuList,
            numberOfThreads=threads
        )
        self.launchProtocol(protMc)

        # --------- CTF estimation with Gctf ------
        protGctf = self.newProtocol(
            gctf.protocols.ProtTsGctf,
            inputTiltSeries=protMc.outputTiltSeries,
            gpuList=gpuList,
            numberOfThreads=threads,
        )
        self.launchProtocol(protGctf)

        # # -------- Basic alignment and reconstruction with IMOD ------
        # protImodAuto = self.newProtocol(
        #     imod.protocols.ProtImodAuto3D,
        #     inputTiltSeries=protGctf.outputTiltSeries,
        #     excludeList=1,
        #     rotationAngle=90,
        #     zWidth=400,
        #     useRaptor=True,
        #     markersDiameter=20,
        #     markersNumber=20
        # )
        # self.launchProtocol(protImodAuto)


class TestTomoAssignAlignment(BaseTest):
    """This class check if the protocol assign alignments for subtomograms works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runPreviousProtocols(self):
        protImport = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)
        protImport2 = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport2)
        from xmipp2.protocols import Xmipp2ProtMLTomo
        protMltomo = self.newProtocol(Xmipp2ProtMLTomo,
                                      inputVolumes=protImport.outputSubTomograms,
                                      randomInitialization=True,
                                      numberOfReferences=1,
                                      numberOfIters=3,
                                      angularSampling=30)
        self.launchProtocol(protMltomo)
        self.assertIsNotNone(protMltomo.outputSubtomograms,
                         "There was a problem with SetOfSubtomogram output")
        self.assertIsNotNone(protMltomo.outputClassesSubtomo,
                         "There was a problem with SetOfSubtomogram output")
        return protImport2, protMltomo

    def _assignAlignment(self):
        protImport2, protMltomo = self._runPreviousProtocols()
        assign = self.newProtocol(tomo.protocols.ProtAlignmentAssignSubtomo,
                                 inputSubtomos=protImport2.outputSubTomograms,
                                 inputAlignment=protMltomo.outputSubtomograms)
        self.launchProtocol(assign)
        self.assertIsNotNone(assign.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return assign

    def test_assignAlignment(self):
        assign = self._assignAlignment()
        self.assertTrue(getattr(assign, 'outputSubtomograms'))
        self.assertTrue(assign.outputSubtomograms.getFirstItem().hasTransform())
        return assign


class TestTomoAssignTomo2Subtomo(BaseTest):
    """This class check if the protocol assign tomograms to subtomograms works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')
        cls.setOfTomograms = cls.dataset.getFile('*.em')

    def _runPreviousProtocols(self):
        protImportTomo = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                        filesPath=self.setOfTomograms,
                                        filesPattern='',
                                        acquisitionAngleMax=40,
                                        acquisitionAngleMin=-40,
                                        samplingRate=1.35)
        self.launchProtocol(protImportTomo)
        protImportSubtomo = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                              filesPath=self.setOfSubtomograms,
                                              samplingRate=5)
        self.launchProtocol(protImportSubtomo)
        return protImportTomo, protImportSubtomo

    def _assignTomos2subtomos(self):
        protImportTomo, protImportSubtomo = self._runPreviousProtocols()
        tomo2subtomo = self.newProtocol(tomo.protocols.ProtAssignTomo2Subtomo,
                                 inputSubtomos=protImportSubtomo.outputSubTomograms,
                                 inputTomos=protImportTomo.outputTomograms)
        self.launchProtocol(tomo2subtomo)
        self.assertIsNotNone(tomo2subtomo.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return tomo2subtomo

    def test_assignTomos2subtomos(self):
        tomo2subtomo = self._assignTomos2subtomos()
        self.assertTrue(getattr(tomo2subtomo, 'outputSubtomograms'))
        self.assertFalse(tomo2subtomo.outputSubtomograms.getFirstItem().getVolName())
        return tomo2subtomo


class TestTomoSplitEvenOdd(BaseTest):
    """This class check if the protocol split even/odd tomos or subtomos works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfTomograms = cls.dataset.getFile('*.em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _splitTomoSet(self):
        protImportTomo = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                        filesPath=self.setOfTomograms,
                                        filesPattern='',
                                        acquisitionAngleMax=40,
                                        acquisitionAngleMin=-40,
                                        samplingRate=1.35)
        self.launchProtocol(protImportTomo)
        split = self.newProtocol(tomo.protocols.ProtSplitEvenOddTomoSet,
                                 inputSet=protImportTomo.outputTomograms)
        self.launchProtocol(split)
        self.assertIsNotNone(split.outputset_even,
                             "There was a problem with even tomograms output")
        self.assertIsNotNone(split.outputset_odd,
                             "There was a problem with odd tomograms output")
        return split

    def _splitSubtomoSet(self):
        protImportSubtomo = self.newProtocol(tomo.protocols.ProtImportSubTomograms,
                                              filesPath=self.setOfSubtomograms,
                                              samplingRate=5)
        self.launchProtocol(protImportSubtomo)
        split = self.newProtocol(tomo.protocols.ProtSplitEvenOddTomoSet,
                                 inputSet=protImportSubtomo.outputSubTomograms)
        self.launchProtocol(split)
        self.assertIsNotNone(split.outputset_even,
                             "There was a problem with even subtomograms output")
        self.assertIsNotNone(split.outputset_odd,
                             "There was a problem with odd subtomograms output")
        return split

    def test_splitTomos(self):
        split = self._splitTomoSet()
        self.assertTrue(getattr(split, 'outputset_even'))
        self.assertTrue(getattr(split, 'outputset_odd'))
        self.assertTrue(split.outputset_even.getFirstItem().getObjId() % 2 == 0,
                        "Elements id in even set of tomos are not even")
        self.assertTrue(split.outputset_odd.getFirstItem().getObjId() % 2 != 0,
                        "Elements id in odd set of tomos are not odd")
        self.assertSetSize(split.outputset_even, 1, "Size of even set of tomos is not correct")
        self.assertSetSize(split.outputset_odd, 1, "Size of odd set of tomos is not correct")

        return split

    def test_splitSubtomos(self):
        split = self._splitSubtomoSet()
        self.assertTrue(getattr(split, 'outputset_even'))
        self.assertTrue(getattr(split, 'outputset_odd'))
        self.assertTrue(split.outputset_even.getFirstItem().getObjId() % 2 == 0,
                        "Elements id in even set of subtomos are not even")
        self.assertTrue(split.outputset_odd.getFirstItem().getObjId() % 2 != 0,
                        "Elements id in odd set of subtomos are not odd")
        self.assertSetSize(split.outputset_even, 2, "Size of even set of subtomos is not correct")
        self.assertSetSize(split.outputset_odd, 2, "Size of odd set of subtomos is not correct")

        return split


if __name__ == 'main':
    pass
