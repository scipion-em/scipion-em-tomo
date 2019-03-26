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

import pyworkflow as pw
from pyworkflow.tests import BaseTest, setupTestOutput, DataSet, setupTestProject
from pyworkflow.em import Domain, CTFModel

from tomo.objects import SetOfTiltSeriesM, SetOfTiltSeries
import tomo.protocols


class TestTomoBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        #cls.dataset = DataSet.getDataSet('relion_tutorial')
        #cls.getFile = cls.dataset.getFile

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
            self.assertTrue(e in objects, "%s should be in Domain.getObjects" % e)

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
        cls.dataPath = os.environ.get('SCIPION_TOMO_EMPIAR10164', '')

        if not os.path.exists(cls.dataPath):
            raise Exception("Can not run tomo tests, "
                            "SCIPION_TOMO_EMPIAR10164 variable not defined. ")

    def _runImportTiltSeriesM(self):
        protImport = self.newProtocol(
            tomo.protocols.ProtImportTiltSeries,
            importType=tomo.protocols.ProtImportTiltSeries.IMPORT_TYPE_MOVS,
            filesPath=os.path.join(self.dataPath, 'data', 'frames'),
            filesPattern='{TS}_{TO}_{TA}.mrc',
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

        protMc = self.newProtocol(
            tomo.protocols.ProtMotionCorrTiltSeries,
            binFactor=2.0
        )

        protMc.inputTiltSeriesM.set(protImport.outputTiltSeriesM)
        self.launchProtocol(protMc)


if __name__ == 'main':
    pass