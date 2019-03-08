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
from pyworkflow.em import Domain

from tomo.objects import SetOfTiltSeriesM, SetOfTiltSeries
from tomo.protocols import ProtImportTiltSeries


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

    def test_create_tiltseries(self):
        self._create_tiltseries(SetOfTiltSeries)

    def test_create_tiltseriesM(self):
        self._create_tiltseries(SetOfTiltSeriesM)


class TestTomoImport(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _runImportTiltSeries(self, filesPath, filesPattern,
                             voltage, magnification, samplingRate,
                             importType, dosePerFrame, gainFile=None, darkFile=None,
                             sphericalAberration=2.7, amplitudeContrast=0.1, doseInitial=0):
        return self.newProtocol(ProtImportTiltSeries,
                                importType=importType,
                                filesPath=filesPath,
                                filesPattern=filesPattern,
                                voltage=voltage,
                                magnification=magnification,
                                sphericalAberration=sphericalAberration,
                                amplitudeContrast=amplitudeContrast,
                                samplingRate=samplingRate,
                                doseInitial=doseInitial,
                                dosePerFrame=dosePerFrame,
                                gainFile=gainFile,
                                darkFile=darkFile)

    def test_import_tiltseries(self):
        dataPath = os.environ.get('SCIPION_TOMO_EMPIAR10164', '')

        if not os.path.exists(dataPath):
            raise Exception("Can not run tomo tests, "
                            "SCIPION_TOMO_EMPIAR10164 variable not defined. ")

        filesPath = os.path.join(dataPath, 'data', 'frames')
        protImport = self._runImportTiltSeries(filesPath,
                                               '{TS}_{TO}_{TA}.mrc',
                                               300, 105000, 1.35,
                                               ProtImportTiltSeries.IMPORT_TYPE_MOVS,
                                               0.3)
        self.launchProtocol(protImport)
        output = getattr(protImport, 'outputTiltSeriesM', None)
        self.assertFalse(output is None)


if __name__ == 'main':
    pass