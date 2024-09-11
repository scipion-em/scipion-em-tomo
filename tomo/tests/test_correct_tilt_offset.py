# **************************************************************************
# *
# * Authors:     Alberto García Mena (alberto.garcia@cnb.csic.es)
# *# *
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
from . import DataSet

from tomo.protocols.protocol_correct_tilt_offset import ProtCorrectTiltOffset
from tomo.protocols.protocol_ts_import import ProtImportTs
from tomo.protocols.protocol_import_tomomasks import ProtImportTomomasks
import os



class TestCorrectTiltOffset(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.getFile = cls.dataset.getFile('tutorialData')

    def testCorrectTiltOffSet(self):
        minAngle = -55
        maxAngle = 65
        protImport = self.newProtocol(ProtImportTs,
                                    filesPath=self.getFile,
                                    filesPattern='BB{TS}.st',
                                    minAngle=minAngle,
                                    maxAngle=maxAngle,
                                    stepAngle=2,
                                    voltage=300,
                                    magnification=105000,
                                    sphericalAberration=2.7,
                                    amplitudeContrast=0.1,
                                    samplingRate=1.35,
                                    doseInitial=0,
                                    dosePerFrame=0.3,
                                    tiltAxisAngle=84.7
        )
        self.launchProtocol(protImport)



        self.assertFalse(protImport.outputTiltSeries is None)
        self.assertEqual(protImport.outputTiltSeries.getSize(), 2)
        

        tiltOffset=10
        protCorrectOffset = self.newProtocol(ProtCorrectTiltOffset,
                                             inputTiltSeries=getattr(protImport, 'outputTiltSeries'),
                                             tiltOffset=tiltOffset)
        self.launchProtocol(protCorrectOffset)

        outTsSet = getattr(protCorrectOffset, protCorrectOffset._possibleOutputs.tiltSeries.name)
        acqInfo = outTsSet.getAcquisition()
        self.assertEqual(acqInfo.getAngleMin(), minAngle+tiltOffset)
        self.assertEqual(acqInfo.getAngleMax(), maxAngle+tiltOffset)



