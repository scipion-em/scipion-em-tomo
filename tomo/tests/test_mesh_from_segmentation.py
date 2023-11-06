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
import time
from pyworkflow.tests import BaseTest, setupTestProject
from . import DataSet
import pwem.protocols as emprot
from tomo.protocols.protocol_mesh_from_segmentation import ProtMeshFromSegmentation, OUTPUTATTRIBUTE
from tomo.protocols.protocol_import_tomograms import ProtImportTomograms, OUTPUT_NAME
from tomo.protocols.protocol_import_tomomasks import ProtImportTomomasks
from pyworkflow.plugin import Domain
import pyworkflow.protocol as pwprot


class TestMeshFromSegmentation(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tomoMask')
        #cls.pattern = '*.mrc'


    def testCreateMesh(self):
        # Import tomograms
        protTomograms = self.newProtocol(ProtImportTomograms,
                                         objLabel='Import tomograms',
                                         filesPath='/home/vilas/software/data/tests/tomoMask/',#self.dataset.getPath(),
                                         filesPattern='*.mrc',
                                         samplingRate=20)
        self.launchProtocol(protTomograms)
        outputTomograms = getattr(protTomograms, OUTPUT_NAME, None)
        self.assertIsNotNone(outputTomograms, "There was a problem with tomogram output")


        # Import tomoMask
        protTomoMask = self.newProtocol(ProtImportTomomasks,
                                           objLabel='Import tomoMaks',
                                           inputTomos=outputTomograms,
                                           filesPath='/home/vilas/software/data/tests/tomoMask/',#self.dataset.getPath(),,
                                           filesPattern='vTomo*_flt.mrc')
        self.launchProtocol(protTomoMask)
        self.assertIsNotNone(protTomoMask.outputTomoMasks, 'tomoMasks not imported')

        # Create Mesh
        protMesh = self.newProtocol(ProtMeshFromSegmentation,
                                       objLabel='mesh fromsegmentation',
                                       tomoMasks=protTomoMask.outputTomoMasks,
                                       lowLimit=0.1,
                                       highLimit=1.0,
                                       density=0.05)

        self.launchProtocol(protMesh)
        self.assertIsNotNone(protMesh.outputMeshes, 'mesh from segmentation failed')
