# **************************************************************************
# *
# * Authors:    Estrella Fernandez Gimenez [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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

from tomo.objects import SetOfCoordinates3D, Coordinate3D
from pwem.protocols import EMProtocol
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import (ProtImportCoordinates3D,
                            ProtImportTomograms, TomoProtFitEllipsoid)
from pyworkflow.utils import weakImport


class TestTomoProtFitVesicles(BaseTest):
    """ This class check if the protocol fit vesicles works properly."""
    tomograms = None
    coordinates = None
    meshes = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')

    def __init__(self, methodName):
        super().__init__(methodName)

    def getTomograms(self):
        if self.tomograms is None:
            protImportTomogram = self.newProtocol(ProtImportTomograms,
                                                  filesPath=self.tomogram,
                                                  samplingRate=5)
            self.launchProtocol(protImportTomogram)
            self.assertIsNotNone(protImportTomogram.Tomograms,
                                 "There was a problem with tomogram output")

            self.tomograms = protImportTomogram.Tomograms
        return self.tomograms

    def getCoordinates(self):

        if self.coordinates is None:

            tomograms = self.getTomograms()
            protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                       filesPath=self.coords3D,
                                                       importTomograms=tomograms,
                                                       boxSize=32,
                                                       samplingRate=5)
            self.launchProtocol(protImportCoordinates3d)

            outputCoords = protImportCoordinates3d.outputCoordinates

            # Customize the coordinates

            newCoords =[
                (0, 0, -50, '1', "overview_wbp"),
                (-2, -2, -48, '1', "overview_wbp"),
                (2, 2, -48, '1', "overview_wbp"),
                (-2, 2, -48, '1', "overview_wbp"),
                (2, -2, -48, '1', "overview_wbp"),
                (-50, 0, 0, '1', "overview_wbp"),
                (0, -50, 0, '1', "overview_wbp"),
                (50, 0, 0, '1', "overview_wbp"),
                (0, 50, 0, '1', "overview_wbp"),
                (0, 0, 50, '1', "overview_wbp"),

            ]
            self._replaceCoordinates( outputCoords, newCoords, protImportCoordinates3d)
            self.coordinates= outputCoords

        return self.coordinates

    @classmethod
    def _replaceCoordinates(cls, coordset:SetOfCoordinates3D, newCoords, prot:EMProtocol):
        """ Replaces the coordinates of a set with the coords list

        :param coordset: Scipion's Coordinate set to be replaced
        :param newCoords: list of tuples with x, y, z, groupId, tsId
        :param prot: Protocol holding the coordinates set
        """

        coordset.load()
        coordset.clear()

        for newCoordTuple in newCoords:

            newCord = Coordinate3D()
            x,y,z,group, tsId = newCoordTuple
            newCord._x.set(x)
            newCord._y.set(y)
            newCord._z.set(z)
            newCord.setGroupId(group)
            newCord.setTomoId(tsId)
            coordset.append(newCord)

        coordset.write()

        # Save the protocol to persist the size
        prot._store(coordset)

    def getMeshes(self):
        if self.meshes is None:
            protFitVesicles = self.newProtocol(TomoProtFitEllipsoid,
                                               input=self.getCoordinates(),
                                               inputTomos=self.getTomograms())
            self.launchProtocol(protFitVesicles)
            self.assertIsNotNone(protFitVesicles.outputMeshes, "There was a problem with output vesicles (SetOfMeshes)")
            self.assertEqual(protFitVesicles.outputMeshes.getSize(), 292)

            self.meshes = protFitVesicles.outputMeshes

        return self.meshes

    def test_filterByNormals(self):
        self.getMeshes()


        with weakImport("tomoviz"):

            from tomoviz.protocols import TomovizProtFilterbyNormal
            from tomoviz.protocols.protocol_mesh_normal import OUTPUT_NAME

            def filterByNormal(tolerance, expectedSize):


                filterByNormals = self.newProtocol(TomovizProtFilterbyNormal,
                                                   input=self.getCoordinates(),
                                                   inputMeshes=self.getMeshes(),
                                                   tol=tolerance)

                self.launchProtocol(filterByNormals)

                self.assertSetSize(getattr(filterByNormals, OUTPUT_NAME), expectedSize,
                                   "Filter by normal does not work for %d degree tolerance." % tolerance)

            # NOTE: seems that normals are based on closest triangle that might not be 90 degrees to the coordinates.
            filterByNormal(15, 1)

            filterByNormal(75, 5)



