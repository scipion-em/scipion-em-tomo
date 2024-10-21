# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from glob import glob

from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr

from tomo.tests import DataSet, DataSetRe5STA, RE5_STA
from tomo.protocols import (ProtImportTomograms, ProtImportCoordinates3D,
                            ProtExportCoordinates3D)
from tomo.utils import existsPlugin


class TestTomoExportCoordinates3D(BaseTest):
    """This class checks if the protocol to export 3d coordinates works properly. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds1 = DataSet.getDataSet('tomo-em')
        cls.sRate1 = 5.5
        cls.boxSize1 = 32
        cls.ds2 = DataSet.getDataSet(RE5_STA)

    def _importTomograms(self, filesPath, samplingRate):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=filesPath,
                                              filesPattern=None,
                                              samplingRate=samplingRate)
        self.launchProtocol(protImportTomogram)
        self.assertIsNotNone(protImportTomogram.Tomograms,
                             msg='No tomograms were imported.')

        return protImportTomogram.Tomograms

    @classmethod
    def _runImportCoordinates(cls, filesPath, pattern, program,
                              importTomograms=None,
                              coordsSRate=None, boxSize=None):
        print(magentaStr(f"\n==> Importing data - coordinates 3D ({pattern}):"))
        protImportCoordinates3d = cls.newProtocol(ProtImportCoordinates3D,
                                                  importFrom=0,  # auto
                                                  filesPath=filesPath,
                                                  importTomograms=importTomograms,
                                                  filesPattern=pattern,
                                                  boxSize=boxSize,
                                                  samplingRate=coordsSRate)
        protImportCoordinates3d.setObjLabel(f"Import from {program} - {pattern.lstrip('*.')}")
        cls.launchProtocol(protImportCoordinates3d)

        return protImportCoordinates3d

    @classmethod
    def _runExportCoordinates(cls, inputCoords, program, format="txt"):
        print(magentaStr(f"\n==> Testing tomo - export coordinates 3D to {format}:"))
        outputFormat = ProtExportCoordinates3D._getExportChoices().index(format)
        protExportCoordinates3d = cls.newProtocol(ProtExportCoordinates3D,
                                                  outputFormat=outputFormat,
                                                  inputCoordinates=inputCoords)
        protExportCoordinates3d.setObjLabel(f"Export to {program} - {format}")
        cls.launchProtocol(protExportCoordinates3d)

        return protExportCoordinates3d

    def test_exportCoordinates3D(self):
        importedTomos = self._importTomograms(
            self.ds1.getFile('tomo1'), self.sRate1)

        importedTomos2 = self._importTomograms(
            self.ds2.getFile(DataSetRe5STA.tomosDir.name) + "/*.mrc",
            DataSetRe5STA.tomosSRate.value)

        def _runWorkflow(filesPath, pattern, tomos, program, format,
                         pixelSize, boxSize):
            importCoords = self._runImportCoordinates(filesPath, pattern, program,
                                                      importTomograms=tomos,
                                                      coordsSRate=pixelSize,
                                                      boxSize=boxSize)
            exportCoords = self._runExportCoordinates(importCoords.outputCoordinates,
                                                      program, format)
            self.checkOutputCoords(exportCoords, pattern)

        _runWorkflow(filesPath=self.ds1.getPath(), pattern="*.txt",
                     tomos=importedTomos, program="tomo", format="txt",
                     pixelSize=self.sRate1, boxSize=self.boxSize1)

        if existsPlugin("sphire"):
            _runWorkflow(filesPath=self.ds1.getPath(), pattern="*.cbox",
                         tomos=importedTomos, program="cryolo", format="cbox",
                         pixelSize=self.sRate1, boxSize=self.boxSize1)

        if existsPlugin('emantomo'):
            _runWorkflow(filesPath=self.ds1.getPath(), pattern="*.json",
                         tomos=importedTomos, program="EMAN", format="eman",
                         pixelSize=self.sRate1, boxSize=self.boxSize1)

        if existsPlugin('dynamo'):
            _runWorkflow(filesPath=self.ds1.getPath(), pattern="*.tbl",
                         tomos=importedTomos, program="Dynamo", format="dynamo",
                         pixelSize=self.sRate1, boxSize=self.boxSize1)

        if existsPlugin("reliontomo"):
            from reliontomo.protocols import ProtImportCoordinates3DFromStar
            print(magentaStr(f"\n==> Importing data - coordinates 3D (*.star):"))
            starFile = self.ds2.getFile(DataSetRe5STA.coordsPickedWithRe5Star.name)
            importCoords3dFromStar = self.newProtocol(ProtImportCoordinates3DFromStar,
                                                      starFile=starFile,
                                                      inTomos=importedTomos2,
                                                      boxSize=DataSetRe5STA.boxSize.value)
            self.launchProtocol(importCoords3dFromStar)

            exportCoords = self._runExportCoordinates(importCoords3dFromStar.coordinates,
                                                      program="Relion", format="relion")
            self.checkOutputCoords(exportCoords, "*.star")

    def checkOutputCoords(self, protocol, pattern="*.txt"):
        outputFiles = glob(protocol._getExportPath(pattern))
        self.assertIsNotNone(outputFiles)
