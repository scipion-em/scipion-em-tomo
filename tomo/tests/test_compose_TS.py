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
import os
from pyworkflow.tests import BaseTest, setupTestProject, setupTestOutput
from . import DataSet
import pwem.protocols as emprot
from tomo.protocols.protocol_compose_TS import ProtComposeTS
from pyworkflow.plugin import Domain


class TestTomoComposeTS(BaseTest):
    """ This class check if the protocol create the TS stack in streaming
     properly."""

    @classmethod
    def setUpClass(self):
        setupTestProject(self)
        self.inputDataSet = DataSet.getDataSet('tomo-em')
        self.partFolderPath = self.inputDataSet.getFile('ts_tsM_and_mdocs/*')
        self.pattern = '*.mrc'


    def testCompose(self):
        #Import movies
        protMovieImport = self.newProtocol(emprot.ProtImportMovies,
                        objLabel='Import movies (SPA)',
                        importFrom=emprot.ProtImportMovies.IMPORT_FROM_FILES,
                        filesPath=self.partFolderPath,
                        filesPattern=self.pattern,
                        voltage=300)
        self.launchProtocol(protMovieImport)
        self.assertIsNotNone(protMovieImport.outputMovies, 'Movies not imported')

        #Align movies and create a set of micrographs
        xmipp3 = Domain.importFromPlugin('xmipp3.protocols', doRaise=True)
        protAlign = self.newProtocol(xmipp3.XmippProtMovieCorr,
                                      objLabel='Movie Alignment (SPA)',
                                      alignFrame0=1, alignFrameN=0,
                                      useAlignToSum=True,
                                      doLocalAlignment=False)
        protAlign.inputMovies.set(protMovieImport.outputMovies)
        self.launchProtocol(protAlign)
        self.assertIsNotNone(protAlign.outputMicrographs, 'Micrograph not generated')

        #ComposeTS
        protCompose = self.newProtocol(ProtComposeTS,
                                 objLabel='ComposeTS',
                                 inputMicrographs=protAlign.outputMicrographs,
                                 filesPath=self.partFolderPath,
                                 time4NextTilt=20,
                                 time4NextMic=12,
                                 time4NextTS=60)
        self.launchProtocol(protCompose)
        self.assertIsNotNone(protCompose.TiltSeries, 'Set of Tilt Series not generated')

        #xcor prealigment
        imod = Domain.importFromPlugin('imod.protocols', doRaise=True)
        prealigment= self.newProtocol(imod.ProtImodXcorrPrealignment,
                                      computeAlignment=0,
                                      binning=2)
        prealigment.inputSetOfTiltSeries.set(protCompose.TiltSeries)
        self.launchProtocol(prealigment)
        self.assertIsNotNone(prealigment.TiltSeries, 'TiltSeries not galignment')


