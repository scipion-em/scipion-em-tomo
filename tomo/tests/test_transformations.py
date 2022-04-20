# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
# *
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
# * Foundation, Inc., shifts9 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import weakImport

with weakImport("xmipptomo"):
    from xmipptomo.protocols.protocol_phantom_subtomo import XmippProtPhantomSubtomo

    # Define the class inside so only the test is available if xmipptomo is.
    class TestSubtomoTransformations(BaseTest):
        """ This class check if the transformation matrices are managed properly across different software."""

        @classmethod
        def setUpClass(cls):
            setupTestProject(cls)


        def testOnlyShifts(self):
            # Angles =1 cancels the rotation
            self._runPhantomSubtomo(angles=1)

        def testOnlyAngles(self):
            # Shifts =1 cancels the shifts
            self._runPhantomSubtomo(shifts=1)

        def testAll(self):

            self._runPhantomSubtomo()

        def _addAveragers(self, inputProt:XmippProtPhantomSubtomo):
            """ Add all the averagers available to be run with the subtomo set"""

            # Add reliontomo average method
            with weakImport("reliontomo"):
                from reliontomo.protocols import ProtRelionSubTomoReconstructAvg
                relionAve = self.newProtocol(ProtRelionSubTomoReconstructAvg)
                relionAve.inputSubtomos.set(inputProt)
                relionAve.inputSubtomos.setExtended(inputProt._possibleOutputs.outputSubtomograms.name)
                self.launchProtocol(relionAve)

            # Add eman average method
            with weakImport("emantomo"):
                from emantomo.protocols import EmanProtSubTomoAverage
                emanAve = self.newProtocol(EmanProtSubTomoAverage)
                emanAve.inputSetOfSubTomogram.set(inputProt)
                emanAve.inputSetOfSubTomogram.setExtended(inputProt._possibleOutputs.outputSubtomograms.name)
                self.launchProtocol(emanAve)

            # Add xmipp tomo averager. In this case we do not need the weak import since this test will only run if xmipptomo is available
            from xmipptomo.protocols import XmippProtApplyTransformSubtomo
            xmippAve = self.newProtocol(XmippProtApplyTransformSubtomo)
            xmippAve.inputSubtomograms.set(inputProt)
            xmippAve.inputSubtomograms.setExtended(inputProt._possibleOutputs.outputSubtomograms.name)
            self.launchProtocol(xmippAve)

        def _runPhantomSubtomo(self, shifts=5, angles=90):
            """ Creates a a phantom subtomo protocol with shift and rotation by default

            :param shifts: value for the 3 shifts, to cancel shifts use 1
            :param angles: value for the 3 angles, to cancel rotation use 1
            """

            PHANTOM_DESC="""60 60 60 0
    con = 1 0 0 0 8 30 0 180 0
    con = 1 0 0 0 8 30 0 0 0
    con = 1 0 0 0 8 30 0 90 0
    con = 1 0 0 0 8 30 0 270 0
    con = 1 0 0 0 8 30 90 90 0
    con = 1 0 0 0 8 30 90 270 0"""

            shiftsStr = "NO" if shifts==1 else "0-%s" % shifts
            anglesStr = "NO" if angles==1 else "0-%s" % angles
            label = "ST phantoms SHIFT %s, ROT %s" % (shiftsStr, anglesStr)

            protPhantom = self.newProtocol(XmippProtPhantomSubtomo,
                                           objLabel=label,
                                           option=1, # Use phantom description
                                           create=PHANTOM_DESC,
                                           sampling=4,
                                           nsubtomos=10,
                                           rotmin=0,
                                           rotmax=angles,
                                           tiltmin=0,
                                           tiltmax=angles,
                                           psimin=0,
                                           psimax=angles,
                                           xmin=0,
                                           xmax=shifts,
                                           ymin=0,
                                           ymax=shifts,
                                           zmin=0,
                                           zmax=shifts)

            self.launchProtocol(protPhantom)

            # Add the averagers
            self._addAveragers(protPhantom)


    # Testing phantoms in tomograms, going through coordinates, extraction and average
    from xmipptomo.protocols.protocol_phantom_tomo import XmippProtPhantomTomo

    # Define the class inside so only the test is available if xmipptomo is.
    class TestTomoTransformations(BaseTest):
        """ This class check if the transformation matrices from the 3D coordinates are managed properly
         across different software and along the pipeline."""

        @classmethod
        def setUpClass(cls):
            setupTestProject(cls)

        def testNoAngles(self):
            # angles =1 cancels the rotation
            self._runPhantomSubtomo(angles=1)

        def testAngles(self):

            self._runPhantomSubtomo()

        def _addAveragers(self, inputProt:XmippProtPhantomSubtomo):
            """ Add all the averagers available to be run with the subtomo set"""

            # Add reliontomo average method
            with weakImport("reliontomo"):
                from reliontomo.protocols import ProtRelionSubTomoReconstructAvg
                relionAve = self.newProtocol(ProtRelionSubTomoReconstructAvg)
                relionAve.inputSubtomos.set(inputProt)
                relionAve.inputSubtomos.setExtended(inputProt._possibleOutputs.outputSubtomograms.name)
                self.launchProtocol(relionAve)

            # Add eman average method
            with weakImport("emantomo"):
                from emantomo.protocols import EmanProtSubTomoAverage
                emanAve = self.newProtocol(EmanProtSubTomoAverage)
                emanAve.inputSetOfSubTomogram.set(inputProt)
                emanAve.inputSetOfSubTomogram.setExtended(inputProt._possibleOutputs.outputSubtomograms.name)
                self.launchProtocol(emanAve)

            # Add xmipp tomo averager. In this case we do not need the weak import since this test will only run if xmipptomo is available
            from xmipptomo.protocols import XmippProtApplyTransformSubtomo
            xmippAve = self.newProtocol(XmippProtApplyTransformSubtomo)
            xmippAve.inputSubtomograms.set(inputProt)
            xmippAve.inputSubtomograms.setExtended(inputProt._possibleOutputs.outputSubtomograms.name)
            self.launchProtocol(xmippAve)

        def _runPhantomSubtomo(self, angles=90):
            """ Creates a a phantom subtomo protocol with rotation by default

            :param angles: value for the 3 angles, to cancel rotation use 1
            """

            anglesStr = "NO" if angles==1 else "0-%s" % angles
            label = "TOMO phantoms ROT %s" % anglesStr

            protPhantom = self.newProtocol(XmippProtPhantomTomo,
                                           objLabel=label,
                                           dimensions="100 100 50",
                                           sampling=4,
                                           nparticles=10,
                                           ntomos=3,
                                           rotmin=0,
                                           rotmax=angles,
                                           tiltmin=0,
                                           tiltmax=angles,
                                           psimin=0,
                                           psimax=angles)

            self.launchProtocol(protPhantom)

            # Add the averagers
            #self._addAveragers(protPhantom)
