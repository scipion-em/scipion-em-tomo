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
from tomo.protocols import ProtTomoMisalignTiltSeries


with weakImport("xmipptomo"):
    from xmipptomo.protocols.protocol_phantom_subtomo import XmippProtPhantomSubtomo, OutputPhantomSubtomos


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
            self._runPhantomSubtomo(expected=(51,14,71,4,4,1))

        def _runPhantomSubtomo(self, shifts=5, angles=90, expected=None):
            """ Creates a phantom subtomo protocol with shift and rotation by default

            :param shifts: value for the 3 shifts, to cancel shifts use 1
            :param angles: value for the 3 angles, to cancel rotation use 1
            :param expected: List of Tuple or just a tuple where the tuple are expected values for rot, tilt, psi, x, y, z
            """

            PHANTOM_DESC = """60 60 60 0
    con = 1 0 0 0 8 30 0 180 0
    con = 1 0 0 0 8 30 0 0 0
    con = 1 0 0 0 8 30 0 90 0
    con = 1 0 0 0 8 30 0 270 0
    con = 1 0 0 0 8 30 90 90 0
    con = 1 0 0 0 8 30 90 270 0"""

            # This is an asymmetric phantom. For more simple approach use the description above
            PHANTOM_DESC = """60 60 60 0
    cub = 1 0 -4 0 24 5 5 0 0 0
    ell = 1 0 0 -2 5 8 5 0 0 0
    con = 1 0 0 0 6 24 0 0 0"""

            shiftsStr = "NO" if shifts == 1 else "0-%s" % shifts
            anglesStr = "NO" if angles == 1 else "0-%s" % angles
            label4Averagers = "SHIFT %s, ROT %s" % (shiftsStr, anglesStr)
            label = "ST phantoms %s" % label4Averagers

            protPhantom = self.newProtocol(XmippProtPhantomSubtomo,
                                           objLabel=label,
                                           option=1,  # Use phantom description
                                           create=PHANTOM_DESC,
                                           sampling=4,
                                           nsubtomos=10,
                                           randomseed=True, # Force random function to return the same values.
                                           rotate=True,
                                           uniformAngularDistribution=False,
                                           rotmin=0,
                                           rotmax=angles,
                                           tiltmin=0,
                                           tiltmax=angles,
                                           psimin=0,
                                           psimax=angles,
                                           applyShift=True,
                                           xmin=0,
                                           xmax=shifts,
                                           ymin=0,
                                           ymax=shifts,
                                           zmin=0,
                                           zmax=shifts)

            self.launchProtocol(protPhantom)


            output = getattr(protPhantom, OutputPhantomSubtomos.outputSubtomograms.name)
            if expected:
                if not isinstance(expected, list):
                    expected = [expected]

                for index, expectedvalues in enumerate(expected):

                    rot, tilt, psi, x, y, z = expectedvalues

                    item = output[index+1]

                    self.assertEquals(rot,item.phantom_rot.get(), "Rot value not expected fot item %s" % index )
                    self.assertEquals(tilt, item.phantom_tilt.get(), "Tilt value not expected fot item %s" % index)
                    self.assertEquals(psi, item.phantom_psi.get(), "Psi value not expected fot item %s" % index)
                    self.assertEquals(x, item.phantom_shiftX.get(), "Shift X value not expected fot item %s" % index)
                    self.assertEquals(y, item.phantom_shiftY.get(), "Shift Y value not expected fot item %s" % index)
                    self.assertEquals(z, item.phantom_shiftZ.get(), "Shift Z value not expected fot item %s" % index)


            # Add the averagers
            addAveragers(self, protPhantom, XmippProtPhantomSubtomo._possibleOutputs.outputSubtomograms.name, label4Averagers)


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
            self._runPhantomTomo("NO ANGLES", rot=1, tilt=1, psi=1)

        def testAngles(self):
            self._runPhantomTomo("ALL ANGLES")

        def testOnlyRot(self):
            self._runPhantomTomo("ONLY ROT", tilt=1, psi=1)

        def testOnlyTilt(self):
            self._runPhantomTomo("ONLY TILT", rot=1, psi=1)

        def testOnlyPsi(self):
            self._runPhantomTomo("ONLY PSI", rot=1, tilt=1)

        def _runPhantomTomo(self, label, rot=90, tilt=90, psi=90):
            """ Creates a phantom tomo protocol with rotation by default

            :param angles: value for the 3 angles, to cancel rotation use 1
            """

            protPhantom = self.newProtocol(XmippProtPhantomTomo,
                                           objLabel=label,
                                           dimensions="200 200 100",
                                           sampling=4,
                                           nparticles=10,
                                           ntomos=3,
                                           rotmin=0,
                                           rotmax=rot,
                                           tiltmin=0,
                                           tiltmax=tilt,
                                           psimin=0,
                                           psimax=psi)

            self.launchProtocol(protPhantom)

            with weakImport("emantomo"):
                # Add eman extraction
                from emantomo.protocols.protocol_extraction_from_tomo import EmanProtTomoExtraction
                stExtraction = self.newProtocol(EmanProtTomoExtraction,
                                                inputCoordinates=protPhantom.coordinates3D,
                                                boxSize=protPhantom.coordinates3D.getBoxSize(),
                                                doInvert=True
                                                )

                self.launchProtocol(stExtraction)

                # Add the averagers
                addAveragers(self, stExtraction, EmanProtTomoExtraction._possibleOutputs.subtomograms.name, label=label)


    # Define the class inside so only the test is available if xmipptomo is.
    class TestTiltSeriesTransformations(BaseTest):
        """ This class generates a phantom based workflow with contrlolled
        tilt series information from a phantom tomogram."""

        @classmethod
        def setUpClass(cls):
            setupTestProject(cls)

        def testNoAngles(self):
            # angles =1 cancels the rotation
            self._runPhantomTomo("NO ANGLES", rot=1, tilt=1, psi=1)

        def _runPhantomTomo(self, label, rot=90, tilt=90, psi=90):
            """ Creates a phantom tomo protocol with rotation by default

            :param angles: value for the 3 angles, to cancel rotation use 1
            """

            protPhantom = self.newProtocol(XmippProtPhantomTomo,
                                           objLabel=label,
                                           dimensions="200 200 100",
                                           sampling=4,
                                           nparticles=10,
                                           ntomos=3,
                                           rotmin=0,
                                           rotmax=rot,
                                           tiltmin=0,
                                           tiltmax=tilt,
                                           psimin=0,
                                           psimax=psi,
                                           addNoise=False)

            self.launchProtocol(protPhantom)

            with weakImport("imod"):
                # Project the tomogram
                from imod.protocols.protocol_tomoProjection import ProtImodTomoProjection
                tomoProjection = self.newProtocol(ProtImodTomoProjection,
                                                  inputSetOfTomograms=protPhantom.tomograms)

                self.launchProtocol(tomoProjection)

                self.addMisaligner("SHIFT X", tomoProjection.TiltSeries, shiftXNoiseToggle=True, a6param=3)
                self.addMisaligner("SHIFT Y", tomoProjection.TiltSeries, shiftYNoiseToggle=True, b6param=3)
                self.addMisaligner("JUST ANGLES", tomoProjection.TiltSeries, angleNoiseToggle=True, c6param=3)
                self.addMisaligner("XY & A ", tomoProjection.TiltSeries,
                                   shiftXNoiseToggle=True, a6param=3,
                                   shiftYNoiseToggle=True, b6param=3,
                                   angleNoiseToggle=True, c6param=3
                                   )

        def addMisaligner(self, label, inputTs, **kwargs):
            """ Adds a misalinger with the label and parameters passes

            :param label: label for the protocol
            :param kwargs: params to pass to the misalign protocol
            """



            missAligner = self.newProtocol(ProtTomoMisalignTiltSeries,
                                           objLabel=label,
                                           inputSetOfTiltSeries=inputTs,
                                           applyMatrix=True,
                                           addInverseMatrix=True,
                                           **kwargs)
            self.launchProtocol(missAligner)


def addAveragers(test, inputProt, outputName, label):
    """ Add all the averagers available to be run with the subtomo set"""

    # Add reliontomo average method
    with weakImport("reliontomo"):
        from reliontomo.protocols import ProtRelionSubTomoReconstructAvg
        relionAve = test.newProtocol(ProtRelionSubTomoReconstructAvg,
                                     objLabel="Relion rec. - %s" % label)
        relionAve.inputSubtomos.set(inputProt)
        relionAve.inputSubtomos.setExtended(outputName)
        test.launchProtocol(relionAve)

    # Add eman average method
    with weakImport("emantomo"):
        from emantomo.protocols import EmanProtSubTomoAverage
        emanAve = test.newProtocol(EmanProtSubTomoAverage,
                                   msWedge=0,
                                   objLabel = "Eman rec. - %s" % label)
        emanAve.inputSetOfSubTomogram.set(inputProt)
        emanAve.inputSetOfSubTomogram.setExtended(outputName)
        test.launchProtocol(emanAve)

    # Add xmipp tomo averager. In this case we do not need the weak import since this test will only run if xmipptomo is available
    from xmipptomo.protocols import XmippProtApplyTransformSubtomo
    xmippAve = test.newProtocol(XmippProtApplyTransformSubtomo,
                                objLabel="Xmipp rec. - %s" % label)
    xmippAve.inputSubtomograms.set(inputProt)
    xmippAve.inputSubtomograms.setExtended(outputName)
    test.launchProtocol(xmippAve)
