# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es) [1]
# *
# * [1] Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
import tempfile

from pyworkflow.tests import BaseTest
from tomo.constants import SCIPION
from tomo.objects import (SetOfTiltSeriesCoordinates, TiltSeriesCoordinate,
                          SetOfSubTomograms, SetOfTomograms, Tomogram,
                          SetOfCoordinates3D, Coordinate3D, SubTomogram,
                          SetOfTiltSeries, TiltSeries, TiltImage, LandmarkModel,
                          CTFTomo)

TS_1 = "TS_1"
TS_2 = "TS_2"

MYTSID = "MYTSID"

SAMPLING_RATE = 10

Z_VALUE = 30

Y_VALUE = 20

X_VALUE = 10


class TestTomoModel(BaseTest):
    """ This class check if the Object model behaves as expected"""

    @classmethod
    def setUpClass(cls):
        cls.setupTestOutput()

    def test_set_of_tilt_series_coordinates(self):
        """ Tests the SetOfTiltSeries model"""

        tsCoord = TiltSeriesCoordinate()

        tsCoord.setX(X_VALUE)
        tsCoord.setY(Y_VALUE)
        tsCoord.setZ(Z_VALUE)

        self.assertEqual(X_VALUE, tsCoord.getX(), "X is wrong")
        self.assertEqual(Y_VALUE, tsCoord.getY(), "Y is wrong")
        self.assertEqual(Z_VALUE, tsCoord.getZ(), "Z is wrong")

        posAngs = tsCoord.getPosition()
        self.assertEqual(X_VALUE, posAngs[0], "X is wrong")
        self.assertEqual(Y_VALUE, posAngs[1], "Y is wrong")
        self.assertEqual(Z_VALUE, posAngs[2], "Z is wrong")

        posPixels = tsCoord.getPosition(sampling_rate=SAMPLING_RATE)
        self.assertEqual(X_VALUE / SAMPLING_RATE, posPixels[0], "X is wrong in pixels")
        self.assertEqual(Y_VALUE / SAMPLING_RATE, posPixels[1], "Y is wrong in pixels")
        self.assertEqual(Z_VALUE / SAMPLING_RATE, posPixels[2], "Z is wrong in pixels")

        tsCoord.setPosition(X_VALUE / SAMPLING_RATE, Y_VALUE / SAMPLING_RATE,
                            Z_VALUE / SAMPLING_RATE, sampling_rate=SAMPLING_RATE)
        posAngs = tsCoord.getPosition()
        self.assertEqual(X_VALUE, posAngs[0], "X is wrong when set with sampling rate")
        self.assertEqual(Y_VALUE, posAngs[1], "Y is wrong when set with sampling rate")
        self.assertEqual(Z_VALUE, posAngs[2], "Z is wrong when set with sampling rate")

        tsCoord.setTsId(MYTSID)
        self.assertEqual(MYTSID, tsCoord.getTsId(), "TSID is wrong")

        set = SetOfTiltSeriesCoordinates.create(self.outputPath)
        set.append(tsCoord)

        for tsCoord in set:
            posAngs = tsCoord.getPosition()
            self.assertEqual(X_VALUE, posAngs[0], "X is wrong when set with sampling rate")
            self.assertEqual(Y_VALUE, posAngs[1], "Y is wrong when set with sampling rate")
            self.assertEqual(Z_VALUE, posAngs[2], "Z is wrong when set with sampling rate")
            self.assertEqual(MYTSID, tsCoord.getTsId(), "TSID is wrong")

    def test_set_of_subtomograms(self):
        """ Tests the SetOfSubtomograms model"""

        # Create tomograms set
        tomos = SetOfTomograms.create(self.outputPath)

        tomo1 = Tomogram()
        tomo1.setTsId(TS_1)
        tomos.append(tomo1)

        tomo2 = Tomogram()
        tomo2.setTsId("TS_2")
        tomos.append(tomo2)

        # Create Coordinates set with just one coordinate pointing to just one tomogram
        coords = SetOfCoordinates3D.create(self.outputPath)

        coord1 = Coordinate3D()
        coord1.setX(0, SCIPION)
        coord1.setY(0, SCIPION)
        coord1.setZ(0, SCIPION)
        coord1.setTomoId(TS_1)
        coords.append(coord1)

        coords.setPrecedents(tomos)

        coords.write()

        # Test precedents involved
        self.assertEqual(len(coords.getPrecedentsInvolved()), 1, "getPrecedentsInvolved does not seem to work")

        self.assertEqual(coords.getTSIds(), [TS_1], "SetOfCoordinates.getTSIds not working.")

        subtomos = SetOfSubTomograms.create(self.outputPath)
        subtomos.setCoordinates3D(coords)

        subtomo = SubTomogram()
        subtomo.setCoordinate3D(coord1)
        subtomos.append(subtomo)

        subtomos.write()

        self.assertEqual(len(subtomos.getTomograms()), 1, "getTomograms does not return the right element count")

        tomo = subtomos.getTomogram(subtomo)
        self.assertEqual(tomo.getTsId(), TS_1, "Recovered tomogram from SOST.getTomogram() does not work.")

    def test_set_of_ctfTomo_series(self):
        # test enabled is cloned
        newCtf = CTFTomo()

        newCtf.setEnabled(False)

        clonedCTF = newCtf.clone()
        self.assertFalse(clonedCTF.isEnabled(), "Enabled not cloned for the tomo CTF")

    def test_set_of_tilt_series(self):
        """ Tests the SetOfTiltSeries model"""

        # Create the set
        tiltseries = SetOfTiltSeries.create(self.outputPath)
        tiltseries.setAnglesCount(3)

        def addTiltSerie(tsId):
            ts = TiltSeries()
            ts.setTsId(tsId)

            # We need to append the tilt series before adding tilt images
            tiltseries.append(ts)

            # Add tilt images
            ti = TiltImage()
            ti.setTiltAngle(-3)
            ti.setTsId(tsId)
            ts.append(ti)

            # Update properties
            tiltseries.update(ts)

            # This should not persist properties at all
            ts.write(properties=True)

        # Add tilt series
        addTiltSerie(TS_1)
        addTiltSerie(TS_2)

        # Persist the whole set, this should persist properties of the set
        tiltseries.write()

        tiltSerieFromFile = SetOfTiltSeries(filename=tiltseries.getFileName())
        tiltSerieFromFile.loadAllProperties()

        self.assertEqual(tiltSerieFromFile.getAnglesCount(), 1, "SetOfTiltSeries.getAnglesCount is wrong.")
        self.assertEqual(tiltSerieFromFile.getTSIds(), [TS_1, TS_2], "SetOfTiltSeries.getTSIds not working.")

        # Accessing through [] and mapper sync

        ts1 = tiltseries[1]
        ti1 = ts1.getFirstItem()
        self.assertEqual(ti1.getTsId(), TS_1, "Tilt image wrong for the tilt series 1")

        ts2 = tiltseries[2]
        ti2 = ts2.getFirstItem()
        self.assertEqual(ti2.getTsId(), TS_2, "Tilt image wrong for the tilt series 2")

        # test enabled is cloned
        newTi= TiltImage()

        newTi.setEnabled(False)

        clonedTi = newTi.clone()
        self.assertFalse(clonedTi.isEnabled(), "Enabled not cloned for the tilt image")

    def test_landmarks(self):
        """ Test the Landmark model"""

        temp_name = next(tempfile._get_candidate_names())
        # Verify the count
        lm = LandmarkModel(fileName=temp_name)

        str(lm)  # Should not fail

        self.assertEqual(0, lm.getCount(), "Count not initialized to 0")

        # Chain id 4
        lm.addLandmark(1, 2, 3, 4, None, None)
        self.assertEqual(1, lm.getCount(), "Count not increased.")

        # same chain id
        lm.addLandmark(1, 2, 3, 4, None, None)
        self.assertEqual(1, lm.getCount(), "Count has been increased?")

        # Add a new chainId 5
        lm.addLandmark(1, 2, 3, 5, None, None)
        self.assertEqual(2, lm.getCount(), "Count not increased when not empty.")

        str(lm)
