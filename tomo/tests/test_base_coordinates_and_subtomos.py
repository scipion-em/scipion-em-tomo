# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from os.path import exists

from flatbuffers.builder import np

from pwem.objects import Transform
from pyworkflow.tests import BaseTest
from tomo.constants import TR_SCIPION
from tomo.objects import SetOfSubTomograms


class TestUtilsSetsOfCoordsAndSubtomos(BaseTest):
    """Test class that contains auxiliary methods to exhaustively check the results of operations carried out
    with SetOf3DCoordinates or SetOfSubTomograms"""

    def check3dTransformMatrix(self, outMatrix, orientedParticles=False):
        """Checks the shape and coarsely the contents of the transformation matrix provided.

        :param outMatrix: transformation matrix of a subtomogram or coordinate.
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be an
        eye matrix
        (False) or not (True).
        """
        transfMatrixShape = (4, 4)
        self.assertIsNotNone(outMatrix)
        if type(outMatrix) != np.ndarray:
            outMatrix = np.array(outMatrix)
        self.assertIsNotNone(outMatrix)
        self.assertTrue(outMatrix.shape, transfMatrixShape)
        if orientedParticles:
            self.assertFalse(np.array_equal(outMatrix, np.eye(4)))
        else:
            self.assertTrue(np.array_equal(outMatrix, np.eye(4)))

    def checkShiftsScaling(self, inTransform, outTransform, scaleFactor):
        """Check if the shifts were scaled properly in the case of subtomogram extraction from a different
        tomo source than the tomograms in which the picking was carried out.

        :param inTransform: Transform object or transformation matrix attribute of a coordinate before the extraction.
        :param outTransform: Transform object or transformation matrix attribute of the corresponding extracted
        subtomogram.
        :param scaleFactor: Scale factor expected between the shifts corresponding to the outTransform respecting the
        inputTransform.
        """
        sx, sy, sz = self.getShiftsFromTransformOrMatrix(inTransform)
        osx, osy, osz = self.getShiftsFromTransformOrMatrix(outTransform)
        for inShift, outShift in zip([sx, sy, sz], [osx, osy, osz]):
            self.assertEqual(outShift, scaleFactor * inShift)

    @staticmethod
    def getShiftsFromTransformOrMatrix(inTransform):
        """Returns the shifts from a Transformation object or a transformation matrix.

        :param inTransform: it can be a Transformation object or a transformation matrix"""
        if isinstance(inTransform, Transform):
            return inTransform.getShifts()
        else:
            return inTransform[0, 3], inTransform[1, 3], inTransform[2, 3]

    @staticmethod
    def getMinAndMaxCoordValuesFromSet(inSet):
        """Get the extreme values of the coordinates from an introduced set. The output is a numpy array 
         with the following values [x_min, x_max, y_min, y_max, z_min, z_max]. It's very useful to compare 
         sets of subtomograms or coordinates. The coordinates are expected to be referred to the center of
         the tomogram (Scipion-like)
         
         :param inSet: it can be a SetOf3DCoordinates or a SetOfSubTomograms."""
        if type(inSet) == SetOfSubTomograms:
            inSet = inSet.getCoordinates3D()

        dataDict = inSet.aggregate(['MAX'], '_tomoId', ['_x', '_y', '_z'])
        xcoords, ycoords, zcoords = zip(*[(d['_x'], d['_y'], d['_z']) for d in dataDict])
        return np.array([min(xcoords),
                         max(xcoords),
                         min(ycoords),
                         max(ycoords),
                         min(zcoords),
                         max(zcoords)])

    def compareCoordinates(self, set1, set2, scaleFactor=1):
        """Determines if two sets of coordinates or subtomograms are the same by getting and comparing their
        extreme coordinate values.

        :param set1: it can be a SetOf3DCoordinates or a SetOfSubTomograms.
        :param set2: a set of the same type as set1.
        :param scaleFactor: 1 by default. Used to indicate the scaling factor of the set2 respecting the set1."""
        extremes1 = self.getMinAndMaxCoordValuesFromSet(set1)
        extremes2 = self.getMinAndMaxCoordValuesFromSet(set2)
        self.assertTrue(np.array_equal(extremes1, extremes2 * scaleFactor))

    def areInTomogram(self, inSet, tomogram):
        """Checks if the extreme coordinates that corresponds to the provided set are contained in the range
        [0, tomogramDimension] for x, y and z, considering the coordinates referred to the center of the tomogram
        (Scipion-like) as the method getMinAndMaxCoordValuesFromSet used to get the coordinate extreme values operates
        directly on the sqlite file corresponding to the inSet.

         :param inSet: it can be a SetOf3DCoordinates or a SetOfSubTomograms.
         :param tomogram: a Tomogram object."""
        xt, yt, zt = tomogram.getDim()
        extremes = self.getMinAndMaxCoordValuesFromSet(inSet)
        xc_max = extremes[1]
        yc_max = extremes[3]
        zc_max = extremes[5]
        self.assertTrue(xc_max < xt / 2 and yc_max < yt / 2 and zc_max < zt / 2)


class TestUtilsExtractCoords(TestUtilsSetsOfCoordsAndSubtomos):
    """Test class that contains auxiliary methods to exhaustively check the results of a coordinates
    extraction protocol"""

    def checkExtracted3dCoordinates(self, inSet, outCoords, expectedSetSize=-1, expectedBoxSize=-1,
                                    tomoId=None, expectedSRate=-1, convention=TR_SCIPION, orientedParticles=False):
        if type(inSet) == SetOfSubTomograms:
            inSet = inSet.getCoordinates3D()
        inCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(inSet)
        outCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(outCoords)
        coordScaleFactor = inSet.getSamplingRate() / outCoords.getSamplingRate()
        shiftsScaleFactor = 1 / coordScaleFactor
        # Check the coordinate extremes
        self.assertTrue(np.array_equal(outCoordsExtremes, coordScaleFactor * inCoordsExtremes))
        # Check the sampling rate
        self.assertEqual(outCoords.getSamplingRate(), expectedSRate)
        # Check the set size
        self.assertSetSize(outCoords, expectedSetSize)
        # Check the box size
        self.assertEqual(outCoords.getBoxSize(), expectedBoxSize)
        # Other checks per coordinate
        for inElement, outCoord in zip(inSet, outCoords):
            # Check the transformation matrices and shifts
            inSetTrMatrix = inElement.getMatrix(convention=convention)
            outCoordTrMatrix = outCoord.getMatrix(convention=convention)
            super().check3dTransformMatrix(outCoordTrMatrix, orientedParticles=orientedParticles)
            super().checkShiftsScaling(inSetTrMatrix, outCoordTrMatrix, shiftsScaleFactor)
            # Check the tomoId
            self.assertEqual(outCoord.getTomoId(), tomoId)


class TestUtilsExtractSubtomos(TestUtilsSetsOfCoordsAndSubtomos):
    """Test class that contains auxiliary methods to exhaustively check the results of a sutomogram
    extraction protocol"""

    def checkExtractedSubtomos(self, inCoords, outSubtomos, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                               convention=TR_SCIPION, orientedParticles=False, expectedTomoId=None):
        """Checks exhaustively the subtomograms generated after having carried out a subtomogram extraction

        :param inCoords: SetOf3DCoordinates introduced for the subtomo extraction.
        :param outSubtomos: the resulting SetOfSubTomograms.
        :param expectedSetSize: expected box size to compare.
        :param expectedSRate: expected sampling rate to compare.
        :param expectedBoxSize: expected box size to compare.
        :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be and eye matrix
        (False) or not (True).
        :param expectedTomoId: expected tomoId to compare.
        :"""
        scaleFactor = outSubtomos.getSamplingRate() / inCoords.getSamplingRate()
        # Check the critical properties of the set
        self.assertSetSize(outSubtomos, expectedSetSize)
        self.assertEqual(outSubtomos.getSamplingRate(), expectedSRate)
        self.assertEqual(outSubtomos.getDimensions(), (expectedBoxSize, expectedBoxSize, expectedBoxSize))
        self.assertTrue(outSubtomos.hasCoordinates3D())
        # Check the subtomograms that compose the set
        for subtomo in outSubtomos:
            subtomoTr = subtomo.getTransform(convention=convention)
            subtomoMatrix = subtomoTr.getMatrix()
            coordinate = subtomo.getCoordinate3D()
            coordTr = coordinate._eulerMatrix
            coordMatrix = coordinate.getMatrix(convention=convention)
            self.assertTrue(exists(subtomo.getFileName()))
            self.assertEqual(subtomo.getSamplingRate(), expectedSRate)
            # The shifts in the subtomograms transformation matrix should have been scaled properly
            super().checkShiftsScaling(coordTr, subtomoTr, scaleFactor)
            # Imported coordinates were picked using PySeg, so they must have an orientation
            super().check3dTransformMatrix(subtomoMatrix, orientedParticles=orientedParticles)
            # Also, at this point the transformation matrix should be the same as the coordinate matrix as the angles
            # have not been refined yet
            super().check3dTransformMatrix(coordMatrix, orientedParticles=orientedParticles)
            self.assertTrue(np.array_equal(subtomoMatrix, coordMatrix))
            # Check the tomoId
            self.assertEqual(coordinate.getTomoId(), expectedTomoId)

        # Check that the coordinates remain the same (the scaling is only applied to the shifts of the
        # transformation matrix, while the coordinates are only scaled in the coordinates extraction protocol
        # from the plugin scipion-em-tomo
        currentCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(outSubtomos)
        unbinnedCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(inCoords)
        self.assertTrue(np.array_equal(currentCoordsExtremes, unbinnedCoordsExtremes))


class TestUtilsAverageOfSubtomos(BaseTest):
    """Test class that contains auxiliary methods to exhaustively check the results of an average of
    subtomograms protocol"""

    def checkAverage(self, avg, expectedSRate=-1, expectedBoxSize=-1, hasHalves=True):
        testBoxSize = (expectedBoxSize, expectedBoxSize, expectedBoxSize)
        self.assertTrue(exists(avg.getFileName()), "Average %s does not exists" % avg.getFileName())
        self.assertTrue(avg.getFileName().endswith(".mrc"))
        # The imported coordinates correspond to a binned 2 tomogram
        self.assertEqual(avg.getSamplingRate(), expectedSRate)
        self.assertEqual(avg.getDimensions(), testBoxSize)
        # Check the halves
        if hasHalves:
            self.assertTrue(avg.hasHalfMaps(), "Halves not registered")
            half1, half2 = avg.getHalfMaps().split(',')
            self.assertTrue(exists(half1), msg="Average 1st half %s does not exists" % half1)
            self.assertTrue(exists(half2), msg="Average 2nd half %s does not exists" % half2)

