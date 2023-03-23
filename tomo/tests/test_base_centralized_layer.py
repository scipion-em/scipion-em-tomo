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
from tomo.constants import TR_SCIPION, SCIPION
from tomo.objects import SetOfSubTomograms


class TestBaseCentralizedLayer(BaseTest):

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

    def checkSetGeneralProps(self, inSet, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=0):
        # Check the set
        self.assertSetSize(inSet, expectedSetSize)
        self.assertEqual(inSet.getSamplingRate(), expectedSRate)
        if expectedBoxSize:
            self.assertEqual(inSet.getBoxSize(), expectedBoxSize)

    def checkExtracted3dCoordinates(self, inSet, outCoords, expectedSetSize=-1, expectedBoxSize=-1,
                                    expectedSRate=-1, convention=TR_SCIPION, orientedParticles=False):
        """Checks the results of a coordinate extraction protocol.

        :param inSet: input set from which the coordinates were extracted. It can be a SetOf3DCoordinates or a
        SetOfSubTomograms.
        :param outCoords: the resulting SetOf3DCoordinates after the coordinate extraction.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py.
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be
        and eye matrix (False) or not (True)."""
        if type(inSet) == SetOfSubTomograms:
            inSet = inSet.getCoordinates3D()
        # First, check the set size, sampling rate, and box size
        self.checkSetGeneralProps(outCoords,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate,
                                  expectedBoxSize=expectedBoxSize)
        # Check the coordinate extremes
        inCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(inSet)
        outCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(outCoords)
        coordScaleFactor = inSet.getSamplingRate() / outCoords.getSamplingRate()
        shiftsScaleFactor = 1 / coordScaleFactor
        self.assertTrue(np.array_equal(outCoordsExtremes, coordScaleFactor * inCoordsExtremes))
        # Other checks per coordinate
        for inElement, outCoord in zip(inSet, outCoords):
            # Check the transformation matrices and shifts
            inSetTrMatrix = inElement.getMatrix(convention=convention)
            outCoordTrMatrix = outCoord.getMatrix(convention=convention)
            self.check3dTransformMatrix(outCoordTrMatrix, orientedParticles=orientedParticles)
            self.checkShiftsScaling(inSetTrMatrix, outCoordTrMatrix, shiftsScaleFactor)
            # Check the tomoId
            self.assertEqual(outCoord.getTomoId(), inElement.getTomoId())

    def checkAverage(self, avg, expectedSRate=-1, expectedBoxSize=-1, hasHalves=True):
        """Checks the main properties of an average subtomogram, which can be the result of an average of subtomograms,
        an initial model or refinement of subtomograms.

        :param avg: AverageSubtomogram.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param hasHalves: True by default. Used to indicate if the average is expected to have halves associated."""
        testBoxSize = (expectedBoxSize, expectedBoxSize, expectedBoxSize)
        self.assertTrue(exists(avg.getFileName()), "Average %s does not exists" % avg.getFileName())
        self.assertTrue(avg.getFileName().endswith(".mrc"))
        # The imported coordinates correspond to a binned 2 tomogram
        self.assertEqual(avg.getSamplingRate(), expectedSRate)
        self.assertEqual(avg.getDimensions(), testBoxSize)
        # Check the halves
        if hasHalves:
            self.assertTrue(avg.hasHalfMaps(), "Halves not registered.")
            half1, half2 = avg.getHalfMaps().split(',')
            self.assertTrue(exists(half1), msg="Average 1st half %s does not exists" % half1)
            self.assertTrue(exists(half2), msg="Average 2nd half %s does not exists" % half2)

    def checkImportedSubtomograms(self, subtomograms, expectedSetSize=-1, expectedBoxSize=-1, expectedSRate=-1, tomograms=None):
        """Checks the main properties of an imported set subtomograms.

        :param subtomograms: the resulting SetOfSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param tomograms: SetOfTomograms. If provided, additional features regarding how each subtomogram is referred
        to the corresponding tomogram will be carried out. If not provided, it will be considered that the subtomograms
        are not referred to any tomograms (as in some subtomogram importing cases).
        """
        # Check the set
        self.checkSetGeneralProps(subtomograms, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        if tomograms:
            for tomo in tomograms:
                tomoName = tomo.getFileName()
                tomoObjId = tomo.getObjId()
                tomoOrigin = tomo.getOrigin().getMatrix()
                tomoId = tomo.getTsId()
                for subtomo in subtomograms.iterSubtomos(volume=tomo):
                    self.checkAverage(subtomo, expectedSRate=expectedSRate, expectedBoxSize=expectedBoxSize, hasHalves=False)
                    self.assertEqual(subtomo.getVolName(), tomoName)
                    self.assertEqual(subtomo.getVolId(), tomoObjId)
                    self.assertTrue(np.array_equal(subtomo.getOrigin().getMatrix(), tomoOrigin))
                    self.assertEqual(subtomo.getCoordinate3D().getTomoId(), tomoId)
        else:
            for subtomo in subtomograms:
                self.checkAverage(subtomo, expectedSRate=expectedSRate, expectedBoxSize=expectedBoxSize, hasHalves=False)
                self.assertFalse(subtomo.hasCoordinate3D())

    def checkExtractedSubtomos(self, inCoords, outSubtomos, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                               convention=TR_SCIPION, orientedParticles=False):
        """Checks exhaustively the subtomograms generated after having carried out a subtomogram extraction

        :param inCoords: SetOf3DCoordinates introduced for the subtomo extraction.
        :param outSubtomos: the resulting SetOfSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be and eye matrix
        (False) or not (True)."""
        scaleFactor = inCoords.getSamplingRate() / outSubtomos.getSamplingRate()
        # Check the critical properties of the set
        # First, check the set size, sampling rate, and box size
        self.checkSetGeneralProps(outSubtomos,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate)
        self.assertEqual(outSubtomos.getDimensions(), (expectedBoxSize, expectedBoxSize, expectedBoxSize))
        self.assertTrue(outSubtomos.hasCoordinates3D())
        # Check that the coordinates remain the same (the scaling is only applied to the shifts of the
        # transformation matrix, while the coordinates are only scaled in the coordinates extraction protocol
        # from the plugin scipion-em-tomo)
        currentCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(outSubtomos)
        unbinnedCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(inCoords)
        self.assertTrue(np.array_equal(currentCoordsExtremes, unbinnedCoordsExtremes))
        # Check the subtomograms that compose the set
        for incoord, outSubtomo in zip(inCoords, outSubtomos):
            subtomoTr = outSubtomo.getTransform(convention=convention)
            subtomoMatrix = subtomoTr.getMatrix()
            coordinate = outSubtomo.getCoordinate3D()
            coordTr = coordinate._eulerMatrix
            coordMatrix = coordinate.getMatrix(convention=convention)
            self.assertTrue(exists(outSubtomo.getFileName()))
            self.assertEqual(outSubtomo.getSamplingRate(), expectedSRate)
            # The shifts in the subtomograms transformation matrix should have been scaled properly
            self.checkShiftsScaling(coordTr, subtomoTr, scaleFactor)
            # Imported coordinates were picked using PySeg, so they must have an orientation
            self.check3dTransformMatrix(subtomoMatrix, orientedParticles=orientedParticles)
            # Check the tomoId
            self.assertEqual(coordinate.getTomoId(), incoord.getTomoId())

    def checkRefinedSubtomograms(self, inSubtomos, outSubtomos, expectedSetSize=-1, expectedBoxSize=-1, expectedSRate=-1,
                                 convention=TR_SCIPION, orientedParticles=False, angTol=0.05, shiftTol=1):
        angTolMat = np.ones([3, 3]) * angTol
        shiftTolMat = np.ones(3) * shiftTol
        # Check the set
        self.checkSetGeneralProps(outSubtomos, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        for inSubtomo, outSubtomo in zip(inSubtomos, outSubtomos):
            # Check the subtomogram main properties
            self.checkAverage(outSubtomo, expectedSRate=expectedSRate, expectedBoxSize=expectedBoxSize, hasHalves=False)
            # Check the transformation matrix
            inSubtomoMat = inSubtomo.getTransform(convention=convention).getMatrix()
            outSubtomoMat = outSubtomo.getTransform(convention=convention).getMatrix()
            self.check3dTransformMatrix(outSubtomoMat, orientedParticles=orientedParticles)
            # The input and output matrices should be different
            diffMatrix = np.absolute(outSubtomoMat - inSubtomoMat)
            diffAngularPart = diffMatrix[:3, :3]
            diffShiftPart = diffMatrix[3, :-1]
            self.assertTrue(np.any(np.absolute(diffAngularPart - angTolMat) > 0))
            self.assertTrue(np.any(np.absolute(diffShiftPart - shiftTolMat) > 0))
            # Check that the input and output particles match (We're comparing tow sets, so the convention doesn't
            # matter as long as the coordinates of both sets are retrieved using the same convention)
            inCoord = inSubtomo.getCoordinate3D()
            outCoord = outSubtomo.getCoordinate3D()
            self.assertEqual(inCoord.getPosition(SCIPION), outCoord.getPosition(SCIPION))
            self.assertEqual(inCoord.getTomoId(), outCoord.getTomoId())
