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
from os.path import exists, isabs, islink
from typing import List, Union
import numpy as np
from pwem import ALIGN_NONE
from pwem.objects import Transform
from pyworkflow.tests import BaseTest
from tomo.constants import TR_SCIPION, SCIPION
from tomo.objects import SetOfSubTomograms, SetOfCoordinates3D, Coordinate3D, Tomogram, CTFTomoSeries, SetOfTiltSeries, \
    TomoAcquisition, SetOfTiltSeriesM, TiltSeries, TiltImage, SetOfTomograms


class TestBaseCentralizedLayer(BaseTest):

    def checkSetGeneralProps(self, inSet, expectedSetSize: int, expectedSRate: float, streamState: int = 2):
        """
        :param inSet: A set of Scipion Tomo objects.
        :param expectedSetSize: expected set site to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param streamState: expected stream state, being 2 (default) a stream that is closed.
        """
        self.assertSetSize(inSet, expectedSetSize)
        self.assertAlmostEqual(inSet.getSamplingRate(), expectedSRate, delta=0.001)
        self.assertEqual(inSet.getStreamState(), streamState)

    # TILT SERIES ######################################################################################################
    def checkTiltSeriesM(self, inTsMSet: SetOfTiltSeriesM,  expectedSetSize: int, expectedSRate: float,
                         expectedDimensions: List[int] = None,
                         testAcqObj: TomoAcquisition = None,
                         anglesCount: Union[dict, int] = None,
                         checkIds: bool = False):
        """
        :param inTsMSet: SetOfTiltSeries.
        :param expectedSetSize: expected set site to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param expectedDimensions: list containing the expected X,Y, in pixels, and no. images, to check. A dict of
        structure {key --> tsId: value: expectedDimensions} is also admitted in the case of heterogeneous sets, e.g.
        TS with different number of tilt images.
        :param testAcqObj: TomoAcquisition object generated to test the acquisition associated to the set of tilt
        series Movies.
        :param anglesCount: expected number of tilt images in each tilt series that compose the introduced set. A dict
        of structure {key --> tsId: value: number of tilt images} is also accepted if the set is heterogeneous.
        :param checkIds: boolean used to perform id checking: objId must start with 1 --> i +1. When using {TO} the
        objId ends up with that value. For movies, it is 0, 39, 40
        """
        # CHECK THE SET ------------------------------------------------------------------------------------------------
        self.checkSetGeneralProps(inTsMSet, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        if anglesCount:
            self.checkAnglesCount(inTsMSet, anglesCount)
        if testAcqObj:
            self.checkTomoAcquisition(testAcqObj, inTsMSet.getAcquisition())
        # CHECK THE TILT SERIES ----------------------------------------------------------------------------------------
        for tsM in inTsMSet:
            tsId = tsM.getTsId()
            if testAcqObj:
                self.checkTomoAcquisition(testAcqObj, tsM.getAcquisition(), tsId=tsId)
            # Check the expected dimensions
            if expectedDimensions:
                x, y, z = tsM.getDimensions()
                if type(expectedDimensions) is dict:
                    self.assertEqual([x, y, z], expectedDimensions[tsId])
                else:
                    self.assertEqual([x, y, z], expectedDimensions)
            # Check the angles count
            if anglesCount:
                self.checkAnglesCount(tsM, anglesCount)
            # Check the filenames
            for i, ti in enumerate(tsM):
                self.assertFalse(isabs(ti.getFileName()), "Tilt image file %s is not absolute!. Should be relative.")
                self.assertTrue(islink(ti.getFileName()), "Tilt series file %s is not a link." % ti.getFileName())
                # objId must start with 1 --> i +1
                # When using {TO} the objId ends up with that value. For movies, it is 0, 39, 40
                if checkIds:
                    self.assertEqual(i + 1, ti.getObjId(), "Tilt image Movie objId is incorrect")

    def checkTiltSeries(self, inTsSet: SetOfTiltSeries, expectedSetSize: int, expectedSRate: float,
                        imported: bool = False,
                        expectedDimensions: Union[List[int], dict] = None,
                        testSetAcqObj: TomoAcquisition = None,
                        testAcqObj: Union[dict, TomoAcquisition] = None,
                        alignment: str = ALIGN_NONE,
                        isPhaseFlipped: bool = False,
                        isAmplitudeCorrected: bool = False,
                        hasAlignment: bool = False,
                        hasOddEven: bool = False,
                        anglesCount: Union[dict, int] = None,
                        anglesCountSet: int = None,
                        hasCtfCorrected: bool = False,
                        isInterpolated: bool = False,
                        excludedViewsDict: dict = None):
        """
        :param inTsSet: SetOfTiltSeries.
        :param expectedSetSize: expected set site to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param imported: boolean used to indicate if the introduced set correspond to an imported set, which means, in
        terms of tilt image transformation matrix, not to have it.
        :param expectedDimensions: list containing the expected X,Y, in pixels, and no. images, to check. A dict of
        structure {key --> tsId: value: expectedDimensions} is also admitted in the case of heterogeneous sets, e.g.
        TS with different number of tilt images.
        :param testSetAcqObj: TomoAcquisition object generated to test the acquisition associated to the set of tilt
        series. It may not be the same as testAcqObj, as in the case of heterogeneous sets of tilt series.
        :param testAcqObj: TomoAcquisition object generated to test the acquisition associated to the tilt series. A
        dict of structure {key --> tsId: value: TomoAcquisition object} is also accepted if the set is heterogeneous.
        :param alignment: Alignment type expected (see scipion-em > pwem > constants.py).
        :param isPhaseFlipped: False by default. Used to check if the tilt series were phase flipped.
        :param isAmplitudeCorrected: False by default. The same as before, but for the amplitude correction.
        :param hasAlignment: False by default. Used to indicate if the tilt series are expected to have been aligned
        (thus, the transformation matrix is not the identity) or not, and then the tilt series attribute _hasAlignment
        will be expected to be False and the tilt images will be expected not to have any transformation matrix
        associated.
        :param hasOddEven: False by default. Used to indicate if the set of tilt series are expected to have even/odd
        tilt series associated (generated with the even/odd angles or frames).
        :param anglesCount: expected number of tilt images in each tilt series that compose the introduced set. A dict
        of structure {key --> tsId: value: number of tilt images} is also accepted if the set is heterogeneous.
        :param anglesCountSet: only for heterogeneous sets. Number of angles expected in the set descrioption
        (Properties table).
        :param hasCtfCorrected: False by default. Used to indicate if the tilt series have the CTF corrected.
        :param isInterpolated: False by default. Used to indicate if the tilt series have an alignment applied (True)
        or not (False). Also, if interpolated, the expected transformation matrix will be expected to be the identity.
        :param excludedViewsDict: a dict of structure {key --> tsId: value: list of the indices of images that are
        expected to have been excluded}. An excluded view means that its corresponding attribure _objEnabled is set
        to False and the expected transformation matrix is the identity.
        """
        # TODO: check if attribute hasCtfCorrected makes sense here or if it's inherited from SPA and then does not.
        # CHECK THE SET ------------------------------------------------------------------------------------------------
        self.checkSetGeneralProps(inTsSet, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        if testSetAcqObj:
            self.checkTomoAcquisition(testSetAcqObj, inTsSet.getAcquisition())
        self.assertEqual(inTsSet.hasAlignment(), hasAlignment)
        self.assertEqual(inTsSet.getAlignment(), alignment)
        self.assertEqual(inTsSet.isPhaseFlipped(), isPhaseFlipped)
        self.assertEqual(inTsSet.isAmplitudeCorrected(), isAmplitudeCorrected)
        self.assertEqual(inTsSet.interpolated(), isInterpolated)
        self.assertEqual(inTsSet.hasOddEven(), hasOddEven)
        if anglesCountSet:
            self.checkAnglesCount(inTsSet, anglesCountSet)
        self.assertEqual(inTsSet.ctfCorrected(), hasCtfCorrected)

        # CHECK THE TILT SERIES ----------------------------------------------------------------------------------------
        for ts in inTsSet:
            tsId = ts.getTsId()
            # Check the dimensions
            if expectedDimensions:
                x, y, z = ts.getDimensions()
                if type(expectedDimensions) is dict:
                    self.assertEqual([x, y, z], expectedDimensions[tsId])
                else:
                    self.assertEqual([x, y, z], expectedDimensions)
            # Check the acquisition
            if testAcqObj:
                tsAcq = ts.getAcquisition()
                self.checkTomoAcquisition(testAcqObj, tsAcq, tsId=tsId)
            # Check the angles count
            if anglesCount:
                self.checkAnglesCount(ts, anglesCount, tsId=tsId)
            # Sampling rate
            self.assertAlmostEqual(ts.getSamplingRate(), expectedSRate, delta=0.001)
            # Alignment
            self.assertEqual(ts.hasAlignment(), hasAlignment)
            self.assertEqual(ts.getAlignment(), alignment)
            # Interpolated
            self.assertEqual(ts.interpolated(), isInterpolated)

            # CHECK THE TILT IMAGES ------------------------------------------------------------------------------------
            for ind, ti in enumerate(ts):
                # # Acquisition
                # if testAcqObj:
                #     self.checkTomoAcquisition(testAcqObj, ti.getAcquisition(), tsId=tsId)
                # Excluded view
                if excludedViewsDict:
                    isExcludedView = True if ind in excludedViewsDict[tsId] else False
                else:
                    isExcludedView = False
                self.checkObjectEnabled(ti, isExcludedView, tsId, ind)
                # Odd/Even
                if hasOddEven:
                    self.assertTrue(exists(ti.getEven()))
                    self.assertTrue(exists(ti.getOdd()))
                # Alignment matrix
                self.checkTiTransformMatrix(ti,
                                            isImported=imported,
                                            isInterpolated=isInterpolated,
                                            is2d=True,
                                            isExcludedView=isExcludedView)
                # Filename
                self.assertTrue(exists(ti.getFileName()))
                # Sampling rate
                self.assertAlmostEqual(ti.getSamplingRate(), expectedSRate, delta=0.001)

    def checkAnglesCount(self,
                         inSet: Union[SetOfTiltSeries, TiltSeries],
                         anglesCount: Union[int, dict],
                         tsId: str = None):
        """
        :param inSet: it may be a set of tilt series or a single tilt series.
        :param anglesCount: expected number of tilt images in each tilt series that compose the introduced set. A dict
        of structure {key --> tsId: value: number of tilt images} is also accepted if the set is heterogeneous.
        :param tsId: tilt series identifier. Used to get the corresponding testAcq in case it's a dict.
        """
        tsAnglesCount = anglesCount[tsId] if type(anglesCount) is dict else anglesCount
        self.assertEqual(inSet.getAnglesCount(), tsAnglesCount)

    def checkTiTransformMatrix(self, ti: TiltImage,
                               isImported: bool = False,
                               isInterpolated: bool = False,
                               isExcludedView: bool = False,
                               is2d: bool = False):
        """Checks the shape and coarsely the contents of the transformation matrix provided. Expected behavior is:
        -> If interpolated, no matrix is associated.
        -> Else:
            -> If excluded view, the expected transformation matrix is the Identity.
            -> Else: the transformation matrix exists, and it is not the identity.

        :param ti: current tilt image.
        :param isImported: boolean used to indicate if the introduced set correspond to an imported set, which means, in
        terms of tilt image transformation matrix, not to have it.
        :param is2d: False by default. Used to indicate if the expected transformation matrix should be of size 3 x 3
        (True) or 4 x 4 (False).
        :param isInterpolated: if True, the expected transformation matrix is the Identity
        :param isExcludedView: it behaves the same as param isInterpolated
        """
        if isImported or isInterpolated:
            self.assertIsNone(ti.getTransform())
        else:
            size = 3 if is2d else 4
            transfMatrixShape = (size, size)
            identityMatrix = np.eye(size)
            outMatrix = ti.getTransform().getMatrix()
            self.assertIsNotNone(outMatrix)
            if type(outMatrix) is not np.ndarray:
                outMatrix = np.array(outMatrix)
            self.assertIsNotNone(outMatrix)
            self.assertEqual(outMatrix.shape, transfMatrixShape)
            if isExcludedView:
                self.assertTrue(np.array_equal(outMatrix, identityMatrix))
            else:
                self.assertFalse(np.array_equal(outMatrix, identityMatrix))

    # TOMO ACQUISITION #################################################################################################
    def checkTomoAcquisition(self, testAcq, currentAcq, tsId=None):
        """It compares two TomoAcquisition objects, considering the following attributes:
        * Magnification.
        * Voltage.
        * Spherical aberration.
        * Initial dose.
        * Dose per frame.
        * Accumulated dose.
        * Tilt axis angle.
        * Min tilt angle.
        * Max tilt angle.
        * Tilt angle step.

        :param testAcq: reference TomoAcquisition. It can also be a dict of structure {key --> tsId: value:
        TomoAcquisition object} is also accepted if the set is heterogeneous.
        :param currentAcq: TomoAcquisition to be tested.
        :param tsId: tilt series identifier. Used to get the corresponding testAcq in case it's a dict.
        """
        testAcq = testAcq[tsId] if type(testAcq) is dict else testAcq
        self.assertAlmostEqual(testAcq.getMagnification(), currentAcq.getMagnification(), delta=1)
        self.assertAlmostEqual(testAcq.getVoltage(), currentAcq.getVoltage(), delta=1)
        self.assertAlmostEqual(testAcq.getSphericalAberration(), currentAcq.getSphericalAberration(), delta=0.01)
        self.assertAlmostEqual(testAcq.getDoseInitial(), currentAcq.getDoseInitial(), delta=0.01)
        self.assertAlmostEqual(testAcq.getDosePerFrame(), currentAcq.getDosePerFrame(), delta=0.01)
        self.assertAlmostEqual(testAcq.getAccumDose(), currentAcq.getAccumDose(), delta=0.01)
        self.assertAlmostEqual(testAcq.getTiltAxisAngle(), currentAcq.getTiltAxisAngle(), delta=0.5)
        self.assertAlmostEqual(testAcq.getAngleMin(), currentAcq.getAngleMin(), delta=0.01)
        self.assertAlmostEqual(testAcq.getAngleMax(), currentAcq.getAngleMax(), delta=0.01)
        self.assertAlmostEqual(testAcq.getStep(), currentAcq.getStep(), delta=0.1)

    # CTF ##############################################################################################################
    def checkCTFs(self, inCtfSet, expectedSetSize=-1, expectedPsdFile=None, excludedViewsDict=None, streamState=2):
        """
        :param inCtfSet: A SetOfCtfTomoSeries.
        :param expectedSetSize: expected set site to check.
        :param expectedPsdFile: boolean used to indicate if the psd file is expected to be present or not in the
        metadata stored for each CTFTomoSeries that composes SetOfCtfTomoSeries tested. In case it is expected to be
        present for some of the CTFTomoSeries, a dict of structure {key --> tsId: value --> boolean} is also admitted
        :param excludedViewsDict: a dictionary of structure {key --> tsId: value --> list of the indices of images
        that are expected to have been excluded}. An excluded view in a CTFTomo means that its attribute _objEnabled is
        set to False and the expected values for defocusU, defocusV, and DefocusAngle will be -999, -1,
        and -999, correspondingly, as what's done by the CTFModel method called setWrongDefocus.
        :param streamState: expected stream state, being 2 (default) a stream that is closed.
        """
        self.assertSetSize(inCtfSet, expectedSetSize)
        self.assertEqual(inCtfSet.getStreamState(), streamState)
        expectPsdFile = False
        for ctf in inCtfSet:
            tsId = ctf.getTsId()
            for ind, ctfi in enumerate(ctf):
                # Excluded view
                if excludedViewsDict:
                    isExcludedView = True if ind in excludedViewsDict[tsId] else False
                else:
                    isExcludedView = False
                # Expected psd file
                if type(expectedPsdFile) is dict:
                    expectPsdFile = expectedPsdFile[tsId]
                elif type(expectedPsdFile) is bool:
                    expectPsdFile = expectedPsdFile
                self.checkObjectEnabled(ctfi, isExcludedView, tsId, ind)
                self.checkCtfTomo(ctfi, isExcludedView, expectPsdFile)
        # TODO: Check if the CTFs could be checked more exhaustively

    # TOMOGRAMS ########################################################################################################
    def checkTomograms(self, inTomoSet: SetOfTomograms, expectedSetSize: int, expectedSRate: float,
                       expectedDimensions: Union[List[int], dict] = None,
                       hasOddEven: bool = False,
                       expectedOriginShifts: List[float] = None,
                       ctfCorrected: bool = False,
                       hasHalves: bool = False):
        """
        :param inTomoSet: SetOfSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param expectedDimensions: list containing the expected X,Y, and Z dimensions, in pixels, to check. A dict of
        structure {key --> tsId: value: expectedDimensions} is also admitted in the case of heterogeneous sets, e.g.
        TS with different number of tilt images.
        :param hasOddEven: False by default. Used to indicate if the set of tomograms is expected to have even/odd
        halves.
        :param expectedOriginShifts: list containing the expected shifts of the tomogram center in the X, Y, and Z
        directions.
        :param ctfCorrected: False by default: Used to indicate if the tomograms have the CTF corrected or not.
        :param hasHalves: False by default. Used to indicate if there should be halves associated to each tomogram that
        compose the set. If True, it will be checked if the corresponding halves files exist.
        """
        checkMsgPattern = 'Expected and resulting %s are different.'
        checkSizeMsg = checkMsgPattern % 'dimensions'
        checkSRateMsg = checkMsgPattern % 'sampling rate'
        checkOriginMsg = checkMsgPattern % 'origin shifts'

        # Check the set
        self.checkSetGeneralProps(inTomoSet, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        if hasOddEven:
            self.assertEqual(inTomoSet.hasOddEven, hasOddEven)
        if ctfCorrected is not None:
            self.assertEqual(inTomoSet.ctfCorrected(), ctfCorrected)
        # Check the set elements main attributes
        for tomo in inTomoSet:
            # Check if the filename exists
            self.assertTrue(exists(tomo.getFileName()))
            # Check the sampling rate
            self.assertAlmostEqual(tomo.getSamplingRate(), expectedSRate, delta=1e-3, msg=checkSRateMsg)
            # Check the ctf correction
            if ctfCorrected is not None:
                self.assertEqual(inTomoSet.ctfCorrected(), ctfCorrected)
            # At least, check that the TsId was not lost in metadata generation or copying methods
            self.assertIsNotNone(getattr(tomo, Tomogram.TS_ID_FIELD, None),
                                 msg=f'Tomogram {tomo.getFileName()}\ndoes not have attribute {Tomogram.TS_ID_FIELD} '
                                     f'or it is empty.')
            # Check the dimensions
            if expectedDimensions:
                x, y, z = tomo.getDimensions()
                if type(expectedDimensions) is dict:
                    self.assertEqual([x, y, z], expectedDimensions[tomo.getTsId()], msg=checkSizeMsg)
                else:
                    self.assertEqual([x, y, z], expectedDimensions, msg=checkSizeMsg)

            # Check the origin
            if expectedOriginShifts:
                x, y, z = tomo.getOrigin().getShifts()
                for i, j in zip([x, y, z], expectedOriginShifts):
                    self.assertAlmostEqual(i, j, delta=0.5, msg=checkOriginMsg)
            # Check the halves
            if hasHalves:
                tsId = tomo.getTsId()
                self.assertTrue(tomo.hasHalfMaps(), "Halves not registered.")
                half1, half2 = tomo.getHalfMaps().split(',')
                self.assertTrue(exists(half1), msg="Tomo %s 1st half %s does not exists" % (tsId, half1))
                self.assertTrue(exists(half2), msg="Tomo %s 2nd half %s does not exists" % (tsId, half2))

    # COORDINATES ######################################################################################################
    def checkCoordinates(self, outCoords, expectedSetSize=-1, expectedBoxSize=-1, expectedSRate=-1,
                         orientedParticles=False):
        """Checks the general properties of a SetOfCoordinates3D.

        :param outCoords: SetOf3DCoordinates.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be
        and eye matrix (False) or not (True).
        """

        # First, check the set size, sampling rate, and box size
        self.checkCoordsOrPartsSetGeneralProps(outCoords,
                                               expectedSetSize=expectedSetSize,
                                               expectedSRate=expectedSRate,
                                               expectedBoxSize=expectedBoxSize)
        for tomo in outCoords.getPrecedents():
            for coord in outCoords.iterCoordinates(volume=tomo):
                self.checkTransformMatrix(coord.getMatrix(), alignment=orientedParticles)
                self.assertEqual(coord.getTomoId(), tomo.getTsId())

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
        self.checkCoordsOrPartsSetGeneralProps(outCoords,
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
            self.checkTransformMatrix(outCoordTrMatrix, alignment=orientedParticles)
            self.checkShiftsScaling(inSetTrMatrix, outCoordTrMatrix, shiftsScaleFactor)
            # Check the tomoId
            self.assertEqual(outCoord.getTomoId(), inElement.getTomoId())

    # SUBTOMOGRAMS #####################################################################################################
    def checkImportedSubtomograms(self, subtomograms, expectedSetSize=-1, expectedBoxSize=-1, expectedSRate=-1,
                                  tomograms=None):
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
        self.checkCoordsOrPartsSetGeneralProps(subtomograms, expectedSetSize=expectedSetSize,
                                               expectedSRate=expectedSRate)
        if tomograms:
            for tomo in tomograms:
                tomoName = tomo.getFileName()
                tomoObjId = tomo.getObjId()
                tomoOrigin = tomo.getOrigin().getMatrix()
                tomoId = tomo.getTsId()
                for subtomo in subtomograms.iterSubtomos(volume=tomo):
                    self.checkAverage(subtomo, expectedSRate=expectedSRate, expectedBoxSize=expectedBoxSize,
                                      hasHalves=False)
                    self.assertEqual(subtomo.getVolName(), tomoName)
                    self.assertEqual(subtomo.getVolId(), tomoObjId)
                    self.assertTrue(np.array_equal(subtomo.getOrigin().getMatrix(), tomoOrigin))
                    self.assertEqual(subtomo.getCoordinate3D().getTomoId(), tomoId)
        else:
            for subtomo in subtomograms:
                self.checkAverage(subtomo, expectedSRate=expectedSRate, expectedBoxSize=expectedBoxSize,
                                  hasHalves=False)
                self.assertFalse(subtomo.hasCoordinate3D())

    def checkExtractedSubtomos(self, inCoords, outSubtomos, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                               convention=TR_SCIPION, orientedParticles=False, expectedExtension='.mrc',
                               isStack2d=False, noImgs=None):
        """Checks exhaustively the subtomograms generated after having carried out a subtomogram extraction

        :param inCoords: SetOf3DCoordinates introduced for the subtomo extraction.
        :param outSubtomos: the resulting SetOfSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be an
        eye matrix (False) or not (True).
        SED TO CHECK EACH SUBTOMOGRAM OF THE SET
        :param expectedExtension: '.mrc' by default. Used to indicate the expected extension of the generated particle,
        :param isStack2d: False by default. Used to indicate if the subtomogram is a stack of 2d particles or not.
        :param noImgs: Number of images of the tilt series. Only applies if isStack2d to check the dimensions, that
        are expected to be expectedBoxSize x expectedBoxSize x n <= noImgs (there may be less 2d particles than the
        noImgs because of image removal, such as in Relion 5 when using the dose threshold when extracting the subtomos.
        """
        scaleFactor = inCoords.getSamplingRate() / outSubtomos.getSamplingRate()
        # Check the critical properties of the set
        # First, check the set size, sampling rate, and box size
        self.checkCoordsOrPartsSetGeneralProps(outSubtomos,
                                               expectedSetSize=expectedSetSize,
                                               expectedSRate=expectedSRate)
        if expectedBoxSize:
            self.getSubtomosDims(outSubtomos,
                                 is2dStack=isStack2d,
                                 expectedBoxSize=expectedBoxSize,
                                 noImgs=noImgs)
        self.assertTrue(outSubtomos.hasCoordinates3D())
        # Check that the coordinates remain the same (the scaling is only applied to the shifts of the
        # transformation matrix, while the coordinates are only scaled in the coordinates extraction protocol
        # from the plugin scipion-em-tomo)
        currentCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(outSubtomos)
        unbinnedCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(inCoords)
        self.assertTrue(np.array_equal(currentCoordsExtremes, unbinnedCoordsExtremes))
        # Check the subtomograms that compose the set
        presentTsIds = inCoords.getUniqueValues(Coordinate3D.TOMO_ID_ATTR)
        presentPrecedents = [tomo.clone() for tomo in inCoords.getPrecedents() if tomo.getTsId() in presentTsIds]
        for tomo in presentPrecedents:
            for incoord, outSubtomo in zip(inCoords.iterCoordinates(volume=tomo),
                                           outSubtomos.iterSubtomos(volume=tomo)):
                self.checkSubtomogram(outSubtomo,
                                      expectedSRate=expectedSRate,
                                      expectedBoxSize=expectedBoxSize,
                                      isStack2d=isStack2d,
                                      noImgs=noImgs,
                                      expectedExtension=expectedExtension)
                subtomoTr = outSubtomo.getTransform(convention=convention)
                subtomoMatrix = subtomoTr.getMatrix()
                coordinate = outSubtomo.getCoordinate3D()
                coordTr = coordinate._eulerMatrix
                self.assertTrue(exists(outSubtomo.getFileName()))
                # The shifts in the subtomograms transformation matrix should have been scaled properly
                self.checkShiftsScaling(coordTr, subtomoTr, scaleFactor)
                # Imported coordinates were picked using PySeg, so they must have an orientation
                self.checkTransformMatrix(subtomoMatrix, alignment=orientedParticles)
                # Check the tomoId
                self.assertEqual(coordinate.getTomoId(), incoord.getTomoId())

    def checkRefinedSubtomograms(self, inSubtomos, outSubtomos, expectedSetSize=-1, expectedBoxSize=-1,
                                 expectedSRate=-1, convention=TR_SCIPION, orientedParticles=False, angTol=0.05,
                                 shiftTol=1, expectedExtension='mrc', isStack2d=False, noImgs=None):
        """Checks exhaustively the subtomograms generated after having carried out a subtomogram refinement

        :param inSubtomos: SetOfSubTomograms introduced for the subtomo refinement.
        :param outSubtomos: the resulting SetOfSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be
        and eye matrix.
        (False) or not (True).
        USED TO CHECK EACH SUBTOMOGRAM OF THE SET
        :param expectedExtension: '.mrc' by default. Used to indicate the expected extension of the generated particle,
        :param isStack2d: False by default. Used to indicate if the subtomogram is a stack of 2d particles or not.
        :param noImgs: Number of images of the tilt series. Only applies if isStack2d to check the dimensions, that
        are expected to be expectedBoxSize x expectedBoxSize x n <= noImgs (there may be less 2d particles than the
        noImgs because of image removal, such as in Relion 5 when using the dose threshold when extracting the subtomos.
        TOLERANCE PARAMETERS:
        :param angTol: angular tolerance, in degrees. Used to compare the input and output subtomograms transformation
        matrices.
        :param shiftTol: shift tolerance, in pixels. Used to compare the input and output subtomograms transformation
        matrices.
        """
        angTolMat = np.ones([3, 3]) * angTol
        shiftTolMat = np.ones(3) * shiftTol
        # Check the set
        self.checkCoordsOrPartsSetGeneralProps(outSubtomos,
                                               expectedSetSize=expectedSetSize,
                                               expectedSRate=expectedSRate)
        if expectedBoxSize:
            self.getSubtomosDims(outSubtomos,
                                 is2dStack=isStack2d,
                                 expectedBoxSize=expectedBoxSize,
                                 noImgs=noImgs)
        for inSubtomo, outSubtomo in zip(inSubtomos, outSubtomos):
            # Check the subtomogram main properties
            self.checkSubtomogram(outSubtomo,
                                  expectedSRate=expectedSRate,
                                  expectedBoxSize=expectedBoxSize,
                                  isStack2d=isStack2d,
                                  noImgs=noImgs,
                                  expectedExtension=expectedExtension)
            # Check the transformation matrix
            inSubtomoMat = inSubtomo.getTransform(convention=convention).getMatrix()
            outSubtomoMat = outSubtomo.getTransform(convention=convention).getMatrix()
            self.checkTransformMatrix(outSubtomoMat, alignment=orientedParticles)
            # The input and output matrices should be different
            diffMatrix = np.absolute(outSubtomoMat - inSubtomoMat)
            diffAngularPart = diffMatrix[:3, :3]
            diffShiftPart = diffMatrix[3, :-1]
            self.assertTrue(np.any(np.absolute(diffAngularPart - angTolMat) > 0))
            self.assertTrue(np.any(np.absolute(diffShiftPart - shiftTolMat) > 0))
            # Check that the input and output particles match (We're comparing two sets, so the convention doesn't
            # matter as long as the coordinates of both sets are retrieved using the same convention)
            inCoord = inSubtomo.getCoordinate3D()
            outCoord = outSubtomo.getCoordinate3D()
            self.assertEqual(inCoord.getPosition(SCIPION), outCoord.getPosition(SCIPION))
            self.assertEqual(inCoord.getTomoId(), outCoord.getTomoId())

    def checkSubtomogram(self, subtomo, expectedSRate=-1, expectedBoxSize=-1, expectedExtension='.mrc',
                         isStack2d=False, noImgs=None):
        """Checks the main properties of a subtomogram.

        :param subtomo: SubTomogram.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param expectedExtension: '.mrc' by default. Used to indicate the expected extension of the generated particle,
        which may vary in some cases, such as in Relion 5 2D particles (.mrcs).
        :param isStack2d: False by default. Used to indicate if the subtomogram is a stack of 2d particles or not.
        :param noImgs: Number of images of the tilt series. Only applies if isStack2d to check the dimensions, that
        are expected to be expectedBoxSize x expectedBoxSize x n <= noImgs (there may be less 2d particles than the
        noImgs because of image removal, such as in Relion 5 when using the dose threshold when extracting the subtomos.
        """
        self.assertTrue(exists(subtomo.getFileName()), "Average %s does not exists" % subtomo.getFileName())
        self.assertTrue(subtomo.getFileName().endswith(expectedExtension))
        self.assertTrue(abs(subtomo.getSamplingRate() - expectedSRate) <= 1e-4)
        self.getSubtomosDims(subtomo,
                             is2dStack=isStack2d,
                             expectedBoxSize=expectedBoxSize,
                             noImgs=noImgs)

    def getSubtomosDims(self, inObj, is2dStack=None, expectedBoxSize=None, noImgs=None):
        """It can be used with a sigle SubTomogram or a set of them."""
        if is2dStack:
            x, y, nIm = inObj.getDimensions()
            self.assertEqual(x, expectedBoxSize)
            self.assertEqual(y, expectedBoxSize)
            self.assertLessEqual(nIm, noImgs)
        else:
            testBoxSize = (expectedBoxSize, expectedBoxSize, expectedBoxSize)
            self.assertEqual(inObj.getDimensions(), testBoxSize)

    # AVERAGE ##########################################################################################################
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
        self.assertTrue(abs(avg.getSamplingRate() - expectedSRate) <= 1e-4)
        self.assertEqual(avg.getDimensions(), testBoxSize)
        # Check the halves
        if hasHalves:
            self.assertTrue(avg.hasHalfMaps(), "Halves not registered.")
            half1, half2 = avg.getHalfMaps().split(',')
            self.assertTrue(exists(half1), msg="Average 1st half %s does not exists" % half1)
            self.assertTrue(exists(half2), msg="Average 2nd half %s does not exists" % half2)

    # SUBTOMOGRAM CLASSES ##############################################################################################
    def checkClasses(self, classesSet, expectedSRate=-1, expectedSetSize=-1, hasRepresentatives=True):
        """Checks exhaustively the classes generated after having carried out a subtomogram classification
        :param classesSet: SetOfClassesSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param hasRepresentatives: flag used to indicate if the set of classes is expected to have representatives or
        not.
        """
        self.checkSetGeneralProps(classesSet, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        self.assertEqual(classesSet.hasRepresentatives(), hasRepresentatives)
        representativeFileList = []
        for subtomoClass in classesSet:
            if hasRepresentatives:
                repFileName = subtomoClass.getRepresentative().getFileName()
                representativeFileList.append(repFileName)
                self.assertTrue(exists(repFileName))
            self.assertEqual(subtomoClass.getSamplingRate(), expectedSRate)

        if hasRepresentatives:
            msg = 'At least one of the representative filenames is repeated, which should not be possible.'
            self.assertEqual(len(set(representativeFileList)), expectedSetSize, msg=msg)

    # COORDINATES AND PARTICLES TEST UTILS #############################################################################
    def checkCoordsOrPartsSetGeneralProps(self, inSet, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=0):
        # Check the set
        self.checkSetGeneralProps(inSet, expectedSetSize=expectedSetSize, expectedSRate=expectedSRate)
        if expectedBoxSize:
            self.assertEqual(inSet.getBoxSize(), expectedBoxSize)

    def checkTransformMatrix(self, outMatrix, alignment=False, is2d=False, isInterpolated=False, isExcludedView=False):
        """Checks the shape and coarsely the contents of the transformation matrix provided.

        :param outMatrix: transformation matrix of a subtomogram or coordinate.
        :param alignment: False by default. Used to specify if the expected transformation matrix should be an
        eye matrix (False) or not (True).
        :param is2d: False by default. Used to indicate if the expected transformation matrix should be of size 3 x 3
        (True) or 4 x 4 (False).
        :param isInterpolated: if True, the expected transformation matrix is the Identity
        :param isExcludedView: it behaves the same as param isInterpolated
        """
        size = 3 if is2d else 4
        transfMatrixShape = (size, size)
        identityMatrix = np.eye(size)
        self.assertIsNotNone(outMatrix)
        if type(outMatrix) is not np.ndarray:
            outMatrix = np.array(outMatrix)
        self.assertIsNotNone(outMatrix)
        self.assertEqual(outMatrix.shape, transfMatrixShape)
        if isInterpolated or isExcludedView:
            self.assertTrue(np.array_equal(outMatrix, identityMatrix))
        else:
            if alignment:
                self.assertFalse(np.array_equal(outMatrix, identityMatrix))
            else:
                self.assertTrue(np.array_equal(outMatrix, identityMatrix))

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
        if type(inSet) is not SetOfCoordinates3D:
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

    def checkObjectEnabled(self, obj, isExcludedView, tsId, ind):
        enb = obj.isEnabled()
        objType = 'CTF' if type(obj) is CTFTomoSeries else 'Tilt image'
        if isExcludedView:
            self.assertFalse(enb, msg='TsId = %s: %s %i was expected not to be Enabled' % (tsId, objType, ind))
        else:
            self.assertTrue(enb, msg='TsId = %s: %s %i was expected to be Enabled' % (tsId, objType, ind))

    def checkCtfTomo(self, ctf, isExcluded, expectPsdFile):
        defocusU = ctf.getDefocusU()
        defocusV = ctf.getDefocusV()
        defocusAngle = ctf.getDefocusAngle()
        defocusVals = [defocusU, defocusV, defocusAngle]
        self.assertGreater(ctf.getAcquisitionOrder(), -1)
        if isExcluded:
            expectedDefocusU = -999
            expectedDefocusV = -1
            expectedDefocusAngle = -999
            self.assertTrue(np.allclose(np.array(defocusVals),
                                        np.array([expectedDefocusU, expectedDefocusV, expectedDefocusAngle])))
        else:
            # Check the defocusU, defocusV and defocus angle values
            for val in defocusVals:
                self.assertGreaterEqual(val, 0)
            # Check the psd file
            if expectPsdFile:
                self.assertTrue(exists(ctf.getPsdFile().rsplit('@', 1)[-1]))  # Expected syntax is index@psdFile
