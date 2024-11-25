# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import math
from enum import Enum
from typing import Tuple

from pyworkflow.tests import DataSet
from tomo.objects import TomoAcquisition

DataSet(name='tomo-em', folder='tomo-em',
        files={
            'tomo1': 'overview_wbp.em',
            'tomo2': 'overview_wbp2.em',
            'tomo3': 'tomo_8_mn.mrc',
            'subtomo': 'basename.hdf',
            'eman_coordinates': 'coordinates3Deman2',
            'etomo': 'tutorialData',
            'empiar': 'EMPIAR-10164',
            'tsMParentFolder': 'ts_tsM_and_mdocs',
            'tsM10Dir': 'ts_tsM_and_mdocs/Tomo_10',
            'tsM31Dir': 'ts_tsM_and_mdocs/Tomo_31',
            'empiarMdocDirOk': 'ts_tsM_and_mdocs/mdocs/realFromEmpiar/complete',
            'empiarMdocDirNoOk': 'ts_tsM_and_mdocs/mdocs/realFromEmpiar/incomplete',
            'realFileNoVoltage1': 'tomo4_delay.st.mdoc',
            'realFileNoVoltage2': 'TS_54.mrc.mdoc',
            'simErrorMdocDir': 'ts_tsM_and_mdocs/mdocs/editedForErrorSimulation',
            'noMaginficationMdoc': 'NoMagnification.mdoc',
            'noSamplingRateMdoc': 'NoPixelSpacing.mdoc',
            'noVoltagenoSRateMdoc': 'NoVoltage_NoPixelSpacing.mdoc',
            'someMissingAnglesMdoc': 'SomeTiltAnglesMissing_1_7_48.mdoc',
            'noDoseMdoc': 'noDose.mdoc',
            'ts1_imod': 'tutorialData/BBa.st',
            'ts2_imod': 'tutorialData/BBb.st'
        })

DataSet(name='reliontomo', folder='reliontomo',
        files={
            'tomo1': '64K_defocus_m2_tomo_10_bin1_WBP_CatBinned1.mrc',
            'tomo2': '64K_defocus_m2_tomo_12_bin1_WBP_CatBinned1.mrc',
        })

DataSet(name='tomoMask', folder='tomoMask',
        files={
            'vTomo1': 'vTomo1_flt.mrc',
            'vTomo2': 'vTomo2_flt.mrc',
            'meshCoordinates': 'meshCoordinates.csv',
        })

########################################################################################################################
EMD_10439 = 'emd_10439'
EMPIAR_10491 = 'empiar_10491'


class DataSetEmd10439(Enum):
    tomogramsDir = 'tomograms'
    tomoEmd10439 = 'tomograms/emd_10439.mrc'
    coords3dStarFile = 'importFromStarFiles/picking_001_parts.star'
    coords3dStarFileWithSRate = 'importFromStarFiles/picking_001_parts_with_sRate.star'
    subtomogramsStarFile = 'importFromStarFiles/class_ap_r_ali_k1_split.star'
    scipionSqlite3dCoords = 'importFromScipionSqlite/coordinates.sqlite'
    scipionSqlite3dCoordsSomeBad = 'importFromScipionSqlite/coordinates_3badCoords.sqlite'
    tomomasksAnnotatedDir = 'tomomasksAnnotated'
    tomoMaskAnnotated = 'tomomasksAnnotated/emd_10439_materials.mrc'
    coords39Sqlite = 'coordinates/coordinates39.sqlite'
    nParticles = 39
    binFactor = 2
    bin2BoxSize = 44
    unbinnedBoxSize = 88
    unbinnedSRate = 13.68
    bin2SRate = 27.36


DataSet(name=EMD_10439, folder=EMD_10439, files={el.name: el.value for el in DataSetEmd10439})

########################################################################################################################
RE4_STA_TUTO = 'relion40_sta_tutorial_data'

TS_01 = 'TS_01'
TS_03 = 'TS_03'
TS_43 = 'TS_43'
TS_45 = 'TS_45'
TS_54 = 'TS_54'
nAnglesDict = {TS_03: 40,
               TS_54: 41}
# Acquisition
voltage = 300
sphericalAb = 2.7
amplitudeContrast = 0.07
magnification = 105000
tiltAxisAngle = 85.3
tiltStep = 3
initialDose = 0
dosePerTiltImgWithTltFile = 3

# Acquisition common parameters
dosePerTiltImg = dosePerTiltImgWithTltFile
testAcq = TomoAcquisition(voltage=voltage,
                          sphericalAberration=sphericalAb,
                          amplitudeContrast=amplitudeContrast,
                          magnification=magnification,
                          doseInitial=initialDose,
                          tiltAxisAngle=tiltAxisAngle,
                          dosePerFrame=dosePerTiltImg,
                          angleMax=60,
                          accumDose=123,  # The tlt files indicate they have the same accum dose, despite they have a
                          # different number of images
                          step=tiltStep)
# Acquisition of TS_03
testAcq03 = testAcq.clone()
testAcq03.setAngleMin(-57)
# Acquisition of TS_54
testAcq54 = testAcq.clone()
testAcq54.setAngleMin(-60)
# Tilt series acq dict
tsAcqDict = {TS_03: testAcq03,
             TS_54: testAcq54}
# Acquisition of TS_03 interpolated
testAcq03Interp = testAcq03.clone()
testAcq03Interp.setTiltAxisAngle(0)
# Acquisition of TS_54 interpolated
testAcq54Interp = testAcq54.clone()
testAcq54Interp.setTiltAxisAngle(0)
# Tilt series interpolated acq dict
tsAcqInterpDict = {TS_03: testAcq03Interp,
                   TS_54: testAcq54Interp}


class DataSetRe4STATuto(Enum):
    voltage = voltage
    sphericalAb = sphericalAb
    amplitudeContrast = amplitudeContrast
    magnification = magnification
    unbinnedPixSize = 1.35
    boxSizeBin4 = 192
    boxSizeBin2 = 256
    croppedBoxSizeBin4 = 96
    croppedBoxSizeBin2 = 128
    tiltAxisAngle = tiltAxisAngle
    initialDose = initialDose
    dosePerTiltImg = 3.05  # Mean dose
    dosePerTiltImgWithTltFile = dosePerTiltImgWithTltFile
    exclusionWordsTs03 = 'output 01 43 45 54'
    exclusionWordsTs54 = 'output 01 43 45 03'
    exclusionWordsTs03ts54 = 'output 01 43 45'
    correctCoordsFormula = 'item._y.set(item._y.get() + 18)'
    symmetry = 'C6'

    # Acquisition
    testAcq01 = testAcq03.clone()
    testAcq03 = testAcq03
    testAcq43 = testAcq54.clone()
    testAcq45 = testAcq54.clone()
    testAcq54 = testAcq54
    tsAcqDict = tsAcqDict
    testAcq03Interp = testAcq03Interp
    testAcq54Interp = testAcq54Interp
    tsAcqInterpDict = tsAcqInterpDict

    # Tilt series: 5 TS with no. tilt-images:
    #   - TS_01 = 40
    #   - TS_03 = 40
    #   - TS_43 = 41
    #   - TS_45 = 41
    #   - TS_54 = 41
    nTs = 5
    tsPath = 'tomograms'
    tsPattern = '*/{TS}.mrc'
    transformPattern = 'TS*/*.xf'
    ctfPattern = '*/*.defocus'
    nAnglesDict = nAnglesDict
    dimsTsBin1Dict = {TS_03: [3710, 3838, 40], TS_54: [3710, 3838, 41]}
    dimsTsInterpBin1Dict = {TS_03: [3838, 3710, 40], TS_54: [3838, 3710, 41]}
    tsDims40 = [3710, 3838, 40]
    tsDims41 = [3710, 3838, 41]

    # Tomograms: 5 tomograms with a thickness of (px):
    #   - TS_01 = 340
    #   - TS_03 = 280
    #   - TS_43 = 300
    #   - TS_45 = 300
    #   - TS_54 = 280
    nTomos = 5
    tomosPath = 'tomograms'
    tomosPattern = '*.mrc'
    sRateBin4 = 5.4
    tomoDimsThk300 = [928, 928, 300]
    tomoDimsThk280 = [928, 928, 280]
    tomoDimsThk340 = [928, 928, 340]

    # Coordinates
    coordsStarSubset = 'input/coords_subset_ts03_ts54.star'
    nCoordsFromTs03 = 99
    nCoordsFromTs54 = 111
    nCoordsTotal = 210

    # Initial model
    initModelRelion = 'initial_model_relion_bin4.mrc'

    # Rec particle
    recParticleBin6 = 'rec_particle_bin6.mrc'
    recParticleBin2 = 'recParticle_bin2.mrc'
    recParticleHalf1Bin2 = 'half1_bin2.mrc'
    recParticleHalf2Bin2 = 'half2_bin2.mrc'

    # Masks
    maskBin4 = 'masks/mask_align_bin4.mrc'
    maskBin2 = 'masks/mask_align_bin2.mrc'
    maskFscBin2 = 'masks/mask_fsc_bin2.mrc'

    # For EMAN testing
    initVolByEman = 'testEman/initialModel0ByEman.mrc'

    # For CISTEM Ctffind import testing
    cistemFilesPath = 'testCtfFind'

    # For Aretomo2 CTF import testing
    aretomoCtfFilesPath = 'testAreTomoCtf'

    @classmethod
    def genTestTsDicts(cls,
                       tsIdList: Tuple = (TS_01, TS_03, TS_43, TS_45,  TS_54),
                       unbinnedXYDims: Tuple = (3710, 3838),
                       binFactor: int = 1,
                       swapXY: bool = False,
                       isImod: bool = False):
        testAcqObjDict = dict()
        expectedDimensionsDict = dict()
        anglesCountDict = dict()
        nImgs40 = 40
        nImgs41 = 41
        expectedDims40 = cls._getExpectedTsDims(unbinnedXYDims,
                                                nImgs=nImgs40,
                                                binFactor=binFactor,
                                                swapXY=swapXY,
                                                isImod=isImod)
        expectedDims41 = cls._getExpectedTsDims(unbinnedXYDims,
                                                nImgs=nImgs41,
                                                binFactor=binFactor,
                                                swapXY=swapXY,
                                                isImod=isImod)
        if TS_01 in tsIdList:
            testAcqObjDict[TS_01] = cls.testAcq01.value
            expectedDimensionsDict[TS_01] = expectedDims40
            anglesCountDict[TS_01] = nImgs40
        if TS_03 in tsIdList:
            testAcqObjDict[TS_03] = cls.testAcq03.value
            expectedDimensionsDict[TS_03] = expectedDims40
            anglesCountDict[TS_03] = nImgs40
        if TS_43 in tsIdList:
            testAcqObjDict[TS_43] = cls.testAcq43.value
            expectedDimensionsDict[TS_43] = expectedDims41
            anglesCountDict[TS_43] = nImgs41
        if TS_45 in tsIdList:
            testAcqObjDict[TS_45] = cls.testAcq45.value
            expectedDimensionsDict[TS_45] = expectedDims41
            anglesCountDict[TS_45] = nImgs41
        if TS_54 in tsIdList:
            testAcqObjDict[TS_54] = cls.testAcq54.value
            expectedDimensionsDict[TS_54] = expectedDims41
            anglesCountDict[TS_54] = nImgs41

        return testAcqObjDict, expectedDimensionsDict, anglesCountDict

    @staticmethod
    def _getExpectedTsDims(unbinnedXYDims: Tuple,
                           nImgs: int,
                           binFactor: int = 1,
                           swapXY: bool = False,
                           isImod: bool = False):
        dims = []
        for iDim in unbinnedXYDims:
            newDim = math.ceil(iDim / binFactor)
            # Imod always generates images with even dimensions, at least for the TS
            if isImod:
                if newDim % 2 != 0:
                    newDim += 1
            dims.append(newDim)

        if swapXY:
            dims.reverse()

        return dims + [nImgs]

    @classmethod
    def genTestTomoDicts(cls, tsIdList: Tuple = (TS_01, TS_03, TS_43, TS_45, TS_54)):
        testAcqObjDict = dict()
        expectedDimensionsDict = dict()
        if TS_01 in tsIdList:
            testAcqObjDict[TS_01] = cls.testAcq01.value
            expectedDimensionsDict[TS_01] = cls.tomoDimsThk340.value
        if TS_03 in tsIdList:
            testAcqObjDict[TS_03] = cls.testAcq03.value
            expectedDimensionsDict[TS_03] = cls.tomoDimsThk280.value
        if TS_43 in tsIdList:
            testAcqObjDict[TS_43] = cls.testAcq43.value
            expectedDimensionsDict[TS_43] = cls.tomoDimsThk300.value
        if TS_45 in tsIdList:
            testAcqObjDict[TS_45] = cls.testAcq45.value
            expectedDimensionsDict[TS_45] = cls.tomoDimsThk300.value
        if TS_54 in tsIdList:
            testAcqObjDict[TS_54] = cls.testAcq54.value
            expectedDimensionsDict[TS_54] = cls.tomoDimsThk280.value

        return testAcqObjDict, expectedDimensionsDict


DataSet(name=RE4_STA_TUTO, folder=RE4_STA_TUTO, files={el.name: el.value for el in DataSetRe4STATuto})

########################################################################################################################
RE5_STA = 'relion50_sta'


class DataSetRe5STA(Enum):  # Extends the enumeration DataSetRe4STATuto
    # Coordinates
    coordsPickedWithRe5Star = 'coordsPickedWithRe5.star'
    nCoords = 1527  # Corresponding to the 2 tomograms imported, TS_03 and TS_54
    coordsSRate = 10
    boxSize = 96
    # Tomogrmas
    tomosDir = 'tomosRecWithRe5'
    nTomos = 2
    binFactor = 6
    tomosSRate = 10
    tomosDims = [540, 540, 270]


DataSet(name=RE5_STA, folder=RE5_STA, files={el.name: el.value for el in DataSetRe5STA})


DataSet(name=EMPIAR_10491, folder=EMPIAR_10491, files={})
