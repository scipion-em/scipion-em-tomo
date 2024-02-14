# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
from enum import Enum

from pyworkflow.tests import DataSet

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

EMD_10439 = 'emd_10439'


class DataSetEmd10439(Enum):
    tomoEmd10439 = 'tomograms/emd_10439.mrc'
    coords3dStarFile = 'importFromStarFiles/picking_001_parts.star'
    coords3dStarFileWithSRate = 'importFromStarFiles/picking_001_parts_with_sRate.star'
    subtomogramsStarFile = 'importFromStarFiles/class_ap_r_ali_k1_split.star'
    scipionSqlite3dCoords = 'importFromScipionSqlite/coordinates.sqlite'
    scipionSqlite3dCoordsSomeBad = 'importFromScipionSqlite/coordinates_3badCoords.sqlite'
    tomoMaskAnnotated = 'tomomasksAnnotated/emd_10439_materials.mrc'
    coords39Sqlite = 'coordinates/coordinates39.sqlite'
    nParticles = 39
    binFactor = 2
    bin2BoxSize = 44
    unbinnedBoxSize = 88
    unbinnedSRate = 13.68
    bin2SRate = 27.36


DataSet(name=EMD_10439, folder=EMD_10439, files={el.name: el.value for el in DataSetEmd10439})


class DataSetRe4STATuto(Enum):
    voltage = 300
    sphericalAb = 2.7
    amplitudeContrast = 0.07
    magnification = 105000
    unbinnedPixSize = 1.35
    croppedBoxSizeBin4 = 96
    croppedBoxSizeBin2 = 128
    tiltAxisAngle = 85.3
    initialDose = 0
    dosePerTiltImg = 3.05
    accumDose = 122
    exclusionWordsTs03 = 'output 01 43 45 54'
    exclusionWordsTs03ts54 = 'output 01 43 45'
    correctCoordsFormula = 'item._y.set(item._y.get() + 18)'
    # Tilt series
    tsPath = 'tomograms'
    tsPattern = '*/{TS}.mrc'
    transformPattern = 'TS*/*.xf'
    ctfPattern = '*/*.defocus'
    # Tomograms
    tomosPath = 'tomograms'
    tomosPattern = '*.mrc'
    sRateBin4 = 5.4
    # Coordinates
    coordsStarSubset = 'input/coords_subset_ts03_ts54.star'
    nCoordsFromTs03 = 99
    nCoordsFromTs54 = 111
    # Masks
    maskBin4 = 'masks/mask_align_bin4.mrc'
    maskBin2 = 'masks/mask_align_bin2.mrc'
    # For EMAN testing
    initVolByEman = 'testEman/initialModel0ByEman.mrc'
    # For CISTEM Ctffind import testing
    cistemFilesPath = 'testCtfFind'
    # For Aretomo2 CTF import testing
    aretomoCtfFilesPath = 'testAreTomoCtf'


RE4_STA_TUTO = 'relion40_sta_tutorial_data'
DataSet(name=RE4_STA_TUTO, folder=RE4_STA_TUTO, files={el.name: el.value for el in DataSetRe4STATuto})
