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

