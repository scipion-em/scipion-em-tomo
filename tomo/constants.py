# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] National Center for Biotechnology, Madrid, Spain
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


# ---------------------------- ORIGIN Conventions ----------------------------

# Origin conversion functions based on different conventions (to move the coordinates
# from convention to bottom left corner)
# Input parameter of lambda functions is a tuple/list with the dimensions of a Tomogram
# This functions return the vector in VOXELS. No decimal places allowed --> (int(vx), int(vy), int(vz))
BOTTOM_LEFT_CORNER = lambda dim: (0, 0, 0)  # Eman, Xmipp, Imod, crYOLO
TOP_LEFT_CORNER = lambda dim: (0, dim[1], 0)
CENTER_GRAVITY = lambda dim: (0.5 * dim[0], 0.5 * dim[1], 0.5 * dim[2])

# This is a bit different as we do not know beforehand where te
# user has set the origin (Coordinate3D methods will use the None
# to retrieve the correct origin)
SCIPION = lambda dim: None

# ----------------------------------------------------------------------------
# Error message when non-matching tomograms are found in protocol import 3d coordinates from Scipion sqlite
ERR_COORDS_FROM_SQLITE_NO_MATCH = 'No matching tomograms were found to the introduced coordinates. Coordinates ' \
                                  'attribute tomoId should be the same as the corresponding tomogram attribute tsId ' \
                                  'or be contained in its name.'

# Import tomomasks error messages
ERR_NO_TOMOMASKS_GEN = 'No tomomasks were generated. Check the output log for more details.'
ERR_NON_MATCHING_TOMOS = 'No matching tomograms were found with the input tomomasks. The match is carried out by ' \
                         'checking if the tsId of each tomogram is contained in the basename of the tomomaks. If ' \
                         'the tomograms do not have that attribute, then the base name is considered.'


# ---------------------------- TRANSFORMATION MATRIX Conventions ----------------------------

# A detailed description of the transformation convention adopted by Scipion and its conversion to other packages is
# available in the help of the function "convertMatrix" defined in the "objects.py" file
TR_SCIPION = None
TR_RELION = "relion"
TR_DYNAMO = "dynamo"
TR_EMAN = "eman"

# Conversion direction constants
SET = "set"
GET = "get"

# --------------------------- Napari variables --------------------------------

def getNaparyEnvName(version):
    return f"napari-{version}"


V0_4_17 = '0.4.17'

NAPARI_DEF_VER = V0_4_17
NAPARI_ACTIVATION_CMD = 'conda activate %s' % getNaparyEnvName(NAPARI_DEF_VER)
NAPARI_ENV_ACTIVATION = 'NAPARI_ENV_ACTIVATION'
