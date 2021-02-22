# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import numpy as np

import pyworkflow.utils as pwutils

def getMeshVolFileName(volId):
    return 'Meshes_Vol%d.txt' % volId

def setOfMeshes2Files(meshes, path):

    def writeFile():
        fnInputCoor = getMeshVolFileName(currentVolId)
        pathInputCoor = pwutils.join(path, fnInputCoor)
        np.savetxt(pathInputCoor, np.asarray(coords), fmt='%d', delimiter=",")

    currentVolId = None
    coords = []
    for coor in meshes.iterCoordinates(orderBy="_volId"):
        if coor.getVolId() != currentVolId:
            if currentVolId != None:
                writeFile()
            currentVolId = coor.getVolId()
            coords = []
        coords.append([coor.getX(), coor.getY(), coor.getZ(), coor.getGroupId()])
    if coords:
        writeFile()
