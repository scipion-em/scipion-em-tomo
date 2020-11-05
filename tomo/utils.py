# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os, re
import numpy as np

import pyworkflow.utils as pwutils


def _getUniqueFileName(pattern, filename, filePaths=None):
 if filePaths is None:
     filePaths = [re.split(r'[$*#?]', pattern)[0]]

 commPath = pwutils.commonPath(filePaths)
 return filename.replace(commPath + "/", "").replace("/", "_")

def _matchFileNames(originalName, importName):
 return os.path.basename(importName) in originalName

def normalFromMatrix(transformation):
    rotation = transformation[:3, :3]
    axis = np.array([0, 0, 1])
    normal = np.linalg.inv(rotation).dot(axis)
    return normal

def extractVesicles(coordinates):
    tomos = coordinates.getPrecedents()
    tomoNames = [pwutils.removeBaseExt(tomo.getFileName()) for tomo in tomos]
    vesicleIds = set([coord._vesicleId.get() for coord in coordinates.iterCoordinates()])
    tomo_vesicles = {tomoField: {'vesicles': [], 'normals': [], 'ids': []}
                     for tomoField in tomoNames}

    for idt, tomo in enumerate(tomos.iterItems()):
        for idv in vesicleIds:
            vesicle = []
            normals = []
            ids = []
            for coord in coordinates.iterCoordinates(volume=tomo):
                if coord._vesicleId == idv:
                    vesicle.append(coord.getPosition())
                    trMat = coord.getMatrix()
                    normals.append(normalFromMatrix(trMat))
                    ids.append(coord.getObjId())
            tomo_vesicles[tomoNames[idt]]['vesicles'].append(np.asarray(vesicle))
            tomo_vesicles[tomoNames[idt]]['normals'].append(np.asarray(normals))
            tomo_vesicles[tomoNames[idt]]['ids'].append(np.asarray(ids))
    return tomo_vesicles