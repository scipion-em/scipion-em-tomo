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
import pyvista as pv
import numpy as np

import pyworkflow.utils as pwutils


def _getUniqueFileName(pattern, filename, filePaths=None):
 if filePaths is None:
     filePaths = [re.split(r'[$*#?]', pattern)[0]]

 commPath = pwutils.commonPath(filePaths)
 return filename.replace(commPath + "/", "").replace("/", "_")

def _matchFileNames(originalName, importName):
 return os.path.basename(importName) in originalName

def rotation_matrix_from_vectors(vec1, vec2):
    # a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    a, b = vec1, vec2
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s != 0:
        tr = np.zeros([4, 4])
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        tr[0:3, 0:3] = rotation_matrix
        tr[-1, -1] = 1
        return tr
    else:
        return np.eye(4)

def delaunayTriangulation(cloud, adjustCloud=True):
    cloud_pv = pv.PolyData(cloud)
    mesh = cloud_pv.delaunay_3d()
    shell = mesh.extract_geometry().triangulate()

    # If points lie outside the shell, they are readjusted to compute their normals
    if adjustCloud:
        shell.subdivide(4, inplace=True)
        areZero = np.where((shell.point_normals == (0, 0, 0)).all(axis=1))
        points = np.asarray(shell.points)
        points[areZero] = np.array((0, 0, 0))
        newCoords = np.asarray([points[np.argmin(np.linalg.norm(points - point, axis=1))]
                                for point in cloud])
        shell = delaunayTriangulation(newCoords, adjustCloud=False)
    return shell

def computeNormals(triangulation):
    triangulation.compute_normals(inplace=True)
    return triangulation.point_normals