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


def delaunayTriangulation(cloud):
    cloud = pv.PolyData(cloud)
    mesh = cloud.delaunay_3d()
    shell = mesh.extract_geometry().triangulate()
    return shell


def computeNormals(triangulation):
    triangulation.compute_normals(inplace=True)
    return triangulation.point_normals


def fit_ellipsoid(x, y, z):
    """ Fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
    and A + B + C = 3 constraint removing one extra parameter (from Matlab function "Fit Ellipsoid"). """

    # OUTPUT:
    # center: ellispoid center coordinates[xc, yc, zc]
    # radii: ellipsoid radii[a, b, c]
    # evecs: the radii directions as columns of the 3x3 matrix
    # v: the 10 parameters describing the ellipsoid algebraically:
    #     Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
    # chi2: residual sum of squared errors(chi^2), in the coordinate frame in which the ellipsoid is a unit sphere

    D = np.array([x*x + y*y - 2*z*z, x*x + z*z - 2*y*y, 2*x*y, 2*x*z, 2*y*z, 2*x, 2*y, 2*z, 1 + 0*x])
    D = D.transpose()

    # Solve the normal system of equations
    d2 = x*x + y*y + z*z  # The RHS of the llsq problem (y's)
    cD = D.conj().transpose()
    a = cD @ D
    b = cD @ d2
    u = np.linalg.lstsq(a, b, rcond=None)[0]  # Solution to the normal equations

    # Find the ellipsoid parameters
    # Convert back to the conventional algebraic form
    v = np.zeros(10)
    v[0] = u[0] + u[1] - 1
    v[1] = u[0] - 2 * u[1] - 1
    v[2] = u[1] - 2 * u[0] - 1
    v[3:10] = u[2:9]

    # Form the algebraic form of the ellipsoid
    A = np.array([[v[0], v[3], v[4], v[6]],
                  [v[3], v[1], v[5], v[7]],
                  [v[4], v[5], v[2], v[8]],
                  [v[6], v[7], v[8], v[9]]])

    # Find the center of the ellipsoid
    center = np.linalg.lstsq(-A[0:3, 0:3], v[6:9], rcond=None)[0]

    # Form the corresponding translation matrix
    T = np.eye(4)
    T[3, 0:3] = center.conj().transpose()

    # Translate to the center
    R = T * A * T.conj().transpose()

    # Solve the eigenproblem
    [evals, evecs] = np.linalg.eig(R[0:3, 0:3] / -R[3, 3])
    radii = np.sqrt(1/abs(evals))
    sgns = np.sign(evals)
    radii = radii * sgns

    # Calculate difference of the fitted points from the actual data normalized by the conic radii
    d = np.array([x - center[0], y - center[1], z - center[2]])  # shift data to origin
    d = d.transpose() @ evecs  # Rotate to cardinal axes of the conic
    d = d.transpose()
    d = np.array([d[:, 0] / radii[0], d[:, 1] / radii[1], d[:, 2] / radii[2]])  # normalize to the conic radii
    chi2 = np.sum(np.abs(1 - np.sum(np.dot((d**2), np.tile(sgns.conj().transpose(), (d.shape[0], 1)).transpose()), 1)))

    if np.abs(v[-1]) > 1e-6:
        v = -v / v[-1]  # Normalize to the more conventional form with constant term = -1
    else:
        v = -np.sign(v[-1]) * v

    return center, radii, v, evecs, chi2


def generatePointCloud(algDesc, v, tomoDim):
    # xgrid = np.linspace(0, tomoDim[0], tomoDim[0], dtype=int)
    # ygrid = np.linspace(0, tomoDim[1], tomoDim[1], dtype=int)
    # zgrid = np.linspace(0, tomoDim[2], tomoDim[2], dtype=int)

    xgrid = np.linspace(0, tomoDim[0], 100, dtype=int)
    ygrid = np.linspace(0, tomoDim[1], 100, dtype=int)
    zgrid = np.linspace(0, tomoDim[2], 100, dtype=int)

    pointCloud = []
    epsilon = 1e-5

    # if 'x*x' in algDesc:
    print('-v: ', v)
    if abs(v[0]) > epsilon:
        print('X2')  # x2 from algDesc for z, y values
        for z in zgrid:
            for y in ygrid:
                A = v[0]
                B = (2*v[3]*y) + (2*v[4]*z) + (2*v[6])
                C = (v[1]*y*y) + (v[2]*z*z) + (2*v[5]*y*z) + (2*v[7]*y) + (2*v[8]*z) + v[9]
                D = B*B - (4*A*C)
                if D > 0:
                    sqrtD = np.sqrt(D)
                    x1 = ((-1)*B + sqrtD) / 2*A
                    x2 = ((-1)*B - sqrtD) / 2*A
                    if x1 > epsilon:
                        pointCloud.append([x1, y, z])
                    if x2 > epsilon:
                        pointCloud.append([x2, y, z])

    # elif 'x' in algDesc:
    elif abs(v[3]) > epsilon or abs(v[4]) > epsilon or abs(v[6]) > epsilon:
        print('X')  # x from algDesc for z, y values
        for z in zgrid:
            for y in ygrid:
                A = (2*v[3]*y) + (2*v[4]*z) + (2*v[6])
                B = v[1]*y*y + v[2]*z*z + 2*v[5]*y*z + 2*v[7]*y + 2*v[8]*z + v[9]
                x = (-1)*B/A
                pointCloud.append([x, y, z])

    # elif 'y*y' in algDesc:
    elif abs(v[1]) > epsilon in algDesc:
        print('Y2')  # y2 from algDesc for z values, x=0
        for z in zgrid:
            A = v[1]
            B = (2*v[5]*z) + (2*v[7])
            C = (v[2]*z*z) + (2*v[8]*z) + v[9]
            D = B*B - (4*A*C)
            if D > 0:
                sqrtD = np.sqrt(D)
                y1 = ((-1)*B + sqrtD) / 2*A
                y2 = ((-1)*B - sqrtD) / 2*A
                if y1 > epsilon:
                    pointCloud.append([0, y1, z])
                if y2 > epsilon:
                    pointCloud.append([0, y2, z])

    # elif 'y' in algDesc:
    elif abs(v[5]) > epsilon or abs(v[7]) > epsilon:
        print('Y')  # y from algDesc for z values, x=0
        for z in zgrid:
            A = (2*v[5]*z) + (2*v[7])
            B = (v[2]*z*z) + (2*v[8]*z) + v[9]
            y = (-1)*B/A
            pointCloud.append([0, y, z])

    # elif 'z*z' in algDesc:
    elif abs(v[2]) > epsilon:
        print('Z2')  # if algDesc with z2 = 0 for z values, x=y=0
        for z in zgrid:
            result = v[2]*z*z + 2*v[8]*z + v[9]
            if result == 0:
                pointCloud.append([0, 0, z])

    # elif 'z' in algDesc:
    elif abs(v[8]) > epsilon:
        print('Z')  # if algDesc with z = 0 for z values, x=y=0
        for z in zgrid:
            result = 2*v[8]*z + v[9]
            if result == 0:
                pointCloud.append([0, 0, z])

    return pointCloud
