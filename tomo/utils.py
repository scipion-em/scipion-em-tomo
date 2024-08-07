# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *              Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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

import os
import re
import importlib
from typing import List, Set

import numpy as np
import math
import logging
logger = logging.getLogger(__name__)

import pyworkflow.utils as pwutils

import tomo.constants as const
from tomo.objects import SetOfCoordinates3D, SetOfSubTomograms, SetOfTiltSeries, Coordinate3D, SubTomogram, TiltSeries, \
    CTFTomoSeries, TiltImage, CTFTomo


def existsPlugin(plugin):
    return importlib.util.find_spec(plugin)


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


def initDictVesicles(coordinates):
    tomos = coordinates.getPrecedents()
    volIds = coordinates.aggregate(["MAX"], "_volId", ["_volId"])
    volIds = [d['_volId'] for d in volIds]
    tomoNames = [pwutils.removeBaseExt(tomos[volId].getFileName()) for volId in volIds]
    dictVesicles = {tomoField: {'vesicles': [], 'normals': [], 'ids': [], 'volId': volIds[idt]}
                    for idt, tomoField in enumerate(tomoNames)}
    return dictVesicles, tomoNames


def extractVesicles(coordinates, dictVesicles, tomoName):
    # tomoId = list(dictVesicles.keys()).index(tomoName) + 1
    tomoId = dictVesicles[tomoName]['volId']
    groupIds = coordinates.aggregate(["MAX"], "_volId", ["_groupId", "_volId"])
    groupIds = [d['_groupId'] for d in groupIds if d['_volId'] == tomoId]
    if not dictVesicles[tomoName]['vesicles']:
        for idv in groupIds:
            vesicle = []
            normals = []
            ids = []
            for coord in coordinates.iterCoordinates(volume=coordinates.getPrecedents()[tomoId]):
                if coord.getGroupId() == idv:
                    vesicle.append(coord.getPosition(const.SCIPION))
                    trMat = coord.getMatrix()
                    normals.append(normalFromMatrix(trMat))
                    ids.append(coord.getObjId())
            dictVesicles[tomoName]['vesicles'].append(np.asarray(vesicle))
            dictVesicles[tomoName]['normals'].append(np.asarray(normals))
            dictVesicles[tomoName]['ids'].append(np.asarray(ids))
    return dictVesicles


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

    D = np.array(
        [x * x + y * y - 2 * z * z, x * x + z * z - 2 * y * y, 2 * x * y, 2 * x * z, 2 * y * z, 2 * x, 2 * y, 2 * z,
         1 + 0 * x])
    D = D.transpose()

    # Solve the normal system of equations
    d2 = x * x + y * y + z * z  # The RHS of the llsq problem (y's)
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
    radii = np.sqrt(1 / abs(evals))
    sgns = np.sign(evals)
    radii = radii * sgns

    # Calculate difference of the fitted points from the actual data normalized by the conic radii
    d = np.array([x - center[0], y - center[1], z - center[2]])  # shift data to origin
    d = d.transpose() @ evecs  # Rotate to cardinal axes of the conic
    d = d.transpose()
    d = np.array([d[:, 0] / radii[0], d[:, 1] / radii[1], d[:, 2] / radii[2]])  # normalize to the conic radii
    chi2 = np.sum(
        np.abs(1 - np.sum(np.dot((d ** 2), np.tile(sgns.conj().transpose(), (d.shape[0], 1)).transpose()), 1)))

    if np.abs(v[-1]) > 1e-6:
        v = -v / v[-1]  # Normalize to the more conventional form with constant term = -1
    else:
        v = -np.sign(v[-1]) * v

    return center, radii, v, evecs, chi2


def generatePointCloud(v, tomoDim):
    ygrid = np.linspace(0, 1, 100, dtype=float)
    zgrid = np.linspace(0, 1, 100, dtype=float)

    pointCloud = []
    epsilon = 1e-6

    # a*x*x + b*y*y + c*z*z + 2*d*x*y + 2*e*x*z + 2*f*y*z + 2*g*x + 2*h*y + 2*i*z + j = 0

    if abs(v[0]) > epsilon:
        v = v / v[0]
        a = 1
        b = v[1]
        c = v[2]
        d = v[3]
        e = v[4]
        f = v[5]
        g = v[6]
        h = v[7]
        i = v[8]
        j = v[9]
        print('X^2')
        for z in zgrid:
            for y in ygrid:
                A = a
                B = (2 * d * y) + (2 * e * z) + (2 * g)
                C = (b * y * y) + (c * z * z) + (2 * f * y * z) + (2 * h * y) + (2 * i * z) + j
                D = B * B - (4 * A * C)
                if D == 0:
                    x = (-1) * B / 2 * A
                    pointCloud.append([int(x * tomoDim[0]), int(y * tomoDim[1]), int(z * tomoDim[2])])
                if D > 0:
                    sqrtD = np.sqrt(D)
                    x1 = ((-1) * B + sqrtD) / 2 * A
                    x2 = ((-1) * B - sqrtD) / 2 * A
                    pointCloud.append([int(x1 * tomoDim[0]), int(y * tomoDim[1]), int(z * tomoDim[2])])
                    pointCloud.append([int(x2 * tomoDim[0]), int(y * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[3]) > epsilon:
        v = v / v[3]
        b = v[1]
        c = v[2]
        d = 1
        e = v[4]
        f = v[5]
        g = v[6]
        h = v[7]
        i = v[8]
        j = v[9]
        print('X')
        for z in zgrid:
            for y in ygrid:
                A = (2 * d * y) + (2 * e * z) + (2 * g)
                B = b * y * y + c * z * z + 2 * f * y * z + 2 * h * y + 2 * i * z + j
                x = (-1) * B / A
                pointCloud.append([int(x * tomoDim[0]), int(y * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[4]) > epsilon:
        v = v / v[4]
        b = v[1]
        c = v[2]
        e = 1
        f = v[5]
        g = v[6]
        h = v[7]
        i = v[8]
        j = v[9]
        print('X')
        for z in zgrid:
            for y in ygrid:
                A = (2 * e * z) + (2 * g)
                B = b * y * y + c * z * z + 2 * f * y * z + 2 * h * y + 2 * i * z + j
                x = (-1) * B / A
                pointCloud.append([int(x * tomoDim[0]), int(y * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[6]) > epsilon:
        v = v / v[6]
        b = v[1]
        c = v[2]
        f = v[5]
        g = 1
        h = v[7]
        i = v[8]
        j = v[9]
        print('X')
        for z in zgrid:
            for y in ygrid:
                A = 2 * g
                B = b * y * y + c * z * z + 2 * f * y * z + 2 * h * y + 2 * i * z + j
                x = (-1) * B / A
                pointCloud.append([int(x * tomoDim[0]), int(y * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[1]) > epsilon:
        v = v / v[1]
        b = 1
        c = v[2]
        f = v[5]
        h = v[7]
        i = v[8]
        j = v[9]
        print('Y^2')
        for z in zgrid:
            A = b
            B = (2 * f * z) + (2 * h)
            C = (c * z * z) + (2 * i * z) + j
            D = B * B - (4 * A * C)
            if D > 0:
                sqrtD = np.sqrt(D)
                y1 = ((-1) * B + sqrtD) / 2 * A
                y2 = ((-1) * B - sqrtD) / 2 * A
                pointCloud.append([0, int(y1 * tomoDim[1]), int(z * tomoDim[2])])
                pointCloud.append([0, int(y2 * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[5]) > epsilon:
        v = v / v[5]
        c = v[2]
        f = 1
        h = v[7]
        i = v[8]
        j = v[9]
        print('Y')
        for z in zgrid:
            A = (2 * f * z) + (2 * h)
            B = (c * z * z) + (2 * i * z) + j
            y = (-1) * B / A
            pointCloud.append([0, int(y * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[7]) > epsilon:
        v = v / v[7]
        c = v[2]
        h = 1
        i = v[8]
        j = v[9]
        print('Y')
        for z in zgrid:
            A = 2 * h
            B = (c * z * z) + (2 * i * z) + j
            y = (-1) * B / A
            pointCloud.append([0, int(y * tomoDim[1]), int(z * tomoDim[2])])

    elif abs(v[2]) > epsilon:
        v = v / v[2]
        c = 1
        i = v[8]
        j = v[9]
        print('Z^2')  # if algDesc with z2 = 0 for z values, x=y=0
        for z in zgrid:
            result = c * z * z + 2 * i * z + j
            if result == 0:
                pointCloud.append([0, 0, int(z * tomoDim[2])])

    elif abs(v[8]) > epsilon:
        v = v / v[8]
        i = 1
        j = v[9]
        print('Z')  # if algDesc with z = 0 for z values, x=y=0
        for z in zgrid:
            result = 2 * i * z + j
            if result == 0:
                pointCloud.append([0, 0, int(z * tomoDim[2])])

    return pointCloud


def isMatchingByTsId(set1, set2):
    return True if getattr(set1.getFirstItem(), _getTsIdLabel(set1), None) and \
                   getattr(set2.getFirstItem(), _getTsIdLabel(set2), None) else False


def _getTsIdLabel(setObject):
    """This attribute is named tsId in all the tomography objects excepting in coordinates or subtomograms (via the
    corresponding coordinate)"""
    setType = type(setObject)
    if setType == SetOfCoordinates3D:
        return Coordinate3D.TOMO_ID_ATTR
    elif setType == SetOfSubTomograms:
        return SubTomogram.VOL_NAME_FIELD
    else:
        return TiltSeries.TS_ID_FIELD


def _recoverObjFromRelations(sourceObj, protocol, stopSearchCallback):
    logger.debug("Retrieving relations for %s." % sourceObj)
    p = protocol.getProject()
    graph = p.getSourceGraph(True)  # Graph with all the relations
    sourceNode = graph.getNode(sourceObj.strId())  # Node corresponding to the source object
    # Climb up in the relations graph until the target condition provided in the callback input is fulfilled
    nodes = sourceNode.getParents()
    while nodes:
        sourceNode = nodes.pop()
        if not sourceNode.isRoot():
            relatedOutput = sourceNode.pointer.get()
            logger.debug("Checking related object: %s" % relatedOutput)
            if stopSearchCallback(relatedOutput):
                return relatedOutput
            else:
                parents = sourceNode.getParents()
                if parents is not None:
                    nodes += parents
    return None


def getNonInterpolatedTsFromRelations(sourceObj, prot):
    def stopSearchCallback(pObj):
        return type(pObj) == SetOfTiltSeries and pObj.getFirstItem().getFirstItem().hasTransform()
    return _recoverObjFromRelations(sourceObj, prot, stopSearchCallback)


def getObjFromRelation(sourceObj, prot, targetObj):
    def stopSearchCallback(pObj):
        return type(pObj) == targetObj
    return _recoverObjFromRelations(sourceObj, prot, stopSearchCallback)


def getRotationAngleAndShiftFromTM(ti):
    """ This method calculates que tilt image in-plane rotation angle and shifts from its associated transformation
    matrix."""

    tm = ti.getTransform().getMatrix()
    cosRotationAngle = tm[0][0]
    sinRotationAngle = tm[1][0]
    rotationAngle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

    shifts = [tm[0][2], tm[0][2]]

    return rotationAngle, shifts


def scaleTrMatrixShifts(inTrMatrix, scaleFactor):
    """In Scipion data model, the shifts are in pixels, so they must be scaled properly when the reference size
    (normally the tomograms from which the coordinates were picked) changes, for example, extracting the particles
    or the coordinates, to another set of tomograms.

    :param inTrMatrix: transformation matrix from which the shifts will be read.
    :param scaleFactor: scale factor that will be used to scale the shifts properly.
    :return: a transformation matrix with the shifts properly scaled."""
    if scaleFactor != 1:  # It can be lower (smaller source) or higher (bigger source) than one
        shifts = np.array([inTrMatrix[0, 3], inTrMatrix[1, 3], inTrMatrix[2, 3]])
        scaledShifts = scaleFactor * shifts
        inTrMatrix[0, 3] = scaledShifts[0]
        inTrMatrix[1, 3] = scaledShifts[1]
        inTrMatrix[2, 3] = scaledShifts[2]
    return inTrMatrix


def getCommonTsAndCtfElements(ts: TiltSeries, ctfTomoSeries: CTFTomoSeries, onlyEnabled: bool = True) -> Set[int]:
    """Given a tilt-series and a CTFTomoSeries, it finds the common elements present and enabled in both sets, and
    returns a list with the corresponding acquisition orders or indices, if acquisition order is not present in the
    CTFTomoSeries introduced (old versions, backwards compatibility). By default, it takes the common active elements,
    but it may take common elements no matter if they're enabled or not by setting the input onlyEnabled to False.
    """
    # Attribute _acqOrder was recently added to CTFTomo, so it will be used to discriminate
    ctfTomoSeries._getMapper()  # Avoid finding closed mappers when combining cached sets of sets (TS, CTF) and
    # calls to getFirstItem(). The second closes the first and so on
    firstCtfTomo = ctfTomoSeries.getFirstItem()
    acqOrder = getattr(firstCtfTomo, CTFTomo.ACQ_ORDER_FIELD, None)
    if acqOrder:
        msgStr = 'acquisition order'
        if onlyEnabled:
            tsAcqOrderSet = {ti.getAcquisitionOrder() for ti in ts if ti.isEnabled()}
            ctfAcqOrderSet = {ctf.getAcquisitionOrder() for ctf in ctfTomoSeries if ctf.isEnabled()}
        else:
            tsAcqOrderSet = {ti.getAcquisitionOrder() for ti in ts}
            ctfAcqOrderSet = {ctf.getAcquisitionOrder() for ctf in ctfTomoSeries}
    else:
        msgStr = 'index'
        if onlyEnabled:
            tsAcqOrderSet = {ti.getIndex() for ti in ts if ti.isEnabled()}
            ctfAcqOrderSet = {ctf.getIndex() for ctf in ctfTomoSeries if ctf.isEnabled()}
        else:
            tsAcqOrderSet = {ti.getIndex() for ti in ts}
            ctfAcqOrderSet = {ctf.getIndex() for ctf in ctfTomoSeries}

    logger.debug(f'getCommonTsAndCtfElements: tsId = {ts.getTsId()}, matching used field is {msgStr}')
    return tsAcqOrderSet & ctfAcqOrderSet
