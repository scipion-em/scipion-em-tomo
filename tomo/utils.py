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
import json
import os
import random
import re
import importlib
from os.path import join, exists
from typing import Set, List, Union, Protocol, Any
import time
import numpy as np
import math
import logging
import pyworkflow.utils as pwutils
import tomo.constants as const
from pwem.objects import Transform
from pyworkflow.utils import cyanStr
from tomo.objects import SetOfTiltSeries, TiltSeries, \
    CTFTomoSeries, CTFTomo, TiltImage, TomoAcquisition, LandmarkModel

logger = logging.getLogger(__name__)


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
        return type(pObj) == SetOfTiltSeries and pObj.hasAlignment()

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


# typing.Protocol declaring that inputs must implement .getTSIds()
class HasGetTsIds(Protocol):

    def getTSIds(self) -> Union[List[Any], Set[Any]]: ...


def getTsIdsIntersection(
        *emSets: HasGetTsIds,
        validateIntersectAndDiff: bool = True,
        allowEmptyIntersect: bool = False) -> Set[str]:
    """Extracts TS IDs from N objects using .getTSIds() and computes their
    intersection and generalized symmetric difference (union - intersection).
    """
    if not emSets:
        return set()

    # Extract IDs from each object via .getTsIds() and convert to set
    sets = [set(obj.getTSIds()) for obj in emSets]

    # Intersection: IDs present in ALL objects
    intersection = set.intersection(*sets)

    # Union: IDs present in AT LEAST ONE object
    union = set.union(*sets)

    # Union - Intersection
    difference = union - intersection

    # Do validation if required
    if validateIntersectAndDiff:
        _validateIntersectAndDiff(intersection, difference, allowEmptyIntersect=allowEmptyIntersect)

    return intersection


def _validateIntersectAndDiff(
        tsIdsIntersec: Set[str],
        tsIdsDiff: Set[str],
        allowEmptyIntersect: bool = False) -> None:
    if len(tsIdsIntersec) <= 0 and not allowEmptyIntersect:
        raise Exception("There isn't any common tsIds among the EM sets introduced.")

    if len(tsIdsDiff) > 0:
        logger.info(cyanStr(f"TsIds not common in the introduced EM sets are: {tsIdsDiff}"))


# STREAMING ############################################################################################
def sleepRandomly(lowTimeRange: float = 1.0,
                  highTimeRange: float = 3.0) -> None:
    """Throttle a streaming poll loop with a small random delay.

    The delay (a) keeps the polling loop from CPU/metadata-server spinning and
    (b) JITTERS the timing so multiple concurrent consumers do not synchronise
    their journal/heartbeat reads. The default range was lowered from 4-10s to
    1-3s for a more responsive stream: it keeps a non-zero spread (preserving the
    de-synchronisation jitter) and a >=1s lower bound that stays at/above the
    journal-read debounce window (Set._STREAM_JOURNAL_REFRESH_DEBOUNCE), so faster
    polling does not thrash the cached journal snapshot.
    """
    time.sleep(random.uniform(lowTimeRange, highTimeRange))


# ---------------------------------------------------------------------------
# Per-tilt-series metadata "sidecar" files
#
# A streaming producer (e.g. ProtComposeTS) writes one small JSON sidecar per
# finished tilt-series next to its <tsId>.ready marker. Downstream consumers
# rebuild the TiltSeries (and its TiltImages) fully in memory from the sidecar,
# so they NEVER open the producer's live SQLite set. This removes the
# cross-process SHARED-read vs producer-EXCLUSIVE-commit contention on the
# single shared `tiltseries.sqlite` (the deadlock under journal_mode=DELETE on
# NFS): the producer writes its DB freely; consumers read finished files only.
#
# The schema is explicit (not a generic getObjDict dump) so it is stable and
# round-trippable: it captures exactly what ProtComposeTS produces for a freshly
# composed, not-yet-aligned TiltSeries. Extend `_ACQ_FIELDS` / the per-image
# fields if a producer needs to publish more.
# ---------------------------------------------------------------------------
TS_META_VERSION = 2  # v2 adds per-tilt-image alignment transform + interpolated flag
# (getter, setter) names on TomoAcquisition that ProtComposeTS populates.
_ACQ_FIELDS = (
    'Voltage', 'Magnification', 'SphericalAberration', 'AmplitudeContrast',
    'DosePerFrame', 'AngleMin', 'AngleMax', 'Step', 'AccumDose', 'TiltAxisAngle',
)


def getTsSidecarPath(streamingDir: str, tsId: str) -> str:
    return join(streamingDir, f'{tsId}{const.TS_META_EXT}')


def tsSidecarExists(streamingDir: str, tsId: str) -> bool:
    return exists(getTsSidecarPath(streamingDir, tsId))


def _acqToDict(acq: TomoAcquisition) -> dict:
    if acq is None:
        return {}
    return {f: getattr(acq, 'get' + f)() for f in _ACQ_FIELDS}


def _dictToAcq(d: dict) -> TomoAcquisition:
    acq = TomoAcquisition()
    for f in _ACQ_FIELDS:
        if d.get(f) is not None:
            getattr(acq, 'set' + f)(d[f])
    return acq


def writeTsSidecar(streamingDir: str, ts: TiltSeries,
                   tiltImages: List[TiltImage]) -> None:
    """Atomically write the metadata sidecar for a composed tilt-series.

    Built entirely from the IN-MEMORY ``ts`` / ``tiltImages`` the producer
    already holds in ``registerOutputs`` — it performs NO database read.
    Written to a temp file and ``os.replace``-d into place so a consumer never
    observes a half-written sidecar. Call this BEFORE touching ``<tsId>.ready``
    so the marker only appears once the sidecar is complete.
    """
    sRate = ts.getSamplingRate()
    data = {
        'version': TS_META_VERSION,
        'tsId': ts.getTsId(),
        'samplingRate': sRate,
        # Alignment/interpolation are needed by downstream consumers that align
        # or track fiducials (e.g. ProtImodFiducialModel): without the per-tilt
        # transforms below, the rebuilt TS would report hasAlignment()==False and
        # no .xf prealignment would be written, breaking autofidseed/beadtrack.
        'interpolated': ts.interpolated() if hasattr(ts, 'interpolated') else False,
        'acquisition': _acqToDict(ts.getAcquisition()),
        'tiltImages': [],
    }
    for ti in tiltImages:
        tiAcq = ti.getAcquisition()
        data['tiltImages'].append({
            'index': ti.getIndex(),
            'fileName': ti.getFileName(),
            'tiltAngle': ti.getTiltAngle(),
            'acquisitionOrder': ti.getAcquisitionOrder(),
            'samplingRate': ti.getSamplingRate(),
            'enabled': ti.isEnabled(),
            'doseInitial': tiAcq.getDoseInitial() if tiAcq else None,
            'accumDose': tiAcq.getAccumDose() if tiAcq else None,
            'oddEven': [ti.getOdd(), ti.getEven()] if ti.hasOddEven() else [],
            # Per-tilt 2D alignment matrix (list-of-lists) so a sidecar-rebuilt TS
            # preserves hasAlignment() and genXfFile can regenerate the .xf.
            'transform': ti.getTransform().getMatrix().tolist() if ti.hasTransform() else None,
        })
    path = getTsSidecarPath(streamingDir, ts.getTsId())
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f)
    os.replace(tmp, path)  # atomic publish of the sidecar


def readTsSidecar(streamingDir: str, tsId: str):
    """Rebuild ``(TiltSeries, [TiltImage])`` fully in memory from the sidecar.

    No SQLite access at all -> no lock contention with the producer. This is the
    drop-in replacement for a consumer's ``fetchNewTs`` + ``loadTiltImgsInMemory``
    on the producer's live set.
    """
    with open(getTsSidecarPath(streamingDir, tsId)) as f:
        data = json.load(f)

    sRate = data.get('samplingRate')
    ts = TiltSeries(tsId=data['tsId'])
    ts.setAcquisition(_dictToAcq(data.get('acquisition', {})))
    if sRate is not None:
        ts.setSamplingRate(sRate)
    if data.get('interpolated'):
        ts.setInterpolated(True)

    tiltImages = []
    for d in data['tiltImages']:
        ti = TiltImage()
        ti.setTsId(data['tsId'])
        ti.setIndex(d['index'])
        ti.setFileName(d['fileName'])
        ti.setTiltAngle(d['tiltAngle'])
        ti.setAcquisitionOrder(d['acquisitionOrder'])
        ti.setSamplingRate(d.get('samplingRate', sRate))
        ti.setEnabled(d.get('enabled', True))
        tiAcq = _dictToAcq(data.get('acquisition', {}))
        if d.get('doseInitial') is not None:
            tiAcq.setDoseInitial(d['doseInitial'])
        if d.get('accumDose') is not None:
            tiAcq.setAccumDose(d['accumDose'])
        ti.setAcquisition(tiAcq)
        if d.get('oddEven'):
            ti.setOddEven(d['oddEven'])
        # Restore the per-tilt alignment transform so the rebuilt TS is
        # equivalent to the producer's DB one (hasAlignment + genXfFile work).
        if d.get('transform') is not None:
            ti.setTransform(Transform(matrix=np.array(d['transform'])))
        tiltImages.append(ti)

    # Mirror TiltSeriesBase.append: the TS is aligned iff its tilt-images carry
    # transforms. Set the flag explicitly because the consumer attaches items via
    # setInMemoryTiltImages (not append, which is what normally sets it).
    ts.setHasAlignment(any(ti.hasTransform() for ti in tiltImages))
    return ts, tiltImages


# ---------------------------------------------------------------------------
# Per-CTF-tomo-series metadata "sidecar" files (the CTF analog of the TS sidecar
# above). A streaming producer publishes one JSON sidecar per finished
# CTFTomoSeries; downstream consumers rebuild the CTFTomoSeries (and its CTFTomos)
# fully in memory from it (see SetOfCTFTomoSeries.fetchNewCtfs), so they NEVER
# open the producer's live SQLite set. Explicit, round-trippable schema (not a
# generic getObjDict dump) capturing the standard per-tilt CTF estimation values.
# ---------------------------------------------------------------------------
CTF_META_VERSION = 1


def getCtfSidecarPath(streamingDir: str, tsId: str) -> str:
    return join(streamingDir, f'{tsId}{const.CTF_META_EXT}')


def ctfSidecarExists(streamingDir: str, tsId: str) -> bool:
    return exists(getCtfSidecarPath(streamingDir, tsId))


def writeCtfSidecar(streamingDir: str, ctfTomoSeries: CTFTomoSeries,
                    ctfTomos: List[CTFTomo]) -> None:
    """Atomically write the metadata sidecar for a CTFTomoSeries.

    Built entirely from the IN-MEMORY ``ctfTomoSeries`` / ``ctfTomos`` the
    producer already holds — it performs NO database read. Written to a temp file
    and ``os.replace``-d into place so a consumer never observes a half-written
    sidecar. Call this BEFORE publishing the tsId to the stream journal so the
    journal id only appears once the sidecar is complete. Mirrors writeTsSidecar.
    """
    # Note: the CTFTomoSeries-level defocus-deviation flags are deliberately not
    # serialized: their getters/setters are stubs in the data model (no real
    # state), so the sidecar only carries genuinely round-trippable data.
    data = {
        'version': CTF_META_VERSION,
        'tsId': ctfTomoSeries.getTsId(),
        'ctfTomos': [],
    }
    for ctf in ctfTomos:
        data['ctfTomos'].append({
            'index': ctf.getIndex(),
            'acquisitionOrder': ctf.getAcquisitionOrder(),
            'enabled': ctf.isEnabled(),
            'defocusU': ctf.getDefocusU(),
            'defocusV': ctf.getDefocusV(),
            'defocusAngle': ctf.getDefocusAngle(),
            'resolution': ctf.getResolution(),
            'fitQuality': ctf.getFitQuality(),
            'phaseShift': ctf.getPhaseShift() if ctf.hasPhaseShift() else None,
        })
    path = getCtfSidecarPath(streamingDir, ctfTomoSeries.getTsId())
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f)
    os.replace(tmp, path)  # atomic publish of the sidecar


def readCtfSidecar(streamingDir: str, tsId: str):
    """Rebuild ``(CTFTomoSeries, [CTFTomo])`` fully in memory from the sidecar.

    No SQLite access at all -> no lock contention with the producer. This is the
    CTF analog of readTsSidecar and the building block of
    SetOfCTFTomoSeries.fetchNewCtfs.
    """
    with open(getCtfSidecarPath(streamingDir, tsId)) as f:
        data = json.load(f)

    cts = CTFTomoSeries(tsId=data['tsId'])
    cts.setTsId(data['tsId'])

    ctfTomos = []
    for d in data['ctfTomos']:
        ctf = CTFTomo()
        ctf.setIndex(d['index'])
        ctf.setAcquisitionOrder(d['acquisitionOrder'])
        ctf.setEnabled(d.get('enabled', True))
        defU, defV, defAngle = d.get('defocusU'), d.get('defocusV'), d.get('defocusAngle')
        if None not in (defU, defV, defAngle):
            ctf.setStandardDefocus(defU, defV, defAngle)
        if d.get('resolution') is not None:
            ctf.setResolution(d['resolution'])
        if d.get('fitQuality') is not None:
            ctf.setFitQuality(d['fitQuality'])
        if d.get('phaseShift') is not None:
            ctf.setPhaseShift(d['phaseShift'])
        ctfTomos.append(ctf)
    return cts, ctfTomos


# ---------------------------------------------------------------------------
# Per-landmark-model metadata "sidecar" files (the LandmarkModel analog of the
# TS/CTF sidecars above). A streaming producer publishes one JSON sidecar per
# finished LandmarkModel; a downstream consumer can rebuild the LandmarkModel
# fully in memory from it WITHOUT opening the producer's live SetOfLandmarkModels
# SQLite. The landmark coordinates themselves are NOT embedded here: they already
# live in the referenced '.sfid' file on shared storage (written by
# LandmarkModel.addLandmark) and are read from it lock-free on demand
# (LandmarkModel.retrieveInfoTable). This sidecar carries only the round-trippable
# object metadata needed to reconstruct the LandmarkModel wrapper.
# ---------------------------------------------------------------------------
LANDMARK_META_VERSION = 1


def getLandmarkSidecarPath(streamingDir: str, tsId: str) -> str:
    return join(streamingDir, f'{tsId}{const.LANDMARK_META_EXT}')


def landmarkSidecarExists(streamingDir: str, tsId: str) -> bool:
    return exists(getLandmarkSidecarPath(streamingDir, tsId))


def writeLandmarkSidecar(streamingDir: str, landmarkModel: LandmarkModel) -> None:
    """Atomically write the metadata sidecar for a LandmarkModel.

    Built entirely from the IN-MEMORY ``landmarkModel`` the producer already
    holds — it performs NO database read. Written to a temp file and
    ``os.replace``-d into place so a consumer never observes a half-written
    sidecar. Call this BEFORE publishing the tsId to the stream journal so the
    journal id only appears once the sidecar is complete. Mirrors
    writeTsSidecar / writeCtfSidecar.

    The landmark rows are not serialized here (see module note above): they are
    in the referenced ``fileName`` (.sfid) file.
    """
    data = {
        'version': LANDMARK_META_VERSION,
        'tsId': landmarkModel.getTsId(),
        'fileName': landmarkModel.getFileName(),  # .sfid file with the landmark rows
        'modelName': landmarkModel.getModelName(),  # .fid model file
        'size': landmarkModel.getSize(),  # bead diameter (Å)
        'count': landmarkModel.getCount(),  # number of chains/landmarks
        'applyTSTransformation': landmarkModel.applyTSTransformation(),
        'hasResidualInfo': landmarkModel.hasResidualInfo().get(),
    }
    path = getLandmarkSidecarPath(streamingDir, landmarkModel.getTsId())
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f)
    os.replace(tmp, path)  # atomic publish of the sidecar


def readLandmarkSidecar(streamingDir: str, tsId: str) -> LandmarkModel:
    """Rebuild a ``LandmarkModel`` fully in memory from the sidecar.

    No SQLite access at all -> no lock contention with the producer. This is the
    LandmarkModel analog of readTsSidecar / readCtfSidecar. The associated
    tilt-series pointer is intentionally left unset (it is not persisted on the
    item, ``objDoStore=False``); a consumer associates it via the set's
    ``completeLandmarkModel``.
    """
    with open(getLandmarkSidecarPath(streamingDir, tsId)) as f:
        data = json.load(f)

    lm = LandmarkModel(tsId=data['tsId'],
                       fileName=data.get('fileName'),
                       modelName=data.get('modelName'),
                       size=data.get('size'),
                       applyTSTransformation=data.get('applyTSTransformation', True),
                       hasResidualInfo=data.get('hasResidualInfo', False))
    lm.setTsId(data['tsId'])
    lm.setCount(data.get('count', 0))
    return lm
