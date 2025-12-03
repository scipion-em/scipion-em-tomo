# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

## NOTE: All this code was moved from https://github.com/I2PC/scipion-em-xmipptomo/blob/v3.23.03.2/xmipptomo/protocols/protocol_fit_ellipsoid.py

import numpy as np
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import MeshPoint, Ellipsoid, SubTomogram, SetOfCoordinates3D, Coordinate3D, SetOfMeshes
from tomo.utils import fit_ellipsoid, generatePointCloud
import tomo.constants as const


class FitVesicleOutputs(Enum):
    meshes = SetOfMeshes


class TomoProtFitEllipsoid(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfSubtomograms with coordinates assigned or a SetOfCoordinates3D, to a vesicle
    (ellipsoid), defining regions of interest (SetOfMeshes) for each vesicle as output."""

    _label = 'fit vesicles'
    _devStatus = BETA
    _possibleOutputs = FitVesicleOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass="SetOfSubTomograms, SetOfCoordinates3D",
                      label='Subtomograms/Coordinates3D',
                      help='Subtomograms or coordinates3D picked in vesicles. If there are more than one vesicle per '
                           'tomogram, input subtomograms or coordinates should have assigned groupId.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.fitEllipsoidStep)

    # --------------------------- STEPS functions -------------------------------
    def fitEllipsoidStep(self):
        inCoords = self.input.get()
        if type(inCoords) is not SetOfCoordinates3D:  # It may be a SetOfSubTomograms
            inCoords = inCoords.getCoordinates3D()
        precedentsDict = inCoords.getPrecedentsInvolved()
        outSet = self._createSetOfMeshes(inCoords.getPrecedents())
        outSet.copyInfo(inCoords)
        vesicleDict = {}
        for tomoId, tomo in precedentsDict.items():
            self.info("Fitting tomogram %s" % tomo)
            tomoX, tomoY, tomoZ = tomo.getDim()
            groupIds = inCoords.getUniqueValues(Coordinate3D.GROUP_ID_ATTR, where='%s="%s"' % (Coordinate3D.TOMO_ID_ATTR, tomoId))
            for groupId in groupIds:
                vesicleDict['%s_%s' % (tomoId, groupId)] = {'x': [], 'y': [], 'z': [], 'tomo': tomo}
            for coord in inCoords.iterCoordinates(volume=tomo):
                key = '%s_%s' % (tomoId, coord.getGroupId())
                vesicleDict[key]['x'].append(float(coord.getX(const.BOTTOM_LEFT_CORNER)) / tomoX)
                vesicleDict[key]['y'].append(float(coord.getY(const.BOTTOM_LEFT_CORNER)) / tomoY)
                vesicleDict[key]['z'].append(float(coord.getZ(const.BOTTOM_LEFT_CORNER)) / tomoZ)

        for vesicleInd, coordDict in enumerate(vesicleDict.values()):
            x = coordDict['x']
            y = coordDict['y']
            z = coordDict['z']
            tomo = coordDict['tomo']
            [center, radii, v, _, chi2] = fit_ellipsoid(np.array(x), np.array(y), np.array(z))
            algDesc = '%f*x*x + %f*y*y + %f*z*z + 2*%f*x*y + 2*%f*x*z + 2*%f*y*z + 2*%f*x + 2*%f*y + 2*%f*z + %f ' \
                      '= 0' % (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])
            adjEllipsoid = Ellipsoid()
            adjEllipsoid.setAlgebraicDesc(algDesc)
            adjEllipsoid.setCenter(str(center))
            adjEllipsoid.setRadii(str(radii))
            self.debug(algDesc)
            self.debug('Chi2: %s' % chi2)
            pointCloud = generatePointCloud(v, tomo.getDim())
            if not pointCloud:
                raise Exception("It does not seem like any output is produced!")

            for point in pointCloud:
                meshPoint = MeshPoint()
                meshPoint.setVolume(tomo)
                meshPoint.setX(point[0], const.BOTTOM_LEFT_CORNER)
                meshPoint.setY(point[1], const.BOTTOM_LEFT_CORNER)
                meshPoint.setZ(point[2], const.BOTTOM_LEFT_CORNER)
                meshPoint.setGroupId(vesicleInd)
                meshPoint.setDescription(adjEllipsoid)
                meshPoint.setVolumeName(tomo.getTsId())
                outSet.append(meshPoint)

        outSet.setNumberOfMeshes(len(vesicleDict.keys()))

        self._defineOutputs(**{self._possibleOutputs.meshes.name: outSet})
        self._defineSourceRelation(self.input.get(), outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        input = self.input.get()
        if input.getSize() < 9:
            validateMsgs.append('At least 9 subtomograms/coordinates are required to fit a unique ellipsoid')
        if self._getInputisSubtomo(self.input.get().getFirstItem()):
            if not input.getFirstItem().hasCoordinate3D():
                validateMsgs.append('Subtomograms should have coordinates')
        return validateMsgs

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output vesicles not ready yet.")
        else:
            summary.append("%i subtomograms/coordinates3D\n%i vesicles adjusted" %
                           (self.input.get().getSize(), self.meshes.getNumberOfMeshes()))
        return summary

    def _methods(self):
        if not self.isFinished():
            return ["Protocol not finished yet"]
        else:
            return ["Fit an ellipsoid and compute a 3D set of points (mesh) for each of the %d vesicles adjusted." %
                    self.outputMeshes.getNumberOfMeshes()]

    # --------------------------- UTILS functions --------------------------------------------
    @staticmethod
    def _getInputisSubtomo(item):
        if isinstance(item, SubTomogram):
            return True
        else:
            return False

    def _getCoor(self, item):
        if self._getInputisSubtomo(item):
            coor = item.getCoordinate3D()
        else:
            coor = item
        return coor


# To keep backwards compatibility in existing projects
class XmippProtFitEllipsoid(TomoProtFitEllipsoid):
    @classmethod
    def isDisabled(cls):
        return True
