# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Oier Lauzirika Zarrabeitia
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
import numpy as np
import matplotlib.pyplot as plt

from pyworkflow.viewer import DESKTOP_TKINTER
from pyworkflow.protocol.params import StringParam

from pwem.viewers.viewer_base import EmProtocolViewer

import tomo.constants as const
from tomo.objects import SetOfCoordinates3D, Coordinate3D
from tomo.protocols import ProtGeometricPicking

class GeometricPickingViewer(EmProtocolViewer):
    """ Viewer for oriented coordinates"""
    _label = 'Oriented coodinate viewer'
    _targets = [ProtGeometricPicking]
    _environments = [DESKTOP_TKINTER]

    def _defineParams(self, form):
        form.addSection(label='Viewer')
        form.addParam('view', StringParam, label='Show coordinates')

    def _getVisualizeDict(self):
        return {'view': self._showCoordinates}

    def _showCoordinates(self, e):
        OUTPUT_NAME = 'Coordinates'
        coordinates: SetOfCoordinates3D = getattr(self.protocol, OUTPUT_NAME)
        tsId: str = self.view.get()
        
        positions = []
        normals = []
        query = "%s='%s'" % (Coordinate3D.TOMO_ID_ATTR, tsId)
        for coordinate in coordinates.iterItems(where=query):
            matrix = coordinate.getMatrix()
            normal = matrix[2,:3]
            position = coordinate.getPosition(const.SCIPION)
            positions.append(position)
            normals.append(normal)
            
        positions = np.array(positions)
        normals = np.array(normals)
        normals *= coordinates.getBoxSize() / 2
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.quiver(
            positions[:,0], positions[:,1], positions[:,2],
            normals[:,0], normals[:,1], normals[:,2]
        )
        
        lim = np.max(abs(positions))
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        
        return [fig]
