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

import pyworkflow.utils as pwutils
from pwem.objects.data_tiltpairs import Angles
from pwem.emlib.metadata import (MetaData, MDL_XCOOR, MDL_YCOOR, MDL_ZCOOR)

class TomoImport:

    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = protocol.getCopyOrLink()

    def importCoordinates3D(self, fileName, addCoordinate):
        from tomo.objects import Coordinate3D
        if pwutils.exists(fileName):
            ext = pwutils.getExt(fileName)

        if ext == ".txt":
            md = MetaData()
            md.readPlain(fileName, "xcoor ycoor zcoor")
            for objId in md:
                x = md.getValue(MDL_XCOOR, objId)
                y = md.getValue(MDL_YCOOR, objId)
                z = md.getValue(MDL_ZCOOR, objId)
                coord = Coordinate3D()
                coord.setPosition(x, y, z)
                addCoordinate(coord)

        else:
            raise Exception('Unknown extension "%s" to import Eman coordinates' % ext)