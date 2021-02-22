# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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


import os

import pyworkflow.wizard as pwizard

import pwem.objects as emobj
import pwem.convert as emconv

from tomo import protocols as tomoprot


class ImportOriginTomoWizard(pwizard.Wizard):
    _targets = [(tomoprot.ProtImportTomograms, ['x', 'y', 'z'])]

    def show(self, form, *params):
        protocol = form.protocol
        filesPath = protocol.filesPath.get()
        filesPattern = protocol.filesPattern.get()
        if filesPattern:
            fullPattern = os.path.join(filesPath, filesPattern)
        else:
            fullPattern = filesPath

        sampling = protocol.samplingRate.get()
        for fileName, fileId in protocol.iterFiles():
            inputVol = emobj.Volume()
            inputVol.setFileName(fileName)
            if ((str(fullPattern)).endswith('mrc') or
                    (str(fullPattern)).endswith('map')):
                ccp4header = emconv.Ccp4Header(fileName, readHeader=True)
                x, y, z = ccp4header.getOrigin(changeSign=True)  # In Angstroms
            else:
                x, y, z = self._halfOriginCoordinates(inputVol, sampling)

            form.setVar('x', x)
            form.setVar('y', y)
            form.setVar('z', z)

    @classmethod
    def _halfOriginCoordinates(cls, volume, sampling):
        xdim, ydim, zdim = volume.getDim()
        if zdim > 1:
            zdim = zdim / 2.
        x = xdim / 2. * sampling
        y = ydim / 2. * sampling
        z = zdim * sampling
        return x, y, z