# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import NEW
from pwem.protocols import EMProtocol

from tomo.constants import BOTTOM_LEFT_CORNER, TR_DYNAMO
from tomo.utils import existsPlugin

EXPORT_TO_TXT = 'txt'
EXPORT_TO_STAR = 'relion'
EXPORT_TO_EMAN = 'eman'
EXPORT_TO_DYNAMO = 'dynamo'
EXPORT_TO_CBOX = 'cbox'


class ProtExportCoordinates3D(EMProtocol):
    """ Export 3D subtomogram coordinates to be used outside Scipion. """

    _label = 'export 3D coordinates'
    _devStatus = NEW

    @staticmethod
    def _getExportChoices():
        """ Return a list of possible choices for export.
        """
        exportChoices = [EXPORT_TO_TXT]
        if existsPlugin('reliontomo'):
            exportChoices.append(EXPORT_TO_STAR)
        if existsPlugin('emantomo'):
            exportChoices.append(EXPORT_TO_EMAN)
        if existsPlugin('dynamo'):
            exportChoices.append(EXPORT_TO_DYNAMO)
        if existsPlugin('sphire'):
            exportChoices.append(EXPORT_TO_CBOX)

        return exportChoices

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        exportChoices = self._getExportChoices()

        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      important=True,
                      label="Input 3D coordinates")
        form.addParam('outputFormat', params.EnumParam,
                      choices=exportChoices, default=0,
                      label='Export to',
                      help='Select the output format')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        objId = self.inputCoordinates.get().getObjId()
        format = self._getExportChoices()[self.outputFormat.get()]
        self._insertFunctionStep(self.exportCoordsStep, format, objId)

    # --------------------------- STEPS functions -----------------------------
    def exportCoordsStep(self, format, coordsId):
        inputCoords = self.inputCoordinates.get()
        tomoIds = inputCoords.getTSIds()
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())

        if format == EXPORT_TO_TXT:
            def _writeFunc(coord, f):
                x, y, z = map(int, coord.getPosition(BOTTOM_LEFT_CORNER))
                f.write(f"{x} {y} {z}\n")

            self._writeTxt(inputCoords, "txt", _writeFunc)

        elif format == EXPORT_TO_STAR:
            from reliontomo.convert import writeSetOfCoordinates
            writeSetOfCoordinates(inputCoords,
                                  self._getExportPath("coords.star"),
                                  tomoIds,
                                  sRate=inputCoords.getSamplingRate(),
                                  coordsScale=1)

        elif format == EXPORT_TO_EMAN:
            from emantomo.convert import setCoords3D2Jsons
            json_files = [self._getExportPath(f"{tsId}_info.json") for tsId in tomoIds]
            setCoords3D2Jsons(json_files, inputCoords)

        elif format == EXPORT_TO_DYNAMO:
            from dynamo.convert import matrix2eulerAngles

            def _writeFunc(coord, f):
                x, y, z = coord.getPosition(BOTTOM_LEFT_CORNER)
                # Get alignment information
                #if coord.hasTransform():  # FIXME
                #    tdrot, tilt, narot, shiftx, shifty, shiftz = matrix2eulerAngles(coord.getMatrix())
                #else:
                tdrot, tilt, narot, shiftx, shifty, shiftz = 0, 0, 0, 0, 0, 0
                f.write(f"{coord.getObjId()} 1 1 {shiftx} {shifty} {shiftz} "
                        f"{tdrot} {tilt} {narot} 0 0 0 1 0 0 0 0 0 0 0 0 1 0 "
                        f"{x} {y} {z} 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")

            self._writeTxt(inputCoords, "tbl", _writeFunc)

        elif format == EXPORT_TO_CBOX:
            from sphire.convert import writeSetOfCoordinates3D
            writeSetOfCoordinates3D(self._getExportPath(), inputCoords)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        return validateMsgs

    def _summary(self):
        summary = []

        summary.append(f"Output is written to: \n"
                       f"{os.path.abspath(self._getExportPath())}\n")

        return summary
    
    # --------------------------- UTILS functions -----------------------------
    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)

    def _writeTxt(self, inputCoords, ext="txt", writeCoord=None):
        """ Iterate over coords by tomoId and write output. """
        f = None
        lastTomoId = None
        for coord in inputCoords.iterCoordinates(orderBy="_tomoId"):
            tomoId = coord.getTomoId()
            if tomoId != lastTomoId:
                if f:  # we need to close previous opened file
                    f.close()
                f = open(self._getExportPath(f"{tomoId}.{ext}"), "w")
                lastTomoId = tomoId
            writeCoord(coord, f)
        if f:
            f.close()
