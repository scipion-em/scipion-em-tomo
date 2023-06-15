# coding=utf-8
# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
from os.path import basename

from imod import utils
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pyworkflow.plugin import Domain

from ..objects import SetOfCoordinates3D, SetOfTiltSeriesCoordinates, TiltSeriesCoordinate
from .protocol_base import ProtTomoImportFiles
from ..convert import TomoImport, EmTableCoordImport
from ..utils import existsPlugin

import tomo.constants as const

IMPORT_FROM_TXT = 'txt'
IMPORT_FROM_IMOD = 'imod'
OUTPUT_TS_COORDINATES_NAME = "TiltSeriesCoordinates"


class ProtImportTiltSeriesCoordinates(ProtTomoImportFiles):
    """Protocol to import a set of tilt-series coordinates 3D"""
    _outputClassName = 'SetOfCoordinates3D'
    _label = 'import tilt-series coordinates'
    _devStatus = BETA

    def _getImportChoices(self):
        """ Return a list of possible choices from which the import can be done.
        """
        importChoices = [IMPORT_FROM_TXT]
        if existsPlugin('imod'):
            importChoices.append(IMPORT_FROM_IMOD)
        return importChoices

    def _defineAcquisitionParams(self, form):
        pass

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

        self.OUTPUT_PREFIX = "outputCoordinates"

        self.TiltSeriesCoordinates = None

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

    def _insertAllSteps(self):
        self.fileList = [file for file, _ in self.iterFiles()]

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('importCoordinatesStep',
                                     ts.getObjId())

    # --------------------------- STEPS functions -----------------------------

    def importCoordinatesStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        self.getOutputSetOfTiltSeriesCoordinates(self.inputSetOfTiltSeries.get())

        for coordFilePath in self.fileList:
            if tsId in coordFilePath:
                break

        coordList, xDim, yDim = utils.format3DCoordinatesList(coordFilePath)

        for element in coordList:
            newCoord3D = TiltSeriesCoordinate()
            newCoord3D.setTsId(ts.getTsId())
            newCoord3D.setPosition(element[0] - (xDim / 2),
                                   element[1] - (yDim / 2),
                                   element[2],
                                   sampling_rate=ts.getSamplingRate())

            self.TiltSeriesCoordinates.append(newCoord3D)
        self.TiltSeriesCoordinates.write()
        self._store()

    # ------------------ UTILS functions --------------------------------------
    def getImportFrom(self):
        importFrom = self._getImportChoices()[self.importFrom.get()]

        return importFrom

    def getOutputSetOfTiltSeriesCoordinates(self, setOfTiltSeries=None):

        if self.TiltSeriesCoordinates:
            self.TiltSeriesCoordinates.enableAppend()

        else:
            outputSetOfCoordinates3D = SetOfTiltSeriesCoordinates.create(self._getPath(),
                                                                         suffix='tsCoords')

            outputSetOfCoordinates3D.setSetOfTiltSeries(setOfTiltSeries)
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TS_COORDINATES_NAME: outputSetOfCoordinates3D})
            self._defineSourceRelation(setOfTiltSeries, outputSetOfCoordinates3D)

        return self.TiltSeriesCoordinates


    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods

    def _validate(self):
        validateMsgs = []

        match = False

        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()

            for coordsFilePath, _ in self.iterFiles():
                if tsId in coordsFilePath:
                    match = True
                    break

            if not match:
                validateMsgs.append("No coordinates file found for tilt-series %s: image file is %s and have not "
                                    "found its exact match." % (tsId, ts.getFileName()))

            match = False

        return validateMsgs

    def _warnings(self):
        pass
        warnings = []

        if not existsPlugin('imod'):
            warnings.append('Plugin *scipion-em-imod* has not being installed. Please, install the Plugin to '
                            'import IMOD related formats (currently supported formats: ".txt"). Otherwise, the '
                            'protocol may have unexpected outputs if Eman files are attempted to be imported.\n')

        return warnings

