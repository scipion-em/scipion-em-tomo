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

from pyworkflow import BETA
import pyworkflow.utils as pwutils

import pyworkflow.protocol.params as params
from pyworkflow.plugin import Domain

from ..objects import SetOfCoordinates3D
from .protocol_base import ProtTomoImportFiles
from ..convert import TomoImport, EmTableCoordImport
from ..utils import existsPlugin

import tomo.constants as const

IMPORT_FROM_TXT = 'txt'
IMPORT_FROM_IMOD = 'imod'


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

    def  _defineAcquisitionParams(self, form):
        pass

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)
        self.OUTPUT_PREFIX = "outputCoordinates"

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)

    def _insertAllSteps(self):
        self._insertFunctionStep('importCoordinatesStep')

    # --------------------------- STEPS functions -----------------------------

    def importCoordinatesStep(self, samplingRate):
        pass

    # ------------------ UTILS functions --------------------------------------
    def getImportFrom(self):
        importFrom = self._getImportChoices()[self.importFrom.get()]

        return importFrom

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods

    def _validate(self):
        errors = []

        return errors

    def _warnings(self):
        warnings = []
        tomoFiles = [pwutils.removeBaseExt(file) for file in self.importTomograms.get().getFiles()]
        coordFiles = [pwutils.removeBaseExt(file) for file, _ in self.iterFiles()]
        numberMatches = len(set(tomoFiles) & set(coordFiles))
        if not existsPlugin('imod'):
            warnings.append('Plugin *scipion-em-imod* has not being installed. Please, install the Plugin to '
                            'import IMOD related formats (currently supported formats: ".txt"). Otherwise, the '
                            'protocol may have unexpected outputs if Eman files are attempted to be imported.\n')
        if numberMatches < max(len(tomoFiles), len(coordFiles)):
            warnings.append("Couldn't find a correspondence between all cordinate and tomogram files. "
                            "Association is performed in terms of the file name of the Tomograms and the coordinates. "
                            "(without the extension). For example, if a Tomogram file is named Tomo_1.mrc, the coordinate "
                            "file to be associated to it should be named Tomo_1.ext (being 'ext' any valid extension "
                            "- '.txt', '.tbl', '.json').\n")
            # mismatches_coords = set(coordFiles).difference(tomoFiles)
            # if mismatches_coords:
            #     warnings.append("The following coordinate files will not be associated to any Tomogram "
            #                     "(name without extension):")
            #     for file in mismatches_coords:
            #         warnings.append("\t%s" % file)
            #     warnings.append("\n")
            # mismatches_tomos = set(tomoFiles).difference(coordFiles)
            # if mismatches_tomos:
            #     warnings.append("The following Tomogram files will not be associated to any coordinates "
            #                     "(name without extension):")
            #     for file in mismatches_tomos:
            #         warnings.append("\t%s" % file)
            #     warnings.append("\n")
        return warnings

