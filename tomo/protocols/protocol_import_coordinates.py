# coding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team  (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
from pyworkflow.protocol import LEVEL_ADVANCED
from pyworkflow.utils import replaceBaseExt, removeBaseExt

from ..objects import SetOfCoordinates3D
from .protocol_base import ProtTomoImportFiles
from ..convert import TomoImport, EmTableCoordImport
from ..utils import existsPlugin

import tomo.constants as const

IMPORT_FROM_AUTO = 'auto'
IMPORT_FROM_TXT = 'txt'
IMPORT_FROM_EMAN = 'eman'
IMPORT_FROM_DYNAMO = 'dynamo'
IMPORT_FROM_CBOX = 'cbox'


class ProtImportCoordinates3D(ProtTomoImportFiles):

    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfCoordinates3D'
    _label = 'import coordinates 3D'

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        importChoices = [IMPORT_FROM_AUTO, IMPORT_FROM_TXT, IMPORT_FROM_CBOX]
        if existsPlugin('emantomo'):
            importChoices.append(IMPORT_FROM_EMAN)
        if existsPlugin('dynamo'):
            importChoices.append(IMPORT_FROM_DYNAMO)
        return importChoices

    def _getDefaultChoice(self):
        return 0

    def __init__(self, **args):
        super().__init__(**args)
        self.OUTPUT_PREFIX = "outputCoordinates"
        self.scaleFactor = 1.0
        self.tomoSRate = None

    def _defineParams(self, form):
        super()._defineImportParams(form)
        form.addParam('samplingRate', params.FloatParam,
                      label='Coordinates sampling rate [Å/pix] (opt.)',
                      allowsNull=True,
                      help="If empty, the coordinates' sampling rate will be considered to be the same as the "
                           "tomograms'.\n"
                           "*IMPORTANT*: If a value is provided, the ratio of both tomograms and coordinates sampling "
                           "rate will be used to scale the coordinates properly to the tomograms introduced.")

        form.addParam('boxSize', params.IntParam,
                      label='Box Size [pix]',
                      default=20,
                      help='It will be re-scaled to the tomogram size considering the coordinates and tomograms '
                           'ratio between their corresponding sampling rates.')

        form.addParam('importTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input tomograms',
                      help='Select the tomograms to which the coordinates should be referred to.\n'
                           'The file names of the tomogram and coordinate files must be the same.')

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.importCoordinatesStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tomoSRate = self.importTomograms.get().getSamplingRate()
        coordsSRate = self.samplingRate.get()
        if coordsSRate:
            self.scaleFactor = coordsSRate / tomoSRate
        self.tomoSRate = tomoSRate

    def importCoordinatesStep(self):
        importTomograms = self.importTomograms.get()
        coordsSet = SetOfCoordinates3D.create(self.getPath(), template='coordinates%s.sqlite')
        coordsSet.setSamplingRate(self.tomoSRate)
        coordsSet.setPrecedents(importTomograms)
        coordsSet.setBoxSize(self.boxSize.get() * self.scaleFactor)

        ci = self.getImportClass()
        for tomo in importTomograms.iterItems():
            tomoName = removeBaseExt(tomo.getFileName())
            for coordFile, fileId in self.iterFiles():
                fileName = removeBaseExt(coordFile)
                if tomo is not None and tomoName == fileName:
                    # Parse the coordinates in the given format for this micrograph
                    if self.getImportFrom() in [IMPORT_FROM_EMAN, IMPORT_FROM_TXT, IMPORT_FROM_CBOX]:
                        def addCoordinate(coord, x, y, z):
                            coord.setVolume(tomo.clone())

                            x = x * self.scaleFactor
                            y = y * self.scaleFactor
                            z = z * self.scaleFactor

                            coord.setPosition(x, y, z, const.BOTTOM_LEFT_CORNER)
                            coordsSet.append(coord)
                        ci.importCoordinates3D(coordFile, addCoordinate)
                    elif self.getImportFrom() == IMPORT_FROM_DYNAMO:
                        ci(coordFile, coordsSet, tomo.clone(), scaleFactor=self.scaleFactor)

        args = {self.OUTPUT_PREFIX: coordsSet}
        self._defineOutputs(**args)
        self._defineSourceRelation(self.importTomograms, coordsSet)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return self.hasAttribute(self.OUTPUT_PREFIX)

    def _getCoordsMessage(self):
        return "Coordinates %s" % self.getObjectTag(self.OUTPUT_PREFIX)

    def _summary(self):
        summary = []
        if self._hasOutput():
            summary.append("%s imported from:\n%s" % (self._getCoordsMessage(), self.getPattern()))
            tomoSRate = self.importTomograms.get().getSamplingRate()
            coordsSRate = self.samplingRate.get()
            if coordsSRate:
                scaleFactor = coordsSRate / tomoSRate
                summary.append("*Coordinates were scaled by a factor of %.2f.*" % scaleFactor)
        return summary

    def _methods(self):
        methods = []
        if self._hasOutput():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getCoordsMessage(), self.samplingRate.get()),)
        return methods

    def _getVolumeFileName(self, fileName, extension=None):
        if extension is not None:
            baseFileName = "import_" + str(basename(fileName)).split(".")[0] + ".%s" % extension
        else:
            baseFileName = "import_" + str(basename(fileName)).split(":")[0]

        return self._getExtraPath(baseFileName)

    def _validate(self):
        errors = []
        try:
            next(self.iterFiles())
        except StopIteration:
            errors.append('No files matching the pattern %s were found.' % self.getPattern())
        else:
            tomoFiles = [pwutils.removeBaseExt(file) for file in self.importTomograms.get().getFiles()]
            coordFiles = [pwutils.removeBaseExt(file) for file, _ in self.iterFiles()]
            numberMatches = len(set(tomoFiles) & set(coordFiles))
            if numberMatches == 0:
                errors.append("Cannot relate tomogram and coordinate files. In order to stablish a "
                              "relation, the filename of the corresponding tomogram and coordinate "
                              "files must be equal.")
        return errors

    def _warnings(self):
        warnings = []
        tomoFiles = [pwutils.removeBaseExt(file) for file in self.importTomograms.get().getFiles()]
        coordFiles = [pwutils.removeBaseExt(file) for file, _ in self.iterFiles()]
        numberMatches = len(set(tomoFiles) & set(coordFiles))
        if not existsPlugin('emantomo'):
            warnings.append('Plugin *scipion-em-emantomo* has not being installed. Please, install the Plugin to '
                            'import Eman related formats (currently supported formats: ".json"). Otherwise, the protocol '
                            'may have unexpected outputs if Eman files are attempted to be imported.\n')
        if not existsPlugin('dynamo'):
            warnings.append('Plugin *scipion-em-dynamo* has not being installed. Please, install the Plugin to '
                            'import Dynamo related formats (currently supported formats: ".tbl"). Otherwise, the protocol '
                            'may have unexpected outputs if Dynamo files are attempted to be imported.\n')
        if numberMatches < max(len(tomoFiles), len(coordFiles)):
            warnings.append("Couldn't find a correspondence between all cordinate and tomogram files. "
                            "Association is performed in terms of the file name of the Tomograms and the coordinates. "
                            "(without the extension). For example, if a Tomogram file is named Tomo_1.mrc, the coordinate "
                            "file to be associated to it should be named Tomo_1.ext (being 'ext' any valid extension "
                            "- '.txt', '.tbl', '.json').\n")
            mismatches_coords = set(coordFiles).difference(tomoFiles)
            if mismatches_coords:
                warnings.append("The following coordinate files will not be associated to any Tomogram "
                                "(name without extension):")
                for file in mismatches_coords:
                    warnings.append("\t%s" % file)
                warnings.append("\n")
            mismatches_tomos = set(tomoFiles).difference(coordFiles)
            if mismatches_tomos:
                warnings.append("The following Tomogram files will not be associated to any coordinates "
                                "(name without extension):")
                for file in mismatches_tomos:
                    warnings.append("\t%s" % file)
                warnings.append("\n")
        return warnings

    # ------------------ UTILS functions --------------------------------------
    def getImportFrom(self):
        importFrom = self._getImportChoices()[self.importFrom.get()]
        if importFrom == IMPORT_FROM_AUTO:
            importFrom = self.getFormat()
        return importFrom

    def getFormat(self):
        for coordFile, _ in self.iterFiles():
            if coordFile.endswith('.txt'):
                return IMPORT_FROM_TXT
            elif coordFile.endswith('.json') and existsPlugin('emantomo'):
                return IMPORT_FROM_EMAN
            elif coordFile.endswith('.tbl') and existsPlugin('dynamo'):
                return IMPORT_FROM_DYNAMO
            elif coordFile.endswith('.cbox'):
                return IMPORT_FROM_CBOX
        return -1

    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        importFrom = self.getImportFrom()

        if importFrom == IMPORT_FROM_EMAN:
            EmanImport = Domain.importFromPlugin('emantomo.convert', 'EmanTomoImport',
                                                 errorMsg='Eman is needed to import .json or '
                                                          '.box files',
                                                 doRaise=True)
            return EmanImport(self, None)

        elif importFrom == IMPORT_FROM_DYNAMO:
            readDynCoord = Domain.importFromPlugin("dynamo.convert.convert", "readDynCoord")
            return readDynCoord

        elif importFrom == IMPORT_FROM_CBOX:
            return EmTableCoordImport("cryolo", "CoordinateX", "CoordinateY", "CoordinateZ")

        elif importFrom == IMPORT_FROM_TXT:
            return TomoImport(self)

        else:
            self.importFilePath = ''
            return None
