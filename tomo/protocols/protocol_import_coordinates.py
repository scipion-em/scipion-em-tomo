# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
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
from ..convert import TomoImport
from ..utils import existsPlugin

import tomo.constants as const


class ProtImportCoordinates3D(ProtTomoImportFiles):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfCoordinates3D'
    _label = 'import set of coordinates 3D'
    _devStatus = BETA

    IMPORT_FROM_AUTO = 'auto'
    IMPORT_FROM_TXT = 'txt'
    IMPORT_FROM_EMAN = 'eman'
    IMPORT_FROM_DYNAMO = 'dynamo'

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        importChoices = ['auto', 'txt']
        if existsPlugin('emantomo'):
            importChoices.append('eman')
        if existsPlugin('dynamo'):
            importChoices.append('dynamo')
        return importChoices

    def _getDefaultChoice(self):
        return 0

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)
        self.OUTPUT_PREFIX = "outputCoordinates"


    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)
        form.addParam('boxSize', params.IntParam, label='Box size')

        form.addParam('importTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input tomograms',
                      help='Select the tomograms/tomogram for which you '
                            'want to import coordinates. The file names of the tomogram and '
                            'coordinate files must be the same.')

    def _insertAllSteps(self):
        self._insertFunctionStep('importCoordinatesStep',
                                 self.samplingRate.get())

    # --------------------------- STEPS functions -----------------------------

    def importCoordinatesStep(self, samplingRate):
        importTomograms = self.importTomograms.get()
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(importTomograms, suffix)
        coordsSet.setBoxSize(self.boxSize.get())
        coordsSet.setSamplingRate(samplingRate)
        coordsSet.setPrecedents(importTomograms)
        ci = self.getImportClass()
        for tomo in importTomograms.iterItems():
            tomoName = basename(os.path.splitext(tomo.getFileName())[0])
            for coordFile, fileId in self.iterFiles():
                fileName = basename(os.path.splitext(coordFile)[0])
                if tomo is not None and tomoName == fileName:
                    # Parse the coordinates in the given format for this micrograph
                    if self.getImportFrom() == self.IMPORT_FROM_EMAN or self.getImportFrom() == self.IMPORT_FROM_TXT:
                        def addCoordinate(coord, x, y, z):
                            coord.setVolume(tomo.clone())
                            coord.setPosition(x, y, z, const.BOTTOM_LEFT_CORNER)
                            coordsSet.append(coord)
                        ci.importCoordinates3D(coordFile, addCoordinate)
                    elif self.getImportFrom() == self.IMPORT_FROM_DYNAMO:
                        ci(coordFile, coordsSet, tomo.clone())

        args = {}
        args[self.OUTPUT_PREFIX] = coordsSet
        self._defineOutputs(**args)
        self._defineSourceRelation(self.importTomograms, coordsSet)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return self.hasAttribute('outputCoordinates3D')


    def _getCoordsMessage(self):
        return "Coordinates %s" % self.getObjectTag('outputCoordinates3D')

    def _summary(self):
        summary = []
        if self._hasOutput():
            summary.append("%s imported from:\n%s"
                           % (self._getCoordsMessage(), self.getPattern()))

            summary.append(u"Sampling rate: *%0.2f* (Å/px)" %
                           self.samplingRate.get())
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
        if importFrom == self.IMPORT_FROM_AUTO:
            importFrom = self.getFormat()
        return importFrom

    def getFormat(self):
        for coordFile, _ in self.iterFiles():
            if coordFile.endswith('.txt'):
                return self.IMPORT_FROM_TXT
            if coordFile.endswith('.json') and existsPlugin('emantomo'):
                return self.IMPORT_FROM_EMAN
            if coordFile.endswith('.tbl') and existsPlugin('dynamo'):
                return self.IMPORT_FROM_DYNAMO
        return -1

    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        importFrom = self.getImportFrom()

        if importFrom == self.IMPORT_FROM_EMAN:
            EmanImport = Domain.importFromPlugin('emantomo.convert', 'EmanTomoImport',
                                          errorMsg='Eman is needed to import .json or '
                                                   '.box files',
                                          doRaise=True)
            return EmanImport(self, None)

        elif importFrom == self.IMPORT_FROM_DYNAMO:
            readDynCoord = Domain.importFromPlugin("dynamo.convert.convert", "readDynCoord")
            return readDynCoord

        elif importFrom == self.IMPORT_FROM_TXT:
            return TomoImport(self)
        else:
            self.importFilePath = ''
            return None
