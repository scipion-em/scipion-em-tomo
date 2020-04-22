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

import pyworkflow.protocol.params as params
from pyworkflow.plugin import Domain

from ..objects import SetOfCoordinates3D

from .protocol_base import ProtTomoImportFiles


class ProtImportCoordinates3D(ProtTomoImportFiles):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfCoordinates3D'
    _label = 'import set of coordinates 3D'

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_EMAN = 1
    IMPORT_FROM_DYNAMO = 2

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['auto', 'eman', 'dynamo']

    def _getDefaultChoice(self):
        return self.IMPORT_FROM_AUTO

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
                fileName = "import_" + basename(os.path.splitext(coordFile)[0])
                if tomo is not None and tomoName == fileName:
                    # Parse the coordinates in the given format for this micrograph
                    if self.getImportFrom() == self.IMPORT_FROM_EMAN:
                        def addCoordinate(coord):
                            coord.setVolume(tomo.clone())
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

    # ------------------ UTILS functions --------------------------------------
    def getImportFrom(self):
        importFrom = self.importFrom.get()
        if importFrom == self.IMPORT_FROM_AUTO:
            importFrom = self.getFormat()
        return importFrom

    def getFormat(self):
        for coordFile, _ in self.iterFiles():
            if coordFile.endswith('.json') or coordFile.endswith('.txt'):
                return self.IMPORT_FROM_EMAN
            if coordFile.endswith('.tbl'):
                return self.IMPORT_FROM_DYNAMO
        return -1

    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        importFrom = self.getImportFrom()

        if importFrom == self.IMPORT_FROM_EMAN:
            EmanImport = Domain.importFromPlugin('eman2.convert', 'EmanImport',
                                          errorMsg='Eman is needed to import .json or '
                                                   '.box files',
                                          doRaise=True)
            return EmanImport(self, None)

        elif importFrom == self.IMPORT_FROM_DYNAMO:
            readDynCoord = Domain.importFromPlugin("dynamo.convert.convert", "readDynCoord")
            return readDynCoord
        else:
            self.importFilePath = ''
            return None
