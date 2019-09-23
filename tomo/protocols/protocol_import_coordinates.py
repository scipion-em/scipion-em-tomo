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

from os.path import basename

import pyworkflow.protocol.params as params
from pyworkflow.utils import importFromPlugin
import os

from .protocol_base import ProtTomoImportFiles

from tomo.objects import SetOfCoordinates3D


class ProtImportCoordinates3D(ProtTomoImportFiles):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfCoordinates3D'
    _label = 'import set of coordinates 3D'

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_XMIPP = 1
    IMPORT_FROM_RELION = 2
    IMPORT_FROM_EMAN = 3
    IMPORT_FROM_DOGPICKER = 4

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['eman']

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)


    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)
        form.addParam('boxSize', params.IntParam, label='Box size')

        form.addParam('importTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input tomograms',
                      help='Select the tomograms/tomogram for which you '
                            'want to import coordinates.')

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
        coordsSet.setVolumes(importTomograms)
        ci = self.getImportClass()
        for tomo in importTomograms.iterItems():
            tomoName = os.path.splitext(tomo.getFileName())[0]
            tomoName = basename(tomoName)
            for coordFile, fileId in self.iterFiles():
                fileName = os.path.splitext(coordFile)[0]
                fileName = "import_" + basename(fileName)
                if tomo is not None and tomoName == fileName:
                    def addCoordinate(coord):
                        coord.setVolume(tomo.clone())
                        coordsSet.append(coord)

                    # Parse the coordinates in the given format for this micrograph
                    ci.importCoordinates3D(coordFile, addCoordinate)

        self._defineOutputs(outputCoordinates=coordsSet)
        self._defineSourceRelation(self.importTomograms, coordsSet)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return self.hasAttribute('outputCoordinates')


    def _getCoordsMessage(self):
        return "Coordinates %s" % self.getObjectTag('outputCoordinates')

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
            baseFileName="import_" + basename(fileName).split(".")[0] + ".%s"%extension
        else:
            baseFileName="import_" + basename(fileName).split(":")[0]

        return self._getExtraPath(baseFileName)

    # ------------------ UTILS functions --------------------------------------
    def getImportFrom(self):
        importFrom = self.importFrom.get()
        if importFrom == self.IMPORT_FROM_AUTO:
            importFrom = self.IMPORT_FROM_EMAN
        return importFrom

    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        filesPath = self.filesPath.get()
        importFrom = self.getImportFrom()

        if importFrom == self.IMPORT_FROM_XMIPP:
            XmippImport = importFromPlugin('xmipp3.convert', 'XmippImport',
                                           'Xmipp is needed to import .xmd files',
                                           doRaise=True)
            return XmippImport(self, filesPath)

        elif importFrom == self.IMPORT_FROM_RELION:
            RelionImport = importFromPlugin('relion.convert', 'RelionImport',
                                            errorMsg='Relion is needed to import .star files',
                                            doRaise=True)
            return RelionImport(self, filesPath)

        elif importFrom == self.IMPORT_FROM_EMAN:
            EmanImport = importFromPlugin('eman2.convert', 'EmanImport',
                                          errorMsg='Eman is needed to import .json or '
                                                   '.box files',
                                          doRaise=True)
            return EmanImport(self, None)

        elif importFrom == self.IMPORT_FROM_DOGPICKER:
            DogpickerImport = importFromPlugin('appion.convert', 'DogpickerImport',
                                               errorMsg='appion plugin is needed to import '
                                                        'dogpicker files',
                                               doRaise=True)
            return DogpickerImport(self)
        else:
            self.importFilePath = ''
            return None



