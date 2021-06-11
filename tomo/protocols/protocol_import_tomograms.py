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


import re
from os.path import abspath, basename

from pwem.emlib.image import ImageHandler
from pwem.objects import Transform

from pyworkflow import BETA
import pyworkflow.utils as pwutils
from pyworkflow.utils.path import createAbsLink
import pyworkflow.protocol.params as params


from .protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from ..objects import Tomogram
from ..utils import _getUniqueFileName


class ProtImportTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfTomograms'
    _label = 'import tomograms'
    _devStatus = BETA

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)
        ProtTomoImportAcquisition._defineParams(self, form)
        form.addSection('Origin Info')
        form.addParam('setOrigCoord', params.BooleanParam,
                      condition='importFrom == IMPORT_FROM_FILES',
                      label="Set origin of coordinates",
                      help="Option YES:\nA new volume will be created with "
                           "the "
                           "given ORIGIN of coordinates. This ORIGIN will be "
                           "set in the map file header.\nThe ORIGIN of "
                           "coordinates will be placed at the center of the "
                           "whole volume if you select n(x)/2, n(y)/2, "
                           "n(z)/2 as "
                           "x, y, z coordinates (n(x), n(y), n(z) are the "
                           "dimensions of the whole volume). However, "
                           "selecting "
                           "0, 0, 0 as x, y, z coordinates, the volume will be "
                           "placed at the upper right-hand corner.\n\n"
                           "Option NO:\nThe ORIGIN of coordinates will be "
                           "placed at the center of the whole volume ("
                           "coordinates n(x)/2, n(y)/2, n(z)/2 by default). "
                           "This "
                           "ORIGIN will NOT be set in the map file header.\n\n"
                           "WARNING: In case you want to process "
                           "the volume with programs requiring a specific "
                           "symmetry regarding the origin of coordinates, "
                           "for example the protocol extract unit "
                           "cell, check carefully that the coordinates of the "
                           "origin preserve the symmetry of the whole volume. "
                           "This is particularly relevant for loading "
                           "fragments/subunits of the whole volume.\n",
                      default=False)
        line = form.addLine('Offset',
                            help="A wizard will suggest you possible "
                                 "coordinates for the ORIGIN. In MRC volume "
                                 "files, the ORIGIN coordinates will be "
                                 "obtained from the file header.\n "
                                 "In case you prefer set your own ORIGIN "
                                 "coordinates, write them here. You have to "
                                 "provide the map center coordinates in "
                                 "Angstroms (pixels x sampling).\n",
                            condition='setOrigCoord')
        # line.addParam would produce a nicer looking form
        # but them the wizard icon is drawn outside the visible
        # window. Until this bug is fixed form is a better option
        form.addParam('x', params.FloatParam, condition='setOrigCoord',
                      label="x", help="offset along x axis (Angstroms)")
        form.addParam('y', params.FloatParam, condition='setOrigCoord',
                      label="y", help="offset along y axis (Angstroms)")
        form.addParam('z', params.FloatParam, condition='setOrigCoord',
                      label="z", help="offset along z axis (Angstroms)")

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        """
        return ['eman']

    def _insertAllSteps(self):
        self._insertFunctionStep('importTomogramsStep',
                                 self.getPattern(),
                                 self.samplingRate.get())

    # --------------------------- STEPS functions -----------------------------

    def importTomogramsStep(
            self,
            pattern,
            samplingRate
    ):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)

        # Create a Volume template object
        tomo = Tomogram()
        tomo.setSamplingRate(samplingRate)

        imgh = ImageHandler()

        tomoSet = self._createSetOfTomograms()
        tomoSet.setSamplingRate(samplingRate)

        self._parseAcquisitionData()
        for fileName, fileId in self.iterFiles():
            x, y, z, n = imgh.getDimensions(fileName)
            if fileName.endswith('.mrc') or fileName.endswith('.map'):
                fileName += ':mrc'
                if z == 1 and n != 1:
                    zDim = n
                    n = 1
                else:
                    zDim = z
            else:
                zDim = z

            origin = Transform()

            if self.setOrigCoord.get():
                origin.setShiftsTuple(self._getOrigCoord())
            else:
                origin.setShifts(x / -2. * samplingRate,
                                 y / -2. * samplingRate,
                                 zDim / -2. * samplingRate)

            tomo.setOrigin(origin)  # read origin from form

            newFileName = _getUniqueFileName(self.getPattern(), fileName.split(':')[0])

            # newFileName = abspath(self._getVolumeFileName(newFileName))

            if fileName.endswith(':mrc'):
                fileName = fileName[:-4]
            createAbsLink(fileName, abspath(self._getExtraPath(newFileName)))
            if n == 1:
                tomo.cleanObjId()
                tomo.setFileName(self._getExtraPath(newFileName))
                tomo.setAcquisition(self._extractAcquisitionParameters(fileName))
                tomoSet.append(tomo)
            else:
                for index in range(1, n+1):
                    tomo.cleanObjId()
                    tomo.setLocation(index, self._getExtraPath(newFileName))
                    tomo.setAcquisition(self._extractAcquisitionParameters(fileName))
                    tomoSet.append(tomo)

        self._defineOutputs(outputTomograms=tomoSet)

    # --------------------------- UTILS functions ------------------------------
    def _getOrigCoord(self):
        return -1. * self.x.get(), -1. * self.y.get(), -1. * self.z.get()

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return self.hasAttribute('outputTomograms')

    def _getTomMessage(self):
        return "Tomograms %s" % self.getObjectTag('outputTomograms')

    def _summary(self):

        try:
            summary = []
            if self._hasOutput():
                summary.append("%s imported from:\n%s"
                               % (self._getTomMessage(), self.getPattern()))

                if self.samplingRate.get():
                    summary.append(u"Sampling rate: *%0.2f* (Å/px)" % self.samplingRate.get())

                outputTomograms = getattr(self, 'outputTomograms')

                ProtTomoImportAcquisition._summary(self, summary, outputTomograms)

                x, y, z = self.outputTomograms.getFirstItem().getShiftsFromOrigin()
                summary.append(u"Tomograms Origin (x,y,z):\n"
                               u"    x: *%0.2f* (Å/px)\n"
                               u"    y: *%0.2f* (Å/px)\n"
                               u"    z: *%0.2f* (Å/px)" % (x, y, z))

        except Exception as e:
            print(e)

        return summary

    def _methods(self):
        methods = []
        if self._hasOutput():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getTomMessage(), self.samplingRate.get()),)
        return methods

    def _getVolumeFileName(self, fileName, extension=None):
        if extension is not None:
            baseFileName = "import_" + str(basename(fileName)).split(".")[0] + ".%s"%extension
        else:
            baseFileName = "import_" + str(basename(fileName)).split(":")[0]

        return self._getExtraPath(baseFileName)

    def _validate(self):
        errors = []
        try:
            next(self.iterFiles())
        except StopIteration:
            errors.append('No files matching the pattern %s were found.' % self.getPattern())
        return errors


