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
import pyworkflow.utils as pwutils
from pyworkflow.utils.path import createAbsLink


from .protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from ..objects import Tomogram


class ProtImportTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfTomograms'
    _label = 'import tomograms'

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)
        ProtTomoImportAcquisition._defineParams(self, form)

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        """
        return ['eman2']

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

            origin.setShifts(x / -2. * samplingRate,
                        y / -2. * samplingRate,
                        zDim / -2. * samplingRate)

            tomo.setOrigin(origin)  # read origin from form

            newFileName = self._getUniqueFileName(fileName)

            newFileName = abspath(self._getVolumeFileName(newFileName))

            if fileName.endswith(':mrc'):
                fileName = fileName[:-4]
            createAbsLink(fileName, newFileName)
            if n == 1:
                tomo.cleanObjId()
                tomo.setFileName(newFileName)
                tomo.setAcquisition(self._extractAcquisitionParameters(fileName))
                tomoSet.append(tomo)
            else:
                for index in range(1, n+1):
                    tomo.cleanObjId()
                    tomo.setLocation(index, newFileName)
                    tomo.setAcquisition(self._extractAcquisitionParameters(fileName))
                    tomoSet.append(tomo)

        self._defineOutputs(outputTomograms=tomoSet)

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

    def _getUniqueFileName(self, filename, filePaths=None):
        if filePaths is None:
            filePaths = [re.split(r'[$*#?]', self.getPattern())[0]]

        commPath = pwutils.commonPath(filePaths)
        return filename.replace(commPath + "/", "").replace("/", "_")

