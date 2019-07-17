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

from os.path import abspath, basename

from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Transform
from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils.path import createAbsLink

from .protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from tomo.objects import SubTomogram


class ProtImportSubTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfSubTomograms'
    _label = 'import subtomograms'

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)

        form.addParam('importCoordinates', PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      allowsNull=True,
                      label='Input coordinates 3D',
                      help='Select the coordinates for which the '
                            'subtomograms were extracted.')

        ProtTomoImportAcquisition._defineParams(self, form)


    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['eman2']


    def _insertAllSteps(self):
        self._insertFunctionStep('importSubTomogramsStep',
                                 self.getPattern(),
                                 self.samplingRate.get())

    # --------------------------- STEPS functions -----------------------------

    def importSubTomogramsStep(
            self,
            pattern,
            samplingRate
    ):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)

        # Create a Volume template object
        subtomo = SubTomogram()
        subtomo.setSamplingRate(samplingRate)

        imgh = ImageHandler()

        subtomoSet = self._createSetOfSubTomograms()
        subtomoSet.setSamplingRate(samplingRate)

        if self.importCoordinates.get():
            self.coordDict = []
            for coord3DSet in self.importCoordinates.get().iterCoordinates():
                self.coordDict.append(coord3DSet.clone())
            subtomoSet.setCoordinates3D(self.importCoordinates)

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

            origin.setShifts(x/-2. * samplingRate,
                        y/-2. * samplingRate,
                        zDim/-2. * samplingRate)

            subtomo.setOrigin(origin)  # read origin from form

            newFileName = abspath(self._getVolumeFileName(fileName))

            if fileName.endswith(':mrc'):
                fileName = fileName[:-4]
            createAbsLink(fileName, newFileName)

            if n == 1:
                subtomo.cleanObjId()
                subtomo.setFileName(newFileName)
                if self.hasAttribute('coordDict'):
                    subtomo.setCoordinate3D(self.coordDict.pop(0))
                subtomo.setAcquisition(self._extractAcquisitionParameters(fileName))
                subtomoSet.append(subtomo)
            else:
                for index in range(1, n+1):
                    subtomo.cleanObjId()
                    subtomo.setLocation(index, newFileName)
                    if self.hasAttribute('coordDict'):
                        subtomo.setCoordinate3D(self.coordDict.pop(0))
                    subtomo.setAcquisition(self._extractAcquisitionParameters(fileName))
                    subtomoSet.append(subtomo)

        if subtomoSet.getSize() > 1:
            self._defineOutputs(outputSubTomograms=subtomoSet)
        else:
            self._defineOutputs(outputSubTomogram=subtomo)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return (self.hasAttribute('outputSubTomogram')
                or self.hasAttribute('outputSubTomograms'))

    def _getSubTomMessage(self):
        if self.hasAttribute('outputSubTomogram'):
            return "SubTomogram %s" % self.getObjectTag('outputSubTomogram')
        else:
            return "SubTomograms %s" % self.getObjectTag('outputSubTomograms')

    def _summary(self):
        summary = []
        if self._hasOutput():
            summary.append("%s imported from:\n%s"
                           % (self._getSubTomMessage(), self.getPattern()))

            if self.samplingRate.get():
                summary.append(u"Sampling rate: *%0.2f* (Å/px)" %
                               self.samplingRate.get())
            if self.hasAttribute('outputSubTomogram'):
                outputSubTomograms = [getattr(self, 'outputSubTomogram')]
            else:
                outputSubTomograms = getattr(self, 'outputSubTomograms')

            ProtTomoImportAcquisition._summary(self, summary, outputSubTomograms)

        return summary

    def _methods(self):
        methods = []
        if self._hasOutput():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getSubTomMessage(), self.samplingRate.get()))
        return methods

    def _getVolumeFileName(self, fileName, extension=None):
        if extension is not None:
            baseFileName="import_" + basename(fileName).split(".")[0] + ".%s"%extension
        else:
            baseFileName="import_" + basename(fileName).split(":")[0]

        return self._getExtraPath(baseFileName)


