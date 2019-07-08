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
from pyworkflow.protocol.params import PointerParam, FloatParam, PathParam, EnumParam
from pyworkflow.utils.path import createAbsLink

from .protocol_base import ProtTomoImportFiles
from tomo.objects import SubTomogram


class ProtImportSubTomograms(ProtTomoImportFiles):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfSubTomograms'
    _label = 'import subtomograms'
    MANUAL_IMPORT = 0
    FROM_FILE_IMPORT = 1

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)
        importAcquisitionChoices = self._getAcquisitionImportChoices()

        form.addParam('importCoordinates', PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      allowsNull=True,
                      label='Input coordinates 3D',
                      help='Select the coordinates for which the '
                            'subtomograms were extracted.')

        form.addSection(label='Acquisition Info')

        form.addParam('importAcquisitionFrom', EnumParam,
                      choices=importAcquisitionChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('acquisitionData', PathParam,
                      label="Acquisition parameters file",
                      help="File with the acquisition paramenters for every "
                           "subtomogram to import. File must be in .txt format. The file must contain a row per file to be imported "
                           "and have the following parameters in order: \n"
                           "\n"
                           "'File_name AcquisitionAngleMin AcquisitionAngleMax Step AngleAxis1 AngleAxis2' \n"
                           "\n"
                           "An example would be: \n"
                           "subtomo1.em -40 40 3 25 30 \n"
                           "subtomo2.em -45 50 2 15 15 \n",
                      condition="importAcquisitionFrom == %d" % self.FROM_FILE_IMPORT)

        form.addParam('acquisitionAngleMax', FloatParam,
                      allowsNull=True,
                      default=90,
                      label='Acquisition angle max',
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      help='Enter the positive limit of the acquisition angle')

        form.addParam('acquisitionAngleMin', FloatParam,
                      allowsNull=True,
                      default=-90,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Acquisition angle min',
                      help='Enter the negative limit of the acquisition angle')

        form.addParam('step', FloatParam,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Step',
                      help='Enter the step size for the import')

        form.addParam('angleAxis1', FloatParam,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Angle axis 1',
                      help='Enter the angle axis 1')

        form.addParam('angleAxis2', FloatParam,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Angle Axis 2',
                      help='Enter the angle axis 2')


    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['eman2']

    def _getAcquisitionImportChoices(self):
        return ['Manual', 'From file']

    def _insertAllSteps(self):
        self._insertFunctionStep('importSubTomogramsStep',
                                 self.getPattern(),
                                 self.samplingRate.get(),
                                 self.importAcquisitionFrom.get(),
                                 self.acquisitionAngleMax.get(),
                                 self.acquisitionAngleMin.get(),
                                 self.step.get(),
                                 self.angleAxis1.get(),
                                 self.angleAxis2.get(),
                                 self.acquisitionData.get())

    # --------------------------- STEPS functions -----------------------------

    def importSubTomogramsStep(
            self,
            pattern,
            samplingRate,
            importAcquisitionFrom,
            acquisitionAngleMax,
            acquisitionAngleMin,
            step,
            angleAxis1,
            angleAxis2,
            acquisitionDataPath,
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

        if importAcquisitionFrom == self.FROM_FILE_IMPORT:
            acquisitionParams = self._parseAcquisitionData(acquisitionDataPath)

        for fileName, fileId in self.iterFiles():
            onlyName = fileName.split('/')[-1]

            if importAcquisitionFrom == self.MANUAL_IMPORT:
                subtomo.setAcquisitionAngleMax(acquisitionAngleMax)
                subtomo.setAcquisitionAngleMin(acquisitionAngleMin)
                subtomo.setStep(step)
                subtomo.setAngleAxis1(angleAxis1)
                subtomo.setAngleAxis2(angleAxis2)
            else:
                try:
                    subtomo.setAcquisitionAngleMin(acquisitionParams[onlyName]['acquisitionAngleMin'])
                    subtomo.setAcquisitionAngleMax(acquisitionParams[onlyName]['acquisitionAngleMax'])
                    subtomo.setStep(acquisitionParams[onlyName]['step'])
                    subtomo.setAngleAxis1(acquisitionParams[onlyName]['angleAxis1'])
                    subtomo.setAngleAxis2(acquisitionParams[onlyName]['angleAxis2'])
                except:
                    raise Exception('Acqusition data file missing parameters')

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
                subtomoSet.append(subtomo)
            else:
                for index in range(1, n+1):
                    subtomo.cleanObjId()
                    subtomo.setLocation(index, newFileName)
                    subtomoSet.append(subtomo)
                subtomoSet.setCoordinates3D(self.importCoordinates)

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

            if self.importAcquisitionFrom.get() == self.MANUAL_IMPORT:
                summary.append(u"Acquisition angle max: *%0.2f*" % self.acquisitionAngleMax.get())
                summary.append(u"Acquisition angle min: *%0.2f*" % self.acquisitionAngleMin.get())
                if self.step.get():
                    summary.append(u"Step: *%d*" % self.step.get())
                if self.angleAxis1.get():
                    summary.append(u"Angle axis 1: *%0.2f*" % self.angleAxis1.get())
                if self.angleAxis2.get():
                    summary.append(u"Angle axis 2: *%0.2f*" % self.angleAxis2.get())
            else:
                for key, outputSubTomogram in enumerate(outputSubTomograms, 1):
                    summary.append(u"File %d" % key)
                    summary.append(u"Acquisition angle max: *%0.2f*" % outputSubTomogram.getAcquisitionAngleMax())
                    summary.append(u"Acquisition angle min: *%0.2f*" % outputSubTomogram.getAcquisitionAngleMin())
                    if outputSubTomogram.getStep():
                        summary.append(u"Step: *%d*" % outputSubTomogram.getStep())
                    if outputSubTomogram.getAngleAxis1():
                        summary.append(u"Angle axis 1: *%0.2f*" % outputSubTomogram.getAngleAxis1())
                    if outputSubTomogram.getAngleAxis2():
                        summary.append(u"Angle axis 2: *%0.2f*" % outputSubTomogram.getAngleAxis2())

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

    def _parseAcquisitionData(self, pathToData):
        params = open(pathToData, "r")
        lines = params.readlines()
        parameters = {}
        for l in lines:
            words = l.split()
            try:
                parameters.update({words[0]: {
                    'acquisitionAngleMin': float(words[1]),
                    'acquisitionAngleMax': float(words[2]),
                    'step': int(words[3]),
                    'angleAxis1': float(words[4]),
                    'angleAxis2': float(words[5])
                }})
            except:
                raise Exception('Wrong acquisition data file format')
        return parameters
