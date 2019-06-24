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
import inspect

from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Transform
from pyworkflow.utils.path import createAbsLink
from pyworkflow.protocol.params import FloatParam, EnumParam, PathParam


from .protocol_base import ProtTomoImportFiles
from tomo.objects import Tomogram


class ProtImportTomograms(ProtTomoImportFiles):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfTomograms'
    _label = 'import tomograms'
    MANUAL_IMPORT = 0
    FROM_FILE_IMPORT = 1

    def __init__(self, **args):
        ProtTomoImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        ProtTomoImportFiles._defineParams(self, form)

        importAcquisitionChoices = self._getAcquisitionImportChoices()

        form.addSection(label='Acquisition Info')

        form.addParam('importAcquisitionFrom', EnumParam,
                      choices=importAcquisitionChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('acquisitionData', PathParam,
                      label="Acquisition parameters file",
                      help="File with the acquisition paramenters for every "
                            "tomogram to import. File must be in .txt format. The file must contain a row per file to be imported "
                           "and have the following parameters in order: \n"
                           "\n"
                           "'File_name AcquisitionAngleMin AcquisitionAngleMax Step AngleAxis1 AngleAxis2' \n"
                           "\n"
                           "An example would be: \n"
                           "tomo1.em -40 40 3 25 30 \n"
                           "tomo2.em -45 50 2 15 15 \n",
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
                      help='Enter the angle axis angle 1')

        form.addParam('angleAxis2', FloatParam,
                      allowsNull=True,
                      condition="importAcquisitionFrom == %d" % self.MANUAL_IMPORT,
                      label='Angle Axis 2',
                      help='Enter the angle axis 2')


    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        """
        return ['eman2']

    def _getAcquisitionImportChoices(self):
        return ['Manual', 'From file']

    def _insertAllSteps(self):
        self._insertFunctionStep('importTomogramsStep',
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

    def importTomogramsStep(
            self,
            pattern,
            samplingRate,
            importAcquisitionFrom,
            acquistionAngleMax,
            acquistionAngleMin,
            step,
            angleAxis1,
            angleAxis2,
            acquistionDataPath,
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

        if importAcquisitionFrom == self.FROM_FILE_IMPORT:
            acquisitionParams = self._parseAcquisitionData(acquistionDataPath)
        for fileName, fileId in self.iterFiles():
            onlyName = fileName.split('/')[-1]

            if importAcquisitionFrom == self.MANUAL_IMPORT:
                tomo.setAcquisitionAngleMax(acquistionAngleMax)
                tomo.setAcquisitionAngleMin(acquistionAngleMin)
                tomo.setStep(step)
                tomo.setAngleAxis1(angleAxis1)
                tomo.setAngleAxis2(angleAxis2)
            else:
                try:
                    tomo.setAcquisitionAngleMin(acquisitionParams[onlyName]['acquisitionAngleMin'])
                    tomo.setAcquisitionAngleMax(acquisitionParams[onlyName]['acquisitionAngleMax'])
                    tomo.setStep(acquisitionParams[onlyName]['step'])
                    tomo.setAngleAxis1(acquisitionParams[onlyName]['angleAxis1'])
                    tomo.setAngleAxis2(acquisitionParams[onlyName]['angleAxis2'])
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

            tomo.setOrigin(origin)  # read origin from form

            newFileName = abspath(self._getVolumeFileName(fileName))

            if fileName.endswith(':mrc'):
                fileName = fileName[:-4]
            createAbsLink(fileName, newFileName)
            if n == 1:
                tomo.cleanObjId()
                tomo.setFileName(newFileName)
                tomoSet.append(tomo)
            else:
                for index in range(1, n+1):
                    tomo.cleanObjId()
                    tomo.setLocation(index, newFileName)
                    tomoSet.append(tomo)

        if tomoSet.getSize() > 1:
            self._defineOutputs(outputTomograms=tomoSet)
        else:
            self._defineOutputs(outputTomogram=tomo)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return (self.hasAttribute('outputTomogram')
                or self.hasAttribute('outputTomograms'))

    def _getTomMessage(self):
        if self.hasAttribute('outputTomogram'):
            return "Tomogram %s" % self.getObjectTag('outputTomogram')
        else:
            return "Tomograms %s" % self.getObjectTag('outputTomograms')

    def _summary(self):
        summary = []
        if self._hasOutput():
            summary.append("%s imported from:\n%s"
                           % (self._getTomMessage(), self.getPattern()))
            summary.append(u"Sampling rate: *%0.2f* (Å/px)" %
                           self.samplingRate.get())

            if self.hasAttribute('outputTomogram'):
                outputTomograms = [getattr(self, 'outputTomogram')]
            else:
                outputTomograms = getattr(self, 'outputTomograms')
            i = 1
            for outputTomogram in outputTomograms:
                summary.append(u"File %d" % i)
                summary.append(u"Acquisition angle max: *%0.2f*" % outputTomogram.getAcquisitionAngleMax())
                summary.append(u"Acquisition angle min: *%0.2f*" % outputTomogram.getAcquisitionAngleMin())
                summary.append(u"Step: *%d*" % outputTomogram.getStep())
                summary.append(u"Angle axis 1: *%0.2f*" % outputTomogram.getAngleAxis1())
                summary.append(u"Angle axis 2: *%0.2f*" % outputTomogram.getAngleAxis2())
                i += 1

        return summary

    def _methods(self):
        methods = []
        if self._hasOutput():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getTomMessage(), self.samplingRate.get()),)
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
            except: raise Exception('Wrong acquisition data file format')
        return parameters
