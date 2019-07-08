# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow as pw
import pyworkflow.em as pwem
from pyworkflow.protocol.params import PointerParam, EnumParam, PathParam, FloatParam
from pyworkflow.mapper.sqlite_db import SqliteDb
from pyworkflow.utils.properties import Message

import tomo.objects


class ProtTomoBase:
    def __createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pw.utils.cleanPath(setFn)

        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj

    def _createSetOfTiltSeriesM(self, suffix=''):
            return self.__createSet(tomo.objects.SetOfTiltSeriesM,
                                    'tiltseriesM%s.sqlite', suffix)

    def _createSetOfTiltSeries(self, suffix=''):
            self._ouputSuffix = ''
            return self.__createSet(tomo.objects.SetOfTiltSeries,
                                    'tiltseries%s.sqlite', suffix)

    def _createSetOfCoordinates3D(self, volSet, suffix=''):
        coord3DSet = self.__createSet(tomo.objects.SetOfCoordinates3D,
                                      'coordinates%s.sqlite', suffix,
                                      indexes=['_volId'])
        coord3DSet.setVolumes(volSet)
        return coord3DSet

    def _createSetOfTomograms(self, suffix=''):
        return self.__createSet(tomo.objects.SetOfTomograms,
                                'tomograms%s.sqlite', suffix)

    def _createSetOfSubTomograms(self, suffix=''):
        return self.__createSet(tomo.objects.SetOfSubTomograms,
                                'subtomograms%s.sqlite', suffix)

    def _createSetOfClassesSubTomograms(self, subTomograms, suffix=''):
        classes = self.__createSet(pwem.SetOfClasses3D,
                                   'classes3D%s.sqlite', suffix)
        classes.setImages(subTomograms)

        return classes

    def _getOutputSuffix(self, cls):
        """ Get the name to be used for a new output.
        For example: output3DCoordinates7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(cls):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1 # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else '' # empty if not output


class ProtTomoReconstruct(pwem.EMProtocol, ProtTomoBase):
    """ Base class for Tomogram reconstruction protocols. """
    pass

class ProtTomoPicking(pwem.ProtImport, ProtTomoBase):

    OUTPUT_PREFIX = 'output3DCoordinates'

    """ Base class for Tomogram boxing protocols. """
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputTomogram', PointerParam, label="Input Tomogram", important=True,
                      pointerClass='Tomogram',
                      help='Select the Tomogram to be used during picking.')


class ProtTomoImportFiles(pwem.ProtImportFiles, ProtTomoBase):

    def _defineParams(self, form):
        self._defineImportParams(form)

        self._defineAcquisitionParams(form)

    def _defineImportParams(self, form):
        """ Override to add options related to the different types
        of import that are allowed by each protocol.
        """
        importChoices = self._getImportChoices()

        form.addSection(label='Import')
        if len(importChoices) > 1: # not only from files
            form.addParam('importFrom', pwem.params.EnumParam,
                          choices=importChoices, default=self._getDefaultChoice(),
                          label='Import from',
                          help='Select the type of import.')
        else:
            form.addHidden('importFrom', pwem.params.EnumParam,
                          choices=importChoices, default=self.IMPORT_FROM_FILES,
                          label='Import from',
                          help='Select the type of import.')

        form.addParam('filesPath', pwem.params.PathParam,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/Particles/data/day??_micrographs/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/Particles/data/day*_micrographs/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/Particles/data/day#_micrographs/\n"
                           "'#' represents one digit that will be used as "
                           "micrograph ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', pwem.params.StringParam,
                      label='Pattern',
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addParam('copyFiles', pwem.params.BooleanParam, default=False,
                      expertLevel=pwem.params.LEVEL_ADVANCED,
                      label="Copy files?",
                      help="By default the files are not copied into the "
                           "project to avoid data duplication and to save "
                           "disk space. Instead of copying, symbolic links are "
                           "created pointing to original files. This approach "
                           "has the drawback that if the project is moved to "
                           "another computer, the links need to be restored.")

    def _defineAcquisitionParams(self, form):
        """ Override to add options related to acquisition info.
        """
        form.addParam('samplingRate', pwem.params.FloatParam,
                      label=Message.LABEL_SAMP_RATE)

    def _validate(self):
        pass

class ProtTomoImportAcquisition:

    MANUAL_IMPORT = 0
    FROM_FILE_IMPORT = 1

    def _defineParams(self, form):

        """ Override to add options related to acquisition info.
        """

        importAcquisitionChoices = ['Manual', 'From file']

        form.addSection(label='Acquisition Info')

        form.addParam('importAcquisitionFrom', EnumParam,
                      choices=importAcquisitionChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('acquisitionData', PathParam,
                      label="Acquisition parameters file",
                      help="File with the acquisition parameters for every "
                           "subtomogram to import. File must be in plain format. The file must contain a row per file to be imported "
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

    def _extractAcquisitionParameters(self, object, fileName):

        if self.importAcquisitionFrom.get() == self.MANUAL_IMPORT:
            object.setAcquisitionAngleMax(self.acquisitionAngleMax.get())
            object.setAcquisitionAngleMin(self.acquisitionAngleMin.get())
            object.setStep(self.step.get())
            object.setAngleAxis1(self.angleAxis1.get())
            object.setAngleAxis2(self.angleAxis2.get())
        else:
            try:
                onlyName = fileName.split('/')[-1]
                acquisitionParamsFile = self._parseAcquisitionData(self.acquisitionData.get())
                acquisitionParams = acquisitionParamsFile[onlyName]
                object.setAcquisitionAngleMin(acquisitionParams['acquisitionAngleMin'])
                object.setAcquisitionAngleMax(acquisitionParams['acquisitionAngleMax'])
                object.setStep(acquisitionParams['step'])
                object.setAngleAxis1(acquisitionParams['angleAxis1'])
                object.setAngleAxis2(acquisitionParams['angleAxis2'])
            except:
                raise Exception('Acquisition data file missing parameters')

    def _summary(self, summary, object):
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
            for key, outputSubTomogram in enumerate(object, 1):
                summary.append(u"File %d" % key)
                summary.append(u"Acquisition angle max: *%0.2f*" % outputSubTomogram.getAcquisitionAngleMax())
                summary.append(u"Acquisition angle min: *%0.2f*" % outputSubTomogram.getAcquisitionAngleMin())
                if outputSubTomogram.getStep():
                    summary.append(u"Step: *%d*" % outputSubTomogram.getStep())
                if outputSubTomogram.getAngleAxis1():
                    summary.append(u"Angle axis 1: *%0.2f*" % outputSubTomogram.getAngleAxis1())
                if outputSubTomogram.getAngleAxis2():
                    summary.append(u"Angle axis 2: *%0.2f*" % outputSubTomogram.getAngleAxis2())






