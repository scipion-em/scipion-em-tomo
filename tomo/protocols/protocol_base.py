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
from pyworkflow.protocol.params import (PointerParam, EnumParam, PathParam,
                                        FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.mapper.sqlite_db import SqliteDb
from pyworkflow.utils.properties import Message
from pwem.protocols import ProtImport, EMProtocol, ProtImportFiles

import tomo.objects


class ProtTomoBase:
    def _createSet(self, SetClass, template, suffix, **kwargs):
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
        return self._createSet(tomo.objects.SetOfTiltSeriesM,
                               'tiltseriesM%s.sqlite', suffix)

    def _createSetOfTiltSeries(self, suffix=''):
        self._ouputSuffix = ''
        return self._createSet(tomo.objects.SetOfTiltSeries,
                               'tiltseries%s.sqlite', suffix)

    def _createSetOfCoordinates3D(self, volSet, suffix=''):
        coord3DSet = self._createSet(tomo.objects.SetOfCoordinates3D,
                                     'coordinates%s.sqlite', suffix,
                                     indexes=['_volId'])
        coord3DSet.setPrecedents(volSet)
        return coord3DSet

    def _createSetOfTomograms(self, suffix=''):
        return self._createSet(tomo.objects.SetOfTomograms,
                               'tomograms%s.sqlite', suffix)

    def _createSetOfSubTomograms(self, suffix=''):
        return self._createSet(tomo.objects.SetOfSubTomograms,
                               'subtomograms%s.sqlite', suffix)

    def _createSetOfAverageSubTomograms(self, suffix=''):
        return self._createSet(tomo.objects.SetOfAverageSubTomograms,
                                  'avgSubtomograms%s.sqlite', suffix)

    def _createSetOfClassesSubTomograms(self, subTomograms, suffix=''):
        classes = self._createSet(tomo.objects.SetOfClassesSubTomograms,
                                  'subtomogramClasses%s.sqlite', suffix)
        classes.setImages(subTomograms)

        return classes

    def _createSetOfLandmarkModels(self, suffix=''):
        return self._createSet(tomo.objects.SetOfLandmarkModels, 'setOfLandmarks%s.sqlite', suffix)

    def _createSetOfMeshes(self, volSet, suffix=''):
        meshSet = self._createSet(tomo.objects.SetOfMeshes,
                                  'meshes%s.sqlite', suffix)
        meshSet.setPrecedents(volSet)
        return meshSet

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


class ProtTomoPicking(ProtImport, ProtTomoBase):

    OUTPUT_PREFIX = 'output3DCoordinates'

    """ Base class for Tomogram boxing protocols. """
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputTomograms', PointerParam, label="Input Tomograms", important=True,
                      pointerClass='SetOfTomograms',
                      help='Select the Tomogram to be used during picking.')

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append("Output 3D Coordinates not ready yet.")

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:* \n %s " % (key, output.getSummary()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

class ProtTomoImportFiles(ProtImportFiles, ProtTomoBase):

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
            form.addParam('importFrom', EnumParam,
                          choices=importChoices, default=self._getDefaultChoice(),
                          label='Import from',
                          help='Select the type of import.')
        else:
            form.addHidden('importFrom', EnumParam,
                          choices=importChoices, default=self.IMPORT_FROM_FILES,
                          label='Import from',
                          help='Select the type of import.')

        form.addParam('filesPath', PathParam,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/Tomograms/data/day??_tomograms/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/Tomograms/data/day*_tomograms/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/Tomograms/data/day#_tomograms/\n"
                           "'#' represents one digit that will be used as "
                           "tomogram ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', StringParam,
                      label='Pattern',
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addParam('copyFiles', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
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
        form.addParam('samplingRate', FloatParam,
                      label=Message.LABEL_SAMP_RATE)

    def _validate(self):
        pass


class ProtTomoSubtomogramAveraging(EMProtocol, ProtTomoBase):
    """ Base class for subtomogram averaging protocols. """
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
                      choices=importAcquisitionChoices,
                      default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.')

        form.addParam('acquisitionData', PathParam,
                      label="Acquisition parameters file",
                      help="File with the acquisition parameters for every "
                           "subtomogram to import. File must be in plain format."
                           " The file must contain a row per file to be imported "
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

    def _parseAcquisitionData(self):
        if self.importAcquisitionFrom.get() == self.MANUAL_IMPORT:
            self.acquisitionParameters = {
                        'angleMin': self.acquisitionAngleMin.get(),
                        'angleMax': self.acquisitionAngleMax.get(),
                        'step': self.step.get(),
                        'angleAxis1': self.angleAxis1.get(),
                        'angleAxis2': self.angleAxis2.get()
                    }
        else:
            params = open(self.acquisitionData.get(), "r")
            self.acquisitionParameters = {}
            for line in params.readlines():
                param = line.split()
                try:
                    self.acquisitionParameters.update({param[0]: {
                        'angleMin': float(param[1]),
                        'angleMax': float(param[2]),
                        'step': int(param[3]),
                        'angleAxis1': float(param[4]),
                        'angleAxis2': float(param[5])
                    }})
                except Exception as e:
                    print('Wrong acquisition data file format', e)

    def _extractAcquisitionParameters(self, fileName):
        if self.importAcquisitionFrom.get() == self.FROM_FILE_IMPORT:
            onlyName = fileName.split('/')[-1]
            acquisitionParams = self.acquisitionParameters[onlyName]
        else:
            acquisitionParams = self.acquisitionParameters

        return tomo.objects.TomoAcquisition(**acquisitionParams)

    def _summary(self, summary, setOfObject):
        for obj in setOfObject:
            if obj.hasAcquisition():
                summary.append(u"File: %s" % obj.getFileName())
                summary.append(u"Acquisition angle max: *%0.2f*" % obj.getAcquisition().getAngleMax())

                summary.append(u"Acquisition angle min: *%0.2f*" % obj.getAcquisition().getAngleMin())
                if obj.getAcquisition().getStep():
                    summary.append(u"Step: *%d*" % obj.getAcquisition().getStep())
                if obj.getAcquisition().getAngleAxis1():
                    summary.append(u"Angle axis 1: *%0.2f*" % obj.getAcquisition().getAngleAxis1())
                if obj.getAcquisition().getAngleAxis2():
                    summary.append(u"Angle axis 2: *%0.2f*" % obj.getAcquisition().getAngleAxis2())
