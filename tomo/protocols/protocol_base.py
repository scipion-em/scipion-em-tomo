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
from pyworkflow.protocol.params import PointerParam
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





