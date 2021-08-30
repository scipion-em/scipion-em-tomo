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

import os
import re
from glob import glob
import time
from datetime import timedelta, datetime
from collections import OrderedDict
from os.path import join, exists
from pathlib import PureWindowsPath
from statistics import mean

import numpy as np
from sqlite3 import OperationalError

import pyworkflow as pw
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.objects import Transform
from pyworkflow.object import Integer
from pyworkflow.utils import getParentFolder, removeBaseExt, yellowStr
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtImport

from tomo.convert import getAnglesFromHeader, getAnglesFromMdoc, getAnglesFromTlt
from tomo.objects import TomoAcquisition

from .protocol_base import ProtTomoBase


class ProtImportTsBase(ProtImport, ProtTomoBase):
    """ Base class for Tilt-Series and Tilt-SeriesMovies import protocols.
    """
    IMPORT_FROM_FILES = 0

    # How to handle the input files into the project
    IMPORT_COPY_FILES = 0
    IMPORT_LINK_ABS = 1
    IMPORT_LINK_REL = 2

    ANGLES_FROM_FILENAME = 'Filename'
    ANGLES_FROM_HEADER = 'Header'
    ANGLES_FROM_MDOC = 'Mdoc'
    ANGLES_FROM_TLT = 'Tlt'
    ANGLES_FROM_RANGE = 'Range'

    NOT_MDOC_GUI_COND = 'filesPattern is None or (filesPattern is not None and ".mdoc" not in filesPattern)'
    MDOC_DATA_SOURCE = False

    acquisitions = None
    sRates = None
    accumDoses = None
    incomingDose = None
    meanDosesPerFrame = None

    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.skippedMdocs = Integer()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Root directory of the tilt-series (or movies) files.")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      help="It determines if the tilt series are going to be imported using the mdoc file or the tilt "
                           "series files. To import from the mdoc files, the word '.mdoc' must appear in the pattern, "
                           "if not, a tilt series pattern is expected. In the first case, the angular and acquisition "
                           "data are directly read from the corresponding mdoc file, while in the second it is read "
                           "the base name of the matching files, according to the pattern introduced.\n\n"
                           "*IMPORTING WITH MDOC FILES*\n\n"
                           "For *tilt series movies*, a mdoc per tilt series movies is expected. "
                           "The corresponding movie file/s must be located in the same "
                           "path as the mdoc file. The tilt series id will be the base name of the mdoc files, "
                           "so the names of the mdoc files must be different, even if they're located in "
                           "different paths.\n\n"
                           "For *tilt series*, the only difference is that a stack .mrcs file is expected for each "
                           "mdoc, which means, per each tilt series desired to be imported.\n\n"
                           "*IMPORTING WITH A PATTERN OF THE TILT SERIES FILE NAMES*\n\n"
                           "The pattern can contain standard wildcards such as *, ?, etc.\n\n"
                           "It should also contains the following special tags:\n"
                           "   *{TS}*: tilt series identifier, which can be any UNIQUE part of the path. This must be "
                           "an alpha-numeric sequence (avoid symbols as -) that can not start with a number.\n"
                           "   *{TO}*: acquisition order, an integer value (important for dose).\n"
                           "   *{TA}*: tilt angle, a positive or negative float value.\n\n"
                           "Examples:\n\n"
                           "To import a set of image stacks (tilt-series or tilt-series movies) as: \n"
                           "TiltSeries_a_001_0.0.mrc\n"
                           "TiltSeries_a_002_3.0.mrc\n"
                           "TiltSeries_a_003_-3.0.mrc\n"
                           "...\n"
                           "TiltSeries_b_001_0.0.mrc\n"
                           "TiltSeries_b_002_3.0.mrc\n"
                           "TiltSeries_b_003_-3.0.mrc\n"
                           "...\n"
                           "The pattern TiltSeries_{TS}_{TO}_{TA}.mrc will identify:\n"
                           "{TS} as a, b, ...\n"
                           "{TO} as 001, 002, 003, ...\n"
                           "{TA} as 0.0, 3.0, -3.0, ...\n")
        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the path should not have",
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('mdocInfo', params.LabelParam,
                      condition='not (%s)' % self.NOT_MDOC_GUI_COND,
                      label='Acquisition values provided below will override the mdoc corresponding values',
                      important=True)
        self._defineAngleParam(form)
        form.addParam('importAction', params.EnumParam,
                      default=self.IMPORT_LINK_REL,
                      choices=['Copy files',
                               'Absolute symlink',
                               'Relative symlink'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Import action on files",
                      help="This parameters determine how the project will deal "
                           "with imported files. It can be: \n"
                           "*Copy files*: Input files will be copied into your "
                           "project. (this will duplicate the raw data)."
                           "*Absolute symlink*: Create symbolic links to the "
                           "absolute path of the files."
                           "*Relative symlink*: Create symbolic links as "
                           "relative path from the protocol run folder. ")

        self._defineAcquisitionParams(form)

        form.addSection('Streaming')

        form.addParam('dataStreaming', params.BooleanParam, default=False,
                      label="Process data in streaming?",
                      help="Select this option if you want import data as it is "
                           "generated and process on the fly by next protocols. "
                           "In this case the protocol will keep running to check "
                           "new files and will update the output Set, which can "
                           "be used right away by next steps.")

        form.addParam('timeout', params.IntParam, default=43200,
                      condition='dataStreaming',
                      label="Timeout (secs)",
                      help="Interval of time (in seconds) after which, if no new file "
                           "is detected, the protocol will end. When finished, "
                           "the output Set will be closed and no more data will be "
                           "added to it. \n"
                           "Note 1:  The default value is  high (12 hours) to avoid "
                           "the protocol finishes during the acq of the "
                           "microscope. You can also stop it from right click and press "
                           "STOP_STREAMING.\n"
                           "Note 2: If you're using individual frames when importing "
                           "movies, the timeout won't be refreshed until a whole "
                           "movie is stacked.")

        form.addParam('fileTimeout', params.IntParam, default=30,
                      condition='dataStreaming',
                      label="File timeout (secs)",
                      help="Interval of time (in seconds) after which, if a file has "
                           "not changed, we consider it as a new file. \n")

        self._defineBlacklistParams(form)

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        pass

    def _defineAcquisitionParams(self, form):
        """ Define acq parameters, it can be overridden
        by subclasses to change what parameters to include.
        """
        group = form.addGroup('Acquisition info - override mdoc values if provided')
        group.addParam('voltage', params.FloatParam,
                       label=Message.LABEL_VOLTAGE,
                       allowsNull=True,
                       help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', params.FloatParam, default=2.7,
                       label=Message.LABEL_SPH_ABERRATION,
                       expertLevel=params.LEVEL_ADVANCED,
                       help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
                       label=Message.LABEL_AMPLITUDE,
                       expertLevel=params.LEVEL_ADVANCED,
                       help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', params.IntParam,
                       label=Message.LABEL_MAGNI_RATE,
                       allowsNull=True,
                       help=Message.TEXT_MAGNI_RATE)
        group.addParam('samplingRate', params.FloatParam,
                       label=Message.LABEL_SAMP_RATE,
                       allowsNull=True,
                       help=Message.TEXT_SAMP_RATE)
        group.addParam('tiltAxisAngle', params.FloatParam,
                       label='Tilt axis angle (deg.)',
                       allowsNull=True)
        line = group.addLine('Dose (electrons/sq.Å)',
                             help="Initial accumulated dose (usually 0) and "
                                  "dose per tilt image (electrons/sq.Å). ")
        line.addParam('doseInitial', params.FloatParam,
                      default=0,
                      label='Initial dose')
        line.addParam('dosePerFrame', params.FloatParam,
                      allowsNull=True,
                      label='Dose per tilt image')

        return group

    def _defineBlacklistParams(self, form):
        """ Override to add options related to blacklist info.
        """
        pass

    # -------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('importStep')

    # -------------------------- STEPS functions -------------------------------
    def importStep(self):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        accumDoseList = []
        incomingDoseList = []
        samplingRate = self.samplingRate.get()
        counter = 0

        if not self.MDOC_DATA_SOURCE:
            self.info("Using glob pattern: '%s'" % self._globPattern)
            self.info("Using regex pattern: '%s'" % self._regexPattern)

        outputSet = getattr(self, self._outputName, None)

        if outputSet is None:
            createSetFunc = getattr(self, self._createOutputName)
            outputSet = createSetFunc()
        elif outputSet.getSize() > 0:
            outputSet.loadAllProperties()

            for ts in outputSet:
                self._existingTs.add(ts.getTsId())

        self._fillAcquisitionInfo(outputSet)
        tsClass = outputSet.ITEM_TYPE
        tiClass = tsClass.ITEM_TYPE

        finished = False
        lastDetectedChange = datetime.now()

        # Ignore the timeout variables if we are not really in streaming mode
        if self.dataStreaming:
            timeout = timedelta(seconds=self.timeout.get())
            fileTimeout = timedelta(seconds=self.fileTimeout.get())
        else:
            timeout = timedelta(seconds=5)
            fileTimeout = timedelta(seconds=5)

        while not finished:
            time.sleep(3)  # wait 3 seconds before check for new files
            someNew = False  # Check if some new TS has been found
            someAdded = False  # Check if some new were added
            # incompleteTs = False  # Check if there are incomplete TS

            matchingFiles = self.getMatchingFiles(fileTimeOut=fileTimeout)

            if self._existingTs:
                outputSet.enableAppend()

            for ts, tiltSeriesList in matchingFiles.items():
                someNew = True
                tsObj = tsClass(tsId=ts)
                # Form value has higher priority than the mdoc values
                samplingRate = float(samplingRate if samplingRate else self.sRates[ts])

                origin = Transform()
                tsObj.setOrigin(origin)
                tsObj.setAnglesCount(len(tiltSeriesList))

                # we need this to set mapper before adding any item
                outputSet.append(tsObj)

                if self.MDOC_DATA_SOURCE:
                    accumDoseList = self.accumDoses[ts]
                    # tsObj.getAcquisition().setAccumDose(accumDoseList[-1])
                    incomingDoseList = self.incomingDose[ts]
                    counter = 0

                tiltSeriesObjList = []
                toList = [ti[1] for ti in tiltSeriesList]
                # Add tilt images to the tiltSeries
                for f, to, ta in tiltSeriesList:
                    try:
                        # Link/move to extra
                        if type(f) == tuple:
                            imageFile = f[1]
                        else:
                            imageFile = f

                        finalDestination = self._getExtraPath(os.path.basename(imageFile))
                        self.copyOrLink(imageFile, finalDestination)

                        if type(f) == tuple:
                            f = f[0], finalDestination
                        else:
                            f = finalDestination

                        ti = tiClass(location=f,
                                     acquisitionOrder=to,
                                     tiltAngle=ta)

                        ti.setAcquisition(tsObj.getAcquisition().clone())
                        if self.MDOC_DATA_SOURCE:
                            dosePerFrame = incomingDoseList[counter]
                            accumDose = accumDoseList[counter]
                        else:
                            dosePerFrame = self.dosePerFrame.get()
                            accumDose = self.dosePerFrame.get() * int(to if min(toList) == 1 else (int(to) + 1))

                        # Incoming dose in current ti
                        ti.getAcquisition().setDosePerFrame(dosePerFrame)
                        # Accumulated dose in current ti
                        ti.getAcquisition().setAccumDose(accumDose)
                        tiltSeriesObjList.append(ti)
                        counter += 1

                    except OperationalError as e:

                        raise Exception("%s is an invalid for the {TS} field, it must be an alpha-numeric sequence "
                                        "(avoid symbols as -) that can not start with a number." % ts)

                # Sort tilt image metadata if importing tilt series
                if not self._isImportingTsMovies():
                    tiltSeriesObjList.sort(key=lambda x: x.getTiltAngle(), reverse=False)

                for ti in tiltSeriesObjList:
                    tsObj.append(ti)

                tsObjFirstItem = tsObj.getFirstItem()
                origin.setShifts(-tsObjFirstItem.getXDim() / 2 * samplingRate,
                                 -tsObjFirstItem.getYDim() / 2 * samplingRate,
                                 0)

                if self.MDOC_DATA_SOURCE:
                    tsObj.getAcquisition().setAccumDose(accumDoseList[-1])
                    # Tilt series object dose per frame has been updated each time the tilt image dose per frame has
                    # been updated before, so the mean value is used to be the reference in the acquisition of the
                    # whole tilt series movie
                    tsObj.getAcquisition().setDosePerFrame(mean(incomingDoseList))
                else:
                    tsObj.getAcquisition().setDosePerFrame(self.dosePerFrame.get())
                    tsObj.getAcquisition().setAccumDose(self.dosePerFrame.get() * len(tiltSeriesList))

                outputSet.update(tsObj)  # update items and size info
                self._existingTs.add(ts)
                someAdded = True

            if someAdded:
                self.debug('Updating output...')
                outputSet.updateDim()
                self._updateOutputSet(self._outputName, outputSet,
                                      state=outputSet.STREAM_OPEN)
                self.debug('Update Done.')

            self.debug('Checking if finished...someNew: %s' % someNew)

            now = datetime.now()

            if not someNew:
                # If there are no new detected files, we should check the
                # inactivity time elapsed (from last event to now) and
                # if it is greater than the defined timeout, we conclude
                # the import and close the output set
                # Another option is to check if the protocol have some
                # special stop condition, this can be used to manually stop
                # some protocols such as import movies
                finished = (now - lastDetectedChange > timeout
                            or self.streamingHasFinished())
                self.debug("Checking if finished:")
                self.debug("   Now - Last Change: %s"
                           % pwutils.prettyDelta(now - lastDetectedChange))
            else:
                # If we have detected some files, we should update
                # the timestamp of the last event
                lastDetectedChange = now
                finished = not self.isInStreaming()

            self.debug("Finished: %s" % finished)

        # Close the output set
        self._updateOutputSet(self._outputName, outputSet,
                              state=outputSet.STREAM_CLOSED)

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []

        if self.getOutputsSize():
            for key, output in self.iterOutputAttributes():
                summary.append("Imported tilt-series %s from: %s" % (
                    "movies" if output.getLastName() == "outputTiltSeriesM" else "",
                    self.filesPath.get()))
                summary.append("Using pattern: %s" % self.filesPattern.get())
                summary.append(u"Sampling rate: *%0.2f* (Å/px)" % output.getSamplingRate())
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        if self.skippedMdocs.get():
            summary.append('*%i* mdoc files were skipped --> check the Output Log tab for more details.'
                           % self.skippedMdocs.get())
        return summary

    def _validate(self):
        self._initialize()
        try:
            matching = self.getMatchingFiles(isValidation=True)

        except Exception as e:
            errorStr = str(e)
            return [errorStr]

        if not matching:
            return ["There are no files matching the pattern %s"
                    % self._globPattern]

        self._firstMatch = list(matching.items())[0]
        self._tiltAngleList = self._getSortedAngles(self._firstMatch[1])

        errMsg = []
        errorMsgAngles = self._validateAngles()
        if errorMsgAngles:
            errMsg.append(errorMsgAngles)
        # In the mdoc case, voltage, magnification and sampling rate are optional inputs. In the user introduces
        # one of these values, it will be considered more prior than the corresponding value read from the mdoc file
        if not self.MDOC_DATA_SOURCE:
            if not self.voltage.get():
                errMsg.append('Voltage should be a float')
            if not self.magnification.get():
                errMsg.append('Magnification should be an integer')
            if not self.samplingRate.get():
                errMsg.append('Sampling rate should be a float')
            if not self.dosePerFrame.get():
                errMsg.append('Dose per frame should be a float')
            if not self.tiltAxisAngle.get():
                errMsg.append('Tilt axis angle should be a float')

        return errMsg

    def _validateAngles(self):
        """ Function to be implemented in subclass to validate
        the angles range.
        """
        return []

    # -------------------------- BASE methods to be overridden -----------------
    # def _getImportChoices(self):
    #     """ Return a list of possible choices
    #     from which the import can be done.
    #     (usually packages form such as: xmipp3, eman2, relion...etc.
    #     """
    #     return ['files']

    # -------------------------- UTILS functions -------------------------------
    def _initialize(self):
        """ Initialize some internal variables such as:
        - patterns: Expand the pattern using environ vars or username
            and also replacing special character # by digit matching.
        - outputs: Output variable names
        """
        self.MDOC_DATA_SOURCE = 'mdoc' in self.filesPattern.get()
        if not self.MDOC_DATA_SOURCE:
            path = self.filesPath.get('').strip()
            pattern = self.filesPattern.get('').strip()
            self._pattern = os.path.join(path, pattern) if pattern else path

            def _replace(p, ts, to, ta):
                p = p.replace('{TS}', ts)
                p = p.replace('{TO}', to)
                p = p.replace('{TA}', ta)
                return p

            self._regexPattern = _replace(self._pattern.replace('*', '(.*)'),
                                          '(?P<TS>.*)',
                                          '(?P<TO>\d+)',
                                          '(?P<TA>[+-]?\d+(\.\d+)?)')
            self._regex = re.compile(self._regexPattern)
            self._globPattern = _replace(self._pattern, '*', '*', '*')

        # Set output names depending on the import type (either movies or images)
        self._outputName = 'outputTiltSeries'
        self._createOutputName = '_createSetOfTiltSeries'

        # Keep track of which existing tilt-series has already been found
        self._existingTs = set()

    def _anglesInPattern(self):
        """ This function should be called after a call to _initialize"""
        return '{TA}' in self._pattern and '{TO}' in self._pattern

    def _getMatchingFilesFromMdoc(self, isValidation):
        """If the list of files provided by the user is a list of mdoc files, then the tilt series movies
        are built from them, following the considerations listed below:
            - For each mdoc file, it and the corresponding movie files must be in the same directory.
            - The tilt series id will be the base name of the mdoc file, by default, so the mdocs must have different
              base name. If another name is desired, the user can introduce the name structure (see advanced parameter)
            """
        fpath = self.filesPath.get()
        mdocList = glob(join(fpath, self.filesPattern.get()))
        hasDoseList = []
        if not mdocList:
            raise Exception('No mdoc files were found in the introduced path:\n%s' % fpath)

        matchingFiles = OrderedDict()
        self.acquisitions = OrderedDict()
        self.sRates = OrderedDict()
        self.accumDoses = OrderedDict()
        self.incomingDose = OrderedDict()
        warningHeadMsg = yellowStr('The following mdoc files were skipped. See details below:\n\n')
        warningDetailedMsg = []
        skippedMdocs = 0

        for mdoc in mdocList:
            # Note: voltage, magnification and sampling rate values are the ones introduced by the user in the
            # protocol's form. Otherwise, the corresponding values considered will be the ones read from the mdoc.
            # This is because because you can't trust mdoc (often dose is not calibrated in serialem, so you get 0;
            # pixel size might be binned as mdoc comes from a binned record not movie and  there are no Cs and amp
            # contrast fields in mdoc)
            mdocObj = MDoc(mdoc, voltage=self.voltage.get() if self.voltage.get() else None,
                           magnification=self.magnification.get() if self.magnification.get() else None,
                           samplingRate=self.samplingRate.get() if self.samplingRate.get() else None,
                           doseProvidedByUser=self.dosePerFrame.get() if self.dosePerFrame.get() else None,
                           tiltAngleProvidedByUser=self.tiltAxisAngle.get() if self.tiltAxisAngle.get() else None)
            validationError = mdocObj.read(isImportingTsMovies=self._isImportingTsMovies())
            hasDoseList.append(mdocObj.mdocHasDose)
            if validationError:
                warningHeadMsg += yellowStr('\t- %s\n' % mdoc)
                warningDetailedMsg.append(validationError)
                skippedMdocs += 1
                # validationErrors.append(validationError)
                # Continue parsing the remaining mdoc files to provide a fully detailed error message
                continue
            acquisition = self._genTsAcquisitionFromMdoc(mdocObj)
            tsId = mdocObj.getTsId()
            fileOrderAngleList = []
            accumulatedDoseList = []
            incomingDoseList = []
            for tiltMetadata in mdocObj.getTiltsMetadata():
                fileOrderAngleList.append((
                    tiltMetadata.getAngleMovieFile(),                    # Filename
                    '{:03d}'.format(tiltMetadata.getAcqOrder()),         # Acquisition order
                    tiltMetadata.getTiltAngle(),                         # Tilt angle
                ))
                accumulatedDoseList.append(tiltMetadata.getAccumDose())
                incomingDoseList.append(tiltMetadata.getIncomingDose())

            # self._getTsIdFromMdocData(fileList)
            matchingFiles[tsId] = fileOrderAngleList
            self.acquisitions[tsId] = acquisition
            self.sRates[tsId] = mdocObj.getSamplingRate()
            self.accumDoses[tsId] = accumulatedDoseList
            self.incomingDose[tsId] = incomingDoseList

        if isValidation:
            if matchingFiles:
                if warningDetailedMsg:
                    self.skippedMdocs.set(skippedMdocs)
                    self._store(self.skippedMdocs)
                    print(warningHeadMsg + ' '.join(warningDetailedMsg))
                return matchingFiles
            else:
                # If the only info missing is the dose related data, it's suggested to introduce it manually
                if not any(hasDoseList):
                    raise Exception('*The dose was not possible to be obtained from any of the provided mdoc files.\n'
                                    'Please check the data of your mdoc files or introduce a dose value in the '
                                    'protocol form.\n\n'
                                    'Dose related mdoc labels are:\n\n'
                                    '- ExposureDose or\n'
                                    '- FrameDosesAndNumber or\n'
                                    '- DoseRate and ExposureTime or\n'
                                    '- MinMaxMean and CountsPerElectron'
                                    )
                else:
                    raise Exception('*All the mdoc files introduced present validation errors.*\n\n%s' %
                                    (warningHeadMsg + ' '.join(warningDetailedMsg)))
        else:
            return matchingFiles

    def _isImportingTsMovies(self):
        return True if type(self) is ProtImportTsMovies else False

    def _genTsAcquisitionFromMdoc(self, mdocObj):
        acq = TomoAcquisition(voltage=mdocObj.getVoltage(),
                              sphericalAberration=self.sphericalAberration.get(),
                              amplitudeContrast=self.amplitudeContrast.get(),
                              magnification=mdocObj.getMagnification(),
                              tiltAxisAngle=mdocObj.getTiltAxisAngle()
                              )

        if hasattr(self, 'doseInitial'):  # This field is only present in the form for TsM import
            acq.setDoseInitial(self.doseInitial.get())

        return acq

    def _excludeByWords(self, files):
        exclusionWords = self.exclusionWords.get()

        if exclusionWords is None:
            return files

        exclusionWordList = exclusionWords.split()

        allowedFiles = []

        for file in files:
            if any(bannedWord in file for bannedWord in exclusionWordList):
                print("%s excluded. Contains any of %s" % (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles

    def getMatchingFiles(self, fileTimeOut=None, isValidation=False):
        """ Return an ordered dict with TiltSeries found in the files as key
        and a list of all tilt images of that series as value.
        """
        if self.MDOC_DATA_SOURCE:
            return self._getMatchingFilesFromMdoc(isValidation=isValidation)
        else:
            return self._getMatchingFilesFromRegExPattern(fileTimeOut=fileTimeOut)

    def _getMatchingFilesFromRegExPattern(self, fileTimeOut):
        filePaths = glob(self._globPattern)

        filePaths = self._excludeByWords(filePaths)

        filePaths.sort(key=lambda fn: os.path.getmtime(fn))
        fileTimeOut = fileTimeOut or timedelta(seconds=5)

        matchingFiles = OrderedDict()

        def _getTsId(match):
            """ Retrieve the TiltSerie ID from the matching object.
            We need to have tsId that starts with character, so
            let's add a prefix if it is not the case
            """
            tsId = match.group('TS')
            return tsId

        def _addOne(fileList, file, match):
            """ Add one file matching to the list. """
            fileList.append((file, int(match.group('TO')), float(match.group('TA'))))

        def _addMany(fileList, file, match):
            """ Add many 'files' (when angles in header or mdoc) to the list. """
            anglesFrom = self.getEnumText('anglesFrom')

            if anglesFrom == self.ANGLES_FROM_HEADER:
                angles = getAnglesFromHeader(file)
            elif anglesFrom == self.ANGLES_FROM_MDOC:
                mdocFn = os.path.splitext(file)[0] + '.mdoc'
                if not os.path.exists(mdocFn):
                    raise Exception("Missing angles file: %s" % mdocFn)
                angles = getAnglesFromMdoc(mdocFn)
            elif anglesFrom == self.ANGLES_FROM_TLT:
                tltFn = os.path.splitext(file)[0] + '.tlt'
                if not os.path.exists(tltFn):
                    raise Exception("Missing angles file: %s" % tltFn)
                angles = getAnglesFromTlt(tltFn)
            elif anglesFrom == self.ANGLES_FROM_RANGE:
                angles = self._getTiltAngleRange()
            else:
                raise Exception('Invalid angles option: %s' % anglesFrom)

            for i, a in enumerate(angles):
                fileList.append(((i + 1, file), i + 1, a))

        addFunc = _addOne if self._anglesInPattern() else _addMany

        # Handle special case of just one TiltSeries, to avoid
        # the user the need to specify {TS}
        if len(filePaths) == 1 and not self.isInStreaming():
            f = filePaths[0]
            ts = pwutils.removeBaseExt(f)  # Base name without extension
            matchingFiles[ts] = []
            _addMany(matchingFiles[ts], f, None)
        else:
            for f in filePaths:
                if self.fileModified(f, fileTimeOut):
                    continue
                m = self._regex.match(f)
                if m is not None:
                    ts = _getTsId(m)
                    # Only report files of new tilt-series
                    if ts not in self._existingTs:
                        if ts not in matchingFiles:
                            matchingFiles[ts] = []
                        addFunc(matchingFiles[ts], f, m)

        return matchingFiles

    def _getCopyOrLink(self):
        """ Returns a function to copy or link files based on user selected option"""

        if self.importAction.get() == self.IMPORT_COPY_FILES:
            return pw.utils.copyFile
        elif self.importAction.get() == self.IMPORT_LINK_REL:
            return pw.utils.createLink
        else:
            return pw.utils.createAbsLink

    def copyOrLink(self, source, destination):
        """ Calls the copy or link method chosen by the user in importAction option"""
        func = self._getCopyOrLink()
        func(source, destination)

    def fileModified(self, fileName, fileTimeout):
        """ Check if the fileName modification time is less
        than a given timeout.
        Params:
            fileName: input filename that will be checked.
            fileTimeout: timeout
        """
        self.debug('Checking file: %s' % fileName)
        mTime = datetime.fromtimestamp(os.path.getmtime(fileName))
        delta = datetime.now() - mTime
        self.debug('   Modification time: %s' % pw.utils.prettyTime(mTime))
        self.debug('   Delta: %s' % pw.utils.prettyDelta(delta))

        return delta < fileTimeout

    def isBlacklisted(self, fileName):
        """ Overwrite in subclasses """
        return False

    @classmethod
    def worksInStreaming(cls):
        # Import protocols always work in streaming
        return True

    def _fillAcquisitionInfo(self, inputTs):
        # TODO Acquisition is historically expected as one per set of tilt series movies,
        # so the first one is hte one used, at least for now
        if self.MDOC_DATA_SOURCE:
            firstTsId = list(self.acquisitions.keys())[0]
            acq = self.acquisitions[firstTsId]
            sRate = self.sRates[firstTsId]
            inputTs.setAcquisition(acq)
            inputTs.setSamplingRate(sRate)
        else:
            inputTs.setSamplingRate(self.samplingRate.get())
            acq = inputTs.getAcquisition()
            acq.setVoltage(self.voltage.get())
            acq.setSphericalAberration(self.sphericalAberration.get())
            acq.setAmplitudeContrast(self.amplitudeContrast.get())
            acq.setMagnification(self.magnification.get())

    def _getTiltAngleRange(self):
        """ Return the list with all expected tilt angles. """
        offset = 1
        if self.minAngle.get() > self.maxAngle.get():
            offset *= -1
        return np.arange(self.minAngle.get(),
                         self.maxAngle.get() + offset,  # also include last angle
                         self.stepAngle.get())

    @staticmethod
    def _getSortedAngles(tiltSeriesList):
        """ Return the sorted angles from a given tiltSeriesList. """
        return sorted(item[2] for item in tiltSeriesList)

    def _sameTiltAngleRange(self, tiltAngleRange, tiltSeriesList):
        # allow some tolerance when comparing tilt-angles
        return np.allclose(tiltAngleRange,
                           self._getSortedAngles(tiltSeriesList),
                           atol=0.1)

    # --------------- Streaming special functions -----------------------
    def _getStopStreamingFilename(self):
        return self._getExtraPath("STOP_STREAMING.TXT")

    def getActions(self):
        """ This method will allow that the 'Stop import' action to appears
        in the GUI when the user right-click in the protocol import box.
        It will allow a user to manually stop the streaming.
        """
        # Only allow to stop if running and in streaming mode
        if self.dataStreaming and self.isRunning():
            return [('STOP STREAMING', self.stopImport)]
        else:
            return []

    def stopImport(self):
        """ Since the actual protocol that is running is in a different
        process that the one that this method will be invoked from the GUI,
        we will use a simple mechanism to place an special file to stop
        the streaming.
        """
        # Just place an special file into the run folder
        f = open(self._getStopStreamingFilename(), 'w')
        f.close()

    def streamingHasFinished(self):
        return (not self.isInStreaming() or
                os.path.exists(self._getStopStreamingFilename()))

    def isInStreaming(self):
        return self.dataStreaming.get()


class ProtImportTs(ProtImportTsBase):
    """Protocol to import tilt series."""
    _label = 'import tilt-series'
    _devStatus = pw.BETA

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        group = form.addGroup('Tilt info',
                              condition=self.NOT_MDOC_GUI_COND)

        group.addParam('anglesFrom', params.EnumParam,
                       default=0,
                       choices=[self.ANGLES_FROM_RANGE,
                                self.ANGLES_FROM_HEADER,
                                # self.ANGLES_FROM_MDOC,
                                self.ANGLES_FROM_TLT],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Import angles from',
                       help="Choose how the tilt angles will be inferred. "
                            "It can be taken from a range: Min, Max, Step "
                            "or from the image header, or from an"
                            "mdoc or tlt file (should have the SAME file name "
                            "but with the .mdoc or .tlt extension at the end).")

        line = group.addLine('Tilt angular range',
                             condition='anglesFrom==0',  # ANGLES_FROM_RANGE
                             help="Specify the tilting angular range. "
                                  "Depending on the collection schema, the "
                                  "order of the acquisition does not need to "
                                  "be the same order of the angular range. ")
        line.addParam('minAngle', params.FloatParam, default=-60, label='min')
        line.addParam('maxAngle', params.FloatParam, default=60, label='max')
        line.addParam('stepAngle', params.FloatParam, default=3, label='step')

    def _validateAngles(self):
        if not self.MDOC_DATA_SOURCE:
            ts, tiltSeriesList = self._firstMatch
            i, fileName = tiltSeriesList[0][0]
            x, y, z, n = ImageHandler().getDimensions(fileName)
            nImages = max(z, n)  # Just handle ambiguity with mrc format
            nAngles = len(self._tiltAngleList)
            if nAngles != nImages:
                return ['Tilt-series %s stack has different number of images (%d) '
                        'than the expected number of tilt angles (%d). '
                        % (fileName, nImages, nAngles)]
        else:
            return []


class ProtImportTsMovies(ProtImportTsBase):
    """Protocol to import tilt series movies."""
    _label = 'import tilt-series movies'
    _devStatus = pw.BETA

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        group = form.addGroup('Tilt info',
                              condition=self.NOT_MDOC_GUI_COND)

        group.addParam('anglesFrom', params.EnumParam,
                       default=0,
                       choices=[self.ANGLES_FROM_FILENAME],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Import angles from',
                       help='Angles will be parsed from the filename pattern.'
                            'The special token {TA} should be specified as part'
                            'of the pattern, that will be used to match the '
                            'value of the angle for each TiltSeriesMovie.')

    def _defineAcquisitionParams(self, form):
        """ Add movie specific options to the acquisition section. """
        group = ProtImportTsBase._defineAcquisitionParams(self, form)
        group.addParam('gainFile', params.FileParam,
                       label='Gain image',
                       help='A gain reference related to a set of movies '
                            'for gain correction')
        group.addParam('darkFile', params.FileParam,
                       label='Dark image',
                       help='A dark image related to a set of movies')
        return group

    def _fillAcquisitionInfo(self, inputTs):
        ProtImportTsBase._fillAcquisitionInfo(self, inputTs)
        inputTs.setGain(self.gainFile.get())
        inputTs.setDark(self.darkFile.get())

        if not self.MDOC_DATA_SOURCE:
            acq = inputTs.getAcquisition()
            acq.setDoseInitial(self.doseInitial.get())
            acq.setDosePerFrame(self.dosePerFrame.get())

    def _initialize(self):
        ProtImportTsBase._initialize(self)
        self._outputName += 'M'
        self._createOutputName += 'M'

    def _validateAngles(self):
        """ Function to be implemented in subclass to validate
        the angles range.
        """
        if not self.MDOC_DATA_SOURCE and self.getEnumText('anglesFrom') == self.ANGLES_FROM_FILENAME:
            if not self._anglesInPattern():
                return ['When importing movies, {TA} and {TO} should be in the '
                        'files pattern.']
        else:
            return []


class MDoc:

    def __init__(self, fileName, voltage=None, magnification=None, samplingRate=None,
                 doseProvidedByUser=None, tiltAngleProvidedByUser=None):

        self._mdocFileName = fileName
        self._tsId = None
        # Dose related attributes
        self.doseProvidedByUser = doseProvidedByUser
        self.mdocHasDose = False
        # Acquisition general attributes
        self._voltage = voltage
        self._magnification = magnification
        self._samplingRate = samplingRate
        self._tiltAxisAngle = tiltAngleProvidedByUser
        # Acquisition specific attributes (per angle)
        self._tiltsMetadata = []

    @staticmethod
    def normalizeTSId(rawTSId):
        """ Normalizes the name of a TS to prevent sqlite errors, it ends up as a table in a set"""
        # remove paths and extension
        normTSID = removeBaseExt(rawTSId)

        # Avoid dots case: TS_234.mrc.mdoc
        normTSID = normTSID.split(".")[0]

        if normTSID[0].isdigit():
            normTSID = "TS_" + normTSID

        return normTSID

    def read(self, isImportingTsMovies=True, ignoreFilesValidation=False):
        validateTSFromMdocErrMsg = ''
        tsFile = None
        mdoc = self._mdocFileName
        headerDict, zSlices = self._parseMdoc()

        # Get acquisition general info
        self._getAcquisitionInfoFromMdoc(headerDict, zSlices[0])
        self._tsId = self.normalizeTSId(mdoc)
        parentFolder = getParentFolder(mdoc)
        if not isImportingTsMovies:
            # Some mdoc files point to an .st file stored in the ImageFile header line
            tsFile = join(parentFolder, headerDict.get("ImageFile", None))
            if not os.path.exists(tsFile):
                tsFile = join(parentFolder, self._tsId + '.mrcs')
            if not os.path.exists(tsFile):
                tsFile = join(parentFolder, self._tsId + '.mrc')
            if not os.path.exists(tsFile):
                tsFile = join(parentFolder, self._tsId + '.st')
            validateTSFromMdocErrMsg = self._validateTSFromMdoc(mdoc, tsFile)

        # Get acquisition specific (per angle) info
        self._getSlicesData(zSlices, tsFile)

        # Check Mdoc info read
        validateMdocContentsErrorMsgList = self._validateMdocInfoRead(
            ignoreFilesValidation=ignoreFilesValidation or not isImportingTsMovies)

        # Check all the possible errors found
        exceptionMsg = ''
        if validateTSFromMdocErrMsg:
            exceptionMsg += validateTSFromMdocErrMsg
        if validateMdocContentsErrorMsgList:
            exceptionMsg += ' '.join(validateMdocContentsErrorMsgList)

        # If ts images we are assuming slices follow the angle order
        if not isImportingTsMovies and not exceptionMsg:
            self._tiltsMetadata.sort(key=lambda x: float(x.getTiltAngle()),
                                     reverse=False)
            for index, tiMd in enumerate(self._tiltsMetadata):
                tiMd.setAngleMovieFile((index+1, tiMd.getAngleMovieFile()))

        return exceptionMsg

    def _parseMdoc(self):
        """
        Parse the mdoc file and return a list with a dict key=value for each
        of the [Zvalue = X] sections and a dictionary for the first lines global variables.

        :return: dictionary (header), list of dictionaries (Z slices)
        """
        headerDict = {}
        headerParsed = False
        zvalueList = []

        with open(self._mdocFileName) as f:
            for line in f:
                if line.startswith('[ZValue'):
                    # We have found a new z value
                    headerParsed = True
                    zvalue = int(line.split(']')[0].split('=')[1])
                    if zvalue != len(zvalueList):
                        raise Exception("Unexpected ZValue = %d" % zvalue)
                    zvalueDict = {}
                    zvalueList.append(zvalueDict)
                elif line.startswith('[T') and not self.getTiltAxisAngle():
                    strLine = line.strip().replace(' ', '').lower()
                    pattern = 'tiltaxisangle='
                    if pattern in strLine:
                        # Example of the most common syntax (after having checked multiple mdocs from EMPIAR)
                        # [T =     Tilt axis angle = 90.1, binning = 1  spot = 9  camera = 0]
                        tiltAxisAngle = strLine.split('tiltaxisangle=')[1].split(',')[0]
                        # Check if it's a string which represents a float or not
                        if tiltAxisAngle.replace('.', '', 1).isdigit():
                            self._tiltAxisAngle = float(tiltAxisAngle)
                elif line.strip():
                    key, value = line.split('=')
                    if not headerParsed:
                        headerDict[key.strip()] = value.strip()
                    if zvalueList:
                        zvalueDict[key.strip()] = value.strip()

        return headerDict, zvalueList

    def _getAcquisitionInfoFromMdoc(self, headerDict, firstSlice):
        """Acquisition data is read from to data sources (from higher to lower priority):
            - From the values introduced in the form by the user (introduced as attributes of the MDoc object)
            - From the first ZSlice data.
            - From the file header data."""
        if not self.getVoltage():
            VOLTAGE = 'Voltage'
            self._voltage = firstSlice.get(VOLTAGE, headerDict.get(VOLTAGE, self._voltage))
        if not self.getMagnification():
            MAGNIFICATION = 'Magnification'
            self._magnification = firstSlice.get(MAGNIFICATION, headerDict.get(MAGNIFICATION, self._magnification))
        if not self.getSamplingRate():
            PIXEL_SPACING = 'PixelSpacing'
            self._samplingRate = firstSlice.get(PIXEL_SPACING, headerDict.get(PIXEL_SPACING, self._samplingRate))

    def _getSlicesData(self, zSlices, tsFile):
        parentFolder = getParentFolder(self._mdocFileName)
        accumulatedDose = 0
        for counter, zSlice in enumerate(zSlices):
            if self.doseProvidedByUser:
                incomingDose = self.doseProvidedByUser
            else:
                incomingDose = self._getDoseFromMdoc(zSlice, self.getSamplingRate())
            accumulatedDose += incomingDose
            self._tiltsMetadata.append(TiltMetadata(
                angle=zSlice.get('TiltAngle', None),
                angleFile=self._getAngleMovieFileName(parentFolder, zSlice, tsFile),
                acqOrder=counter+1,
                accumDose=accumulatedDose,
                incomingDose=incomingDose
            ))
        if round(accumulatedDose) > 0:  # round is used to make the condition more robust, for cases like 0.0000000001
            self.mdocHasDose = True


    @staticmethod
    def _getAngleMovieFileName(parentFolder, zSlice, tsFile):
        if tsFile:
            return tsFile
        else:
            # PureWindowsPath pathlib is ised to make possible deal with different path separators, like \\
            return join(parentFolder, PureWindowsPath(zSlice['SubFramePath']).parts[-1])

    @staticmethod
    def _getDoseFromMdoc(zSlice, pixelSize):
        """It calculates the accumulated dose on the frames represented by zSlice, and add it to the
        previous accumulated dose"""

        EXPOSURE_DOSE = 'ExposureDose'  # Dose on specimen during camera exposure in electrons/sq. Angstrom
        FRAME_DOSES_AND_NUMBERS = 'FrameDosesAndNumber'  # Dose per frame in electrons per square Angstrom followed
        # by number of frames at that dose
        DOSE_RATE = 'DoseRate'  # Dose rate to the camera, in electrons per unbinned pixel per second
        EXPOSURE_TIME = 'ExposureTime'  # Image exposure time
        MIN_MAX_MEAN = 'MinMaxMean'  # Minimum, maximum, and mean value for this image
        COUNTS_PER_ELECTRON = 'CountsPerElectron'
        DIVIDED_BY_TWO = 'DividedBy2'

        def _keysInDict(listOfKeys):
            return all([key in zSlice.keys() for key in listOfKeys])

        # Different ways of calculating the dose, ordered by priority considering the possible variability between
        # different mdoc files
        newDose = 0

        if pixelSize:
            pixelSize = float(pixelSize)
        else:
            pixelSize = 1  # This case cover the possibility of no sampling rate in form nor mdoc, avoiding the error
            # execution before the whole information is gathered and the exception is raised

        # Directly from field ExposureDose
        if EXPOSURE_DOSE in zSlice:
            expDoseVal = zSlice[EXPOSURE_DOSE]
            if expDoseVal:
                newDose = float(expDoseVal)

        # Directly from field FrameDosesAndNumbers
        if not newDose and FRAME_DOSES_AND_NUMBERS in zSlice:
            frameDoseAndNums = zSlice[FRAME_DOSES_AND_NUMBERS]
            if frameDoseAndNums:
                data = frameDoseAndNums.split()  # Get the mean from a string like '0 6'
                dosePerFrame = data[0]
                nFrames = data[1]
                newDose = float(dosePerFrame) * float(nFrames)

        # Calculated from fields DoseRate and ExposureTime
        if not newDose and _keysInDict([DOSE_RATE, EXPOSURE_TIME]):
            doseRate = zSlice[DOSE_RATE]
            expTime = zSlice[EXPOSURE_TIME]
            if doseRate and expTime:
                newDose = float(doseRate) * float(expTime) / pixelSize ** 2

        # Calculated from fields MinMaxMean, PixelSpacing and CountsPerElectron
        if not newDose and _keysInDict([MIN_MAX_MEAN, COUNTS_PER_ELECTRON]):
            minMaxMean = zSlice[MIN_MAX_MEAN]
            counts = zSlice[COUNTS_PER_ELECTRON]
            if all([minMaxMean, pixelSize, counts]):
                meanVal = minMaxMean.split()[-1]  # Get the mean from a string like '-42 2441 51.7968'
                newDose = (float(meanVal) / float(counts)) / pixelSize ** 2

        # # Calculated as in Grigorieff paper --> https://doi.org/10.7554/eLife.06980.001
        # if not newDose and _keysInDict([MIN_MAX_MEAN, COUNTS_PER_ELECTRON, EXPOSURE_TIME, DIVIDED_BY_TWO]):
        #     minMaxMean = zSlice[MIN_MAX_MEAN]
        #     counts = zSlice[COUNTS_PER_ELECTRON]
        #     expTime = zSlice[EXPOSURE_TIME]
        #     divByTwo = zSlice[DIVIDED_BY_TWO]
        #     if all([minMaxMean, counts, expTime]):
        #         meanVal = minMaxMean.split()[-1]  # Get the mean from a string like '-42 2441 51.7968'
        #         newDose = _getDivByTwoFactor(divByTwo) * float(meanVal) / (float(counts) * float(expTime))

        divByTwo = 0
        if _keysInDict([DIVIDED_BY_TWO]):
            divByTwo = int(zSlice[DIVIDED_BY_TWO])
        divByTwoFactor = 2 if divByTwo else 1

        return newDose * divByTwoFactor

    @staticmethod
    def _validateTSFromMdoc(mdoc, tsFile):
        errMsg = ''
        if not exists(tsFile):
            errMsg = '\nMdoc --> %s\nExpected tilt series file not found \n%s' % (mdoc, tsFile)

        return errMsg

    def _validateMdocInfoRead(self, ignoreFilesValidation=False):
        validateMdocContentsErrorMsgList = []
        msg = ['\n*Data not found in file*\n%s:\n' % self._mdocFileName]
        missingFiles = []
        missingAnglesIndices = []
        for i, tiltMetadata in enumerate(self._tiltsMetadata):
            # Check the angles
            if not tiltMetadata.getTiltAngle():
                missingAnglesIndices.append(str(i))
            # Check the angle stack files read
            if not ignoreFilesValidation:
                # Ignore the files validation is sometimes used for test purposes
                file = tiltMetadata.getAngleMovieFile()
                if not exists(file):
                    missingFiles.append(file)

        if not self._voltage:
            msg.append('*Voltage*\n')
        if not self._magnification:
            msg.append('*Magnification*\n')
        if not self._samplingRate:
            msg.append('*PixelSpacing*\n')
        if not self.mdocHasDose:
            msg.append('Not able to get the *dose* with the data read.\n'
                       'Dose related mdoc labels are:\n\n'
                       '- ExposureDose or\n'
                       '- FrameDosesAndNumber or\n'
                       '- DoseRate and ExposureTime or\n'
                       '- MinMaxMean and CountsPerElectron'
                       )
        if not self.getTiltAxisAngle():
            msg.append('*RotationAngle (tilt axis angle)*')
        if missingAnglesIndices:
            msg.append('*TiltAngle*: %s\n' % ' '.join(missingAnglesIndices))
        if missingFiles:
            msg.append('*Missing files*:\n%s\n' % ' '.join(missingFiles))
        if len(msg) > 1:
            validateMdocContentsErrorMsgList.append(' '.join(msg))

        return validateMdocContentsErrorMsgList

    def getFileName(self):
        return self._mdocFileName

    def getTsId(self):
        return self._tsId

    def getVoltage(self):
        return self._voltage

    def getMagnification(self):
        return self._magnification

    def getSamplingRate(self):
        return self._samplingRate

    def getTiltsMetadata(self):
        return self._tiltsMetadata

    def getTiltAxisAngle(self):
        return self._tiltAxisAngle


class TiltMetadata:

    def __init__(self, angle=None, angleFile=None, acqOrder=None, accumDose=None, incomingDose=None):
        self._angle = angle
        self._angleFile = angleFile
        self._acqOrder = acqOrder
        self._accumDose = accumDose
        self._incomingDose = incomingDose

    def setTiltAngle(self, tiltAngle):
        self._angle = tiltAngle

    def setAngleMovieFile(self, angleMovieFile):
        self._angleFile = angleMovieFile

    def setAcqOorder(self, order):
        self._acqOrder = order

    def setAccumDose(self, accumDose):
        self._accumDose = accumDose

    def setIncomingDose(self, incDose):
        self._incomingDose = incDose

    def getTiltAngle(self):
        return self._angle

    def getAngleMovieFile(self):
        return self._angleFile

    def getAcqOrder(self):
        return self._acqOrder

    def getAccumDose(self):
        return self._accumDose

    def getAcqOrder(self):
        return self._acqOrder

    def getIncomingDose(self):
        return self._incomingDose

