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
import logging
import os
import re
from glob import glob
from datetime import datetime
from collections import OrderedDict
from os.path import join
from statistics import mean
import numpy as np
from sqlite3 import OperationalError
import pyworkflow as pw
import pyworkflow.protocol
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
import tomo.objects
from pwem.objects import Transform
from pyworkflow.object import Integer
from pyworkflow.utils.properties import Message
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtImport
from tomo.convert import getAnglesFromHeader, getAnglesFromMdoc, getAnglesAndDosesFromTlt
from tomo.convert.mdoc import normalizeTSId, MDoc
from tomo.objects import TomoAcquisition, SetOfTiltSeries, SetOfTiltSeriesM
from .protocol_base import ProtTomoBase

logger = logging.getLogger(__name__)


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
    ANGLES_FROM_TLT = 'Tlt file'
    ANGLES_FROM_RANGE = 'Range'

    NOT_MDOC_GUI_COND = ('filesPattern is None or ' +
                         '(filesPattern is not None and ".mdoc" ' +
                         'not in filesPattern)')

    MDOC_DATA_SOURCE = False

    acquisitions = None
    sRates = None
    accumDoses = None
    incomingDose = None
    meanDosesPerFrame = None

    OUTPUT_NAME = 'outputTiltSeries'
    _possibleOutputs = {OUTPUT_NAME: SetOfTiltSeries}

    def __init__(self, **args):
        ProtImport.__init__(self, **args)
        ProtTomoBase.__init__(self)
        self.skippedMdocs = Integer()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Root directory of the tilt-series "
                           "(or movies) files.")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      help="This determines if the tilt series / movies are going to "
                           "be imported using the mdoc files or the tilt "
                           "series files. To import from the mdoc files, "
                           "the word '.mdoc' must appear in the pattern, "
                           "if not, a tilt series pattern is expected. "
                           "In the first case, the angular and acquisition "
                           "data are read from the corresponding "
                           "mdoc files, while in the second case they are read from "
                           "the name of the matching files.\n\n"
                           "*IMPORTING WITH MDOC FILES*\n\n"
                           "For *tilt series movies*, ONE mdoc per tilt series "
                           "is expected. The movie files must be located in "
                           "the same folder as the mdocs. The filenames will be "
                           "fetched from the _SubFramePath_ value in mdoc. \n"
                           "Example pattern: _TS*.mdoc_\n\n"
                           "For *tilt series*, ONE _mrcs_ stack should match ONE "
                           "mdoc file per each tilt series. To import unstacked "
                           "images use the filename pattern (see below) instead of mdoc.\n\n"
                           "*IMPORTING WITH A FILENAME PATTERN (tilt series and "
                           "movies)*\n\nThe pattern can contain wildcards such "
                           "as *, ?, etc. It should also contain the following "
                           "special tags:\n\n"
                           "   *{TS}*: tilt series identifier, which can be "
                           "any UNIQUE part of the path. This must be "
                           "an alpha-numeric sequence (avoid dash (-) symbol) "
                           "and can not start with a number.\n"
                           "   *{TO}*: acquisition order, an integer value "
                           "(important for dose information).\n"
                           "   *{TA}*: tilt angle, a positive or negative "
                           "float value.\n\n"
                           "Example:\n\n"
                           "To import a set of images (tilt-series "
                           "or tilt-series movies) like: \n"
                           "TiltSeries_a_001_0.0.mrc\n"
                           "TiltSeries_a_002_3.0.mrc\n"
                           "TiltSeries_a_003_-3.0.mrc\n"
                           "...\n"
                           "TiltSeries_b_001_0.0.mrc\n"
                           "TiltSeries_b_002_3.0.mrc\n"
                           "TiltSeries_b_003_-3.0.mrc\n"
                           "...\n"
                           "Use pattern TiltSeries_{TS}_{TO}_{TA}.mrc, which will "
                           "identify:\n"
                           "{TS} as a, b, ...\n"
                           "{TO} as 001, 002, 003, ...\n"
                           "{TA} as 0.0, 3.0, -3.0, ...\n")
        form.addParam('isTomo5', params.BooleanParam,
                      label='Tomo5 mdoc?',
                      default=False,
                      condition='not (%s)' % self.NOT_MDOC_GUI_COND,
                      help="If these mdocs were generated by the Tomography 5 software, check this box to ensure that the tilt axis angle is converted properly: -1 * TiltAxisAngle - 90")        
        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the path "
                           "should not have",
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('mdocInfo', params.LabelParam,
                      condition='not (%s)' % self.NOT_MDOC_GUI_COND,
                      label='Acquisition values provided below will override '
                            'the corresponding mdoc values',
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
                      help="This parameters determine how the project will "
                           "deal with imported files. It can be: \n"
                           "*Copy files*: Input files will be copied into "
                           "your project. (this will duplicate the raw data)."
                           "*Absolute symlink*: Create symbolic links to the "
                           "absolute path of the files."
                           "*Relative symlink*: Create symbolic links as "
                           "relative path from the protocol run folder. ")

        self._defineAcquisitionParams(form)

        self._defineBlacklistParams(form)

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        pass

    def _defineAcquisitionParams(self, form):
        """ Define acq parameters, it can be overridden
        by subclasses to change what parameters to include.
        """
        group = form.addGroup('Acquisition info - '
                              'override mdoc values if provided')
        group.addParam('voltage', params.FloatParam,
                       label=Message.LABEL_VOLTAGE,
                       allowsNull=True,
                       help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', params.FloatParam, default=2.7,
                       label=Message.LABEL_SPH_ABERRATION,
                       help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
                       label=Message.LABEL_AMPLITUDE,
                       help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', params.IntParam, default=50000,
                       label=Message.LABEL_MAGNI_RATE,
                       expertLevel=params.LEVEL_ADVANCED,
                       allowsNull=True,
                       help=Message.TEXT_MAGNI_RATE)
        group.addParam('samplingRate', params.FloatParam,
                       label=Message.LABEL_SAMP_RATE,
                       allowsNull=True,
                       help=Message.TEXT_SAMP_RATE)
        group.addParam('tiltAxisAngle', params.FloatParam,
                       label='Tilt axis angle (deg.)',
                       allowsNull=True,
                       help="The rotation angle is the angle from the vertical "
                            "to the axis of tilting, where counterclockwise is "
                            "positive.\n See "
                            "https://bio3d.colorado.edu/imod/doc/tomoguide.html#UnknownAxisAngle")
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

    # -------------------------- INSERT functions -----------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.importStep)

    # -------------------------- STEPS functions ------------------------------
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

        outputSet = getattr(self, self.OUTPUT_NAME, None)

        if outputSet is None:
            createSetFunc = getattr(self, self._createOutputName)
            outputSet = createSetFunc()
        elif outputSet.getSize() > 0:
            outputSet.loadAllProperties()

            for ts in outputSet:
                self._existingTs.add(ts.getTsId())

        self._fillAcquisitionInfo(outputSet)
        setAcq = outputSet.getAcquisition()
        tsClass = outputSet.ITEM_TYPE
        tiClass = tsClass.ITEM_TYPE

        if self._existingTs:
            outputSet.enableAppend()

        # Get files that matches the pattern via mdoc or not
        matchingFiles = self.getMatchingFiles()

        # Go through all of them
        for ts, tiltSeriesList in matchingFiles.items():

            self.info("Tilt series found: %s" % ts)
            someNew = True
            tsObj = tsClass(tsId=ts)
            # Form value has higher priority than the mdoc values
            samplingRate = \
                float(samplingRate if samplingRate else self.sRates[ts])

            origin = Transform()
            tsObj.setOrigin(origin)
            tsObj.setAnglesCount(len(tiltSeriesList))
            self.setItemExtraAttributes(tsObj)

            # we need this to set mapper before adding any item
            outputSet.append(tsObj)

            if self.MDOC_DATA_SOURCE:
                accumDoseList = self.accumDoses[ts]
                incomingDoseList = self.incomingDose[ts]
                counter = 0

            tiltSeriesObjList = []

            # Acquisition exists for both the tilt images and the tilt series, but some values will be different when
            # referred to the whole TS than when referred to a tilt image (such as the accumDose). Angle max, angle mix,
            # and step will be generated here as it will be the same for both the TS and the tilt images.
            tiltAngles = sorted([float(tiData[2]) for tiData in tiltSeriesList])
            tsAcq = tsObj.getAcquisition().clone()
            maxTilt = tiltAngles[-1]
            minTilt = tiltAngles[0]
            step = round(mean([tiltAngles[i + 1] - tiltAngles[i] for i in range(len(tiltAngles) - 1)]))

            setAcq.setAngleMin(minTilt)
            setAcq.setAngleMax(maxTilt)
            setAcq.setStep(step)

            tsAcq.setAngleMin(minTilt)
            tsAcq.setAngleMax(maxTilt)
            tsAcq.setStep(step)

            # Add each tilt images to the tiltSeries
            for f, to, ta, accDose in tiltSeriesList:
                try:
                    # Link/move to extra
                    imageFile = f[1] if type(f) is tuple else f

                    # Double underscore is used in EMAN to determine set type e.g. phase flipped particles. We replace
                    # it by a single underscore to avoid possible problems if the user uses EMAN
                    finalDestination = self._getExtraPath(os.path.basename(imageFile))
                    finalDestination = finalDestination.replace('__', '_')
                    self.copyOrLink(imageFile, finalDestination)

                    f = (f[0], finalDestination) if type(f) is tuple else finalDestination

                    ti = tiClass(location=f,
                                 acquisitionOrder=to,
                                 tiltAngle=ta)

                    tiAcq = tsAcq.clone()
                    # Calculate the dose
                    if self.MDOC_DATA_SOURCE:
                        dosePerFrame = incomingDoseList[counter]
                        accDose = accumDoseList[counter]
                        initialDose = self.doseInitial.get() if counter == 0 else sum(incomingDoseList[:counter + 1])
                    else:
                        dosePerFrame = self.dosePerFrame.get()
                        if not accDose:
                            accDose = to * dosePerFrame
                        initialDose = self.doseInitial.get() if to == 1 else accDose - dosePerFrame
                    # Initial dose in current ti
                    tiAcq.setDoseInitial(initialDose)
                    # Incoming dose in current ti
                    tiAcq.setDosePerFrame(dosePerFrame)
                    # Accumulated dose in current ti
                    tiAcq.setAccumDose(accDose)
                    ti.setAcquisition(tiAcq)
                    tiltSeriesObjList.append(ti)
                    counter += 1

                except OperationalError:
                    raise Exception("%s is an invalid {TS} tag. "
                                    "It must be an alpha-numeric sequence "
                                    "(avoid symbols like -) that can not "
                                    "start with a number." % ts)

            # Sort tilt image metadata if importing tilt series
            if not self._isImportingTsMovies():
                tiltSeriesObjList.sort(key=lambda x: x.getTiltAngle(),
                                       reverse=False)

            for ti in tiltSeriesObjList:
                tsObj.append(ti)

            tsObjFirstItem = tsObj.getFirstItem()
            origin.setShifts(-tsObjFirstItem.getXDim() / 2 * samplingRate,
                             -tsObjFirstItem.getYDim() / 2 * samplingRate,
                             0)

            if self.MDOC_DATA_SOURCE:
                tsAcq.setAccumDose(accumDoseList[-1])
                # Tilt series object dose per frame has been updated each
                # time the tilt image dose per frame has
                # been updated before, so the mean value is used to be the
                # reference in the acquisition of the
                # whole tilt series movie
                meanDosePerFrame = mean(incomingDoseList)
                tsAcq.setDosePerFrame(meanDosePerFrame)
                setAcq.setDosePerFrame(meanDosePerFrame)
            else:
                accumDose = self.dosePerFrame.get() * len(tiltSeriesList)
                dosePerFrame = self.dosePerFrame.get()
                tiltAxisAngle = self.tiltAxisAngle.get()

                tsAcq.setDosePerFrame(dosePerFrame)
                tsAcq.setAccumDose(accumDose)
                tsAcq.setTiltAxisAngle(tiltAxisAngle)

                setAcq.setDosePerFrame(dosePerFrame)
                setAcq.setAccumDose(accumDose)
                setAcq.setTiltAxisAngle(tiltAxisAngle)

            tsObj.setAcquisition(tsAcq)
            outputSet.update(tsObj)  # update items and size info

            if self._isImportingTsMovies():
                dim = tsObjFirstItem.getDim()
                framesRange = [1, dim[2], 1]
                outputSet.setFramesRange(framesRange)

            self._existingTs.add(ts)
            someAdded = True

        if someAdded:
            self.debug('Updating output...')
            self._updateOutputSet(self.OUTPUT_NAME, outputSet,
                                  state=outputSet.STREAM_OPEN)
            self.debug('Update Done.')

        self.debug('Checking if finished...someNew: %s' % someNew)

        # Close the output set
        self._updateOutputSet(self.OUTPUT_NAME, outputSet,
                              state=outputSet.STREAM_CLOSED)

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []

        if self.getOutputsSize():
            for key, output in self.iterOutputAttributes():
                summary.append("Imported tilt-series %s from: %s" % (
                    "movies"
                    if output.getLastName() == "outputTiltSeriesM" else "",
                    self.filesPath.get()))
                summary.append("Using pattern: %s" % self.filesPattern.get())
                summary.append(u"Sampling rate: *%0.2f* (Å/px)"
                               % output.getSamplingRate())
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        if self.skippedMdocs.get():
            summary.append('*%i* mdoc files were skipped --> '
                           'check the Output Log tab for more details.'
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
        # In the mdoc case, voltage, magnification and sampling
        # rate are optional inputs. In the user introduces
        # one of these values, it will be considered more
        # prior than the corresponding value read from the mdoc file
        if not self.MDOC_DATA_SOURCE:
            if not self.voltage.get():
                errMsg.append('Voltage should be a float')
            if not self.magnification.get():
                errMsg.append('Magnification should be an integer')
            if not self.samplingRate.get():
                errMsg.append('Sampling rate should be a float')
            if not self.dosePerFrame.get():
                errMsg.append('Dose per tilt should be a float')
            if self.tiltAxisAngle.get() is None:
                errMsg.append('Tilt axis angle should be a float')

        return errMsg

    def _validateAngles(self):
        """ Function to be implemented in subclass to validate
        the angles range.
        """
        return None

    # -------------------------- BASE methods to be overridden ----------------
    # def _getImportChoices(self):
    #     """ Return a list of possible choices
    #     from which the import can be done.
    #     (usually packages form such as: xmipp3, eman2, relion...etc.
    #     """
    #     return ['files']

    # -------------------------- UTILS functions ------------------------------
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

            # Handle the special characters from the pattern before compiling the regexp that will be used
            # later for the file matching with the results obtained from the glob module file search
            specialCharDict = {
                '*': '(.*)',
                '[': '\\[',  # Escape the brackets to make the regex compilation identify them as bracket characters
                ']': '\\]'
            }
            inStr = self._pattern
            for specChar, newVal in specialCharDict.items():
                inStr = inStr.replace(specChar, newVal)
            self._regexPattern = _replace(inStr,
                                          r'(?P<TS>.*)',  # regex pattern for TS
                                          r'(?P<TO>\d+)',  # regex pattern for TO
                                          r'(?P<TA>[+-]?\d+(\.\d+)?)')  # regex pattern for TA
            self._regex = re.compile(self._regexPattern)
            self._globPattern = _replace(self._pattern, '*', '*', '*')
            # Glob module does not handle well the brackets (it does not list them)
            self._globPattern = self._globPattern.replace('[', '*').replace(']', '*')

            # Set output names depending on the import type
        # (either movies or images)
        self._createOutputName = '_createSetOfTiltSeries'

        # Keep track of which existing tilt-series has already been found
        self._existingTs = set()

    def _anglesInPattern(self):
        """ This function should be called after a call to _initialize"""
        return '{TA}' in self._pattern and '{TO}' in self._pattern

    def _getMatchingFilesFromMdoc(self, isValidation):
        """If the list of files provided by the user is
           a list of mdoc files, then the tilt series movies
        are built from them, following the considerations listed below:
            - For each mdoc file, it and the corresponding movie files
              must be in the same directory.
            - The tilt series id will be the base name of the mdoc file,
              by default, so the mdocs must have different
              base name. If another name is desired, the user can introduce
              the name structure (see advanced parameter)
            """
        fpath = self.filesPath.get()
        mdocList = glob(join(fpath, self.filesPattern.get()))  # Get matching files by the introduced file pattern
        mdocList = self._excludeByWords(mdocList)  # Check for exclusion words
        hasDoseList = []
        if not mdocList:
            raise Exception(f'There are no mdoc files matching the pattern '
                            f'{join(fpath, self.filesPattern.get())}')

        matchingFiles = OrderedDict()
        self.acquisitions = OrderedDict()
        self.sRates = OrderedDict()
        self.accumDoses = OrderedDict()
        self.incomingDose = OrderedDict()
        warningHeadMsg = 'The following mdoc files were skipped:\n'
        warningDetailedMsg = []
        skippedMdocs = 0

        for mdoc in mdocList:
            # Note: voltage, magnification and sampling rate values are the
            # ones introduced by the user in the protocol's form.
            # Otherwise, the corresponding values considered will be the ones
            # read from the mdoc.
            # This is because you can't trust mdoc
            # (often dose is not calibrated in serialEM, so you get 0;
            # pixel size might be binned as mdoc comes from a binned record
            # not movie and there are no Cs and amp contrast fields in mdoc)
            mdocObj = MDoc(
                mdoc,
                voltage=self.voltage.get() if self.voltage.get() else None,
                magnification=(self.magnification.get()
                               if self.magnification.get()
                               else None),
                samplingRate=(self.samplingRate.get()
                              if self.samplingRate.get()
                              else None),
                doseProvidedByUser=(self.dosePerFrame.get()
                                    if self.dosePerFrame.get()
                                    else None),
                tiltAngleProvidedByUser=(self.tiltAxisAngle.get()
                                         if self.tiltAxisAngle.get()
                                         else None))
            validationError = mdocObj.read(
                isImportingTsMovies=self._isImportingTsMovies())
            hasDoseList.append(mdocObj.mdocHasDose)
            if validationError:
                warningHeadMsg += '    %s\n' % mdoc
                warningDetailedMsg.append(validationError)
                skippedMdocs += 1
                # validationErrors.append(validationError)
                # Continue parsing the remaining mdoc files to
                # provide a fully detailed error message
                continue
            acquisition = self._genTsAcquisitionFromMdoc(mdocObj)
            tsId = mdocObj.getTsId()
            fileOrderAngleList = []
            accumulatedDoseList = []
            incomingDoseList = []
            for tiltMetadata in mdocObj.getTiltsMetadata():
                fileOrderAngleList.append((
                    tiltMetadata.getAngleMovieFile(),  # Filename
                    '{:03d}'.format(tiltMetadata.getAcqOrder()),  # Acquisition order
                    tiltMetadata.getTiltAngle(),  # Tilt angle
                    tiltMetadata.getAccumDose()  # Accumulated dose
                ))
                accumulatedDoseList.append(tiltMetadata.getAccumDose())
                incomingDoseList.append(tiltMetadata.getIncomingDose())

            # We make Tomo5 TiltAxisAngle match the SerialEM convention, see
            # https://groups.google.com/g/eman2/c/piux0GdX_7M/m/ScQNh9TNAQAJ
            if self.isTomo5 and not self.tiltAxisAngle.get():
                acquisition.setTiltAxisAngle(-1 * acquisition.getTiltAxisAngle() - 90)

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
                    self.info(warningHeadMsg + ' '.join(warningDetailedMsg))
                return matchingFiles
            else:
                raise Exception(warningHeadMsg + ' '.join(warningDetailedMsg))
        else:
            return matchingFiles

    def _isImportingTsMovies(self):
        return True if type(self) is ProtImportTsMovies else False

    def _genTsAcquisitionFromMdoc(self, mdocObj):
        acq = \
            TomoAcquisition(voltage=mdocObj.getVoltage(),
                            sphericalAberration=self.sphericalAberration.get(),
                            amplitudeContrast=self.amplitudeContrast.get(),
                            magnification=mdocObj.getMagnification(),
                            tiltAxisAngle=mdocObj.getTiltAxisAngle()
                            )

        # This field is only present in the form for TsM import
        if hasattr(self, 'doseInitial'):
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
                logger.info("%s excluded. Contains any of %s" %
                            (file, exclusionWords))
                continue
            allowedFiles.append(file)

        return allowedFiles

    def getMatchingFiles(self, isValidation=False):
        """ Return an ordered dict with TiltSeries found in the files as key
        and a list of all tilt images of that series as value.
        """
        if self.MDOC_DATA_SOURCE:
            return self._getMatchingFilesFromMdoc(isValidation=isValidation)
        else:
            return self._getMatchingFilesFromRegExPattern()

    def _getMatchingFilesFromRegExPattern(self):
        filePaths = glob(self._globPattern)
        filePaths = self._excludeByWords(filePaths)

        filePaths.sort(key=lambda fn: os.path.getmtime(fn))

        matchingFiles = OrderedDict()

        def _getTsId(match):
            """ Retrieve the TiltSerie ID from the matching object.
            We need to have tsId that starts with character, so
            let's add a prefix if it is not the case
            """
            tsId = match.group('TS')
            return normalizeTSId(tsId)

        def _addOne(fileList, file, match):
            """ Add one file matching to the list. """

            order = int(match.group('TO'))
            dose = float(self.dosePerFrame.get()) * order
            angle = float(match.group('TA'))
            fileList.append((file, order, angle, dose))

        def _addMany(fileList, file, match):
            """ Add many 'files' (when angles in header or mdoc)
                to the list. """
            anglesFrom = self.getEnumText('anglesFrom')
            doses = []
            tiltorders = []

            if anglesFrom == self.ANGLES_FROM_HEADER:
                angles = getAnglesFromHeader(file)
            elif anglesFrom == self.ANGLES_FROM_MDOC:
                mdocFn = os.path.splitext(file)[0] + '.mdoc'
                if not os.path.exists(mdocFn):
                    raise Exception("Missing angles file: %s" % mdocFn)
                angles = getAnglesFromMdoc(mdocFn)
            elif anglesFrom == self.ANGLES_FROM_TLT:
                angles, doses, tiltorders = self.getFromTlt(file)
            elif anglesFrom == self.ANGLES_FROM_RANGE:
                angles = self._getTiltAngleRange()
            else:
                raise Exception('Invalid angles option: %s' % anglesFrom)

            for i, a in enumerate(angles):
                order = i + 1 if not tiltorders else tiltorders[i]
                dose = doses[i] if doses else int(self.dosePerFrame.get()) * order
                fileList.append(((i + 1, file), order, a, dose))

        addFunc = _addOne if self._anglesInPattern() else _addMany

        # Handle special case of just one TiltSeries, to avoid
        # the user the need to specify {TS}
        if len(filePaths) == 1 and not self.isInStreaming():
            f = filePaths[0]
            self.info("Single match: %s" % f)
            ts = pwutils.removeBaseExt(f)  # Base name without extension
            self.info("Raw tilt series id is %s." % ts)
            ts = normalizeTSId(ts)
            self.info("Normalized tilt series id is %s." % ts)
            matchingFiles[ts] = []
            _addMany(matchingFiles[ts], f, None)
        else:
            for f in filePaths:
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
        """ Returns a function to copy or link files
            based on user selected option"""

        if self.importAction.get() == self.IMPORT_COPY_FILES:
            return pw.utils.copyFile
        elif self.importAction.get() == self.IMPORT_LINK_REL:
            return pw.utils.createLink
        else:
            return pw.utils.createAbsLink

    def copyOrLink(self, source, destination):
        """ Calls the copy or link method chosen by
            the user in importAction option"""
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

    def getFromTlt(self, file):
        """ Gets angles and doses from a tlt or rawtlt file """
        tltFn = os.path.splitext(file)[0] + '.tlt'
        rawtltFn = tltFn.replace(".tlt", ".rawtlt")
        if os.path.exists(tltFn):
            angles, doses, tiltorders = getAnglesAndDosesFromTlt(tltFn)
        elif os.path.exists(rawtltFn):
            angles, doses, tiltorders = getAnglesAndDosesFromTlt(rawtltFn)
        else:
            raise Exception("Missing angles file: %s or %s" % (tltFn, rawtltFn))

        return angles, doses, tiltorders

    @classmethod
    def worksInStreaming(cls):
        # Streaming is done through SPA + tilt series composer
        return False

    def _fillAcquisitionInfo(self, inputTs):
        # TODO Acquisition is historically expected as
        # one per set of tilt series movies,
        # so the first one is the one used, at least for now
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
        return np.arange(
            self.minAngle.get(),
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

    def setItemExtraAttributes(self, tsObj):
        """
        To set extra possible attributes of the TiltSerie or TiltSerieM that may have been defined in the form
        :param tsObj:
        :return: None
        """
        pass


class ProtImportTs(ProtImportTsBase):
    """Protocol to import tilt series."""
    _label = 'import tilt-series'
    _devStatus = pw.BETA

    def _defineAngleParam(self, form: pyworkflow.protocol.Form):
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
                            "They can be taken from a range (Min, Max, Step) "
                            "or from the image header, or from an"
                            "mdoc or tlt file (should have the SAME filename "
                            "but with the .mdoc or .tlt or .rawtlt "
                            "extension at the end). If a tlt or rawtlt file is used, "
                            "it is optional to pass the accumulated dose as second "
                            "column beside each angle separated by space")

        line = group.addLine('Tilt angles range',
                             condition='anglesFrom==0',  # ANGLES_FROM_RANGE
                             help="Specify the tilt angles range. "
                                  "The original acquisition order does not have "
                                  "to match the order of tilt angles.")
        line.addParam('minAngle', params.FloatParam, default=-60, label='min')
        line.addParam('maxAngle', params.FloatParam, default=60, label='max')
        line.addParam('stepAngle', params.FloatParam, default=3, label='step')

        form.addParam('ctfCorrected', params.BooleanParam, default=False,
                      label="Have images been CTF corrected?",
                      help="Select Yes if images have been CTF corrected")

        form.addParam('interpolated', params.BooleanParam, default=False,
                      label="Have images been aligned?",
                      help="Select Yes if images have been rotated/interpolated "
                           "using alignment information.")

    def setItemExtraAttributes(self, tsObj: tomo.objects.TiltSeries):
        """
        Sets ctf corrected parameter and interpolation status.

        :param tsObj: Tilt series instance
        :return: nothing
        """

        tsObj.setInterpolated(self.interpolated.get())
        tsObj.setCtfCorrected(self.ctfCorrected.get())

    def _validateAngles(self):
        if not self.MDOC_DATA_SOURCE:
            # If importing from pattern and getting the tilt data from header,
            # it has to check if IMOD's extract tilts is installed
            if self.getEnumText('anglesFrom') == self.ANGLES_FROM_HEADER:
                from pwem import Domain
                imod = Domain.importFromPlugin('imod')

                if imod is None:
                    return ['Imod plugin is needed to import angles from header.'
                            'Please install it']

            ts, tiltSeriesList = self._firstMatch
            i, fileName = tiltSeriesList[0][0]
            x, y, z, n = ImageHandler().getDimensions(fileName)
            nImages = max(z, n)  # Just handle ambiguity with mrc format
            nAngles = len(self._tiltAngleList)
            if nAngles != nImages:
                return 'Tilt-series %s stack has different number of images ' \
                       '(%d) than the expected number of tilt angles (%d). ' \
                    % (fileName, nImages, nAngles)
        else:
            return None


class ProtImportTsMovies(ProtImportTsBase):
    """Protocol to import tilt series movies."""
    _label = 'import tilt-series movies'
    _devStatus = pw.BETA

    OUTPUT_NAME = ProtImportTsBase.OUTPUT_NAME + "M"
    _possibleOutputs = {OUTPUT_NAME: SetOfTiltSeriesM}

    def _defineAngleParam(self, form):
        """ Used in subclasses to define the option to fetch tilt angles. """
        group = form.addGroup('Tilt info',
                              condition=False)

        group.addHidden('anglesFrom', params.EnumParam,
                        default=0,
                        choices=[self.ANGLES_FROM_FILENAME],
                        display=params.EnumParam.DISPLAY_HLIST,
                        label='Import angles from',
                        help='Angles will be parsed from the filename pattern.'
                             'The special token {TA} should be specified as '
                             ' part of the pattern.'
                        )

    def _defineAcquisitionParams(self, form):
        """ Add movie specific options to the acquisition section. """
        group = ProtImportTsBase._defineAcquisitionParams(self, form)
        group.addParam('gainFile', params.FileParam,
                       label='Gain image',
                       help='A gain reference related to a set of movies '
                            'for gain correction')
        group.addParam('darkFile', params.FileParam,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Dark image',
                       help='A dark image related to a set of movies')
        return group

    def _fillAcquisitionInfo(self, inputTs):
        ProtImportTsBase._fillAcquisitionInfo(self, inputTs)
        inputTs.setGain(self.gainFile.get())
        inputTs.setDark(self.darkFile.get())

    def _initialize(self):
        ProtImportTsBase._initialize(self)
        self._createOutputName += 'M'

    def _validateAngles(self):
        """ Function to be implemented in subclass to validate
        the angles range.
        """
        if not self.MDOC_DATA_SOURCE and \
                self.getEnumText('anglesFrom') == self.ANGLES_FROM_FILENAME:
            if not self._anglesInPattern():
                return 'When importing movies, {TS}, {TA} and {TO} ' \
                       'should be present in the pattern.'
        else:
            return None
