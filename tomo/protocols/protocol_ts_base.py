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


from pwem.protocols import EMProtocol
from pyworkflow.protocol import STEPS_PARALLEL, STATUS_NEW

from ..objects import TiltSeriesDict, Tomogram
from .protocol_base import ProtTomoBase


class ProtTsProcess(EMProtocol, ProtTomoBase):
    """
    Base framework for protocols that process tilt-series and generate new
    tilt-series or tomogram-related outputs.

    AI Generated:

    Tilt-Series Processing Base Framework (ProtTsProcess) — User Manual
        Overview

        ProtTsProcess is the generic processing backbone for tomography
        protocols that take tilt-series as input and execute computations
        either per tilt image, per tilt-series, or both.

        Rather than implementing a specific biological operation,
        this class defines the execution model shared by many tomography
        protocols.

        Typical derived protocols include:

            - tilt-series filtering
            - tilt-series correction
            - CTF estimation
            - tomogram reconstruction

        Its main role is to organize task scheduling, parallel execution,
        streaming-aware processing, and output generation.

        General Processing Model

        The protocol follows a hierarchical execution strategy.

        Processing is divided into three conceptual levels:

            1. Input preparation
            2. Tilt-image processing
            3. Tilt-series finalization

        This allows derived protocols to decide whether the operation
        should be performed:

            - image-by-image
            - once per tilt-series
            - or both

        The class uses parallel execution by default.

        Workflow Initialization

        At execution start, the protocol performs an initialization step.

        It then inserts:

            - one input conversion step
            - one final output-closing step

        Between these two, all processing tasks are inserted dynamically.

        Internally, a TiltSeriesDict object tracks:

            - discovered tilt-series
            - tilt images belonging to each tilt-series
            - completed processing tasks
            - finished tilt-series ready to be written

        This makes the protocol suitable for both standard and
        incremental processing workflows.

        Dynamic Step Insertion

        One of the most important features of ProtTsProcess is that
        processing steps are inserted dynamically as new tilt-series
        become available.

        For each newly discovered tilt-series:

            - optional per-tilt-image steps are inserted
            - one final per-tilt-series step is inserted

        This means subclasses do not need to explicitly manage task
        dependencies.

        The framework automatically ensures that:

            - all tilt-image steps finish first
            - the tilt-series step runs afterward
            - output writing waits until the tilt-series is complete

        Tilt Image vs Tilt-Series Processing

        By default, the protocol inserts one step per tilt image.

        This behavior is controlled by:

            _doInsertTiltImageSteps()

        Default behavior:

            True

        Derived protocols can disable per-image execution if their logic
        naturally works at tilt-series level.

        This is particularly useful for reconstruction algorithms,
        where processing is normally performed on the whole tilt-series
        rather than projection-by-projection.

        Output Generation

        Output generation is incremental.

        Whenever one or more tilt-series are completed:

            - the output set is created if necessary
            - the finished tilt-series are appended
            - metadata is updated
            - the output remains in STREAM_OPEN state

        This allows partially completed results to become available
        during execution.

        Once all tilt-series are done:

            - the final output step is unlocked
            - the output stream is closed

        This streaming-aware design is particularly useful in large
        tomography datasets.

        Subclass Responsibilities

        ProtTsProcess only provides the workflow skeleton.

        Derived protocols are expected to implement the actual logic.

        The most commonly overridden methods are:

            convertInputStep(...)
                Performs input preparation before processing.

            processTiltImageStep(...)
                Performs per-image operations.

            processTiltSeriesStep(...)
                Performs final per-series operations.

        Additional customization is possible through:

            _createOutputSet(...)
            _updateOutputSet(...)
            _getArgs(...)
            _initialize(...)

        This modular design makes the framework highly reusable.

        Practical Interpretation

        Biologically, ProtTsProcess should not be interpreted as a
        scientific protocol itself.

        It does not modify data directly.

        Instead, it provides the computational orchestration required
        by many tomography processing tasks.

        Its practical importance lies in:

            - consistent scheduling
            - robust parallelization
            - correct dependency management
            - safe incremental output writing

        For developers of tomography protocols, it provides a standard
        execution architecture that avoids reimplementing processing
        logistics in every protocol.

        ------------------------------------------------------------

        Tomogram Reconstruction Variant (ProtTomoReconstruct)
        Overview

        ProtTomoReconstruct is a specialization of ProtTsProcess
        designed for tomogram reconstruction.

        Unlike general tilt-series processing protocols, reconstruction
        normally operates on the complete tilt-series as a single unit.

        For this reason:

            _doInsertTiltImageSteps() -> False

        No per-image processing steps are inserted.

        Reconstruction Output

        Instead of generating a new SetOfTiltSeries,
        this subclass produces a SetOfTomograms.

        For each completed tilt-series:

            - a Tomogram object is created
            - the reconstructed file is linked using the expected name
            - the tomogram is appended to the output set

        Output sampling rate is inherited from the input tilt-series.

        If binning is applied:

            output_sampling = input_sampling × bin

        This ensures geometric consistency between the input data and
        the reconstructed tomograms.

        Output Naming

        Each reconstructed tomogram follows the convention:

            <tsId>_tomo.mrc

        This makes the output directly traceable to the original
        tilt-series.

        Final Perspective

        ProtTomoReconstruct adapts the general ProtTsProcess execution
        framework to the specific logic of tomographic reconstruction.

        Its main difference is conceptual:

            - process the tilt-series as a whole
            - generate tomograms instead of tilt-series

        This class therefore serves as the common reconstruction base
        for tomography reconstruction protocols in Scipion.
    """
    stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._initialize()

        inputTsM = self._getInputTs()

        self._ciStepId = self._insertFunctionStep('convertInputStep',
                                                  inputTsM.getObjId())
        self._insertFunctionStep('createOutputStep', wait=True,
                                 prerequisites=[self._ciStepId])
        self._coStep = self._steps[-1]  # get last step

        self._tsDict = TiltSeriesDict(inputTsM, self._getOutputSet(),
                                      newItemsCallback=self._insertNewSteps,
                                      doneItemsCallback=self._updateOutput)
        self._tsDict.update()

    def _insertNewSteps(self, tsIdList):
        """ Insert processing steps for newly discovered tilt-series. """
        for tsId in tsIdList:
            tsSteps = []
            if self._doInsertTiltImageSteps():
                for i, ti in enumerate(self._tsDict.getTiList(tsId)):
                    tiStep = self._insertFunctionStep(
                        'processTiltImageStep', tsId, ti.getObjId(),
                        *self._getArgs(), prerequisites=[self._ciStepId])
                    tsSteps.append(tiStep)

            if not tsSteps:  # no Ti steps used
                tsSteps.append(self._ciStepId)
            tsStepId = self._insertFunctionStep('processTiltSeriesStep', tsId,
                                                prerequisites=tsSteps)
            self._coStep.addPrerequisites(tsStepId)

        self.updateSteps()

    def _stepsCheck(self):
        self._tsDict.update()

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, inputId):
        """ To be implemented in subclasses. """
        pass

    def processTiltImageStep(self, tsId, tiltImageId, *args):
        """ To be implemented in subclasses. """
        pass

    def processTiltSeriesStep(self, tsId):
        """ To be implemented in subclasses. """
        pass

    def _updateOutputSet(self, outputSet, tsIdList):
        """ Update the output set with the finished Tilt-series.
        Params:
            :param tsIdList: list of ids of finished tasks.
        """
        for tsId in tsIdList:
            ts = self._tsDict.getTs(tsId)
            outputSet.append(ts)
            for ti in self._tsDict.getTiList(tsId):
                ts.append(ti)

            outputSet.update(ts)

    def _updateOutput(self, tsIdList):
        """ Update the output set with the finished Tilt-series.
        Params:
            :param tsIdList: list of ids of finished tasks.
        """
        # Flag to check the first time we save output
        self._createOutput = getattr(self, '_createOutput', True)

        outputSet = self._getOutputSet()

        if outputSet is None:
            # Special case just to update the outputSet status
            # but it only makes sense when there is outputSet
            if not tsIdList:
                return
            outputSet = self._createOutputSet()
        else:
            outputSet.enableAppend()
            self._createOutput = False

        # Call the sub-class method to update the output
        self._updateOutputSet(outputSet, tsIdList)
        outputSet.setStreamState(outputSet.STREAM_OPEN)

        if self._createOutput:
            outputSet.updateDim()
            self._defineOutputs(**{self._getOutputName(): outputSet})
            self._defineSourceRelation(self._getInputTsPointer(), outputSet)
            self._createOutput = False
        else:
            outputSet.write()
            self._store(outputSet)

        outputSet.close()

        if self._tsDict.allDone():
            self._coStep.setStatus(STATUS_NEW)

    def createOutputStep(self):

        self._closeOutputSet()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    # --------------------------- UTILS functions ------------------------------
    def _initialize(self):
        """ Allow sub-classes to make some initialization before all steps
        will be inserted in the execution list.
        """
        pass

    def _getInputTsPointer(self):
        return None

    def _getInputTs(self):
        """ Return the tiltSeries input object. """
        return self._getInputTsPointer().get()

    def _getOutputName(self):
        """ Return the output name, by default 'outputTiltSeries'.
        This method can be re-implemented in subclasses that have a
        different output. (e.g outputTomograms).
        """
        return 'outputTiltSeries'

    def _getOutputSet(self):
        return getattr(self, self._getOutputName(), None)

    def _createOutputSet(self, suffix=''):
        """ Method to create the output set.
        By default will a SetOfTiltSeries, but can be re-defined in subclasses.
        """
        outputSet = self._createSetOfTiltSeries(suffix=suffix)
        outputSet.copyInfo(self._getInputTs())
        return outputSet

    def _getArgs(self):
        """ Return a list with parameters that will be passed to the process
        TiltSeries step. It can be redefined by subclasses.
        """
        return []

    def _doInsertTiltImageSteps(self):
        """ Default True, but if return False, the steps for each
        TiltImage will not be inserted. """
        return True


class ProtTomoReconstruct(ProtTsProcess):
    """ Base class for Tomogram reconstruction protocols. """

    def _doInsertTiltImageSteps(self):
        # For reconstruction protocols, usually we don't
        # want one step per each tilt-image
        return False

    def _updateOutputSet(self, outputSet, tsIdList):
        """ Override this method to convert the TiltSeriesM into TiltSeries.
        """
        for tsId in tsIdList:
            t = Tomogram(location=self._getPath(self._getTomoName(tsId)))
            outputSet.append(t)

    def _createOutputSet(self):
        """ Create the output set of Tomograms. """
        outputSet = self._createSetOfTomograms()
        samplingRate = self._getInputTs().getSamplingRate()

        if self.bin > 1:
            samplingRate *= self.bin.get()

        outputSet.setSamplingRate(samplingRate)

        return outputSet

    def _getOutputName(self):
        return 'outputTomograms'

    # --------------------------- UTILS functions ----------------------------
    def _getInputTsPointer(self):
        return self.inputTiltSeries

    def _getTomoName(self, tsId):
        return '%s_tomo.mrc' % tsId
