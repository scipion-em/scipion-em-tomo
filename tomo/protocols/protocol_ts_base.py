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
    Base class for Tilt-Series (images or movies) processing protocols.
    This class should be used by protocols that receive tilt-series as input
    and produce tilt-series as output. This class will contain some common
    functionality about the steps execution (also in streaming) and
    the output generation.
    """
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

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
        outputSet = self._getOutputSet()
        outputSet.setStreamState(outputSet.STREAM_CLOSED)
        outputSet.write()
        self._store(outputSet)

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
