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

import pyworkflow as pw
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params
from pyworkflow.protocol import STEPS_PARALLEL, STATUS_NEW

from pyworkflow.utils.properties import Message

from tomo.objects import TiltSeriesDict
from .protocol_base import ProtTomoBase


class ProtTsProcess(pwem.EMProtocol, ProtTomoBase):
    """
    Base class for Tilt-Series (images or movies) processing protocols.
    This class should be used by protocols that receive tilt-series as input
    and produce tilt-series as output. This class will contain some common
    functionality about the steps execution (also in streaming) and
    the output generation.
    """
    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._initialize()

        inputTs = self._getInputTs()

        self._ciStepId = self._insertFunctionStep('convertInputStep',
                                                inputTs.getObjId())
        self._insertFunctionStep('createOutputStep', wait=True,
                                 prerequisites=[self._ciStepId])
        self._coStep = self._steps[-1]  # get last step

        self._tsDict = TiltSeriesDict(inputTs, self._getOutputTs(),
                                      newItemsCallback=self._insertNewSteps,
                                      doneItemsCallback=self._updateOutput)
        self._tsDict.update()

    def _insertNewSteps(self, tsIdList):
        """ Insert processing steps for newly discovered tilt-series. """
        for tsId in tsIdList:
            tsSteps = []
            for i, ti in enumerate(self._tsDict.getTiList(tsId)):
                tiStep = self._insertFunctionStep('processTiltImageStep',
                                                  tsId, ti.getObjId(),
                                                  *self._getArgs(),
                                                  prerequisites=[self._ciStepId])
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

    def _updateOutputTsSet(self, outputSet, tsIdList):
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

        outputSet = self._getOutputTs()

        if outputSet is None:
            # Special case just to update the outputSet status
            # but it only makes sense when there is outputSet
            if not tsIdList:
                return
            inputTs = self._getInputTs()
            outputSet = self._createSetOfTiltSeries()
            outputSet.copyInfo(inputTs)
        else:
            outputSet.enableAppend()
            self._createOutput = False

        # Call the sub-class method to update the output
        self._updateOutputTsSet(outputSet, tsIdList)

        outputSet.setStreamState(outputSet.STREAM_OPEN)

        if self._createOutput:
            outputSet.updateDim()
            self._defineOutputs(outputTiltSeries=outputSet)
            self._defineSourceRelation(self._getInputTsPointer(), outputSet)
            self._createOutput = False
        else:
            outputSet.write()
            self._store(outputSet)

        outputSet.close()

        if self._tsDict.allDone():
            self._coStep.setStatus(STATUS_NEW)

    def createOutputStep(self):
        outputSet = self._getOutputTs()
        outputSet.setStreamState(outputSet.STREAM_CLOSED)
        outputSet.write
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

    def _getOutputTs(self):
        return getattr(self, 'outputTiltSeries', None)

    def _getArgs(self):
        """ Return a list with parameters that will be passed to the process
        TiltSeries step. It can be redefined by subclasses.
        """
        return []

