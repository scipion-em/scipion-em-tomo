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
import sys

import pyworkflow as pw
import pyworkflow.em as pwem
from pyworkflow.protocol import STEPS_PARALLEL

from tomo.protocols import ProtTsEstimateCTF

from grigoriefflab.protocols.program_ctffind import ProgramCtffind


class ProtTsCtffind(ProtTsEstimateCTF):
    """
    CTF estimation on Tilt-Series using CTFFIND4.
    """
    _label = 'tiltseries ctffind'

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineCtfParamsDict(self):
        ProtTsEstimateCTF._defineCtfParamsDict(self)
        self._ctfProgram = ProgramCtffind(self)

    def _defineProcessParams(self, form):
        ProgramCtffind.defineFormParams(form)

        form.addParallelSection(threads=3, mpi=1)

    # --------------------------- STEPS functions ----------------------------
    def _estimateCtf(self, workingDir, ti):
        try:
            downFactor = self.ctfDownFactor.get()
            micFnMrc = os.path.join(workingDir, self.getTiPrefix(ti) + '.mrc')
            outputLog = os.path.join(workingDir, 'output-log.txt')
            outputPsd = os.path.join(workingDir, self.getPsdName(ti))

            ih = pw.em.ImageHandler()

            if downFactor != 1:
                # Replace extension by 'mrc' because there are some formats
                # that cannot be written (such as dm3)
                ih.scaleFourier(ti, micFnMrc, downFactor)
            else:
                ih.convert(ti, micFnMrc, pw.em.DT_FLOAT)

        except Exception as ex:
            print >> sys.stderr, "Some error happened: %s" % ex
            import traceback
            traceback.print_exc()
        try:
            program, args = self._ctfProgram.getCommand(
                micFn=micFnMrc,
                ctffindOut=outputLog,
                ctffindPSD=outputPsd
            )
            self.runJob(program, args)
            # Move files we want to keep
            pw.utils.moveFile(outputPsd, self._getExtraPath())
            pw.utils.moveFile(outputPsd.replace('.mrc', '.txt'), self._getTmpPath())
        except Exception as ex:
            print >> sys.stderr, "ctffind has failed with micrograph %s" % micFnMrc
            import traceback
            traceback.print_exc()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _getArgs(self):
        """ Return a list with parameters that will be passed to the process
        TiltSeries step. It can be redefined by subclasses.
        """
        return []

    def getPsdName(self, ti):
        return '%s_PSD.mrc' % self.getTiPrefix(ti)

    def getCtf(self, ti):
        """ Parse the CTF object estimated for this Tilt-Image
        """
        psd = self.getPsdName(ti)
        outCtf = self._getTmpPath(psd.replace('.mrc', '.txt'))
        return self._ctfProgram.parseOutputAsCtf(outCtf,
                                                 psdFile=self._getExtraPath(psd))

    def isNewCtffind4(self):
        # This function is needed because it is used in Form params condition
        return ProgramCtffind.isNewCtffind4()

