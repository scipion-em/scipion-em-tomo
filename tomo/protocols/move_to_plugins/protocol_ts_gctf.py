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
import gctf
from gctf.protocols.program_gctf import ProgramGctf


class ProtTsGctf(ProtTsEstimateCTF):
    """
    CTF estimation on Tilt-Series using CTFFIND4.
    """
    _label = 'tiltseries gctf'

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineCtfParamsDict(self):
        ProtTsEstimateCTF._defineCtfParamsDict(self)
        self._gctfProgram = ProgramGctf(self)

    def _defineProcessParams(self, form):
        ProgramGctf.defineFormParams(form)

        form.addParallelSection(threads=3, mpi=1)

    # --------------------------- STEPS functions ----------------------------
    def _estimateCtf(self, workingDir, tiFn, ti):
        try:
            program, args = self._gctfProgram.getCommand(
                scannedPixelSize=self._params['scannedPixelSize']
            )
            args += ' %s/*.mrc' % workingDir
            self.runJob(program, args, env=gctf.Plugin.getEnviron())

            ext = self._gctfProgram.getExt()

            # Move files we want to keep
            micBase = pw.utils.removeBaseExt(tiFn)

            def _getFile(suffix):
                print("File: %s" % os.path.join(workingDir, micBase + suffix))
                return os.path.join(workingDir, micBase + suffix)

            # move output from tmp to extra
            pw.utils.moveFile(_getFile(ext),
                              self._getExtraPath(micBase + '_ctf.mrc'))
            pw.utils.moveFile(_getFile('_gctf.log'), self._getTmpPath())
            pw.utils.moveFile(_getFile('_EPA.log'), self._getTmpPath())

        except Exception as ex:
            print >> sys.stderr, "Some error happened with micrograph %s" % tiFn
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

    def getCtf(self, ti):
        """ Parse the CTF object estimated for this Tilt-Image
        """
        prefix = self.getTiPrefix(ti)
        psd = self._getExtraPath(prefix + '_ctf.mrc')
        outCtf = self._getTmpPath(prefix + '_gctf.log')
        return self._gctfProgram.parseOutputAsCtf(outCtf, psdFile=psd)

