# **************************************************************************
# *
# * Authors:     Alberto García Mena (alberto.garcia@cnb.csic.es) [1]
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

from pwem.protocols.protocol_import.base import ProtImportFiles, ProtImport
import pyworkflow.protocol.params as params

import time
import subprocess


class ProtImportTsStreaming(ProtImport):
    """ class for Tilt-Series import protocols in streaming.
    """

    def __init__(self, **args):
        ProtImport.__init__(self, **args)


    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Root directory of the tilt-series "
                           "(or movies) files.")
        form.addParam('filesPattern', params.StringParam,
                      label='Pattern',
                      help="It determines if the tilt series are going to "
                           "be imported using the mdoc file or the tilt "
                           "series files. To import from the mdoc files, "
                           "the word '.mdoc' must appear in the pattern, "
                           "if not, a tilt series pattern is expected. "
                           "In the first case, the angular and acquisition "
                           "data are directly read from the corresponding "
                           "mdoc file, while in the second it is read "
                           "the base name of the matching files, according to "
                           " the pattern introduced.\n\n"
                           "*IMPORTING WITH MDOC FILES*\n\n"
                           "For *tilt series movies*, a mdoc per tilt series "
                           "movies is expected. "
                           "The corresponding movie file/s must be located in "
                           "the same path as the mdoc file. The tilt series "
                           "id will be the base name of the mdoc files, "
                           "so the names of the mdoc files must be different, "
                           "even if they're located in "
                           "different paths.\n\n"
                           "For *tilt series*, the only difference is that a "
                           "stack .mrcs file is expected for each "
                           "mdoc, which means, per each tilt series desired "
                           "to be imported.\n\n"
                           "*IMPORTING WITH A PATTERN OF THE TILT SERIES FILE "
                           "NAMES*\n\n"
                           "The pattern can contain standard wildcards such "
                           "as *, ?, etc.\n\n"
                           "It should also contains the following special "
                           "tags:\n"
                           "   *{TS}*: tilt series identifier, which can be "
                           "any UNIQUE part of the path. This must be "
                           "an alpha-numeric sequence (avoid symbols as -) "
                           "that can not start with a number.\n"
                           "   *{TO}*: acquisition order, an integer value "
                           "(important for dose).\n"
                           "   *{TA}*: tilt angle, a positive or negative "
                           "float value.\n\n"
                           "Examples:\n\n"
                           "To import a set of image stacks (tilt-series "
                           "or tilt-series movies) as: \n"
                           "TiltSeries_a_001_0.0.mrc\n"
                           "TiltSeries_a_002_3.0.mrc\n"
                           "TiltSeries_a_003_-3.0.mrc\n"
                           "...\n"
                           "TiltSeries_b_001_0.0.mrc\n"
                           "TiltSeries_b_002_3.0.mrc\n"
                           "TiltSeries_b_003_-3.0.mrc\n"
                           "...\n"
                           "The pattern TiltSeries_{TS}_{TO}_{TA}.mrc will "
                           "identify:\n"
                           "{TS} as a, b, ...\n"
                           "{TO} as 001, 002, 003, ...\n"
                           "{TA} as 0.0, 3.0, -3.0, ...\n")
        form.addParam('exclusionWords', params.StringParam,
                      label='Exclusion words:',
                      help="List of words separated by a space that the path "
                           "should not have",
                      expertLevel=params.LEVEL_ADVANCED)


    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=[], wait=True)


    def findTSComplete(self):
        pass


    def matchTSPaths(self):
        pass

    def sortTS(self):
        pass

