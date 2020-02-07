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

import pyworkflow

from .protocol_base import ProtTomoBase, ProtTomoPicking
from .protocol_ts_base import ProtTomoReconstruct
from .protocol_ts_import import ProtImportTsBase, ProtImportTs, ProtImportTsMovies
from .protocol_ts_correct_motion import ProtTsCorrectMotion, ProtTsAverage
from .protocol_ts_estimate_ctf import ProtTsEstimateCTF
from .protocol_import_tomograms import ProtImportTomograms
from .protocol_import_subtomograms import ProtImportSubTomograms
from .protocol_import_coordinates import ProtImportCoordinates3D

# Hack to allow creating subsets of Tomograms and Subtomograms from ShowJ
# We should modify ProtUserSubSet to remove hardcoded methods to create subsets
# This method/class can be pass as parameter to showj
# See https://github.com/I2PC/scipion/issues/2036
emprotocol = pyworkflow.em.protocol.EMProtocol
setattr(emprotocol, "_createSetOfClassesSubTomograms", ProtTomoBase._createSetOfClassesSubTomograms.__func__)
setattr(emprotocol, "_createSetOfSubTomograms", ProtTomoBase._createSetOfSubTomograms.__func__)
setattr(emprotocol, "_createSetOfTomograms", ProtTomoBase._createSetOfTomograms.__func__)
setattr(emprotocol, "_createSet", ProtTomoBase._createSet.__func__)
setattr(emprotocol, "_createSetOfMeshes", ProtTomoBase._createSetOfMeshes.__func__)

# This code extends ImageHandler allowing the use of scaling with splines until this functionality is implemented
# in Scipion
def scaleSplines(inputFn, outputFn, scaleFactor):
    """ Scale an image using splines. """
    import xmippLib
    I = xmippLib.Image(inputFn)
    x, y, z, _ = I.getDimensions()
    I.scale(int(x/scaleFactor), int(y/scaleFactor), int(z/scaleFactor))
    I.write(outputFn)

ih = pyworkflow.em.ImageHandler
setattr(ih, "scaleSplines", staticmethod(scaleSplines))
