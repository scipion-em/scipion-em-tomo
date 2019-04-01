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

import pyworkflow.em as pwem


def writeTiStack(inputTiList, outputStackFn, outputTltFn=None):
    """ Write a given list of tilt images as a single stack
    Params:
        inputTiList: input list of tilted images.
        outputStackFn: output path where to write the stack.
        orderBy: column to sort by, by default tilt angle (ascending)
        excludeList: a list of images to skip (starting at 1)
    Results:
        A new stack file will be created and also a tlt file with tilt-angles
    """
    ih = pwem.ImageHandler()
    i = 0
    f = open(outputTltFn, 'w') if outputStackFn else None

    for ti in inputTiList:
        i += 1
        ih.convert(ti, (i, outputStackFn))
        if f:
            f.write('%f\n' % ti.getTiltAngle())

    if f:
        f.close()


