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
logger = logging.getLogger(__name__)

import numpy
from pwem.emlib.image import ImageHandler
from .convert import *


def writeTiStack(inputTiList, outputStackFn, outputTltFn=None,
                 excludeList=None):
    """ Write a given list of tilt images as a single stack
    Params:
        inputTiList: input list of tilted images.
        outputStackFn: output path where to write the stack.
        orderBy: column to sort by, by default tilt angle (ascending)
        excludeList: a list of indexes of the images to skip (starting at 1)
    Results:
        A new stack file will be created and also a tlt file with tilt-angles
    """
    excludeList = excludeList or []
    ih = ImageHandler()
    j = 0
    f = open(outputTltFn, 'w') if outputStackFn else None

    for i, ti in enumerate(inputTiList):
        if i + 1 not in excludeList:
            j += 1
            ih.convert(ti, (j, outputStackFn))
            if f:
                f.write('%f\n' % ti.getTiltAngle())

    if f:
        f.close()


def getAnglesFromHeader(tsImage):
    """ Extract the tilt angles from the tilt-series stack using the
    IMOD program: extracttilts.
    """
    from pwem import Domain
    imod = Domain.importFromPlugin("imod", "Plugin")
    anglesFn = '/tmp/angles.txt'
    args = '--input %s --output %s' % (tsImage, anglesFn)
    pwutils.runJob(None, imod.getImodCmd('extracttilts'), args)
    angles = []
    with open(anglesFn) as f:
        angles = [float(line) for line in f]
    return angles


def parseMdoc(mdocFn):
    """
    Parse the mdoc file and return a list with a dict key=value for each
    of the [Zvalue = X] sections
    :param mdocFn: Path to the mdoc file
    :return: list of dictonaries
    """
    zvalueList = []

    with open(mdocFn) as f:
        for line in f:
            if line.startswith('[ZValue'):
                # We have found a new Zvalue
                zvalue = int(line.split(']')[0].split('=')[1])
                if zvalue != len(zvalueList):
                    raise Exception("Unexpected ZValue = %d" % zvalue)
                zvalueDict = {}
                zvalueList.append(zvalueDict)
            else:
                if line.strip() and zvalueList:
                    key, value = line.split('=')
                    zvalueDict[key.strip()] = value.strip()

    return zvalueList


def getAnglesFromMdoc(mdocFn):
    """ Return only the angles from the given mdoc file. """
    return [float(d['TiltAngle']) for d in parseMdoc(mdocFn)]


def getAnglesAndDosesFromTlt(tltFn):
    """ Parse the tilt-angles from tlt file. """

    logger.info("Reading %s file for angles and dose." % tltFn)
    angles = []
    doses = []
    orders = []
    with open(tltFn) as f:
        for line in f:
            strippedLine = line.strip()
            if not strippedLine:
                logger.info("Empty line found in %s. Ignoring it." % tltFn)
                continue
            line = strippedLine.split(" ")
            if line:
                angles.append(float(line[0]))
                # If there is a second column, we take it as dose
                if len(line) > 1:
                    doses.append(float(line[1]))
                #If there is a third column, we take it as tilt order
                if len(line) > 2:
                    orders.append(int(line[2]))

    logger.info("Angles found: %s" % angles)

    if doses:
        logger.info("Doses found: %s" % doses)

    if orders:
        logger.info("Tilt order found: %s" % orders)
    elif doses:
        # Calculate tilt order base on dose
        orders = getOrderFromList(doses)
        logger.info("Tilt orders inferred from dose list. %s" % orders)


    return angles, doses, orders

def getOrderFromList(unsortedList):
    """ Return a list with the position of each items in the sorted list
     Example: [ 8, 5, 1 ] --> [3,2,1]

     :param unsortedList: list to be sorted with ideally unique values

    """

    # First we sort items in the unsortedList getting the indices in the original unsortedList
    # Ex: [ 2, 1, 0 ] --> where [1, 5, 8] were in the original unsortedList
    sortedList = numpy.argsort(unsortedList)

    orderList = [0] * len(unsortedList)
    for i in range(len(unsortedList)):
        finalIndex = sortedList[i]
        orderList[finalIndex] = i+1

    return orderList