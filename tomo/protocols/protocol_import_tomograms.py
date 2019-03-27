# coding=utf-8

import pyworkflow.em as pwem
from tomo.objects import Tomogram
from .protocol_base import ProtTomoBase
from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Transform
from os.path import exists, basename, abspath
from pyworkflow.em.convert import Ccp4Header
from pyworkflow.utils.path import createAbsLink


class ProtImportTomograms(pwem.ProtImportVolumes, ProtTomoBase):
    """Protocol to import a set of tomograms to the project"""
    _outputClassName = 'SetOfTomograms'
    _label = 'import tomograms'

    def __init__(self, **args):
        pwem.ProtImportVolumes.__init__(self, **args)

    def _insertAllSteps(self):
        self._insertFunctionStep('importTomogramsStep',
                                 self.getPattern(),
                                 self.samplingRate.get(),
                                 self.setOrigCoord.get())

    # --------------------------- STEPS functions -----------------------------

    def importTomogramsStep(self, pattern, samplingRate, setOrigCoord=False):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)

        # Create a Volume template object
        tomo = Tomogram()
        tomo.setSamplingRate(samplingRate)

        imgh = ImageHandler()

        tomoSet = self._createSetOfTomograms()
        tomoSet.setSamplingRate(samplingRate)

        for fileName, fileId in self.iterFiles():
            x, y, z, n = imgh.getDimensions(fileName)
            if fileName.endswith('.mrc') or fileName.endswith('.map'):
                fileName += ':mrc'
                if (z == 1 and n != 1):
                    zDim = n
                    n = 1
                else:
                    zDim = z
            else:
                zDim = z
            origin = Transform()
            if setOrigCoord:
                origin.setShiftsTuple(self._getOrigCoord())
            else:
                origin.setShifts(x/-2. * samplingRate,
                            y/-2. * samplingRate,
                            zDim/-2. * samplingRate)

            tomo.setOrigin(origin)  # read origin from form

            if self.copyFiles or setOrigCoord:
                newFileName = abspath(self._getVolumeFileName(fileName, "mrc"))
                Ccp4Header.fixFile(fileName, newFileName, origin.getShifts(),
                                   samplingRate, Ccp4Header.ORIGIN)
            else:
                newFileName = abspath(self._getVolumeFileName(fileName))

                if fileName.endswith(':mrc'):
                    fileName = fileName[:-4]
                createAbsLink(fileName, newFileName)
            if n == 1:
                tomo.cleanObjId()
                tomo.setFileName(newFileName)
                tomoSet.append(tomo)
            else:
                for index in range(1, n+1):
                    tomo.cleanObjId()
                    tomo.setLocation(index, newFileName)
                    tomoSet.append(tomo)

        if tomoSet.getSize() > 1:
            self._defineOutputs(outputTomograms=tomoSet)
        else:
            self._defineOutputs(outputTomogram=tomo)


    # --------------------------- INFO functions ------------------------------

    def _getTomMessage(self):
        if self.hasAttribute('outputTomogram'):
            return "Tomogram %s" % self.getObjectTag('outputTomogram')
        else:
            return "Tomograms %s" % self.getObjectTag('outputTomograms')

    def _summary(self):
        summary = []
        if self.hasAttribute('outputTomogram') or \
                self.hasAttribute('outputTomograms'):
            summary.append("%s imported from:\n%s" % (self._getTomMessage(),
                           self.getPattern()))

            summary.append(u"Sampling rate: *%0.2f* (Å/px)" %
                           self.samplingRate.get())
        return summary

    def _methods(self):
        methods = []
        if self.hasAttribute('outputTomogram') or \
                self.hasAttribute('outputTomograms'):
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getTomMessage(), self.samplingRate.get()),)
        return methods

    def _getTomogramFileName(self, fileName, extension=None):
        if extension is not None:
            baseFileName="import_" + basename(fileName).split(".")[0] + ".%s"%extension
        else:
            baseFileName="import_" + basename(fileName).split(":")[0]

        return self._getExtraPath(baseFileName)

    def _getOrigCoord(self):
        return -1.*self.x.get(), -1.*self.y.get(), -1.*self.z.get()