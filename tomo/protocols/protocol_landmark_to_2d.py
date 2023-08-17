# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import os.path

from pwem.emlib.image import ImageHandler
from pyworkflow.protocol.params import PointerParam
from pwem.protocols import EMProtocol
from pwem.objects import SetOfMicrographs, Micrograph, SetOfCoordinates, Coordinate
from pyworkflow.object import Integer, Pointer
from functools import lru_cache
import tomo.objects as tomoObjs


class ProtTomoLandmarksTo2D(EMProtocol):
    """ Converts a set of landmark models into SPA 2d objects. Landmarks will become
    a SetOfCoordinates2D and the tilt series associated will become a micrograph set associated
    to the coordinates. Since SPA methods will not handle micrograph in stacks, tilt series files will
    be unstacked."""

    _label = "Landmarks to 2D coordinates"

    # Possible outputs definition
    OUTPUT_MICS_NAME = "Micrographs"
    OUTPUT_COORDS_NAME = "Coordinates2D"
    _possibleOutputs = {OUTPUT_COORDS_NAME: SetOfCoordinates, OUTPUT_MICS_NAME:SetOfMicrographs}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._micsCache = {} # To cache mics by mic name

    def _defineParams(self, form):

        form.addSection("Input")

        form.addParam('landmarks', PointerParam,
                      pointerClass=tomoObjs.SetOfLandmarkModels,
                      label='Landmarks',
                      help='Set of Landmarks to convert')

    def getLandmarks(self) -> tomoObjs.SetOfLandmarkModels:
        return self.landmarks.get()

    def _insertAllSteps(self):

        self._insertFunctionStep(self.tiltSeriesToMics)
        self._insertFunctionStep(self.landmarksToCoords2D)

    def tiltSeriesToMics(self):
        """ Converts tilt series to micrographs"""
        mics = SetOfMicrographs.create(self.getPath())

        # Get the tilt series of the landmark model
        tsSet = self.getLandmarks().getSetOfTiltSeries()

        # remove extra properties the micrograph set does not define
        acquisition = tsSet.getAcquisition()
        del acquisition._accumDose

        # Previous version of TomoAcquisition had _angleAxis1 and _angleAxis2
        if hasattr(acquisition, "_angleAxis1"):
            del acquisition._angleAxis1
            del acquisition._angleAxis2

        del acquisition._angleMax
        del acquisition._angleMin
        del acquisition._step
        del acquisition._tiltAxisAngle

        mics.copyInfo(tsSet)

        ih = ImageHandler()

        # Iterate the tilt series associated to the landmark
        for ts in tsSet.iterItems():

            # for each tilt image
            for ti in ts:

                # Get the file
                slice, stackFn = ti.getLocation()
                stackFn = os.path.basename(stackFn)
                titlAngleSafe = str(ti.getTiltAngle()).replace(".", "_")
                micFn = self._getExtraPath("S%s_A%s_%s" % (slice,titlAngleSafe,stackFn))
                # Extract it as a single mrcfile
                self.info("Extracting tilt image %s (%s) from %s to %s" % (slice, ti.getTiltAngle(), stackFn, micFn))
                ih.convert(ti, micFn)

                # Create the micrograph
                newMic = Micrograph(location=micFn)
                newMic.setSamplingRate(ti.getSamplingRate())
                newMic.setMicName(self.composeMicName(ti))

                mics.append(newMic)

        mics.write()

        self._defineOutputs(**{self.OUTPUT_MICS_NAME: mics})
        self._defineSourceRelation(tsSet, mics)

    def composeMicName(self, ti:tomoObjs.TiltImage):

        return "%s_A%s" % (ti.getTsId(), ti.getTiltAngle())

    @lru_cache
    def getTiltImageFromLandmark(self, ts:tomoObjs.TiltSeries, slice, tsId ):
        """

        :param ts: TiltSeries
        :param slice:
        :param tsId: TsId (Only use for cache)
        :return:
        """

        self.info("Retrieving tilt image %s@%s" % (slice, tsId))
        return ts[slice].clone()

    @lru_cache
    def getMicFromLandmark(self, ti: tomoObjs.TiltImage):

        micName = self.composeMicName(ti)
        self.info("Retrieving micrograph with micname %s." % micName)

        if True or micName not in self._micsCache:

            # Get the micrographs
            mics = getattr(self, self.OUTPUT_MICS_NAME)

            mic = mics[{'_micName': micName}]

            self._micsCache[micName] = mic.clone()

        return self._micsCache[micName]

    def landmarksToCoords2D(self):
        """ Converts landmarks to 2D coordinates"""
        coords = SetOfCoordinates.create(self.getPath())

        # Get the tilt series of the landmark model
        landmarks = self.getLandmarks()

        firstItem = landmarks.getFirstItem()

        coords.setBoxSize(firstItem.getSize())
        coords.setMicrographs(Pointer(self, extended=self.OUTPUT_MICS_NAME))

        # Iterate the landmark models
        for landmark_model in landmarks.iterItems():

            landmarks.completeLandmarkModel(landmark_model)
            ts = landmark_model.getTiltSeries()
            self.info("Creating 2d coordinates from %s landmark model." % ts.getTsId())

            # for each tilt image
            for landmark_line in landmark_model.retrieveInfoTable():

                newCoord =  Coordinate()
                newCoord.setX(landmark_line[0])
                newCoord.setY(landmark_line[1])
                ti = self.getTiltImageFromLandmark(ts,landmark_line[2], ts.getTsId())
                mic = self.getMicFromLandmark(ti)
                newCoord.setMicrograph(mic)
                newCoord._chainId = Integer(landmark_line[3])

                coords.append(newCoord)

        coords.write()

        self._defineOutputs(**{self.OUTPUT_COORDS_NAME: coords})
        self._defineSourceRelation(landmarks, coords)
