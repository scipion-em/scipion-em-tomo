# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.constants import CENTER_GRAVITY
from tomo.protocols import ProtTomoBase

METADATA_INPUT_COORDINATES = "fiducialCoordinates.xmd"


class ProtTsConvertCoordinates3d(EMProtocol, ProtTomoBase):
    """
    Converts a set of tilt-series 3D coordinates into a standard set of 3D
    coordinates associated with a corresponding set of tomograms.

    AI Generated:

    Convert Tilt-Series Coordinates to Tomogram Coordinates
    (ProtTsConvertCoordinates3d) — User Manual
        Overview

        The ProtTsConvertCoordinates3d protocol transforms a
        SetOfTiltSeriesCoordinates into a SetOfCoordinates3D linked to a
        SetOfTomograms.

        In tomography workflows, 3D coordinates obtained during tilt-series
        alignment usually represent fiducial positions in the coordinate
        system of the tilt-series. However, downstream tomogram-based
        processing requires coordinates expressed directly in the spatial
        reference frame of reconstructed tomograms.

        This protocol performs that conversion.

        Biological Context

        During tilt-series alignment, fiducial markers or reference points
        are commonly tracked to estimate image transformations.

        These coordinates are useful during alignment itself, but later
        processing stages—such as particle picking, subtomogram extraction,
        or spatial annotation—require coordinates attached to the final
        reconstructed tomograms.

        ProtTsConvertCoordinates3d bridges this transition by mapping
        tilt-series coordinates into tomogram-associated coordinates.

        Inputs

        The protocol requires two input objects:

            - a SetOfTiltSeriesCoordinates
            - a SetOfTomograms

        The input coordinates are expected to come from a previous
        tilt-series alignment step.

        The input tomograms provide:

            - the tomogram spatial reference
            - the tomogram sampling rate
            - the tomogram identity used for coordinate assignment

        Matching Strategy

        Matching between coordinates and tomograms is based on tsId.

        Internally, the protocol creates a dictionary where:

            - each tomogram is indexed by its tsId
            - each input coordinate searches for the tomogram with the
              same tsId

        Only coordinates whose tsId exists in the tomogram set are
        converted.

        Coordinates without a corresponding tomogram are ignored.

        Coordinate Conversion

        For each matched coordinate:

            1. A new Coordinate3D object is created.
            2. The coordinate is linked to the corresponding tomogram.
            3. X, Y, and Z values are converted using the tomogram
               sampling rate.
            4. The original score value is preserved.
            5. The tomogram object identifier is stored.

        The conversion rescales coordinates according to:

            coordinate_tomogram = coordinate_tilt_series / sampling_rate

        This step transforms physical coordinate values into the spatial
        grid of the reconstructed tomogram.

        Output Generation

        The output is a SetOfCoordinates3D.

        The output set is initialized with:

            - the sampling rate of the input tomograms
            - the tomogram set as precedents
            - a default box size of 32
            - STREAM_OPEN status during processing

        Converted coordinates are appended progressively.

        Once conversion finishes, the output set is written to disk.

        Streaming Behavior

        After all coordinates have been processed, the protocol executes
        a final closing step.

        During this step:

            - the output stream state is changed to STREAM_CLOSED
            - the output is stored persistently

        This behavior makes the protocol compatible with standard Scipion
        streaming conventions.

        Practical Interpretation

        From a practical tomography perspective, this protocol does not
        perform reconstruction or geometric refinement.

        Its role is purely referential:

            - preserve the original spatial locations
            - express them in the tomogram coordinate frame
            - make them usable by tomogram-based downstream protocols

        This is especially useful when fiducials or alignment landmarks
        need to be propagated from the alignment stage into later
        reconstruction analysis.

        Output Summary

        Once completed, the protocol reports:

            - number of input tilt-series coordinates
            - number of generated tomogram-associated 3D coordinates

        If the protocol has not yet generated outputs, the summary simply
        reports that no output coordinates are available.

        Final Perspective

        ProtTsConvertCoordinates3d is a lightweight but important
        conversion protocol in cryo-electron tomography workflows.

        Its main purpose is to transfer 3D coordinate information from
        tilt-series alignment space into reconstructed tomogram space.

        Although computationally simple, this step is often essential
        before particle picking, landmark propagation, or subtomogram
        extraction workflows.
    """

    _label = 'Tilt-series convert coords3D'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfCoordinates',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeriesCoordinates',
                      important=True,
                      label='Input set of coordinates 3D',
                      help='Set of 3D coordinates indicating the position in space of the fiducials. This set should '
                           'be obtained from the previous alignment step of the tilt-series.')

        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.convertCoordinates)
        self._insertFunctionStep(self.closeOutputSetStep)

    # --------------------------- STEPS functions ----------------------------

    def convertCoordinates(self):
        sotsc3d = self.inputSetOfCoordinates.get()

        sr = self.inputSetOfTomograms.get().getSamplingRate()

        self.getOutputSetOfCoordinates3Ds()
        tomoDict = self.getTomoDict()

        for coor3d in sotsc3d:
            tsId = coor3d.getTsId()

            if tsId in tomoDict.keys():
                tomo = tomoDict[tsId]

                newCoord3D = tomoObj.Coordinate3D()
                newCoord3D.setVolume(tomo)
                newCoord3D.setX(coor3d.getX()/sr, CENTER_GRAVITY)
                newCoord3D.setY(coor3d.getY()/sr, CENTER_GRAVITY)
                newCoord3D.setZ(coor3d.getZ()/sr, CENTER_GRAVITY)
                newCoord3D.setScore(coor3d.getScore())
                
                newCoord3D.setVolId(tomo.getObjId())
                self.outputSetOfCoordinates3D.append(newCoord3D)
                self.outputSetOfCoordinates3D.update(newCoord3D)

        self.outputSetOfCoordinates3D.write()

        self._store()

    def closeOutputSetStep(self):
        self.outputSetOfCoordinates3D.setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------

    def getOutputSetOfCoordinates3Ds(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()

        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.inputSetOfTomograms.get(),
                                                                      suffix='Coords3d')

            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTomograms.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTomograms.get())
            outputSetOfCoordinates3D.setBoxSize(32)

            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTomograms.get(), outputSetOfCoordinates3D)

        return self.outputSetOfCoordinates3D

    def getTomoDict(self):
        tomoDict = {}

        for tomo in self.inputSetOfTomograms.get():
            t = tomo.clone()
            tomoDict[tomo.getTsId()] = t

        return tomoDict

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        if not hasattr(self, 'outputSetOfCoordinates3D'):
            summary.append("No output coordinates generated yet")

        else:
            summary.append("Input tilt-series 3d coordinates: %d\n"
                           "Output 3d coordinates associated to a set of tomograms: %d" %
                           (self.inputSetOfCoordinates.get().getSize(),
                            self.outputSetOfCoordinates3D.getSize()))

        return summary

    def _methods(self):
        methods = []

        if not hasattr(self, 'outputSetOfCoordinates3D'):
            methods.append("No output coordinates generated yet")

        else:
            methods.append("%d 3d coordinates associated to a set of tomograms have been generated from the %d input "
                           "tilt-series 3d coordinates.\n" %
                           (self.outputSetOfCoordinates3D.getSize(),
                            self.inputSetOfCoordinates.get().getSize()))
        return methods
