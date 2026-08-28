# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.constants import NEW
from pwem.protocols import EMProtocol

from tomo.constants import BOTTOM_LEFT_CORNER, TR_DYNAMO
from tomo.utils import existsPlugin

EXPORT_TO_TXT = 'txt'
EXPORT_TO_STAR = 'relion'
EXPORT_TO_EMAN = 'eman'
EXPORT_TO_DYNAMO = 'dynamo'
EXPORT_TO_CBOX = 'cbox'


class ProtExportCoordinates3D(EMProtocol):
    """
    Export 3D Coordinates (ProtExportCoordinates3D) — User Manual

    Overview

    The Export 3D Coordinates protocol exports a set of subtomogram particle coordinates
    from Scipion into external file formats that can be used by other cryo-electron
    tomography software packages.

    Its main purpose is to make 3D particle positions generated or refined inside Scipion
    available for downstream processing, visualization, classification, or subtomogram
    averaging in external tomography environments.

    For a biological user, this protocol acts as an interoperability bridge. It allows
    coordinates identified in one workflow to be reused in other specialized packages
    without manual reformatting.

    Inputs and General Workflow

    The protocol requires a single input:

    - A `SetOfCoordinates3D`, containing particle coordinates associated with one or
      more tomograms.

    During execution, the protocol first identifies the tomograms represented in the
    coordinate set and then creates a clean export directory.

    The coordinates are then written tomogram by tomogram into the selected output format.

    This organization is biologically useful because subtomogram workflows are usually
    performed independently for each tomogram, and maintaining that separation preserves
    experimental traceability.

    Export Formats

    The protocol supports several export formats depending on which tomography plugins
    are available in the current Scipion installation.

    Always available:

    - TXT format

    Optionally available if corresponding plugins are installed:

    - STAR format for Relion Tomography
    - JSON format for EMAN Tomography
    - TBL format for Dynamo
    - CBOX format for SPHIRE

    This dynamic behavior ensures that the export options reflect the actual software
    environment available to the user.

    TXT Export

    In the simplest export mode, the protocol writes plain text coordinate files.

    For each tomogram, a separate file is created containing one coordinate per line:

    - X position
    - Y position
    - Z position

    Coordinates are written using the bottom-left corner convention.

    This format is useful for generic external analysis, quick inspection, scripting,
    or importing into software that accepts simple coordinate lists.

    STAR Export

    When the Relion Tomography plugin is available, the protocol can export coordinates
    into STAR format.

    In this mode:

    - all coordinates are written into a STAR file;
    - tomogram identifiers are preserved;
    - sampling rate information is also included.

    This export is particularly useful when transferring particles from Scipion into
    Relion-based subtomogram averaging workflows.

    EMAN Export

    When the EMAN tomography plugin is installed, the protocol exports coordinates into
    JSON metadata files.

    One JSON file is created for each tomogram.

    This mode is especially useful for workflows involving particle picking,
    subtomogram extraction, or visualization inside EMAN-based environments.

    Dynamo Export

    When the Dynamo plugin is available, the protocol exports coordinates in Dynamo
    table (`.tbl`) format.

    For each particle, the protocol writes:

    - particle identifier,
    - placeholder alignment parameters,
    - particle position.

    At present, alignment parameters are exported as default zero values.

    This means the protocol mainly transfers positional information rather than
    refined particle orientations.

    From a practical perspective, this export is useful when coordinates need to be
    imported into Dynamo for subsequent alignment or subtomogram averaging.

    SPHIRE Export

    When the SPHIRE plugin is installed, the protocol can export coordinates in
    CBOX format.

    This provides compatibility with SPHIRE tomography workflows and allows coordinate
    transfer without additional manual conversion.

    Coordinate Organization by Tomogram

    One important internal feature of the protocol is that coordinates are processed
    ordered by tomogram identifier.

    Whenever the tomogram changes:

    - the current output file is closed;
    - a new tomogram-specific file is created.

    This guarantees that exported coordinates remain naturally grouped according to
    their source tomogram.

    From a biological standpoint, this is essential because particles from different
    tomograms often correspond to different specimens, acquisition conditions, or
    biological states.

    Output Management

    All exported files are written into a dedicated `Export` directory created inside
    the protocol working folder.

    Before writing new results, any previous export directory is removed and recreated.

    This ensures that:

    - old files do not remain mixed with new results;
    - the exported output reflects only the most recent execution.

    After completion, the protocol summary reports the absolute path to the export
    directory so the generated files can be easily located.

    Practical Recommendations

    In routine cryo-electron tomography workflows, this protocol is particularly useful
    when particle coordinates need to be transferred between software environments.

    A few practical considerations are important:

    - Use TXT format for generic coordinate exchange or custom scripting.
    - Use STAR format when preparing downstream processing in Relion Tomography.
    - Use EMAN JSON export when continuing particle analysis in EMAN.
    - Use Dynamo TBL export when coordinates will be refined or averaged in Dynamo.
    - Use CBOX export for SPHIRE-compatible workflows.

    Since coordinates are exported per tomogram, users should verify that tomogram
    identifiers are consistent before exporting large coordinate sets.

    Final Perspective

    The Export 3D Coordinates protocol is a lightweight but highly practical utility
    for tomography workflows.

    Rather than modifying coordinates, it preserves particle spatial information and
    makes it portable across multiple subtomogram analysis platforms.

    In modern cryo-ET pipelines, this protocol plays an important interoperability role,
    allowing coordinate information generated inside Scipion to remain useful throughout
    the broader ecosystem of tomography software.
    """

    _label = 'export 3D coordinates'
    _devStatus = NEW

    @staticmethod
    def _getExportChoices():
        """ Return a list of possible choices for export.
        """
        exportChoices = [EXPORT_TO_TXT]
        if existsPlugin('reliontomo'):
            exportChoices.append(EXPORT_TO_STAR)
        if existsPlugin('emantomo'):
            exportChoices.append(EXPORT_TO_EMAN)
        if existsPlugin('dynamo'):
            exportChoices.append(EXPORT_TO_DYNAMO)
        if existsPlugin('sphire'):
            exportChoices.append(EXPORT_TO_CBOX)

        return exportChoices

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        exportChoices = self._getExportChoices()

        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      important=True,
                      label="Input 3D coordinates")
        form.addParam('outputFormat', params.EnumParam,
                      choices=exportChoices, default=0,
                      label='Export to',
                      help='Select the output format')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        objId = self.inputCoordinates.get().getObjId()
        format = self._getExportChoices()[self.outputFormat.get()]
        self._insertFunctionStep(self.exportCoordsStep, format, objId)

    # --------------------------- STEPS functions -----------------------------
    def exportCoordsStep(self, format, coordsId):
        inputCoords = self.inputCoordinates.get()
        tomoIds = inputCoords.getTSIds()
        pwutils.cleanPath(self._getExportPath())
        pwutils.makePath(self._getExportPath())

        if format == EXPORT_TO_TXT:
            def _writeFunc(coord, f):
                x, y, z = map(int, coord.getPosition(BOTTOM_LEFT_CORNER))
                f.write(f"{x} {y} {z}\n")

            self._writeTxt(inputCoords, "txt", _writeFunc)

        elif format == EXPORT_TO_STAR:
            from reliontomo.convert import writeSetOfCoordinates
            writeSetOfCoordinates(inputCoords,
                                  self._getExportPath("coords.star"),
                                  tomoIds,
                                  sRate=inputCoords.getSamplingRate(),
                                  coordsScale=1)

        elif format == EXPORT_TO_EMAN:
            from emantomo.convert import setCoords3D2Jsons
            json_files = [self._getExportPath(f"{tsId}_info.json") for tsId in tomoIds]
            setCoords3D2Jsons(json_files, inputCoords)

        elif format == EXPORT_TO_DYNAMO:
            from dynamo.convert import matrix2eulerAngles

            def _writeFunc(coord, f):
                x, y, z = coord.getPosition(BOTTOM_LEFT_CORNER)
                # Get alignment information
                #if coord.hasTransform():  # FIXME
                #    tdrot, tilt, narot, shiftx, shifty, shiftz = matrix2eulerAngles(coord.getMatrix())
                #else:
                tdrot, tilt, narot, shiftx, shifty, shiftz = 0, 0, 0, 0, 0, 0
                f.write(f"{coord.getObjId()} 1 1 {shiftx} {shifty} {shiftz} "
                        f"{tdrot} {tilt} {narot} 0 0 0 1 0 0 0 0 0 0 0 0 1 0 "
                        f"{x} {y} {z} 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")

            self._writeTxt(inputCoords, "tbl", _writeFunc)

        elif format == EXPORT_TO_CBOX:
            from sphire.convert import writeSetOfCoordinates3D
            writeSetOfCoordinates3D(self._getExportPath(), inputCoords)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        return validateMsgs

    def _summary(self):
        summary = []

        summary.append(f"Output is written to: \n"
                       f"{os.path.abspath(self._getExportPath())}\n")

        return summary
    
    # --------------------------- UTILS functions -----------------------------
    def _getExportPath(self, *paths):
        return os.path.join(self._getPath('Export'), *paths)

    def _writeTxt(self, inputCoords, ext="txt", writeCoord=None):
        """ Iterate over coords by tomoId and write output. """
        f = None
        lastTomoId = None
        for coord in inputCoords.iterCoordinates(orderBy="_tomoId"):
            tomoId = coord.getTomoId()
            if tomoId != lastTomoId:
                if f:  # we need to close previous opened file
                    f.close()
                f = open(self._getExportPath(f"{tomoId}.{ext}"), "w")
                lastTomoId = tomoId
            writeCoord(coord, f)
        if f:
            f.close()
