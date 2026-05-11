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
import re
from glob import glob
from os.path import join, getmtime

import pyworkflow as pw
from pyworkflow.protocol.params import (PointerParam, EnumParam, PathParam,
                                        FloatParam, StringParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.mapper.sqlite_db import SqliteDb
from pyworkflow.utils.properties import Message
from pwem.protocols import ProtImport, EMProtocol, ProtImportFiles

import tomo.objects
from tomo.constants import TS_LABEL
from tomo.convert.mdoc import normalizeTSId

logger = logging.getLogger(__name__)


class ProtTomoBase:
    def _createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be deleted. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        pw.utils.cleanPath(setFn)

        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj

    def _createSetOfTiltSeriesM(self, suffix='') -> tomo.objects.SetOfTiltSeriesM:
        return self._createSet(tomo.objects.SetOfTiltSeriesM,
                               'tiltseriesM%s.sqlite', suffix)

    def _createSetOfTiltSeries(self, suffix='') -> tomo.objects.SetOfTiltSeries:
        self._ouputSuffix = ''
        return self._createSet(tomo.objects.SetOfTiltSeries,
                               'tiltseries%s.sqlite', suffix)

    def _createSetOfCoordinates3D(self, volSet, suffix='') -> tomo.objects.SetOfCoordinates3D:
        coord3DSet = self._createSet(tomo.objects.SetOfCoordinates3D,
                                     'coordinates%s.sqlite', suffix,
                                     indexes=['_volId'])
        coord3DSet.setPrecedents(volSet)
        return coord3DSet

    def _createSetOfTomograms(self, suffix='') -> tomo.objects.SetOfTomograms:
        return self._createSet(tomo.objects.SetOfTomograms,
                               'tomograms%s.sqlite', suffix)

    def _createSetOfSubTomograms(self, suffix='') -> tomo.objects.SetOfSubTomograms:
        return self._createSet(tomo.objects.SetOfSubTomograms,
                               'subtomograms%s.sqlite', suffix)

    def _createSetOfAverageSubTomograms(self, suffix='') -> tomo.objects.SetOfAverageSubTomograms:
        return self._createSet(tomo.objects.SetOfAverageSubTomograms,
                               'avgSubtomograms%s.sqlite', suffix)

    def _createSetOfClassesSubTomograms(self, subTomograms, suffix='') -> tomo.objects.SetOfClassesSubTomograms:
        classes = self._createSet(tomo.objects.SetOfClassesSubTomograms,
                                  'subtomogramClasses%s.sqlite', suffix)
        classes.setImages(subTomograms)

        return classes

    def _createSetOfLandmarkModels(self, suffix='') -> tomo.objects.SetOfLandmarkModels:
        return self._createSet(tomo.objects.SetOfLandmarkModels, 'setOfLandmarks%s.sqlite', suffix)

    def _createSetOfMeshes(self, volSet, suffix='') -> tomo.objects.SetOfMeshes:
        meshSet = self._createSet(tomo.objects.SetOfMeshes,
                                  'meshes%s.sqlite', suffix)
        meshSet.setPrecedents(volSet)
        return meshSet

    def _getOutputSuffix(self, cls):
        """ Get the name to be used for a new output.
        For example: output3DCoordinates7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(cls):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1  # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter + 1) if maxCounter > 0 else ''  # empty if not output


class ProtTomoPicking(ProtImport, ProtTomoBase):
    OUTPUT_PREFIX = 'output3DCoordinates'

    """ Base class for Tomogram boxing protocols. """

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputTomograms', PointerParam, label="Input Tomograms", important=True,
                      pointerClass='SetOfTomograms',
                      help='Select the Tomogram to be used during picking.')

    def _summary(self):
        summary = []
        if self.isFinished() and self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:*\n%s" % (key, output.getSummary()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary


class ProtTomoImportFiles(ProtImportFiles, ProtTomoBase):
    """Base protocol to import tomography files. There are two modes that should be implemented to import files,
    based on the value of the filesPattern form parameter:

    1. If empty (a single file is imported) or using the classic wildcard patterns.
    2. If the pattern contains the label {TS}, to represent the part of the name desired to be considered as the tsId.

    How to implement an import protocol in Scipion tomo:

    * The method initializeParsing must be called before the steps are generated, directly or as part of an
    initialization method. This method manages all the functionality required to deal with the {TS} pattern. If this
    label is not present, it does nothing, and the files matching must be carried out with the method iterFiles.

    * When processing the data, an if statement must be implemented, asking for a pattern {TS} introduced or not (it's
    the value of the protocol attribute self.regEx, which is filled properly in the execution of the method
    initializeParsing. Then:

        > if self.regEx: the iterator is obtained by calling the method getMatchingFilesFromRegEx
        > else: the iterator is obtained via the method iterFiles.

    * With the iterators defined, the rest of the code to generate the corresponding scipion objects is common to both
    modes.

    A good example is the protocol ProtImportTomograms.
    """

    """
    Imports tomography-related files into the Scipion Tomo framework while providing
    flexible mechanisms to identify, organize, and associate datasets through file
    patterns and tilt-series identifiers (tsId). The protocol is designed to support
    both conventional file import workflows and advanced tomography-oriented import
    strategies where datasets must be linked consistently across processing stages.

    AI Generated:

    Tomography File Import (ProtTomoImportFiles) — User Manual

        Overview

        The ProtTomoImportFiles protocol provides the foundational infrastructure
        for importing tomography data into Scipion. Its primary purpose is to
        standardize how tomography-related files are discovered, identified, and
        associated with specific tilt-series through the use of file patterns and
        metadata extraction strategies.

        In cryo-electron tomography workflows, maintaining a consistent tilt-series
        identifier is critically important because many downstream objects—such as
        tomograms, alignments, CTF estimations, subtomograms, or transformation
        matrices—must remain biologically and computationally linked throughout the
        processing pipeline. This protocol establishes that relationship at the
        import stage.

        Unlike conventional import systems that simply collect files from a directory,
        this protocol introduces tomography-aware parsing logic capable of identifying
        unique tilt-series identifiers directly from filenames. This allows previously
        generated external results, such as IMOD alignments or reconstructed tomograms,
        to be integrated seamlessly into Scipion projects.

        Inputs and General Workflow

        The protocol supports two different import strategies depending on the value
        of the Pattern parameter. The first strategy corresponds to traditional file
        imports using standard wildcard expressions such as '*', '?', or numeric
        placeholders. This mode behaves similarly to classical Scipion SPA import
        protocols and is appropriate when datasets do not require advanced tomography
        associations.

        The second strategy is specifically designed for tomography workflows and
        relies on the special {TS} label embedded in the filename pattern. This label
        defines which region of the filename should be interpreted as the tilt-series
        identifier. During initialization, the protocol converts this pattern into
        regular expressions and glob-compatible patterns that allow automatic matching
        between imported files and existing tomography objects inside the project.

        From a biological and workflow-management perspective, this capability is
        extremely important because it ensures consistency between datasets generated
        externally and the internal tomography structure maintained by Scipion.

        Tilt-Series Identifier Parsing

        One of the central concepts of this protocol is the extraction and normalization
        of tilt-series identifiers. When the {TS} label is present in the import
        pattern, the protocol dynamically builds a regular expression capable of
        isolating the desired substring from each filename.

        This mechanism allows users to adapt the import process to highly heterogeneous
        naming conventions commonly encountered in cryo-ET facilities and collaborative
        projects. For example, a filename such as:

            TS_01_binned.xf

        can generate the tilt-series identifier "TS_01" automatically when using the
        pattern:

            {TS}_binned.xf

        This flexibility becomes essential when importing external alignments,
        reconstruction files, or metadata generated outside Scipion but intended to
        match already imported tilt-series.

        The protocol also normalizes extracted identifiers to ensure compatibility
        across workflows, minimizing issues caused by inconsistent naming conventions.

        Wildcards and Pattern Matching

        The import system supports both standard wildcard expressions and tomography-
        specific placeholders. Standard wildcards allow users to scan directories and
        identify groups of files efficiently, while the tomography-aware {TS} label
        introduces semantic meaning into the import process.

        This dual approach makes the protocol suitable both for simple exploratory
        workflows and for highly structured production pipelines. In practical terms,
        users can import individual files, entire directories, or complex collections
        of tomography outputs generated by external software suites.

        The protocol internally combines glob-based file discovery with regular
        expression matching, ensuring both flexibility and reproducibility when working
        with large tomography datasets.

        File Exclusion and Dataset Filtering

        The protocol includes an exclusion-word filtering mechanism intended to simplify
        dataset curation during import. Users may provide a list of forbidden words,
        and any file containing those substrings will automatically be ignored.

        In biological workflows this functionality becomes useful when directories
        contain intermediate files, temporary outputs, backup reconstructions, or
        partially processed datasets that should not be incorporated into the project.

        This filtering stage helps maintain clean and reproducible datasets while
        reducing the risk of importing unintended files into downstream analyses.

        Acquisition Metadata Management

        The ProtTomoImportAcquisition class complements the import system by handling
        tomography acquisition parameters. These parameters can either be introduced
        manually or extracted from an external metadata file containing acquisition
        information for each imported object.

        The protocol supports biologically relevant acquisition parameters such as
        angular range, angular step size, tilt-axis orientation, microscope voltage,
        spherical aberration, and amplitude contrast. These parameters are fundamental
        for downstream tomographic reconstruction, CTF correction, and subtomogram
        averaging workflows.

        When acquisition parameters are imported from file, the protocol associates
        metadata entries with imported objects through filename matching. This enables
        large tomography datasets to be imported consistently without requiring manual
        metadata entry for each tomogram or subtomogram.

        From a cryo-ET perspective, preserving accurate acquisition metadata is
        essential because errors in angular geometry or microscope parameters can
        propagate into reconstruction artifacts and compromise biological interpretation.

        Streaming and Workflow Integration

        Although this protocol mainly focuses on import functionality, its architecture
        is designed to integrate naturally with Scipion tomography workflows. Imported
        objects preserve the relationships required by downstream processing protocols,
        allowing alignment files, tomograms, CTF estimations, and subtomogram datasets
        to remain synchronized through shared tilt-series identifiers.

        This interoperability is particularly important in modern cryo-ET pipelines,
        where datasets often move between external software environments and Scipion-
        based workflow managers.

        Subtomogram Averaging Context

        The ProtTomoSubtomogramAveraging base class serves as a structural foundation
        for subtomogram averaging protocols within the tomography framework. Although
        the class itself does not yet implement processing logic, it establishes a
        standardized inheritance structure for future averaging protocols.

        In biological workflows, subtomogram averaging is used to improve the signal-
        to-noise ratio of repetitive macromolecular complexes extracted from tomograms.
        The existence of this shared base class helps maintain consistency across
        tomography averaging implementations.

        Outputs and Their Interpretation

        After import, the protocol generates tomography-compatible Scipion objects
        linked to their associated tilt-series identifiers and acquisition metadata.
        These imported datasets can then be used directly in reconstruction, alignment,
        CTF estimation, particle picking, or subtomogram averaging workflows.

        The correctness of the import stage is biologically significant because all
        downstream associations depend on the consistency of the imported identifiers
        and metadata. Incorrect tsId extraction or mismatched acquisition parameters
        may propagate through the entire tomography workflow.

        Practical Recommendations

        In routine cryo-electron tomography workflows, it is generally advisable to
        adopt a consistent filename convention before importing datasets into Scipion.
        Using clear and reproducible tilt-series identifiers greatly simplifies
        downstream data management and minimizes matching errors between related
        tomography objects.

        When importing files generated externally—such as IMOD alignments or external
        reconstructions—the {TS} pattern mechanism provides the most robust solution
        for preserving dataset consistency. Users should verify that extracted tsIds
        match existing tilt-series identifiers already present in the project.

        For large datasets, exclusion filters can help avoid accidental import of
        temporary or irrelevant files. Similarly, importing acquisition metadata from
        external files is often preferable when handling large tomography collections,
        since it improves reproducibility and reduces manual annotation errors.

        Final Perspective

        In cryo-electron tomography, data import is not simply a technical operation
        but a critical organizational step that defines how all downstream datasets
        remain connected throughout the workflow. The ProtTomoImportFiles framework
        provides the flexibility required for heterogeneous tomography environments
        while ensuring that biological datasets remain consistently associated through
        standardized tilt-series identifiers and acquisition metadata management.
    """