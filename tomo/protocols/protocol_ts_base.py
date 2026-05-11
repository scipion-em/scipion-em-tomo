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


from pwem.protocols import EMProtocol
from pyworkflow.protocol import STEPS_PARALLEL, STATUS_NEW

from ..objects import TiltSeriesDict, Tomogram
from .protocol_base import ProtTomoBase


class ProtTsProcess(EMProtocol, ProtTomoBase):
    """
    Base class for Tilt-Series (images or movies) processing protocols.
    This class should be used by protocols that receive tilt-series as input
    and produce tilt-series as output. This class will contain some common
    functionality about the steps execution (also in streaming) and
    the output generation.
    """
    """
    Provides a generic framework for tilt-series processing workflows in cryo-electron
    tomography. This base protocol is designed for methods that receive tilt-series
    as input and generate processed tilt-series or tomographic outputs while supporting
    both parallel and streaming execution.

    AI Generated:

    Tilt-Series Processing Framework (ProtTsProcess & ProtTomoReconstruct) — User Manual
        Overview

        The ProtTsProcess class defines the core infrastructure used for processing
        tilt-series datasets inside tomography workflows. Rather than implementing
        a specific reconstruction or alignment algorithm, this protocol provides
        the execution logic, streaming management, and output handling required
        by tomography processing methods operating on sequential tilt-series data.

        In cryo-electron tomography, tilt-series datasets are often acquired continuously
        and may contain hundreds or thousands of images distributed across multiple
        tomographic series. Efficient processing therefore requires a framework capable
        of dynamically detecting new datasets, launching parallel computations, and
        updating outputs incrementally while acquisition or preprocessing is still
        ongoing.

        The protocol addresses these requirements by organizing the workflow into
        modular execution steps that can process individual tilt images, complete
        tilt-series, or final tomographic reconstructions depending on the needs
        of the derived protocol.

        Inputs and General Workflow

        The framework operates on sets of tilt-series and automatically manages
        their conversion, processing, and output generation. During execution,
        the protocol initializes a dynamic dictionary of tilt-series objects and
        continuously monitors the input dataset for newly available items.

        For each detected tilt-series, the workflow may optionally launch independent
        processing steps for every tilt image before executing a global tilt-series
        processing stage. This hierarchical organization allows developers to implement
        highly parallel tomography pipelines while maintaining synchronization between
        intermediate and final results.

        The framework is particularly suited for streaming environments in which
        tilt-series are progressively generated during microscope acquisition or
        preprocessing stages. As new data become available, the protocol automatically
        inserts the corresponding execution steps and updates the output set without
        requiring manual intervention.

        Parallelization and Streaming Execution

        One of the most important characteristics of this framework is its support
        for parallel and streaming execution. Cryo-ET datasets are computationally
        expensive, especially when processing large numbers of tilt-series or high-resolution
        images. The protocol therefore distributes processing tasks across independent
        execution steps to maximize throughput and scalability.

        The design separates image-level operations from tilt-series-level operations.
        This allows lightweight preprocessing operations to be applied independently
        to each tilt image while preserving a later global reconstruction or refinement
        stage operating on the complete series.

        Streaming support is biologically and operationally important in facility-scale
        cryo-ET pipelines because it enables reconstruction and quality assessment
        to begin before acquisition has fully completed. Early feedback can help
        identify acquisition problems, specimen instability, or alignment artifacts
        before significant microscope time is lost.

        Output Generation and Incremental Updates

        The protocol dynamically maintains output datasets while processing progresses.
        Completed tilt-series are appended incrementally into the output set together
        with their associated tilt images. This strategy minimizes memory overhead
        while allowing downstream protocols to access partially completed results
        during streaming execution.

        The framework automatically manages output creation, append operations,
        metadata propagation, and stream state transitions. Once all tilt-series
        have been processed, the protocol closes the output stream and finalizes
        the dataset.

        From a workflow perspective, this behavior enables seamless integration
        with larger Scipion tomography pipelines in which downstream steps may
        begin consuming reconstructed data immediately after they become available.

        Tomogram Reconstruction Extension

        The ProtTomoReconstruct class extends the generic tilt-series processing
        framework specifically for tomographic reconstruction workflows. Instead
        of producing processed tilt-series, reconstruction protocols derived from
        this class generate tomograms as final outputs.

        Unlike generic tilt-series processing pipelines, reconstruction protocols
        usually operate on the complete tilt-series as a single computational unit.
        For this reason, image-by-image processing steps are disabled by default,
        simplifying execution and reducing scheduling overhead.

        For every processed tilt-series, the protocol generates a corresponding
        tomographic volume stored as an MRC file. The output tomograms inherit
        the sampling information from the original tilt-series, optionally adjusted
        according to the reconstruction binning factor.

        Biologically, the resulting tomograms represent volumetric reconstructions
        of the original specimen and serve as the starting point for downstream
        analyses such as subtomogram averaging, particle picking, segmentation,
        membrane tracing, or spatial organization studies.

        Reconstruction Workflow Design

        The reconstruction-oriented specialization provided by ProtTomoReconstruct
        allows developers to focus exclusively on the mathematical reconstruction
        algorithm while delegating execution management, output synchronization,
        and streaming behavior to the parent framework.

        This separation of responsibilities is particularly valuable for tomography
        software development because reconstruction algorithms may vary widely
        in complexity and computational cost. Weighted backprojection, SIRT, ART,
        compressed sensing, and deep-learning-based reconstruction methods can
        all reuse the same workflow infrastructure while implementing only their
        specific numerical reconstruction logic.

        Biological and Practical Perspective

        From a biological perspective, the framework represented by ProtTsProcess
        and ProtTomoReconstruct forms the operational backbone of many cryo-electron
        tomography pipelines. Although these classes do not directly implement
        biological analysis algorithms, they enable efficient handling of the
        large and complex datasets required for modern structural cell biology.

        Their streaming-oriented architecture is especially important for large-scale
        cryo-ET facilities and automated acquisition environments where throughput,
        robustness, and scalability are critical. By separating tilt-image processing,
        tilt-series management, and tomogram reconstruction into modular execution
        stages, the framework supports flexible and extensible workflows adaptable
        to a broad range of tomography applications.

        Final Perspective

        The ProtTsProcess and ProtTomoReconstruct classes provide a generalized,
        scalable infrastructure for cryo-electron tomography processing pipelines.
        Through dynamic step scheduling, streaming execution, and incremental
        output management, they enable efficient processing of large tilt-series
        datasets while serving as the foundation for advanced tomographic reconstruction
        protocols and downstream structural analysis workflows.
    """