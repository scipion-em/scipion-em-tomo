# scipion-em-tomo

Base Scipion plugin for electron cryo-tomography and subtomogram averaging.

## Current development setup

This plugin is currently under beta-testing and is not ready for production yet.
It can be installed by cloning this repository and executing the command:

.. code-block::

    git clone -b devel https://github.com/scipion-em/scipion-em-tomo.git
    scipion installp -p ./scipion-em-tomo --devel

or by adding the path to the $PYTHONPATH environment variable:

.. code-block::

    export PYTHONPATH=$PYTHONPATH:$DEVEL_HOME/scipion-em-tomo

To check the installation, simply run one of the following tests:

.. code-block::
   
    scipion3 tests tomo.tests.test_transformations.TestTomoTransformations
    scipion3 tests tomo.tests.test_transformations.TestTiltSeriesTransformations
    scipion3 tests tomo.tests.test_transformations.TestSubtomoTransformations
    scipion3 tests tomo.tests.test_tomo_base.TestTomoSubSetsTs
    scipion3 tests tomo.tests.test_tomo_base.TestTomoSubSetsTomograms
    scipion3 tests tomo.tests.test_tomo_base.TestTomoSubSetsSubTomograms
    scipion3 tests tomo.tests.test_tomo_base.TestTomoSplitEvenOdd
    scipion3 tests tomo.tests.test_tomo_base.TestTomoPreprocessing
    scipion3 tests tomo.tests.test_tomo_base.TestTomoCoordinatesOrigin
    scipion3 tests tomo.tests.test_tomo_base.TestTomoBase
    scipion3 tests tomo.tests.test_tomo_base.TestTomoAssignTomo2Subtomo
    scipion3 tests tomo.tests.test_tomo_base.TestTomoAssignAlignment
    scipion3 tests tomo.tests.test_tomo_base.TestParticlesToSubtomograms
    scipion3 tests tomo.tests.test_objects.TestTomoModel
    scipion3 tests tomo.tests.test_import.TestTomoImportTsFromPattern
    scipion3 tests tomo.tests.test_import.TestTomoImportTsFromMdoc
    scipion3 tests tomo.tests.test_import.TestTomoImportTomograms
    scipion3 tests tomo.tests.test_import.TestTomoImportSubTomograms
    scipion3 tests tomo.tests.test_import.TestTomoImportSetOfCoordinates3D
    scipion3 tests tomo.tests.test_import.TestTomoBaseProtocols
    scipion3 tests tomo.tests.test_import.TestImportTomoMasks
    scipion3 tests tomo.tests.test_compose_TS.TestTomoComposeTS


A complete list of tests can also be seen by executing ``scipion tests --grep tomo.tests``

### Prerequisites

* Use scipion v3+
* Use scipion-em-tomo from branch 'devel'
* TestTomoPreprocessing test requires cistem, motioncorr and imod plugins to work properly. If these plugins cannot be found the test will not be executed. For this test to succeed one available GPU is needed, more GPUs can be defined by the var SCIPION_TEST_GPULIST (e.g export SCIPION_TEST_GPULIST='0 1')

### Protocols

* 2D particles to subtomograms
* 2d coordinates to 3d coordinates
* Compose Tilt Serie
* Tilt-series assign alignment
* Tilt-series consensus alignment
* Tilt-series convert coords3D
* assign alignment
* assign tomograms to tomo masks (segmentations)
* assign tomos to subtomos
* average tilt-series movies
* consensus classes subtomo
* ctf validate
* extract 3D coordinates
* import 3D coordinates from scipion
* import set of coordinates 3D
* import subtomograms
* import tilt-series
* import tilt-series movies
* import tomograms
* import tomomasks (segmentations)
* split even/odd tomos/subtomos
* tomograms to micrographs
