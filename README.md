# scipion-em-tomo

Base Scipion plugin for Cryo-Tomography  and Subtomogram averaging.

## Current development setup
This plugin is currently under initial development and it is not ready for production yet. 

In the meantime, it can be used for development, base on Scipion v2.x with plugins. 
 
This tomography plugin can be enabled by cloning this repository and execute the command: 

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-tomo.git
    scipion installp -p ~/scipion-em-tomo --devel


or by adding the path to the $PYTHONPATH environment variable. 

.. code-block::

    # Enable the tomography plugin
    export PYTHONPATH=$PYTHONPATH:$DEVEL_HOME/scipion-em-tomo


To check the installation, simply run one of the following Scipion tests:

.. code-block::
   
       scipion test tomo.tests.test_tomo_base.TestTomoPreprocessing
       scipion test tomo.tests.test_tomo_base.TestTomoImportTs
       scipion test tomo.tests.test_tomo_base.TestTomoImportTomograms
       scipion test tomo.tests.test_tomo_base.TestTomoImportSubTomograms
       scipion test tomo.tests.test_tomo_base.TestTomoImportSetOfCoordinates3D
       scipion test tomo.tests.test_tomo_base.TestTomoBaseProtocols
       scipion test tomo.tests.test_tomo_base.TestTomoBase


 A complete list of tests can also be seen by executing ``scipion test tomo.tests.test_tomo_base``

### Pre-requisites
* Use scipion v3.0.0
* Use scipion-em-tomo from branch 'tomo' (see previous section for a possible setup)
* TestTomoPreprocessing test requires gctf and imod plugins to work properly. If these plugins cannot be found the test will not be executed. For this test to succeed one available GPU is needed, more GPUs can be defined defining the var SCIPION_TEST_GPULIST (e.g export SCIPION_TEST_GPULIST='0 1' )

This test will launch automatically the pre-processing following pipeline:
* import tilt-series (only one)
* run motioncor2
* run gctf
* run imod_auto3d (script ported from J.chichón)

### Protocols

* Import Titl Series
* Import Titl Series Movies
* Import Tomograms
* Import Subtomograms
* Import Coordinates 3D
* Tilt Series Averaging







