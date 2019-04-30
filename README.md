# scipion-em-tomo

Base Scipion plugin for Cryo-Tomography  and Subtomogram averaging.

## Current development setup
This plugin is currently under initial development and it is not ready for production yet. So, there are missing some parts to make it a proper plugin and installable via pip. 

In the meantime, it can be used for development, base on Scipion v2.x with plugins. I have been using the following .bashrc file to load some environment variables. 

```
export DEVEL_HOME=/home/josem/work/development
export SCIPION_HOME=$DEVEL_HOME/scipion-devel

alias scipion='$SCIPION_HOME/scipion --config $SCIPION_HOME/config/scipion-plugins.conf'
alias plugins='scipion python $DEVEL_HOME/scipion-scripts/plugins.py'

. $DEVEL_HOME/scipion-em-plugins/bash-plugins.sh
. $DEVEL_HOME/xmipp-tools/build/xmipp.bashrc

# Enable the tomography plugin
export PYTHONPATH=$PYTHONPATH:$DEVEL_HOME/scipion-em-tomo

# Load IMOD
. /home/josem/installs/imod/imod.bashrc

# Specify test data folder
export SCIPION_TOMO_EMPIAR10164=/data/work_data/empiar-10164
```

The first lines are to load Scipion 2.0 and related plugins in a development environment. 
This tomography plugin should be enabled by adding the path to the $PYTHONPATH environment variable. 

Moreover, for initial tests I have setup the IMOD variable and also the SCIPION_TOMO_EMPIAR10164, as test data. 
This variable should point to empiar-10164 folder containing 'data/frames' subfolder with Tilt Series 01, 03 and 07 from this EMPIAR dataset: 

http://www.ebi.ac.uk/pdbe/emdb/empiar/entry/10164/

## Test of pre-processing protocols

There is a test that will launch automatically the pre-processing following pipeline:
* import tilt-series (only one)
* run motioncor2
* run gctf
* run imod_auto3d (script ported from J.chichón)

### Pre-requisites
* Use scipion v2.0.0 from branch 'tomo'
* Use scipion-em-tomo from branch 'tomo' (see previous section for a possible setup)
* Install IMOD and source it (this defines the IMOD_DIR variable that is required)
* Install 'tomo3d' program from J.J. Fernández (it should be placed into a folder with 'bin', and define variable TOMO3D_DIR pointing to install folder. This may change if we group all JJF programs into a 'package', that was suggested by him)
* This test requires at least one available GPU, more GPUs can be defined defining the var SCIPION_TEST_GPULIST (e.g export SCIPION_TEST_GPULIST='0 1' )

### Launching the test

`scipion test tomo.tests.test_tomo_base.TestTomoPreprocessing` 


