===========
Tomo plugin
===========

Base Scipion plugin for electron cryo-tomography and subtomogram averaging.

This plugin is currently in **BETA** mode.

Protocols
---------

* **import tilt-series**: Protocol to import tilt series.
* **import tilt-series movies**: Protocol to import tilt series movies.
* **import tomograms** : Protocol to import a set of tomograms to the project.
* **import tomomasks (segmentations)**: Protocol to import a set of tomomasks (segmentations) to the project.
* **import 3D coordinates from scipion** : Protocol to import a set of 3d coordinates from Scipion sqlite file.
* **import set of coordinates 3D** : Protocol to import a set of coordinates 3D to the project.
* **import subtomograms** : Protocol to import a set of tomograms to the project.
* **2D particles to subtomograms** : Protocol to create a set of subtomograms from a selected 2D particles.
* **2d coordinates to 3d coordinates** : Turns 2d coordinates into set of 3d coordinates. Works in coordination with 'tomograms to micrographs' protocol.
* **Compose Tilt Serie** : Compose in streaming a set of tilt series based on a sets of micrographs and mdoc files. Three time parameters are abailable for the streaming behaviour: Time for next tilt, Time for next micrograph and Time for next Tilt Serie
* **Tilt-series assign alignment** : Assign the transformation matrices from an input set of tilt-series to a target one.
* **Tilt-series consensus alignment** : Perform a consensus of a set of alignments for the same tilt series. Returns the average alignment matrix of the consensus alignments and its standard deviation of shift and angle.
* **Tilt-series convert coords3D** : Scipion protocol to convert a set of tilt-series coordinates 3d to a set of coordinates 3d associated to a set of tomograms.
* **assign alignment** : Assign the alignment stored in a set of Subtomograms/Coordinates3D to another set. Both sets should have same pixel size (A/px). The Subtomograms/Coordinates3D with the alignment can also be a subset of a bigger set.
* **assign tomograms to tomo masks (segmentations)** : This protocol assign tomograms to tomomasks (segmentations).
* **assign tomos to subtomos** : This protocol assign tomograms to subtomograms that have been imported before without tomograms. Subtomograms should contain the name of the original tomogram in their own file name.
* **average tilt-series movies** : Simple protocol to average TiltSeries movies as basic  motion correction. It is used mainly for testing purposes.
* **consensus classes subtomo** : Compare several SetOfClassesSubTomograms. Return the intersection of the input classes.
* **ctf validate** : Validate a set of CTF tomo series and separate into two sets (good and bad tomo series )
* **extract 3D coordinates** : Extract the coordinates information from a set of subtomograms. This protocol is useful when we want to re-extract the subtomograms (maybe resulting from classification) with the original dimensions. It can be also handy to visualize the resulting subtomograms in their location on the tomograms.
* **split even/odd tomos/subtomos** : Protocol to split set of tomograms or subtomograms in even/odd sets by element id.
* **tomograms to micrographs** : Turns tomograms into set of micrographs to apply SPA picking methods.


**Latest plugin versions**
==========================

If you want to check the latest version and release history go to `CHANGES <https://github.com/scipion-em/scipion-em-tomo/blob/master/CHANGES.txt>`_


**Installing the plugin**
=========================

In order to install the plugin follow these instructions:

.. code-block::

    scipion installp -p scipion-em-tomo


or through the **plugin manager** by launching Scipion and following **Configuration** >> **Plugins**


**To install in development mode**

Clone or download the plugin repository

.. code-block::

    git clone -b devel https://github.com/scipion-em/scipion-em-tomo.git

Install the plugin in developer mode.

.. code-block::

    scipion installp -p local/path/to/scipion-em-tomo --devel


**Configuration variables**
===========================

**NAPARI_ENV_ACTIVATION** (default = conda activate napari-0.4.17):
Command to activate napari environment (used by other tomo plugins).


Buildbot status
---------------

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/tomo_devel.svg


Status production version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/tomo_prod.svg
