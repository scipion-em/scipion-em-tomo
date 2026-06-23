===========
Tomo plugin
===========

.. image:: https://img.shields.io/pypi/v/scipion-em-tomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-tomo
        :alt: PyPI release

Base Scipion plugin for electron cryo-tomography and subtomogram averaging.

Protocols
---------

* **import tilt-series**: Protocol to import tilt series.
* **import tilt-series movies**: Protocol to import tilt series movies.
* **import tomograms** : Protocol to import a set of tomograms to the project.
* **import tomomasks (segmentations)**: Protocol to import a set of tomomasks (segmentations) to the project.
* **import 3D coordinates from scipion** : Protocol to import a set of 3D coordinates from a Scipion sqlite file.
* **import coordinates 3D** : Protocol to import a set of coordinates 3D to the project.
* **import subtomograms** : Protocol to import a set of subtomograms to the project.
* **import tomo CTFs** : Protocol to import the CTF estimation of a set of tilt-series.
* **Compose Tilt Series** : Compose in streaming a set of tilt series based on sets of micrographs and mdoc files. Three time parameters are available for the streaming behaviour: time for next tilt, time for next micrograph and time for next tilt series.
* **correct tilt offset** : Correct the tilt angle offset of a tilt-series (particularly useful for lamellae) by adding a fixed offset to every tilt angle in the metadata.
* **invert tilt angles** : Inverts the physical handedness of the introduced tilt-series by inverting the tilt angles in the metadata associated to each tilt-series.
* **exclude views filter** : Filter a set of aligned tilt-series according to maximum allowed shift after alignment, accumulated dose range, tilt-angle range and minimum number of views.
* **assign excluded views** : Transfer excluded-view annotations (the _enabled state) from a source set of tilt-series to a target one. Tilt-series are matched by tsId and images by acquisition order, so it is robust to re-stacking.
* **Tilt-series assign alignment** : Assign the transformation matrices from an input set of tilt-series to a target one.
* **tilt-series from tomograms** : Get the tilt-series that correspond to a given set of tomograms. Useful for fiducial-less samples, to discard data at tomogram level and keep further processing with the matching tilt-series.
* **misalign tilt-series** : Introduce misalignment in the transformation matrix of a tilt-series (mainly for testing purposes).
* **apply tomomasks to tomograms** : Applies a set of masks to a given set of tomograms. Some operations can be applied to the mask: invert, dilate and apply a gaussian filter.
* **extract 3D coordinates** : Extract the coordinates information from a set of subtomograms. This protocol is useful when we want to re-extract the subtomograms (maybe resulting from classification) with the original dimensions. It can be also handy to visualize the resulting subtomograms in their location on the tomograms.
* **project coordinates** : Project 3D coordinates into a set of landmarks.
* **mask 3d coordinates** : Filter a set of 3D coordinates using a segmentation (tomomask), keeping or excluding the coordinates that fall within a given segmentation label.
* **meshes from tomo mask** : Create meshes based on segmentations or voxel values (TomoMasks).
* **2D particles to subtomograms** : Protocol to create a set of subtomograms from selected 2D particles.
* **split even/odd tomos/subtomos** : Protocol to split a set of tomograms or subtomograms into even/odd sets by element id.
* **export 3D coordinates**: Export 3D subtomogram coordinates to be used outside Scipion.


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

**NAPARI_ENV_ACTIVATION** (default = conda activate napari-0.4.19):
Command to activate napari environment (used by other tomo plugins).


Buildbot status
---------------

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/tomo_devel.svg


Status production version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/tomo_prod.svg
