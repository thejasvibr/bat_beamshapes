.. beamshapes documentation master file, created by
   sphinx-quickstart on Tue Apr 20 19:50:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Beamshapes: generate beamshapes for multiple models
===================================================

Introduction
~~~~~~~~~~~~
The `beamshapes` (currently called `bat_beamshapes`) package implements 
non-trivial sound radiation models (piston in a sphere, piston in a closed circular baffle, etc.)
. Users can generate the beamshape pattern of their model of interest from any of 
these models. Since the beamshape calculation can be computationally intensive, 
the package also comes with a set of pre-calculated examples. 

Who is this useful for? 
~~~~~~~~~~~~~~~~~~~~~~~
Acousticians and bio-acousticians looking to assess model-fits or perform 
parameter estimation on microphone array data. 


.. toctree::
      :maxdepth: 4
      :caption: Examples

      gallery_examples/index.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
