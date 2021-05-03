.. beamshapes documentation master file, created by
   sphinx-quickstart on Tue Apr 20 19:50:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Beamshapes: generate beamshapes for multiple models
===================================================

Introduction
~~~~~~~~~~~~
The `beamshapes` (currently called `bat_beamshapes`) package implements 
non-classical sound radiation models (piston in a sphere, cap in a sphere, 
piston in a closed circular baffle, etc.).

While there are many sound emission models in the literature, there aren't
that many (AFAIK, any?) openly available computational implementations. `bat_beamshapes` 
covers this gap by providing a set of models where you can generate predictions
by providing the relevant parameters through a function call.



Who is this useful for? 
~~~~~~~~~~~~~~~~~~~~~~~
Acousticians and bio-acousticians looking to assess model-fits or perform 
parameter estimation on their sound sources. Check out more on the how
to use this package in the examples page. 


Package under active development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This package is relatively young, and changing dynamically. Things may
break fairly often if you are installing from the GitHub repo directly.
A stable release will hopefully be out by the end of May 2021 or 
beginning of June 2021. Stay tuned!!

.. toctree::
      :maxdepth: 4
      :caption: Examples

      gallery_examples/index.rst


Contributors
~~~~~~~~~~~~
Thejasvi Beleyur (maintainer)
Gaurav Dhariwal

Acknowledgements
~~~~~~~~~~~~~~~~
Many thanks to Tim Mellow for sharing Mathematica code to help with porting to Python.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
