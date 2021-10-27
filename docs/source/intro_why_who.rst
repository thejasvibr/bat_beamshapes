Introduction
============
The :code:`beamshapes` package implements directivity functions (:math:`\frac{D_{\theta}}{D_{0}}`) of various published 
sound radiation models.

What is a directivity function? It's a function that quantifies how sound level changes as you change
the frequency of the emitted sound, and location of the receiver. 

Check out a general introduction to the concepts of the package :doc:`here <general_intro>`.

Why `beamshapes`?
~~~~~~~~~~~~~~~~~
While there are many sound radiation models described in the literature, there aren't
that many (also see `levitate <https://github.com/AppliedAcousticsChalmers/levitate/blob/master/levitate/transducers.py>`_ and `pyroomacoustics <https://pyroomacoustics.readthedocs.io/en/pypi-release/pyroomacoustics.directivities.html>`_) openly available computational implementations of their beamshapes.

Who is this useful for? 
~~~~~~~~~~~~~~~~~~~~~~~
Acousticians and bio-acousticians looking to assess model-fits or perform 
parameter estimation on their sound sources. Check out more on the how
to use this package in the :doc:`examples <gallery_examples/index>`. 


Package installation
~~~~~~~~~~~~~~~~~~~~

`pip installation` : Install the latest stable version with :code:`pip install beamshapes` 

`Local installation` : Install from the `GitHub repo <https://github.com/thejasvibr/bat_beamshapes>`_ directly by cloning and 
following the instructions in the README.


Source models implemented
~~~~~~~~~~~~~~~~~~~~~~~~~
* Point source on a sphere
* Piston in an infinite baffle 
* Oscillating cap of a sphere 
* Piston in a sphere
