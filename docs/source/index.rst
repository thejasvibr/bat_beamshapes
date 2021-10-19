.. beamshapes documentation master file, created by
   sphinx-quickstart on Tue Apr 20 19:50:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Beamshapes: directivity patterns for various sound sources
==========================================================

Introduction
~~~~~~~~~~~~
The `beamshapes` package implements directivity functions (:math:`\frac{D_{\theta}}{D_{0}}`) of various published 
sound radiation models.

What is a directivity function? It's a function that quantifies how sound level changes as you change
the frequency of the emitted sound, and location of the receiver. 

Check out a general introduction to the concepts of the package :doc:`here <general_intro>`.

Why `beamshapes`?
~~~~~~~~~~~~~~~~~
While there are many sound radiation models described in the literature, there aren't
that many (or any?) openly available computational implementations of their beamshapes.

Who is this useful for? 
~~~~~~~~~~~~~~~~~~~~~~~
Acousticians and bio-acousticians looking to assess model-fits or perform 
parameter estimation on their sound sources. Check out more on the how
to use this package in the :doc:`examples <gallery_examples/index>`. 

Source models implemented
~~~~~~~~~~~~~~~~~~~~~~~~~
* Point source on a sphere
* Piston in an infinite baffle 
* Oscillating cap of a sphere 
* Piston in a sphere

Wishlist for future releases 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Piston in a finite closed baffle
* One-sided piston radiator (baffle=cylinder width)
* Rectangular cap of a sphere 
* Piston on a prolate spheroid

Interested?? Contribute, check out the general :doc:`workflow tips here <general_workflow>`.
All of the sound-source models and directivity calculations are from Leo Beranek & Tim Mellow's 
`Acoustics: sound fields and transducers. (Academic Press)` -  check it out for models
that may be of interest or to get an idea of how the wishlist models look like!


Package installation
~~~~~~~~~~~~~~~~~~~~
Install from the `GitHub repo <https://github.com/thejasvibr/bat_beamshapes>`_ directly by cloning and 
following the instructions in the README.


.. toctree::
      :maxdepth: 4
      :caption: Concepts

      ./general_intro.rst

.. toctree::
      :maxdepth: 4
      :caption: Use Cases

      gallery_examples/index.rst
      

Contributors
~~~~~~~~~~~~
Thejasvi Beleyur (maintainer, thejasvib@gmail.com)

Gaurav Dhariwal


Acknowledgements
~~~~~~~~~~~~~~~~

Many thanks to Tim Mellow for sharing Mathematica code to help with porting to Python.
Also thanks to Holger R. Goerlitz for layout feedback (still in progress!) and 
Neetash MR for inspiring the package logo!


.. toctree::
   :maxdepth:1
   :caption: API reference:

   source_models.rst
   misc.rst

.. toctree::
    :maxdepth: 1
    :caption: dev notes

    developer_notes.rst


    
    
    
    



