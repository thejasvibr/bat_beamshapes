.. beamshapes documentation master file, created by
   sphinx-quickstart on Tue Apr 20 19:50:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Beamshapes: directivity patterns for multiple sound sources
===========================================================

Introduction
~~~~~~~~~~~~
The `beamshapes` package implements the 
directivity functions (:math:`\frac{D_{\theta}}{D_{0}}`) of various sound radiation models
(piston in a sphere, cap in a sphere, piston in a closed circular baffle, etc.).


Check out a general introduction to the concepts of the package :doc:`here <general_intro>`.

What gap does this package fill?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
While there are many sound radiation models described in the literature, there aren't
that many (or any?) openly available computational implementations of their beamshapes.

Who is this useful for? 
~~~~~~~~~~~~~~~~~~~~~~~
Acousticians and bio-acousticians looking to assess model-fits or perform 
parameter estimation on their sound sources. Check out more on the how
to use this package in the :doc:`examples <gallery_examples/index>`. 


Models implemented
~~~~~~~~~~~~~~~~~~
* Point source on a sphere
* Piston in an infinite baffle 
* Piston in a sphere
* Oscillating cap of a sphere 

Models under implementation 
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Piston in a finite closed baffle  

Models to be implemented
~~~~~~~~~~~~~~~~~~~~~~~~
* Rectangular cap of a sphere 
* Any others of potential relevance -- suggest!


Package under active development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This package is relatively young, and changing dynamically. Things may
break fairly often if you are installing from the `GitHub repo <https://github.com/thejasvibr/beamshapes>`_ directly.
A stable release will hopefully be out by the end of May 2021 or 
beginning of June 2021. Stay tuned!


.. toctree::
      :maxdepth: 3
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



.. toctree::
   :maxdepth:1
   :caption: API reference:

   source_models.rst
   misc.rst


    
    
    
    



