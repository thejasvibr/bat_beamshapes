.. beamshapes documentation master file, created by
   sphinx-quickstart on Tue Apr 20 19:50:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Beamshapes: directivity patterns for various sound sources
==========================================================


.. include:: ./intro_why_who.rst

.. toctree::
      :maxdepth: 4
      :caption: Concepts

      ./general_intro.rst

.. toctree::
      :maxdepth: 4
      :caption: Use Cases

      gallery_examples/index.rst

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


    
    
    
    



