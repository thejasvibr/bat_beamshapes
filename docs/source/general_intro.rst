Concepts behind :code:`beamshapes`
=============================
Sound sources don't radiate sound uniformly most of the time. Even us humans for instance, we emit more energy while talking to the front than to the back. The exact (non-uniform) pattern in which sound is radiated defines the 'beamshape' of a source. The beamshape is typically a combination of the frequency of emitted sound and the geometry of the vibrating surface and its associated (non-vibrating) surfaces. 

Let's go through the main concepts required to use the `beamshapes` package.

Source models
-------------
The source model refers to the source of sound and its geometric properties, and assumptions of the vibrations etc. According
to the situation in hand, different source models may be physically/biologically relevant! 


There are two main ways to predict how sound will radiate from a source - analytical or numerical modelling. Analytical models
start with equations defining the physics of sound radiation and move on to produce mathematical solutions. Numerical
methods use various numerical algorithms (eg. finite-element method) to computationally simulate sound radiation. :code:`beamshapes`
specifically implements analytical source models with pre-defined solutions. The advantage of using such analytical models is the lower
number of parameters needed to describe and understand the resulting sound radation. 

Below is a brief description of the source models and the parameters relevant to them. For more information on each of the model's please refer to 
the source references. 

Piston in an infinite baffle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
.. image:: ./_static/piston_in_inf_baffle.png.png
	:width: 200

`Piston in an infinite baffle schematic`

Here a rigid circular disk (the 'piston') vibrates back and forth across a hole (of matching size) set in a huge baffle .
The parameters needed to define this model are the wavenumber (`k` - see below for a list of all common abbreviations) and piston radius (`a`).  

Reference : Chp 13, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers. Academic Press.

Point on a sphere
~~~~~~~~~~~~~~~~~
.. image:: ./_static/point_on_sphere.png.png
	:width: 200

`Point on a sphere schematic`

An infinitesimally small portion of a sphere's surface (the 'point') is considered to vibrate. The rest of the sphere does not vibrate.
The parameters needed to define this model are the wavenumber (`k`) and sphere radius (`R`). 

Reference: Chp 12, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

Oscillating cap of a sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ./_static/oscillating_cap_sphere.png.png
	:width: 200
	
`Oscillating cap of a sphere schematic`

The 'sliced' part of a sphere is the 'cap' in this case. The cap moves with an axial velocity :math:`u_{0}`.
Portions of the cap closer to the periphery vibrate less than portions close to the center, this is summarised by the 
relation :math:`u(R,\theta) = u_{0}cos \theta`, where :math:`\theta` is the distended angle from the cap's centre. 

The parameters needed to define this model are wavenumber (`k`), sphere radius (`R`), and the aperture angle of the cap (:math:`\alpha`). 

Reference: Chp 12, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

Piston in a sphere
~~~~~~~~~~~~~~~~~~

.. image:: ./_static/piston_in_a_sphere.png.png
	:width: 200

`Piston in a sphere schematic`

A sphere is sliced, and the 'cap' is discarded. The 'open' portion of the sliced sphere is now replaced with a piston. 
This piston in the sphere vibrates to produce sound. The parameters needed to define this model are wavenumber (`k`), sphere radius (`R`), 
aperture angle of the piston (:math:`\alpha`), and piston radius (`a`). 

Reference: Chp 12, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

Common parameters and abbreviations
-----------------------------------
The inputs will tend to be model-specific, but the common input parameters
to keep in mind are:

    #. `k`, wavenumber. This is :math:`\frac{2\pi}{\lambda}` - this is another wave of defining the frequency of the vibrating part, and the sound produced. :math:`\lambda` is the wavelength of the sound, also defined as :math:`\frac{v_\text{sound}}{\text{frequency}}`. 
    #. `a` : piston radius, wherever applicable
    #. `R` : sphere radius, wherever applicable. 
    #. :math:`\alpha, \theta, \phi`: various angles describing the size of the oscillating portion. These angles are in radians - not degrees!



