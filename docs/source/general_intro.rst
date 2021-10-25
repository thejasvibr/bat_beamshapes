General Intro to `beamshapes`
=============================
Sound sources don't radiate sound uniformly most of the time. Even us humans for instance, we emit more energy while talking to the front than to the back. The exact (non-uniform) pattern in which sound is radiated defines the 'beamshape' of a source. The beamshape is typically a combination of the frequency of emitted sound and the geometry of the vibrating surface and its associated (non-vibrating) surfaces. 

Let's go through the main concepts required to use the `beamshapes` package.

Source models
-------------
A 'beamshape' is tightly bound to the input parameters, and the sound source model behind
it eg. piston in a sphere, cap of a sphere etc. The source model refers to the source
of sound and its geometric properties, and assumptions of the vibrations etc. According
to the situation in hand, different source models may be physically/biologically relevant! 


Piston in an infinite baffle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
.. image:: ./_static/piston_in_inf_baffle.png.png
	:width: 200

`Piston in an infinite baffle schematic`

Here a rigid circular disk (the 'piston') vibrates back and forth between a hole in a huge baffle that matches the disk's size.
The parameters needed to define this model are the wavenumber and piston radius. 

Point on a sphere
~~~~~~~~~~~~~~~~~
.. image:: ./_static/point_on_sphere.png.png
	:width: 200

`Point on a sphere schematic`

An infinitesimally small portion of a sphere's surface (the 'point') is considered to vibrate. The rest of the sphere does not vibrate.
The parameters needed to define this model are the wavenumber and sphere radius. 


Oscillating cap of a sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ./_static/oscillating_cap_sphere.png.png
	:width: 200
	
`Oscillating cap of a sphere schematic`

The 'sliced' part of a sphere is the 'cap' in this case. The cap moves with an axial velocity :math:`u_{0}`.
Portions of the cap closer to the periphery vibrate less than portions close to the center, this is summarised by the 
relation :math:`u(R,\theta) = u_{0}cos \theta`, where :math:`\theta` is the distended angle from the cap's centre. 

The parameters needed to define this model are wavenumber, sphere radius, and the aperture angle of the cap. 

Piston in a sphere
~~~~~~~~~~~~~~~~~~

.. image:: ./_static/piston_in_a_sphere.png.png
	:width: 200

`Piston in a sphere schematic`

A sphere is sliced, and the 'cap' is discarded. The 'open' portion of the sliced sphere is now replaced with a piston. 
This piston in the sphere vibrates to produce sound. The parameters needed to define this model are wavenumber, sphere radius, 
aperture angle of the piston, and piston radius. 

Input Parameters
----------------
The inputs will tend to be model-specific, but the common input parameters
to keep in mind are:

    #. `k`, wavenumber. This is :math:`\frac{2\pi}{\lambda}` - this is another wave of defining the frequency of the vibrating part, and the sound produced. :math:`\lambda`is the wavelength of the sound, also defined as :math:`\frac{v_\text{sound}}{\text{frequency}}`. 
    #. `a` : piston radius, wherever applicable
    #. `R` : sphere radius, wherever applicable. 
    #. :math:`\alpha, \theta, \phi`: various angles describing the size of the oscillating portion. These angles are in radians - not degrees!


Directionality functions
------------------------
Each source model has an associated `directionality` function with it. 
The directionality function returns the :math:`20log_{10}(D_{\theta}/D_{0})`
value for a given set of angles. Some models require a long time to compute 
coefficients, and also output other objects to save time for future runs.

The oscillating cap in a sphere does not produce additional output

.. code-block:: shell

    # .... having chosen certain input parameters and put them into  input_params
    >>> import mpmath
    >>> import beamshapes as beamshapes
    >>> from beamshapes import cap_in_sphere_directionality
    >>> angles = mpmath.linspace(0,pi,10) 
    >>>  _, spherecap_beam = cap_in_sphere_directionality(angles, input_params)

In case the chosen source model's directionality requires intensive calculations
to generate estimates, then the calculated outputs are also returned. 

.. code-block:: shell

    # .... having chosen certain input parameters and put them into input_params 
    >>> import mpmath
    >>> import beamshapes as beamshapes
    >>> from beamshapes import piston_in_sphere_directionality
    >>> angles = mpmath.linspace(0,pi,10) 
    >>> An_out, spherepiston_beam = piston_in_sphere_directionality(angles, input_params)
    # in case you need to calculate more points now - it saves time to do this:
    >>> input_params['An'] = An_out
    >>> new_angles = mpmath.linspace(0,pi,100)
    >>> _, detailed_spherepiston_beam = piston_in_sphere_directionality(new_angles, input_params)




