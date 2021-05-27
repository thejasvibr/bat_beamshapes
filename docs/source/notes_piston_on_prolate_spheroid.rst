Notes for 'Acoustic directivity of rectangular pistons on prolate spheroids'
============================================================================

* A 'prolate spheroid' is obtained when you rotate an ellipse on its major axis - a symmetric egg basically. 
* The coordinate system is defined by :math:`\xi, \eta, \phi` , just like how in spherical coordinates there's :math:`r,\theta,\phi`, except with the ranges (eqn 2):

.. math::

    1 \leq \xi < \infty

   -1 \leq \eta \leq 1

    0 \leq \phi \leq 2 \pi


Defined variables of interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The `Directivity function` is given by eqn. 15
.. math::

    f(\theta,\phi) = \sum_{m=0}^\infty \sum_{l=m}^\infty \frac{\epsilon_{m}S^{(1)}_{ml}(h,cos\theta)}{R^{(4)\prime}_{ml}(h, \xi) N_{ml}}i^{l+1} \times [\tilde I^{s}_{ml} cos m\phi + \tilde I^{c}_{ml} sin m\phi]


* :math:`\epsilon_{m}` has a piecewise nature, where:

.. math::

    \epsilon_{m}=\begin{cases}
          1 \quad &\text{if} \, m = 0 \\
          1 \quad &\text{if} \, m \neq 0 \\
     \end{cases}

* `h`, 'size parameter' :math:`h=kd/2`, where `k` is the wavenumber and `d` is the 'interfocal distance of the generating ellipse'
* :math:`R^{(4)}_{ml}(h, \xi)` : prolate spheroidal radial function of the 4th kind  (eqn.5), where :math:`R^{(4)}_{ml}(h, \xi) = R^{(1)}_{ml}(h, \xi) - iR^{(2)}_{ml}(h, \xi)`. Defined in [2]
* :math:`S^{(1)}_{ml}(h, \eta)`: prolate spheroidal angle function of the 1st kind (eqn. 4). Defined in [2].
* :math:`N_{ml}` : prolate spheroidal angle normalization factor (eqn. 11)
* :math:`A_{ml}, B_{ml}` : unknown coefficients to be estimated -- related to the boundary condition of 'particle velocity at the spheroid-fluid interface'. Defined in [3].
* :math:`\tilde I^{s}_{ml} cos m\phi` (eqn. 12) and :math:`\tilde I^{c}_{ml} sin m\phi` '..define the size, shape, and location of the piston in the spheroisal baffle..'. In detail, 
.. math::

    \tilde I^{s}_{ml} = \int \int_{S_{i}} (\xi^{2}_{0} - \eta^{2})^{1/2}S^{(1)}_{ml}(h, \eta) cos m\phi \:d\eta d\phi

    \tilde I^{c}_{ml} = \int \int_{S_{i}} (\xi^{2}_{0} - \eta^{2})^{1/2}S^{(1)}_{ml}(h, \eta) sin m\phi \:d\eta d\phi

* The `directivity` itself is given by eqn. 16: 

.. math::

    F(\theta, \phi) = f(\theta,\phi)/f(\theta_{0},\phi_{0})

where :math:`\theta_{0},\phi_{0}` 'define the direction of maximum response'


Notes 
~~~~~
* :math:`R^{(1,2,4)}_{ml}(h, \xi)` and :math:`S^{(1)}_{ml}(h, \eta)` are 
* There seem to be something relevant `here <https://docs.scipy.org/doc/scipy/reference/special.html>`_ (Scipy implementations).


References
~~~~~~~~~~
#. Boisvert & Buren 2004, Acoustic directivity of rectangular pistons on prolate spheroids, JASA, 116, 1932 (2004); doi: 10.1121/1.1778840
#. C. Flammer, Spheroidal Wave Functions, Stanford University Press, Stanford, CA, 1957
#. J. E. Boisvert and A. L. Van Buren, ‘‘Acoustic radiation impedance of rectangular pistons on prolate spheroids,’’ J. Acoust. Soc. Am. 111, 867–874 (2002)
