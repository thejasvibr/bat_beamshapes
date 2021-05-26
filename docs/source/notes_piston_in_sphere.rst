Notes for Piston in a Sphere (Beranek & Mellow 2012)
============================================================================

This page specifically highlights some of the discrepancies in code and equations I've noticed.
Given my lack of specialised math, I'm assuming they are the result of typo's and proceeding
with what seems to be the correct versions.

Since the original figures may be the result of typo's how to test the function outputs?
I've resorted to altering only the suspect portions of the model and implementing these in the tests.
This method checks that the rest of the equations/functions are implemented correctly.

Discrepancy 1: eqn. 12.98 (:math:`sin(\theta)` --> :math:`-sin^{2}(\theta)`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Equation 12.98 is 

.. math::

    P^{\prime}_{n}(cos(\theta) = \frac{\partial}{\partial \theta} P_n(cos(\theta)) = - \frac{n(n+1)}{(2n+1)sin(\theta)}(P_{n-1}cos(\theta)-P_{n+1}cos(\theta))

Eqn. 12.98 refers to Appendix II eqn. 65.

However, App.II, eqn. 65 defines the derivative of the Legendre function as:

.. math::
    
    n(n+1)(P_{n+1}(z)-P_{n-1}(z)) = (2n+1)(z^{2}-1)\frac{d}{dz}P_{n}(z)

and moving :math:`\frac{d}{dz}P_{n}(z)` to the LHS gives:

.. math::

    \frac{d}{dz}P_{n}(z) = \frac{n(n+1)(P_{n+1}(z)-P_{n-1}(z))}{(2n+1)(z^{2}-1)}

Now substituting :math:`z = cos(\theta)`, we get:

.. math::
    
    \frac{d}{dz}P_{n}(z) = \frac{n(n+1)}{(2n+1)(cos^{2}(\theta)-1)}(P_{n+1}(cos(\theta))-P_{n-1}(cos(\theta)))

Applying the identity :math:`1 = sin^{2}(\theta) + cos^{2}(\theta)` , and thus substituting :math:`cos^{2}(\theta) - 1 = - sin^{2}(\theta)`
we get:

.. math::

    \frac{\partial}{\partial \theta} P_n(cos(\theta)) = \frac{n(n+1)}{(2n+1)(-sin^{2}(\theta))}(P_{n+1}(cos(\theta))-P_{n-1}(cos(\theta)))

And for better comparability:

.. math::

   \frac{\partial}{\partial \theta} P_n(cos(\theta)) = - \frac{n(n+1)}{(2n+1)(sin^{2}(\theta))}(P_{n-1}(cos(\theta))-P_{n+1}(cos(\theta)))


Here, I therefore think the :math:`sin(\theta)` in eqn. 12.98 term should be :math:`-sin^2{\theta}`. In the original Mathematica code, 
(which was presumably used to generate Fig. 12.23) - the :math:`sin(\theta)` term is present.



Discrepancy 2: coding order of :math:`P^{\prime}_{n}(cos(\theta))` in Mathematica code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solution for :math:`K_{mn}` is given in App.II, eqn. 70. In the case where  :math:`m \neq n`, the solution is:

.. math:: 

    \frac{sin(\alpha)( P_{m}(cos(\alpha))P^{\prime}_{n}(cos(\alpha)) - P_{n}(cos(\alpha))P^{\prime}_{m}(cos(\alpha)))}{m(m+1) - n(n+1)}

When we do the substitutions for :math:`P^{\prime}_{n}(cos(\alpha))` and :math:`P^{\prime}_{m}(cos(\alpha))` the full term will be:

.. math::

    \small
    \frac{sin(\alpha)}{m(m+1) - n(n+1)} \\
    \left( P_{m}(cos(\alpha))\frac{n(n+1)}{(2n+1)(-sin^{2}(\alpha))}(P_{n+1}(cos(\alpha))-P_{n-1}(cos(\alpha))) \\
     - P_{n}(cos(\alpha))\frac{m(m+1)}{(2m+1)(-sin^{2}(\alpha))}(P_{m+1}(cos(\alpha))-P_{m-1}(cos(\alpha))) \right)

TODO 
~~~~
However, in the code what seems to be implemented is actuall


References
~~~~~~~~~~
Boisvert & Buren 2004, Acoustic directivity of rectangular pistons on prolate spheroids, JASA, 116, 1932 (2004); doi: 10.1121/1.1778840
