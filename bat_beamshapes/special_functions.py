#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common special functions to be used across models

"""

from sympy import besselj, bessely, legendre, I, sqrt, pi
from sympy import symbols

p,m,v,n,z = symbols('p m v n z')

# spherical Bessel function
sph_bessel1 = besselj(n+1/2,z)*sqrt(pi/(2*z))
# sphericl Neumann function 
sph_bessel2 = bessely(n+1/2,z)*sqrt(pi/(2*z))
# spherical Hankel function of the second kind 
sph_hankel2 = sph_bessel1 - I*sph_bessel2

# This function can actuall be replaced by the inbuilt 'assoc_legendre'
# pmvz version of legendre function (Appendix II of Beranek & Mello 2012, eqn. 63)
legendre_mvz = ((-1)**m)*((1-z**2)**(m/2))*legendre(v,z).diff((z,m))