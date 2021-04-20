#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common special functions to be used across models

"""

from sympy import besselj, bessely, I, sqrt, pi
from sympy import symbols

n,z = symbols('n z')

# spherical Bessel function
sph_bessel1 = besselj(n+1/2,z)*sqrt(pi/(2*z))
# sphericl Neumann function 
sph_bessel2 = bessely(n+1/2,z)*sqrt(pi/(2*z))
# spherical Hankel function of the second kind 
sph_hankel2 = sph_bessel1 - I*sph_bessel2