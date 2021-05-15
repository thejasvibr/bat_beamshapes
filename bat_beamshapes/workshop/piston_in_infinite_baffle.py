#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Piston in an infinite baffle
============================


Parameters
----------
alpha : 0<float<pi
    The 'gape'/'aperture' of the piston. The half-angle it occupies 
    in radians. 
k : float>0
    Wavenumber, 2*pi/wavelength
a : float>0
    Radius of piston. 

References
----------
Chp 13, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

"""

from sympy import besselj, sin, symbols, lambdify
k,a,theta = symbols('k a theta')
#%%
# eqn. 13.102
kasintheta = k*a*sin(theta)
d_theta_num = 2*besselj(1,kasintheta)
d_theta_denom = kasintheta
d_theta = d_theta_num/d_theta_denom

d_theta_func = lambdify([k,a,theta],d_theta_num)

#%%
def piston_in_infinite_baffle_directionality(angles, param):
    '''
    '''
    kv,av = [ param[each] for each in ['k','a']]
    dtheta_by_dzero = []
    for angle in angles:
        
        dtheta_by_dzero.append(d_theta_func(kv,av,angle)/1.0)
    return None, dtheta_by_dzero


if __name__ == '__main__':
    import math
    paramv = {'k':100.0, 'a':0.1}
    theta_vals = [0, math.pi/4, math.pi/3, math.pi/2]
    
    _, dirns = piston_in_infinite_baffle_directionality(theta_vals, paramv)
    print(dirns)

