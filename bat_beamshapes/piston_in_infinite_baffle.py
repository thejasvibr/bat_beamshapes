#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Piston in an infinite baffle
============================


Parameters
----------
k : float>0
    Wavenumber, 2*pi/wavelength
a : float>0
    Radius of piston. 

References
----------
Chp 13, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

"""

from sympy import besselj, sin, symbols, lambdify, limit
k,a,theta = symbols('k a theta')
#%%
# eqn. 13.102
kasintheta = k*a*sin(theta)
d_theta_num = 2*besselj(1,kasintheta)
d_theta_denom = kasintheta
d_theta = d_theta_num/d_theta_denom

def d_theta_func(kv,av,thetav):
    '''
    The function must be used in the limit version to get valid outputs. 
    '''
    # the limit version 
    subs_dtheta = d_theta.subs({'k':kv,'a':av})
    return lambdify([], limit(subs_dtheta, theta, thetav))()

#d_theta_func = lambdify([k,a,theta],d_theta_num)




#%%
def piston_in_infinite_baffle_directionality(angles, param):
    '''
    '''
    kv,av = [ param[each] for each in ['k','a']]
    dtheta_by_dzero = []
    for angle in angles:
        
        dtheta_by_dzero.append(d_theta_func(kv,av,angle)/1.0)
    return None, dtheta_by_dzero