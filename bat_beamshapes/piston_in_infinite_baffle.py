#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Piston in an infinite baffle
============================

References
----------
Chp 13, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

"""
import numpy as np 
from bat_beamshapes.utilities import dB
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
def piston_in_infinite_baffle_directivity(angles, param):
    '''
    Calculates relative directivity dB (D(theta)/D(0))
    of a piston in an infinite baffle.
    
    
    Parameters
    ----------
    angles : array-like
        Angles at which the directivity is to be calculated in radians. 
    params : dictionary
        Dictionary with at least the following keys:
            k : float>0
                Wavenumber. 
            a : float>0
                Piston radius
            
    Returns 
    -------
    _ : None 
    beamshape : np.array
        Array with dB(D_theta/D_0).
        The number of items is equal to the number of angles.
    
    '''
    kv,av = [ param[each] for each in ['k','a']]
    dtheta_by_dzero = np.zeros(len(angles))
    for i,angle in enumerate(angles):        
        dtheta_by_dzero[i] = dB(d_theta_func(kv,av,angle))
    return None, dtheta_by_dzero