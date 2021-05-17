#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Point source on a sphere
========================
An infinitesimally oscillating on the surface of a sphere. Useful where the 
source is very small in comparison to the sphere.
   
Parameters
----------
k : float > 0 
    WAvenumber
R : float > 0
    Radius of sphere


References 
----------
Chp 12, secn. 12.5, Beranek, L. L., & Mellow, T. (2012). Acoustics: sound fields and transducers.
Academic Press.

"""
import numpy as np 
from sympy import symbols, legendre, lambdify, I, Sum, cos
from bat_beamshapes.special_functions import sph_hankel2
from bat_beamshapes.utilities import dB

k, R, n, theta, NN, kR = symbols('k R n theta NN kR')
costheta = cos(theta)
kR = k*R
#%% eqn. 12.46

numerator_sumterm = I**(n+1)*(2*n+1)**2*legendre(n, costheta)
denominator_sumterm = n*sph_hankel2.subs({'n':n-1,'z':kR}) - (n+1)*sph_hankel2.subs({'n':n+1,'z':kR})
sum_term = Sum(numerator_sumterm/denominator_sumterm, (n, 0, NN))
d_theta = (1/(kR**2))*sum_term

d_theta_func = lambdify([theta,k,R,NN],d_theta)
#%% eqn. 12.47 
onaxis_num_sumterm = I**n*(2*n+1)**2
onaxissum_term = Sum(onaxis_num_sumterm/denominator_sumterm, (n,0,NN))
d_zero = (1/(kR**2))*onaxissum_term

d_zero_func = lambdify([k,R,NN],d_zero)

#%%

def point_source_on_a_sphere_directionality(angles, param, parallel=False):
    '''
        
    Parameters
    ----------
        angles : list/array-like
            The angles for which the directionality is to be calculated. 
            All angles in radians.

        param: dictionary 
            With the following keys and entries:
                k : float > 0 
                    Wavenumber
                R : float > 0
                    Radius of sphere
    
    Returns 
    -------
    _ : None
    dtheta_by_d0 : np.array
        20log10(D_theta/D_0)
             
    '''
    kv, Rv = [param[each] for each in ['k','R']]
    NNv = param.get('NN', int(10+2*kv*Rv)) # trial formula - need to see if it matches w textbook
    
    onaxis_value = abs(d_zero_func(kv,Rv,NNv))
    offaxis_values = []
    for angle in angles:
        offaxis_values.append(abs(d_theta_func(angle, kv, Rv, NNv)))
    dtheta_by_d0 = dB(np.abs(np.array(offaxis_values))/np.abs(onaxis_value))
    return None , dtheta_by_d0

#%% 
if __name__ == '__main__':
    import matplotlib.pyplot as plt    
    import numpy as np
    
    angles = -np.linspace(0,np.pi,100)
    
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.ylim(-40,10);plt.yticks(np.arange(-40,10,10))
    for kv in [100, 50, 30, 10]:
        paramv = {'k':kv,'R':0.1}
        _, dirnlty = point_source_on_a_sphere_directionality(angles,paramv)    
        plt.plot(angles, dirnlty, label='ka: '+str(kv*0.1))
        angles *= -1 
        
    plt.legend(loc=(0.1,0.9))
    
    
