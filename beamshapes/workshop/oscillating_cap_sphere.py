"""
Trying to implement the oscillating cap in a 
rigid sphere described in Chp. 12 of Beranek and Mello 2012


TODO: 
    * Implement the special case of alpha = pi/2, described by equations 
    12.60
"""

from gmpy2 import *
from symengine import * 
import mpmath
#mpmath.mp.dps = 15
from sympy import expand,symbols, Sum,summation, I, cos, sin, legendre, oo,zoo
from sympy import lambdify,pi
from beamshapes.special_functions import sph_hankel2, legendre_mvz

n, z,k,R,alpha,theta = symbols('n z k R alpha theta')

# equation 12.59 
# split the big parenthesis into 3 parts (j/sph_hankel, cos(theta) term and the summation)
d_theta_term1 = I/(2*sph_hankel2.subs({'n':1, 'z':k*R})) 

d_theta_term2_num =  3*(1-cos(alpha)**3)*cos(theta)
d_theta_term2_denom = (sin(alpha)**2)*(sph_hankel2.subs({'n':0, 'z':k*R})-2*sph_hankel2.subs({'n':2, 'z':k*R}))
d_theta_term2 = d_theta_term2_num/d_theta_term2_denom
P_1ncosalpha = legendre_mvz.subs({'m':1, 'v':n,'z':cos(alpha)}).doit()
dtheta_t3_num = (I**(n+1))*((2*n+1)**2)*(sin(alpha)*legendre(n,cos(alpha)) + cos(alpha)*P_1ncosalpha)
dtheta_t3_denom = (n-1)*(n+2)*sin(alpha)*(n*sph_hankel2.subs({'n':n-1, 'z':k*R})-(n+1)*sph_hankel2.subs({'n':n+1, 'z':k*R}))
dtheta_t3_oneterm = (dtheta_t3_num/dtheta_t3_denom)*legendre(n, cos(theta))



#d_theta_term3 = Sum(expand(dtheta_t3_oneterm), (n,2,oo))
#d_theta = 2/(k**2*R**2)*(d_theta_term1 + d_theta_term2 + d_theta_term3)


d_theta_t1_func = lambdify([k,R,alpha,theta], d_theta_term1, 'mpmath')
d_theta_t2_func = lambdify([k,R,alpha,theta], d_theta_term2, 'mpmath')

# %% the problem is here that the summation doesn't happen as 
# 'n' is somehow not treated as an internal integer variable by mpmath
# the alternative is to then to create a function out of it and then 
# perform summation through mpmath.
def dtheta_t3_func(k_v,R_v,alpha_v,theta_v):
    
    version_with_freen = dtheta_t3_oneterm.subs({'k':k_v, 'R':R_v,
                                                 'alpha':alpha_v,
                                                  'theta': theta_v,
                                                  })
    freen_func = lambdify([n],version_with_freen, 'mpmath')
    return mpmath.nsum(freen_func, [2,mpmath.inf])


def d_theta(kv,Rv,alphav,thetav):
    brackets_term1 = d_theta_t1_func(kv,Rv,alphav,thetav)
    brackets_term2 = d_theta_t2_func(kv,Rv,alphav,thetav)
    brackets_term3 = dtheta_t3_func(kv,Rv,alphav,thetav)
    final_d_theta= (2/kv**2*Rv**2)*(brackets_term1+brackets_term2+brackets_term3)
    return final_d_theta

#%% In eqn. 12.61, term 2 differs by the absence of a cos(theta)
d_0_term2_num =  3*(1-cos(alpha)**3)
d_0_term2 = d_0_term2_num/d_theta_term2_denom
d_0_t2_func = lambdify([k,R,alpha,theta], d_0_term2, 'mpmath')


#%% term 3 of eqn. 12.61 doesn't have a Pn(cos(theta)) at the end. 
d_0_t3_oneterm = dtheta_t3_num/dtheta_t3_denom
def d_0_t3_func(k_v,R_v,alpha_v,theta_v):
    
    version_with_freen = d_0_t3_oneterm.subs({'k':k_v, 'R':R_v,
                                                 'alpha':alpha_v,
                                                  'theta': theta_v})
    
    return mpmath.nsum(lambdify([n],version_with_freen, 'mpmath'), [2,mpmath.inf])


def d_zero(kv,Rv,alphav,thetav=0):
    brackets_term1 = d_theta_t1_func(kv,Rv,alphav,thetav)
    brackets_term2 = d_0_t2_func(kv,Rv,alphav,thetav)
    brackets_term3 = d_0_t3_func(kv,Rv,alphav,thetav)
    final_d_0 = (2/kv**2*Rv**2)*(brackets_term1+brackets_term2+brackets_term3)
    return final_d_0







if __name__=='__main__':
    import numpy as np 
    import matplotlib.pyplot as plt
    wavelength = (mpmath.mpf(330.0)/mpmath.mpf(50000))
    ka = 10
    k_v = 2*mpmath.pi/wavelength
    a_v = ka/k_v
    alpha_v = mpmath.pi/3
    R_v = a_v/sin(alpha_v) #mpmath.mpf(0.01)
    
    angles = mpmath.linspace(0,mpmath.pi,50)
    at_angles = mpmath.matrix(1,len(angles))
    for i, angle in enumerate(angles):
        at_angles[i] = d_theta(k_v, R_v, alpha_v, angle)
    on_axis = d_zero(k_v, R_v, alpha_v, angle)
    rato = [20*mpmath.log10(abs(each/on_axis)) for each in at_angles]
    
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(np.array(angles), rato)
    plt.ylim(-40,10);plt.yticks(np.arange(-40,20,10))
    plt.xticks(np.arange(0,2*np.pi,np.pi/6))
    import pandas as pd
    df = pd.read_csv('fig12-17.csv')
    plt.plot(np.radians(df['deg']), df['rel_db_0deg'],'*',
             label='ground-truth, Fig. 12-17')