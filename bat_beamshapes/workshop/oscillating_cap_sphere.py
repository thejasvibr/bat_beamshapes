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
mpmath.mp.dps = 50
from sympy import expand,symbols, Sum,summation, I, cos, sin, legendre, oo
from sympy import lambdify,pi
from bat_beamshapes.special_functions import sph_hankel2, legendre_mvz

n, z,k,R,alpha,theta = symbols('n z k R alpha theta')

# equation 12.59 
# split the big parenthesis into 3 parts (j/sph_hankel, cos(theta) term and the summation)
kR = k*R
d_theta_term1 = I/(2*sph_hankel2.subs({'n':1, 'z':kR})) 

d_theta_term2_num =  3*(1-cos(alpha)**3)*cos(theta)
d_theta_term2_denom = (sin(alpha)**2)*(sph_hankel2.subs({'n':0, 'z':kR})-2*sph_hankel2.subs({'n':2, 'z':kR}))
d_theta_term2 = d_theta_term2_num/d_theta_term2_denom
print('1')
P_1ncosalpha = legendre_mvz.subs({'m':1, 'v':n,'z':cos(theta)})
dtheta_t3_num = (I**(n+1))*((2*n+1)**2)*(sin(alpha)*legendre(n,cos(alpha)) + cos(alpha)*P_1ncosalpha)
dtheta_t3_denom = (n-1)*(n+2)*sin(alpha)*(n*sph_hankel2.subs({'n':n-1, 'z':kR})-(n+1)*sph_hankel2.subs({'n':n+1, 'z':kR}))
dtheta_t3_oneterm = (dtheta_t3_num/dtheta_t3_denom)*legendre(n, cos(theta))



#d_theta_term3 = Sum(expand(dtheta_t3_oneterm), (n,2,oo))
print('2')
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
                                                  'theta': theta_v})
    
    return mpmath.nsum(lambdify([n],version_with_freen, 'mpmath'), [2,mpmath.inf])


#d_theta = lambdify([k,R,alpha,theta], d_theta,'mpmath')

def d_theta(k,R,alpha,theta):
    brackets_term1 = d_theta_t1_func(k,R,alpha,theta)
    brackets_term2 = d_theta_t2_func(k,R,alpha,theta)
    brackets_term3 = dtheta_t3_func(k,R,alpha,theta)
    final_d_theta= (2/k**2*R**2)*(brackets_term1+brackets_term2+brackets_term3)
    return final_d_theta
    

if __name__=='__main__':
    k_v = 2*mpmath.pi/(mpmath.mpf(330.0)/mpmath.mpf(50000))
    R_v = mpmath.mpf(0.01)
    alpha_v = mpmath.pi/10
    theta_v = 0.0
    
    #dtheta_t3_func(k_v, R_v, alpha_v, theta_v)
    print(d_theta(k_v, R_v, alpha_v, theta_v))
    
