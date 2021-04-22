"""
Trying to implement the oscillating cap in a 
rigid sphere described in Chp. 12 of Beranek and Mello 2012
"""

from sympy import symbols, summation, I, cos, sin, legendre, oo
from bat_beamshapes.special_functions import sph_hankel2, legendre_mvz

n,z,k,R,alpha,theta = symbols('n z k R alpha theta')

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
d_theta_term3 = summation(dtheta_t3_oneterm, (n,2,oo))
print('2')
d_theta = 2/(k**2/R**2)*(d_theta_term1 + d_theta_term2 + d_theta_term3)
print('3')