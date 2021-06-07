#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rectangular cap of a sphere
===========================

:math:`\alpha, \beta` are the half-angles defining the width/height of the 
rectangular cap.

:math:`R` is the radius of the sphere. 

"""

import copy
from joblib import Parallel, delayed
import numpy as np
import mpmath
from sympy import symbols,  I, cos, sin, legendre
from sympy import lambdify, atan, tan, sec, sqrt, pi
from sympy import IndexedBase, Sum, Integral, assoc_legendre
import tqdm
import warnings
from beamshapes.special_functions import sph_hankel2, legendre_mvz
from beamshapes.utilities import args_to_mpmath, args_to_str

n, N, z, k, R, alpha, beta, theta = symbols('n N z k R alpha beta theta')
phi = symbols('phi')
#%%
secsq_alpha = sec(alpha)**2
secsq_beta = secsq_alpha.subs(alpha, beta) 
s_pt1 = (atan(tan(alpha)*tan(beta))/ (secsq_alpha + sqrt(secsq_alpha + tan(beta)**2 )) )
s_pt2 = (atan(tan(alpha)*tan(beta))/ (secsq_beta  + sqrt(secsq_beta  + tan(alpha)**2)) )

S = (4*R**2)*(s_pt1+s_pt2)

sval = lambdify([alpha, beta, R], S,'mpmath')

#%%
# First implement the on-axis response. 

subsdict = {'n':1,'alpha':pi/2.5,'beta':pi/2.5}


ion_pt1_1 = (tan(alpha)/(sqrt(cos(phi)**2+tan(alpha)**2)))
ion_pt1_asocleg_term = cos(phi)/(sqrt(cos(phi)**2+tan(alpha)**2))
ion_pt1_2 = assoc_legendre(n,-1,ion_pt1_asocleg_term)
ion_integral_limit1 = atan(tan(beta)/tan(alpha))
ion_pt1 = Integral(ion_pt1_1*ion_pt1_2, (phi,0,ion_integral_limit1))
Ion_pt1 = lambdify([n,alpha,beta], ion_pt1,'sympy')


sqrt_num_pt2 = sqrt(sin(phi)**2+tan(beta)**2)
ion_pt2_1 = tan(beta)/sqrt_num_pt2
ion_pt2_assocleg_term = sin(phi)/sqrt_num_pt2
ion_pt2_2 = assoc_legendre(n,-1,ion_pt2_assocleg_term)
ion2_integral_limit = (pi/2-atan(tan(alpha)/tan(beta)), pi/2+atan(tan(alpha)/tan(beta)))
ion_pt2 = Integral(ion_pt2_1*ion_pt2_2,
                   (phi,ion2_integral_limit[0],ion2_integral_limit[1]))
Ion_pt2 = lambdify([n,alpha,beta], ion_pt2,'sympy')
# print(Ion_pt2(subsdict['n'],subsdict['alpha'],subsdict['beta']).evalf())

sqrt_num_pt3 = sqrt(cos(phi)**2+tan(alpha)**2)
ion_pt3_1 = tan(alpha)/sqrt_num_pt3
ion_pt3_assocleg_term = -cos(phi)/sqrt_num_pt3
ion_pt3_2 = assoc_legendre(n,-1,ion_pt3_assocleg_term)
ion3_integral_limit = (pi-atan(tan(beta)/tan(alpha)), pi)
ion_pt3 = Integral(ion_pt3_1*ion_pt3_2, 
                   (phi,ion3_integral_limit[0],ion3_integral_limit[1]))

Ion_pt3 = lambdify([n,alpha,beta], ion_pt3,'sympy')
# print(Ion_pt3(subsdict['n'],subsdict['alpha'],subsdict['beta']).evalf())

#%%

Ion = ion_pt1 + ion_pt1 + ion_pt3
Ion_func = lambdify([n, alpha, beta], Ion, 'sympy')
print(Ion_func(subsdict['n'],subsdict['alpha'],subsdict['beta']).evalf())


#%%

