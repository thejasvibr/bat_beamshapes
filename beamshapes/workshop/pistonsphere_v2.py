#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
My attempt at RE-building the piston in a sphere. 

Created on Tue May 11 14:18:53 2021

@author: autumn
"""

import copy
from gmpy2 import *
from joblib import Parallel, delayed 
import mpmath
from mpmath import mpf
import numpy as np
from symengine import * 
import sympy
from sympy import symbols, legendre, sin, cos, tan, Sum, I, diff, pi, sqrt
from sympy import Matrix, Piecewise, IndexedBase
from sympy import  lambdify, expand, Integral
from sympy import HadamardProduct as HP
from beamshapes.special_functions import h2_nz
import tqdm
x, alpha, index, k, m,n,p, r1, R, theta, y, z = symbols('x alpha index k m n p r1 R theta,y,z')
jj, nn = symbols('jj nn')

mpmath.mp.dps = dps

from beamshapes.special_functions import sph_hankel2
from beamshapes.utilities import args_to_mpmath, args_to_str

r1 = (R*cos(alpha))/cos(theta)



# eqn 12.106
Pn_costheta = legendre(n, cos(theta))
Pm_costheta = legendre(m, cos(theta))
Pnminus1_costheta = Pn_costheta.subs(n,n-1)
Pnplus1_costheta = Pn_costheta.subs(n,n+1)


h2_n_kr1 = h2_nz.subs({'z':k*r1})
h2_nminus1_kr1 = h2_n_kr1.subs({'n':n-1})
h2_nplut1_kr1 = h2_n_kr1.subs({'n':n+1})

imn_bigbracket = (n*h2_nminus1_kr1-(n+1)*h2_nplut1_kr1)*Pn_costheta*cos(theta)+ n*(n+1)*h2_n_kr1*(Pnminus1_costheta-Pnplus1_costheta)/(k*r1)

Imn = Integral( imn_bigbracket*Pm_costheta*(r1**2/R**2)*tan(theta), (theta, 0, alpha))

Imn_func = lambdify([m,n,alpha,k,R,], Imn, 'mpmath')

# solution to eqn 12.107
prime_pncosalpha = diff(legendre(n, cos(alpha)), alpha)
prime_pmcosalpha = diff(legendre(m, cos(alpha)), alpha)
Pn_cosalpha = legendre(n,cos(alpha))
Pm_cosalpha = legendre(m,cos(alpha))
m_notequal_n = sin(alpha)*(Pm_cosalpha*prime_pncosalpha-Pn_cosalpha*prime_pmcosalpha)/(m*(m+1)-n*(n+1))

Pjj_cosalpha = legendre(jj,cos(alpha))

# use jj instead of 'j' to avoid confusion with imaginary number *j*
m_equal_n = (1+ cos(alpha)*Pm_cosalpha**2 + 2*Sum(Pjj_cosalpha*(Pjj_cosalpha*cos(alpha) - Pjj_cosalpha.subs(jj,jj+1)), (jj, 0,m-1)))/(2*m+1)

Kmn_soln = Piecewise((m_notequal_n, m>n),
                     (m_notequal_n, m<n),
                     (m_equal_n, True))
Kmn_func = lambdify([m,n,alpha], Kmn_soln, 'mpmath')

# eqn 12.108 for bm
L_m = Integral(legendre(m, cos(theta))*(r1**2/R**2)*tan(theta), (theta,0,alpha))
b_m = -I*L_m
bm_func = lambdify([m,alpha,R], b_m, 'mpmath')


# calculate M matrix
h2_nminus1_kR = h2_nz.subs({'n':n-1,'z':k*R})
h2_nplus1_kR = h2_nz.subs({'n':n+1,'z':k*R})
Mmn = (Imn+(n*h2_nminus1_kR - (n+1)*h2_nplus1_kR)*Kmn_soln)/2*n+1
Mmn_func =  lambdify([m,n,k,R,alpha], Mmn, 'mpmath')


def calculate_An(params):
    '''
    '''
    av = params['a']
    kv = params['k']
    Rv = params['R']
    alphav = params['alpha']
    NN = params.get('NN', int(12 + 2*kv*av/mpmath.sin(alphav)))
    
    M_mat = mpmath.matrix(NN,NN)
    b_mat = mpmath.matrix(NN,1)
    for mm in tqdm.trange(NN):
        b_mat[mm] = bm_func(mm, alphav, Rv)
        for nn in range(NN):
            M_mat[mm,nn] =  Mmn_func(mm,nn,kv,Rv,alphav)
    
    An = mpmath.lu_solve(M_mat, b_mat)
    return An
A_n = IndexedBase('A_n')
NN = symbols('NN')

d_zero = -(4/R**2*k**2*sin(alpha)**2)*Sum(A_n[nn]*I**nn, (nn,0,NN))
d_theta = -(4/R**2*k**2*sin(alpha)**2)*Sum(A_n[nn]*I**nn*legendre(nn,cos(theta)), (nn,0,NN))
d_theta_func = lambdify([theta,alpha,k,R,A_n,NN], d_theta,'mpmath')
d_zero_func = lambdify([alpha,k,R,A_n,NN], d_zero,'mpmath')


def calc_dtheta(thetav, alphav, kv,Rv,Anv):
    NNv = len(Anv)-1
    return d_theta_func(thetav, alphav, kv,Rv,Anv,NNv)

def calc_dzero(alphav, kv,Rv,Anv):
    NNv = len(Anv)-1
    return d_zero_func(alphav,kv,Rv,Anv,NNv)


    
if __name__== '__main__':
    dps = 50;
    mpmath.mp.dps = dps
    frequency = mpmath.mpf(50*10**3) # kHz
    vsound = mpmath.mpf(330) # m/s
    wavelength = vsound/frequency
    alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
    k_value = 2*mpmath.pi/(wavelength)
    ka_val = 1
    print(f'Starting piston in sphere for ka={ka_val}')
    ka = mpmath.mpf(ka_val)
    a_value = ka/k_value 
    R_value = a_value/mpmath.sin(alpha_value)  # m
    paramv = {}
    paramv['R'] = R_value
    paramv['alpha'] = alpha_value
    paramv['k'] = k_value
    paramv['a'] = a_value

    an = calculate_An(paramv)
    
    
    direc = []
    for thetav in [0, mpmath.pi/3+mpmath.pi/12, mpmath.pi/4, mpmath.pi]:
        num = abs(calc_dtheta(thetav, paramv['alpha'],
                              paramv['k'],
                              paramv['R'],
                              an))
        denom = abs(calc_dzero(paramv['alpha'],
                              paramv['k'],
                              paramv['R'],
                              an))
        direc.append(20*mpmath.log10(num/denom))
