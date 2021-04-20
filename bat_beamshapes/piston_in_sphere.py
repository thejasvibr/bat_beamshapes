#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TODO: 
    1. IMPLEMENT SOME SORT OF CACHING/MEMOIZATION to prevent recalculation of
    old parameter values. 

Code that calculates piston in a sphere. 
This model is described in Chapter 12 of Beranek & Mellow 2012, 
and parts of it are based on the Mathematica code provided by 
Tim Mellow. 


Parameters
----------
alpha : 0<float<pi
    The 'gape'/'aperture' of the piston. The half-angle it occupies 
    in radians. 
k : float>0
    Wavenumber, 2*pi/wavelength
a : 



References
----------
Beranek, L. L., & Mellow, T. (2012). Acoustics: sound fields and transducers.
Academic Press.

"""
from gmpy2 import *
import mpmath
from mpmath import mpf
import numpy as np
from symengine import * 
import sympy
from sympy import symbols, legendre, sin, cos, tan, summation, I, diff, pi, sqrt
from sympy import Matrix, besselj, bessely, Piecewise
from sympy import  lambdify, integrate, expand,Integral
import tqdm
x, alpha, index, k, m,n,p, r1, R, theta, y, z = symbols('x alpha index k m n p r1 R theta,y,z')
from sympy import N, cse
dps = 50; mpmath.mp.dps = dps # default digit precision set to 50

from special_functions import sph_hankel2

r1 = (R*cos(alpha))/cos(theta)


#%%
# equation 12.106
alternate_hankels = n*sph_hankel2.subs({'n':n-1, 'z':k*r1})-(n+1)*sph_hankel2.subs({'n':n+1, 'z':k*r1})
pt1_postterm = legendre(n,cos(theta))*cos(theta)
Imn_pt1 = alternate_hankels*pt1_postterm

pt2_preterm = n*(n+1)*sph_hankel2.subs({'z':k*r1})
alternate_legendres = (legendre(n-1,cos(theta))-legendre(n+1,cos(theta)))/(k*r1)
Imn_pt2 = pt2_preterm*alternate_legendres

whole_postterm = legendre(m,cos(theta))*(r1**2/R**2)*tan(theta)

Imn_term = (Imn_pt1 + Imn_pt2)*whole_postterm 
Imn = Integral(Imn_term,(theta,0,alpha))
Imn_term_func = lambdify([m,n,k,R,alpha,theta], Imn_term, 'mpmath')
Imn_func = lambdify([m,n,k,R,alpha],Imn,'mpmath') 




#%%
# equation 12.107
Kmn_expr = legendre(n, cos(theta))*legendre(m, cos(theta))*sin(theta) # the integrl of this expression 
# has a solution given in Appendix II, eqn 70
legendre_1stderiv = diff(legendre(n,z),z)
# when m != n
num_legendre_term1 = legendre(m,cos(alpha))*legendre_1stderiv.subs({'z':cos(alpha)})
num_legendre_term2 = legendre(n,cos(alpha))*legendre_1stderiv.subs({'n':m,'z':cos(alpha)})
eqn70_mnoteqn = sin(alpha)*(num_legendre_term1-num_legendre_term2)/(m*(m+1)-n*(n+1))

# when m==n
summn_funcn = legendre(index,cos(alpha))*(legendre(index,cos(alpha))*cos(alpha)-legendre(index+1,cos(alpha)))
# substituting 'j' for 'index' because j is also used for sqrt(-1) all through the book!!

meqn_sumterm = 2*summation(summn_funcn, (index,1,m-1))
eqn70_meqn = (1+ cos(alpha)*legendre(m,cos(alpha))**2 + meqn_sumterm)/(2*m+1)
Kmn = Piecewise((eqn70_mnoteqn,m>n),(eqn70_mnoteqn,m<n), (eqn70_meqn,True), )
Kmn_func = lambdify([m,n,alpha],Kmn,'mpmath')



#%%
# equation 12.108
Lm_expr = legendre(m,cos(theta))*(r1**2/R**2)*tan(theta)
Lm = Integral(Lm_expr, (theta,0,alpha))
Lm_func = lambdify([m, R, alpha], Lm, 'mpmath')


#%%
# b matrix
b = -I*Lm
b_func = lambdify([m,alpha], b,'mpmath') # eqn. 12.104


#%% 
# Setting up the matrices 
# %% 
mmn_hankels = n*sph_hankel2.subs({'n':n-1,'z':k*R})-(n+1)*sph_hankel2.subs({'n':n+1,'z':k*R})
mmn_hankels_func = lambdify([n,k,R], mmn_hankels,'mpmath')


def compute_Mmn(params):
    Nv = 12 + int(2*params['ka']/sin(params['alpha']))
    M_matrix = mpmath.matrix(Nv,Nv)
    for i in tqdm.trange(Nv):
        for j in range(Nv):
            params['m'],params['n'] = i,j
            Imn_value = Imn_func(params['m'], params['n'],params['k'],params['R'],params['alpha'])
            Kmn_value = Kmn_func(params['m'],params['n'],params['alpha'])
            numerator_hankels = mmn_hankels_func(j,params['k'],params['R'])
            numerator = Imn_value+ numerator_hankels*Kmn_value
            denom = 2*params['n']+1
            M_matrix[params['m'],params['n']] = numerator/denom
    return M_matrix

def compute_b(params):
    Nv = 12 + int(2*params['ka']/sin(params['alpha']))
    b_matrix = mpmath.matrix(Nv,1)
    for each_m in range(Nv):
        b_matrix[each_m,:] = b_func(each_m, params['alpha'])
    return b_matrix 

def compute_a(M_mat, b_mat):
    a_matrix = mpmath.inverse(M_mat)*b_mat
    return a_matrix

def d_theta(angle,k_v,R_v,alpha_v,An):
    num = 4 
    N_v = An.rows
    denom  = (k_v**2)*(R_v**2)*mpmath.sin(alpha_v)**2
    part1 = num/denom
    jn_matrix = np.array([1j**f for f in range(N_v)])
    legendre_matrix = mpmath.matrix([legendre(n_v, mpmath.cos(angle)) for n_v in range(N_v)])
    part2_matrix = np.column_stack((An, jn_matrix, legendre_matrix))
    part2 = np.sum(np.apply_along_axis(lambda X: X[0]*X[1]*X[2], 1, part2_matrix))
    
    rel_level = - part1*part2
    return rel_level

def d_zero(k_v,R_v,alpha_v,An):
    num = 4 
    N_v = An.rows
    denom  = (k_v**2)*(R_v**2)*mpmath.sin(alpha_v)**2
    part1 = num/denom
    jn_matrix = np.array([1j**f for f in range(N_v)])
    part2_matrix = np.column_stack((An, jn_matrix))
    part2 = np.sum(np.apply_along_axis(lambda X: X[0]*X[1], 1, part2_matrix))
    rel_level = - part1*part2
    return rel_level

def relative_directionality_db(angle,k_v,R_v,alpha_v,An):
    off_axis = d_theta(angle,k_v,R_v,alpha_v,An)
    on_axis = d_zero(k_v,R_v,alpha_v,An)
    rel_level = 20*mpmath.log10(abs(off_axis/on_axis))
    return rel_level

def check_input_validity(params):
    '''
    '''
    a = params['ka']/params['k']
    if type(a) == mpmath.ctx_mp_python.mpf:
        a_as_expected = a == params['R']*mpmath.sin(params['alpha'])
    else:
        a_as_expected = a == params['R']*np.sin(params['alpha'])
    if not a_as_expected:
        raise ValueError('a != R*sin(alpha), please check the input values')
    
    

def directionality(angles, params):
    '''
    Calculates relative directionality dB (D(theta)/D(0))
    of a piston in a rigid sphere.
    
    
    Parameters
    ----------
    angles : array-like
        Angles at which the directionality is to be calculated in radians. 
    params : dictionary
        Dictionary with at least the following keys:
            k : mpmath.mpf>0
                Wavenumber. 
            ka : mpmath.mpf>0
                Product of wavenumber and radius of piston.
            R : mpmath.mpf>0
                Radius of sphere
            alpha: 0<mpmath.mpf<pi
                Angular aperture of sphere. 
    
    Returns 
    -------
    directionality : list
        List with relative directionalities in 20log10 dB. 
        The number of items is equal to the number of angles.

    Notes
    -----
    While the dictionary entries can be normal floats, it is best 
    to use mpmath.mpf float objects to maintain digit precision. 
    eg. when specifying an angle of 60degrees for the 
    aperture alpha, use ```mpmath.pi/3``` rather than ```np.pi/3```
    for instance.
        
    
    '''
    check_input_validity(params)
    Mmatrix = compute_Mmn(params)
    bmatrix = compute_b(params)
    amatrix = compute_a(Mmatrix, bmatrix)
    
    directionality = []
    for angle_v in angles:
        directionality.append(relative_directionality_db(angle_v,
                                                         params['k'],
                                                         params['R'], 
                                                         params['alpha'],
                                                         amatrix))
    
    return directionality

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    frequency = 50*10**3 # kHz
    vsound = 330 # m/s
    wavelength = vsound/frequency
    alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
    k_value = 2*mpmath.pi/(wavelength)
    ka = 5
    a_value = ka/k_value 
    R_value = a_value/mpmath.sin(alpha_value)  # m
    paramv = {}
    paramv['R'] = R_value
    paramv['alpha'] = alpha_value
    paramv['k'] = k_value
    paramv['ka'] = ka
    
    angles = mpmath.linspace(0,mpmath.pi,100)
    beamshape = directionality(angles, paramv)
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(angles, beamshape)
    plt.ylim(-40,0);plt.yticks(np.arange(-40,10,10))
    plt.xticks(np.arange(0,2*np.pi,np.pi/6))
    
    


