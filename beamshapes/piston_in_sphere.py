#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Piston in a sphere
==================
This source model assumes a piston whose outline matches a cut-out portion of a sphere
exactly. 


See Also 
--------
beamshapes.piston_in_sphere_flint

Notes
-----
This module is an mpmath implementation compatible with multiple operating systems.

The solution for eqn.12.107 (Kmn term) is shown in Appendix II. However, this solution 
has a typo in the 2012 edition. The 'correct' version is used here, which also happens
to match the version used to make the graphs.

For high `ka` values (>3) it can take a few minutes to run. A faster implementation
(but only Linux compatible) is the beamshapes.piston_in_sphere_flint. `piston_in_sphere_flint` 
is better suited for those that need to calculate 100-1000's of beamshapes.


References
----------
Beranek, L. L., & Mellow, T. (2012). Acoustics: sound fields and transducers.
Academic Press.


"""
import copy
from joblib import Parallel, delayed 
import mpmath
from mpmath import mpf
import numpy as np
import platform
import sympy
from sympy import symbols, legendre, sin, cos, tan, Sum, I, diff, pi, sqrt
from sympy import Matrix, Piecewise
from sympy import  lambdify, expand, Integral
from sympy import HadamardProduct as HP
import tqdm
x, alpha, index, k, m,n,p, r1, R, theta, y, z = symbols('x alpha index k m n p r1 R theta,y,z')
legendre_ncosalpha, legendre_mcosalpha = symbols('legendre_ncosalpha legendre_mcosalpha ')

dps = 50;
mpmath.mp.dps = dps

from beamshapes.special_functions import sph_hankel2
from beamshapes.utilities import args_to_mpmath, args_to_str

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
Imn_term_func = lambdify([m,n,k,R,alpha, theta], Imn_term, 'mpmath')
# Imn_func = lambdify([m,n,k,R,alpha],Imn,'mpmath') 

#%%
def Imn_func(mv,nv,kv,Rv,alphav):
    '''
    eqn. 12.106
    The 'gauss-legendre' quadrature method is used here as it provides 
    more accurate output, even with increasing m&n indices. 
    '''
    return mpmath.quad(lambda thetav: Imn_term_func(mv,nv,kv,Rv,alphav, thetav),
                       (0,alphav),
                       method='gauss-legendre')   

#%%
# equation 12.107
Kmn_expr = legendre(n, cos(theta))*legendre(m, cos(theta))*sin(theta) # the integrl of this expression 
# has a solution given in Appendix II, eqn 70
legendre_ncosalpha = legendre(n, cos(alpha))
legendre_mcosalpha = legendre_ncosalpha.subs({'n':m})


# when m != n
num_legendre_term1 = legendre_mcosalpha*diff(legendre_ncosalpha, alpha)
num_legendre_term2 = legendre_ncosalpha*diff(legendre_mcosalpha, alpha)
# Appendix II,70 has a typo in the 2012 edition. This is the correct version here.
eqn70_mnoteqn = sin(alpha)*(num_legendre_term1 - num_legendre_term2)/(n*(n+1)-m*(m+1))

# when m==n
# substituting 'j' for 'index' 
# because j is also used for sqrt(-1) all through the book!!
summn_funcn = legendre(index,cos(alpha))*(legendre(index,cos(alpha))*cos(alpha)-legendre(index+1,cos(alpha)))

meqn_sumterm = 2*Sum(summn_funcn, (index,1,m-1))
eqn70_meqn = (1+ cos(alpha)*legendre(m,cos(alpha))**2 + meqn_sumterm)/(2*m+1)
Kmn = Piecewise((eqn70_mnoteqn,m>n),
                (eqn70_mnoteqn,m<n),
                (eqn70_meqn,True), )
Kmn_func = lambdify([m,n,alpha],Kmn,'mpmath')

#%%
# equation 12.108
Lm_expr = legendre(m,cos(theta))*(r1**2/R**2)*tan(theta)
Lm = Integral(Lm_expr, (theta,0,alpha))
#Lm_func = lambdify([m, R, alpha], Lm, 'mpmath')
Lm_term = lambdify([m,R,alpha,theta], Lm_expr,'mpmath')

def Lm_func(mv,Rv,alphav):
    '''
    eqn. 12.106
    The 'gauss-legendre' quadrature method is used here as it provides 
    more accurate output, even with increasing m&n indices. 
    '''
    return mpmath.quad(lambda thetav: Lm_term(mv,Rv,alphav, thetav),
                       mpmath.linspace(0,alphav,2),
                       method='gauss-legendre')


# mpmath.mp.dps = 300
# errors = []
# nns =  [1,5,50]
# for eachm in nns:
#     out, err= Lm_func(eachm,0.01,mpmath.pi/3)
#     errors.append(err)
# print(errors)
# # what does the Lm term plot like? 
# #%% 
# lm_theta = lambda mv, thetav: Lm_term(mv, paramv['R'],paramv['alpha'], thetav)
# thet = mpmath.linspace(0, mpmath.pi/3,50)
# plt.figure()
# plt.plot(np.degrees(np.float32(thet)), [lm_theta(25, each) for each in thet])

#%%
# b matrix
b = -I*Lm
# b_func = lambdify([m,alpha], b,'mpmath') # eqn. 12.104

# new b_func 
def b_func(mv, Rv, alphav):
    '''
    b = -I*Lm
    '''
    Lm_value = Lm_func(mv,Rv,alphav) # Integral wrt theta
    return -mpmath.j*Lm_value

#%% 
# Setting up the matrices 
# %% 
mmn_hankels = n*sph_hankel2.subs({'n':n-1,'z':k*R})-(n+1)*sph_hankel2.subs({'n':n+1,'z':k*R})
mmn_hankels_func = lambdify([n,k,R], mmn_hankels,'mpmath')

def calc_N(params):
    '''
    Decides the number of terms to run calculations for. 
    
    Keyword Arguments
    -----------------
    trend : 
    
    Returns 
    -------
    num_N : int
    
    '''
    num_N = params.get('trend', int(12+2*params['ka']/sin(params['alpha'])))
    return num_N

### parallel zone 
def calc_one_Mmn_term(**params):
    '''
    '''
    Imn_value = Imn_func(params['m'], params['n'],params['k'],
                         params['R'],params['alpha'])
    Kmn_value = Kmn_func(int(params['m']),int(params['n']),params['alpha'])
    numerator_hankels = mmn_hankels_func(params['n'],params['k'],params['R'])
    numerator = Imn_value+ numerator_hankels*Kmn_value
    denom = 2*params['n']+1
    return numerator/denom

def parallel_calc_one_Mmn_term(**args):
    mpmath_args = args_to_mpmath(**args)
    output = calc_one_Mmn_term(**mpmath_args)
    backto_str = str(output)
    return backto_str

def format_Mmn_to_matrix(string_list):
    '''
    Parameters
    ----------
    string_list : list
        List with string versions of mpmath.mpc numbers. 
    Returns
    -------
    Mmn_matrix: mpmath.matrix
        A matrix with nxn, where n= sqrt(list entries) shape.
    '''
    # convert all list entries to mpc
    mmn_entries = [mpmath.mpmathify(each) for each in string_list]
    # re-format into a matrix
    Nv = int(mpmath.sqrt(len(string_list)))
    Mmn_matrix = mpmath.matrix(Nv,Nv)
    indices = [(row,col) for row in range(Nv) for col in range(Nv)] 
    for (row, col), value in zip(indices, mmn_entries):
        Mmn_matrix[row,col] = value
    return Mmn_matrix

def compute_Mmn_parallel(params):
    '''
    '''
    params['ka'] = mpmath.fmul(params['k'], params['a'])
    Nv = calc_N(params)#12 + int(2*params['ka']/sin(params['alpha']))
    M_matrix = mpmath.matrix(Nv,Nv)

    
    # create multiple paramsets with changing m,n
    multi_paramsets = []
    for i in range(Nv):
        for j in range(Nv):
            this_paramset = copy.deepcopy(params)
            this_paramset['m'] = i
            this_paramset['n'] = j
            multi_paramsets.append(this_paramset)
            
    multi_paramset_str = [args_to_str(**each) for each in multi_paramsets]
    num_cores = int(params.get('num_cores',-1))
    print('miaow miaow')
    if platform.system()=='Windows':
        M_mn_out = Parallel(n_jobs=num_cores)(delayed(parallel_calc_one_Mmn_term, backend='loky')(**inputs) for inputs in tqdm.tqdm(multi_paramset_str))
    else:
        M_mn_out = Parallel(n_jobs=num_cores, backend='multiprocessing')(delayed(parallel_calc_one_Mmn_term)(**inputs) for inputs in tqdm.tqdm(multi_paramset_str))
    M_matrix = format_Mmn_to_matrix(M_mn_out)
    return M_matrix
#%%
#####

def compute_b(params):
    '''
    Keyword Arguments
    -----------------
    m
    R
    alpha
    
    Returns
    -------
    b_matrix : N x 1 
        The size of the matrix is given by the heuristic
        described in ```compute_M```
    See Also
    --------
    compute_M
    compute_a
    
    '''
    params['ka'] = mpmath.fmul(params['k'], params['a'])
    Nv = calc_N(params) #12 + int(2*params['ka']/sin(params['alpha']))
    b_matrix = mpmath.matrix(Nv,1)
    for each_m in range(Nv):
        b_matrix[each_m,:] = b_func(each_m, params['R'], params['alpha'])
    return b_matrix 
#%%
def compute_a(M_mat, b_mat):
    '''
    Keyword Arguments
    -----------------
    M_mat : N x N mpmath.matrix
    b_mat : N x 1 mpmath.matrix
    
    Returns
    -------
    a_matrix : N x 1 
        The size of the matrix is given by the heuristic
        described in ```compute_M```. This is the A_{n}
        used in the summations. 

    See Also
    --------
    compute_M
    compute_b
    '''
    a_matrix = mpmath.lu_solve(M_mat, b_mat)
    return a_matrix
#%%
def d_theta(angle,k_v,R_v,alpha_v,An):
    num = 4 
    N_v = An.rows
    denom  = (k_v**2)*(R_v**2)*mpmath.sin(alpha_v)**2
    part1 = num/denom
    jn_matrix = mpmath.matrix([I**f for f in range(N_v)])
    legendre_matrix = mpmath.matrix([legendre(n_v, mpmath.cos(angle)) for n_v in range(N_v)])
    Anjn = HP(An,jn_matrix).doit()
    part2_matrix = HP(Anjn,legendre_matrix).doit() 
    part2 = sum(part2_matrix)
    rel_level = lambdify([], -part1*part2, 'mpmath')
    return abs(rel_level())

def d_zero(k_v,R_v,alpha_v,An):
    num = 4
    N_v = An.rows
    denom  = (k_v**2)*(R_v**2)*mpmath.sin(alpha_v)**2
    part1 = num/denom
    jn_matrix = mpmath.matrix([I**f for f in range(N_v)])
    part2 = sum(HP(An, jn_matrix).doit())
    rel_level = lambdify([], -part1*part2, 'mpmath')
    return abs(rel_level())

def relative_directivity_db(angle,k_v,R_v,alpha_v,An):
    off_axis = d_theta(angle,k_v,R_v,alpha_v,An)
    on_axis = d_zero(k_v,R_v,alpha_v,An)
    rel_level = 20*mpmath.log10(np.abs(off_axis/on_axis))
    return rel_level

def piston_in_sphere_directivity(angles, params, **kwargs):
    '''
    Calculates relative directivity dB (D(theta)/D(0))
    of a piston in a rigid sphere.
    
    
    Parameters
    ----------
    angles : array-like
        Angles at which the directivity is to be calculated in radians. 
    params : dictionary
        Dictionary with at least the following keys:
            k : mpmath.mpf>0
                Wavenumber. 
            a : mpmath.mpf>0
                Radius of piston. See `Warning` below.
            R : mpmath.mpf>0
                Radius of sphere
            alpha: 0<mpmath.mpf<pi
                Half-angular aperture of piston. 
            num_cores : int, optional
                The number of cores to be used for the Mmn matrix computation. 
                Defaults to using all cores.
    
    Keyword Arguments
    -----------------
    A_n : optional, mpmath.matrix
        If not provided, then the A_n is calculated from scratch, 
        which can take time for certain `ka` values. 

    Returns 
    -------
    A_n : mpmath.matrix
        The solution of the M*b = A equation. 
        A_n contains all the terms used in calculating the directionalities 
        of the piston in a sphere.
    beamshape : np.array
        Array with dB(D_theta/D_0).
        The number of items is equal to the number of angles.

    Notes
    -----
    While the dictionary entries can be normal floats, it is best 
    to use mpmath.mpf float objects to maintain digit precision. 
    eg. when specifying an angle of 60degrees for the 
    aperture alpha, use ```mpmath.pi/3``` rather than ```np.pi/3```
    for instance.
    
    Warning
    -------
    The radius of the piston 'a' needs to follow the 
    relation: a = R*sin(alpha) - please make sure this relationship 
    can hold to have believable results!
    
    '''
    A_n = kwargs.get('A_n', None)
    if A_n is None:
        Mmatrix = compute_Mmn_parallel(params)
        bmatrix = compute_b(params)
        A_n = compute_a(Mmatrix, bmatrix)
        
    directivity = []
    for angle_v in angles:
        rel_dirnlty = relative_directivity_db(angle_v,
                                            params['k'],
                                            params['R'], 
                                            params['alpha'],
                                            A_n)
        directivity.append(rel_dirnlty)
    
    directivity = np.array(directivity,'float32')
    return A_n, directivity
