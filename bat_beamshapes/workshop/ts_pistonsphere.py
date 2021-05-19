"""
Troubleshooting piston in a sphere
==================================
> For ka>3 the calculations show values that deviate from the 
textbook groundtruth by 2-5 dB -- which is not ignorable!!
> I suspected the problem may come from:
    * low mpmath.mp.dps (decimal places) -- changing from 50-400 had no
    effects, which is odd
    * low N (matrix size/number of terms calculcated) -- changing from 
    baseline of 12+f(ka) --> 15+f(ka) had no effect
    * the directivity calculations were done with numpy pre 5th may, 
    and then changed to mpmath backend --> no effect. 

> 2021-06-05: I now suspect the problem lies perhaps with the quadrature 
terms. What if the quadrature is not 'accurate' enough? Here I'll test this idea
    * Some points to support this idea. The default quadrature method behind
    mpmath.quad is the 'tanh-sinh' algorithm. Instead of directly lambdifying 
    the Imn term into a standard mpmath.quad function I 'manually' made a 
    quadrature function for it to manipulate the options. 
    * For dps 200. Using the default 'tanh-sinh' leads to an estimated error 
    of e-203, while using 'gauss-legendre' leads to an estimated error of e-382. 
    Perhaps this is where the error is arising from. I noticed the integration 
    error increases with increasing m,n values. Perhaps this is why for bigger ka's, 
    (ka>3), the predictions get messier than for small ka's? There is at least 
    a connection here. 

Results with Imn - gauss-legendre flavour
-----------------------------------------

* The Imn term is not really the issue -- HOWEVER, using gauss-legendre reduces
run-time by 1/2!!, which is fantastic. 


Quadrture of the Lm term
========================
The other point of errors is the Lm term. What about there, does the exact
quadrature algorithm make a difference here? 

Here too, the choice of quadrature algorithm makes a difference. 'tanh-sinh' gives
higher estimated error terms for any given dps.

|Algorithm  | Error | dps   |
|-----------|-------|-------|
| tanh-sinh | e-74  | 250   |
| tanh-sinh | e-304 | 300   |
| gauss-leg.| e-76  | 250   |
| gauss-leg.| e-316 | 300   |

Broadly-speaking, the level of 'detail' (degree) that is applied during quadrature is proportional 
to the binary/decimal precision set for the code. For whatever reason, I couldn't
bring down the error term down by altering the 'maxdegree' argument of the 
`mpmath.quad` function. Either way, for now the estimated error has been brought
down. Hopefully this should do?


Results with Lm term
--------------------
Apparently the prediction-calculation mismatch doesn't have to do with the Lm term's
integration either!! The same errors appear. 


Is it a `dps` related issue? (2021-05-07 and 2021-05-10)
=========================================


| dps  | mismatch occurs? |
|------|------------------|
| 300  | yes              |
| 500  | yes              |
| 1000 | yes              | 


The mismatch still occurs no matter the actual `dps` - which is odd. 
The bottle-neck must thus be somewhere after the M and b matrix -- is it in 
the `a` term? Perhaps the inverse being calculated is using a fast-but-dirty approach?

Is it the calculation of the `a` matrix? 
========================================
Try lu_solve and qr_solve - one of them is more inaccurate than the other. 

Replace 'mpmath.inverse(M)'*b_matrix with mpmath.lu_solve(M,b) and mpmath.qr_solve
| Function used  | Mismatch occurs? |
|----------------|------------------|
| mpmath.inverse |  yes             |
| mpmath.lu_solve|  yes             |
| mpmath.qr_solve|  yes             |

??? -- the mpmath.mp.dps seems to get weirdly switched back to 50! during/after
the matrix solving -- don't know why!

This is becasue the dps set in the initial part of the module  is not consdiered. 
Setting the dps once more to 100-300 in the '__main__' part of the module
will switch it the the user's  preference



Does setting the D(theta) and D(0) to |D(theta)| and |D(0)| help?
=================================================================
In Tim Mellow's Mathematica code, the values are 'abs-ed' before 
division'  --- NO DIFFERENCE!

Try out mpmath.matrices.linalg -- improve_solution?
===================================================
importing improve_solution isn't so straightforward, perhaps on purpose. 
Needed to check the tests to see hwo to import it. 

>>> from mpmath import *
>>> new_an = mp.improve_solution(Mmn, An, bm, maxsteps=10)

but, the `new_an` had the same residual as the old An, and the function 
didn't make any difference.

Does changing the `N` help?
=============================
Changing the N equation from 12+2*() --> 12+3*() didn't change the plot much. 

There is some relation, and hopefully this will lead me to the source of the 
error. 

The details of the errors to the text-book plots are in 'effect_of_N-dps.csv'.
The patterns were a bit tough to follow through. 

A potential TYPO??
==================
Despite having checked my code and equations multiple times I kept wondering what
was wrong. I now suspect the difference in my results and those in the textbook
are because of a typo in the textbook. 

In eqn. 12.98, the sin(theta) term, should be replaced by a sin(theta)^2 term. 
I think there was an error in the substitution -- this also explains the
deviation with respect to the angle. The on-axis terms are expected to be
pretty similar as at theta around 0, sin(theta) and sin(theta)^2 may be
similar, but as theta increases, the difference increases too!

"""

#%%


import copy
from gmpy2 import *
from joblib import Parallel, delayed 
import mpmath
from mpmath import mpf
import numpy as np
from symengine import * 
import sympy
from sympy import symbols, legendre, sin, cos, tan, Sum, I, diff, pi, sqrt
from sympy import Matrix, Piecewise
from sympy import  lambdify, expand, Integral
from sympy import HadamardProduct as HP
import tqdm
x, alpha, index, k, m,n,p, r1, R, theta, y, z = symbols('x alpha index k m n p r1 R theta,y,z')
dps = 300;
mpmath.mp.dps = dps

from bat_beamshapes.special_functions import sph_hankel2
from bat_beamshapes.utilities import args_to_mpmath, args_to_str

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

# mpmath.mp.dps = 300
# errors = []
# nns =  [1,5,50]
# for eachm in tqdm.tqdm(nns):
#     for eachn in nns:
#         out, err= Imn_func(eachm,eachn, paramv['k'],paramv['R'],paramv['alpha'])
#         errors.append(err)
# print(errors)

# angles = mpmath.linspace(0, paramv['alpha'],100)
# imn_vals = [ (Imn_term_func(20,5,paramv['k'],paramv['R'],paramv['alpha'],ang)) for ang in angles]

# plt.figure()
# plt.plot(np.degrees(np.float32(angles)), [each.real for each in imn_vals], label='real')
# plt.plot(np.degrees(np.float32(angles)), [each.imag for each in imn_vals], label='imag')
# plt.legend()


#%% 
    

#%%
# equation 12.107
Kmn_expr = legendre(n, cos(theta))*legendre(m, cos(theta))*sin(theta) # the integrl of this expression 
# has a solution given in Appendix II, eqn 70
legendre_1stderiv = diff(legendre(n,z),z)

legendre_cosalpha = (n*(n+1)/((-sin(alpha)**2)*(2*n+1)))*(legendre(n+1,cos(alpha))-legendre(n-1, cos(alpha)))

# when m != n
num_legendre_term1 = legendre(m,cos(alpha))*legendre_cosalpha
num_legendre_term2 = legendre(n,cos(alpha))*legendre_cosalpha.subs({'n':m})
eqn70_mnoteqn = sin(alpha)*(num_legendre_term1-num_legendre_term2)/(m*(m+1)-n*(n+1))

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
    rel_level = 20*mpmath.log10(abs(off_axis/on_axis))
    return rel_level

def piston_in_sphere_directivity(angles, params):
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
            ka : mpmath.mpf>0
                Product of wavenumber and radius of piston.
            R : mpmath.mpf>0
                Radius of sphere
            alpha: 0<mpmath.mpf<pi
                Half-angular aperture of piston. 
            num_cores : int, optional
                The number of cores to be used for the Mmn matrix computation. 
                Defaults to using all cores.

    Returns 
    -------
    directivity : list
        List with relative directionalities in  dB. 
        The number of items is equal to the number of angles.

    Notes
    -----
    While the dictionary entries can be normal floats, it is best 
    to use mpmath.mpf float objects to maintain digit precision. 
    eg. when specifying an angle of 60degrees for the 
    aperture alpha, use ```mpmath.pi/3``` rather than ```np.pi/3```
    for instance.
        
    
    '''
    Mmatrix = compute_Mmn_parallel(params)
    bmatrix = compute_b(params)
    amatrix = compute_a(Mmatrix, bmatrix)
    
    # directivity = []
    # for angle_v in angles:
    #     directivity.append(relative_directivity_db(angle_v,
    #                                                      params['k'],
    #                                                      params['R'], 
    #                                                      params['alpha'],
    #                                                      amatrix))
    
    # return directivity
    return amatrix, Mmatrix, bmatrix

if __name__ == '__main__':
    #%% 
    import matplotlib.pyplot as plt
    mpmath.mp.dps = 100
    frequency = mpmath.mpf(25*10**3) # kHz
    vsound = mpmath.mpf(330) # m/s
    wavelength = vsound/frequency
    alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
    k_value = 2*mpmath.pi/(wavelength)
    ka_val = 10
    print(f'Starting piston in sphere for ka={ka_val}')
    ka = mpmath.mpf(ka_val)
    a_value = ka/k_value 
    R_value = a_value/mpmath.sin(alpha_value)  # m
    paramv = {}
    paramv['R'] = R_value
    paramv['alpha'] = alpha_value
    paramv['k'] = k_value
    paramv['a'] = a_value
    paramv['trend'] = int(10+2*paramv['k']*paramv['a']/mpmath.sin(paramv['alpha']))
    

    
    
    import pandas as pd
    # df = pd.read_csv('./ka5_piston_in_sphere.csv')
    df2 = pd.read_csv('../tests/piston_in_sphere_fig12-23.csv')
    ka5 = df2[df2['ka']==ka_val]
    
    #%%
    
    
    mpmath.mp.dps = 25

    angles = mpmath.matrix(np.radians(ka5['angle_deg'])) #mpmath.linspace(0,mpmath.pi,100)
    An, Mmn, bm = piston_in_sphere_directivity(angles, paramv)
    directivity = []
    dzero_value = d_zero(paramv['k'],paramv['R'],
                              paramv['alpha'], An)
    dtheta_values = []
    for angle_v in angles:
        dtheta_values.append(d_theta(angle_v, paramv['k'],paramv['R'],
                                          paramv['alpha'], An))
    
    beamshape = [20*mpmath.log10(abs(each/dzero_value)) for each in dtheta_values]

    beamshape_df = pd.DataFrame(data={'deg':np.degrees(np.float32(angles)),
                                      'd0dthet':np.float32(beamshape)})


    error = beamshape-ka5['relonaxis_db']
    mean_error = np.float(np.mean(np.abs(error)))
    max_error = np.max(np.float32(np.abs(error)))
    matrixsolve_error = mpmath.norm(mpmath.residual(Mmn, An, bm))
    
    #print(f'Decimal precision: {decimalprec} \n Ntrend: {Ntrend}')
    print(f'Matrix solve error : {matrixsolve_error }')
    print(f'Mean abs error: {mean_error}, max error: {max_error}\n')
    

# # %% 
#     an_candidate = mpmath.lu_solve(Mmn,bm)
#     res_mat = mpmath.residual(Mmn, an_candidate, bm)
#     res = mpmath.norm(res_mat)
#     print(f'matrix error is: {np.float64(res)}')

# %%     
    # an_candidate, res = mpmath.qr_solve(Mmn,bm)
    # print(res)
# #%% 
#     an_candidate = mpmath.inverse(Mmn)*bm
#     res_mat = mpmath.residual(Mmn, an_candidate, bm)
#     res = mpmath.norm(res_mat)
#     print(res)
# # %% Improving the solution iteratively 
    # from mpmath import *
    # new_Ax = mp.improve_solution(Mmn, An, bm, maxsteps=3)
    # new_res = mpmath.norm(mpmath.residual(Mmn, new_Ax, bm))
    # print(new_res)

# %%
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(angles, beamshape, '*',label='calculated')
    #plt.plot(angles, beamshape_nonpll, label='serial')
    plt.ylim(-40,0);plt.yticks(np.arange(-40,10,10))
    plt.xticks(np.arange(0,2*np.pi,np.pi/6))
    # load digitised textbook data
    plt.plot(angles, ka5['relonaxis_db'], 'g^', label='actual')
    plt.savefig(f'ka{ka_val}_pistoninasphere.png')
    plt.legend()
    # Also compare the error between prediction and textbook values
#%%
    plt.figure()
    plt.plot(np.degrees(np.float32(angles)), ka5['relonaxis_db'],'*',label='ground truth') # textbook
    plt.plot(np.degrees(np.float32(angles)), beamshape,'-*',label='calculated') # calculated
    plt.plot(np.degrees(np.float32(angles)), beamshape-ka5['relonaxis_db'],'-*',label='error') # relative error
    plt.yticks(np.arange(-36,4,2))
    plt.grid();plt.legend();plt.title('gauss-legendre')
# plt.savefig(f'ka{ka_val}_pistoninasphere_error.png')

    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(np.float32(angles), beamshape-ka5['relonaxis_db'])

