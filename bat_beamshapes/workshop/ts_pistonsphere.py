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
    * the directionality calculations were done with numpy pre 5th may, 
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


TO FOLLOW UP
------------
> The other point of errors is the Lm term. What about there, does the exact
quadrature algorithm make a difference here? 

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
from sympy import symbols, legendre, sin, cos, tan, summation,Sum, I, diff, pi, sqrt
from sympy import Matrix, besselj, bessely, Piecewise
from sympy import  lambdify, integrate, expand,Integral
from sympy import HadamardProduct as HP
import tqdm
x, alpha, index, k, m,n,p, r1, R, theta, y, z = symbols('x alpha index k m n p r1 R theta,y,z')
dps = 100;
mpmath.mp.dps = dps

from bat_beamshapes.special_functions import sph_hankel2
from bat_beamshapes.utilities import args_to_mpmath, args_to_str

r1 = (R*cos(alpha))/cos(theta)

mpmath.quad

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

def Imn_func(mv,nv,kv,Rv,alphav):
    '''
    eqn. 12.106
    The 'gauss-legendre' quadrature method is used here as it provides 
    more accurate output, even with increasing m&n indices. 
    '''
    return mpmath.quad(lambda thetav: Imn_term_func(mv,nv,kv,Rv,alphav, thetav),
                       (0,alphav), method='gauss-legendre')




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

meqn_sumterm = 2*Sum(summn_funcn, (index,1,m-1))
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

def calc_N(params):
    return int(12 + 2*params['ka']/sin(params['alpha']))

def compute_Mmn(params):
    '''
    Obsolete function -- not in use -- here only for historical reference. 
    See computer_Mmn_parallel
    
    Keyword Arguments
    ----------
    alpha : 0<mpmath/float<pi
        half-aperture value in radians
    k : mpmath/float>0
        wavenumber
    a : mpmath/float>0
        Radius of piston 

    Returns
    -------
    M_matrix: mpmath.matrix
        NxN matrix with N defined by the heuristic formula
        12 + 2*ka/sin(alpha)
    
    See Also
    --------
    compute_b
    compute_a
    '''
    params['ka'] = mpmath.fmul(params['k'], params['a'])
    Nv = calc_N(params)# 12 + int(2*params['ka']/sin(params['alpha']))
    M_matrix = mpmath.matrix(Nv,Nv)
    params['R'] = mpmath.fdiv(params['a'], mpmath.sin(params['alpha']))
    
    
    for i in tqdm.trange(Nv):
        for j in range(Nv):
            params['m'],params['n'] = i,j
            Imn_value = Imn_func(params['m'], params['n'],params['k'],
                                 params['R'],params['alpha'])
            Kmn_value = Kmn_func(params['m'],params['n'],params['alpha'])
            numerator_hankels = mmn_hankels_func(j,params['k'],params['R'])
            numerator = Imn_value+ numerator_hankels*Kmn_value
            denom = 2*params['n']+1
            M_matrix[params['m'],params['n']] = numerator/denom
    return M_matrix

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
    params['R'] = mpmath.fdiv(params['a'], mpmath.sin(params['alpha']))
    
    
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

#####

def compute_b(params):
    '''
    Keyword Arguments
    -----------------
    k
    a 
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
        b_matrix[each_m,:] = b_func(each_m, params['alpha'])
    return b_matrix 

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
    a_matrix = mpmath.inverse(M_mat)*b_mat
    return a_matrix

def d_theta(angle,k_v,R_v,alpha_v,An):
    num = 4 
    N_v = An.rows
    denom  = (k_v**2)*(R_v**2)*mpmath.sin(alpha_v)**2
    part1 = num/denom
    #jn_matrix = np.array([1j**f for f in range(N_v)])
    jn_matrix = mpmath.matrix([I**f for f in range(N_v)])
    legendre_matrix = mpmath.matrix([legendre(n_v, mpmath.cos(angle)) for n_v in range(N_v)])
    #part2_matrix = np.column_stack((An, jn_matrix, legendre_matrix))
    Anjn = HP(An,jn_matrix).doit()
    part2_matrix = HP(Anjn,legendre_matrix).doit() 
    #part2 = np.sum(np.apply_along_axis(lambda X: X[0]*X[1]*X[2], 1, part2_matrix))
    part2 = sum(part2_matrix)
    rel_level = lambdify([], - part1*part2, 'mpmath')
    return rel_level()

def d_zero(k_v,R_v,alpha_v,An):
    num = 4 
    N_v = An.rows
    denom  = (k_v**2)*(R_v**2)*mpmath.sin(alpha_v)**2
    part1 = num/denom
    jn_matrix = mpmath.matrix([I**f for f in range(N_v)])
   
    #part2_matrix = np.column_stack((An, jn_matrix))
    #part2 = np.sum(np.apply_along_axis(lambda X: X[0]*X[1], 1, part2_matrix))
    part2 = sum(HP(An, jn_matrix).doit())
    rel_level = lambdify([], - part1*part2, 'mpmath')
    return rel_level()

def relative_directionality_db(angle,k_v,R_v,alpha_v,An):
    off_axis = d_theta(angle,k_v,R_v,alpha_v,An)
    on_axis = d_zero(k_v,R_v,alpha_v,An)
    rel_level = 20*mpmath.log10(abs(off_axis/on_axis))
    return rel_level

def piston_in_sphere_directionality(angles, params, parallel=False):
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
                Half-angular aperture of piston. 
            num_cores : 1>=int, optionals
                The number of cores to be used for the Mmn matrix computation. 
                Defaults to all cores if not specificed.

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
    Mmatrix = compute_Mmn_parallel(params)
    bmatrix = compute_b(params)
    amatrix = compute_a(Mmatrix, bmatrix)
    
    # directionality = []
    # for angle_v in angles:
    #     directionality.append(relative_directionality_db(angle_v,
    #                                                      params['k'],
    #                                                      params['R'], 
    #                                                      params['alpha'],
    #                                                      amatrix))
    
    # return directionality
    return amatrix

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    frequency = mpmath.mpf(50*10**3) # kHz
    vsound = mpmath.mpf(330) # m/s
    wavelength = vsound/frequency
    alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
    k_value = 2*mpmath.pi/(wavelength)
    ka_val = 5
    print(f'Starting piston in sphere for ka={ka_val}')
    ka = mpmath.mpf(ka_val)
    a_value = ka/k_value 
    R_value = a_value/mpmath.sin(alpha_value)  # m
    paramv = {}
    paramv['R'] = R_value
    paramv['alpha'] = alpha_value
    paramv['k'] = k_value
    paramv['a'] = a_value
    
    
    import pandas as pd
    df = pd.read_csv('./ka5_piston_in_sphere.csv')
    df2 = pd.read_csv('../tests/piston_in_sphere_fig12-23.csv')
    ka5 = df2[df2['ka']==ka_val]
    
    
    angles = mpmath.matrix(np.radians(ka5['angle_deg'])) #mpmath.linspace(0,mpmath.pi,100)
    An = piston_in_sphere_directionality(angles, paramv)
    #beamshape_nonpll = piston_in_sphere_directionality(angles, paramv, False)
    directionality = []
    dzero_value = d_zero(paramv['k'],paramv['R'],
                             paramv['alpha'], An)
    dtheta_values = []
    for angle_v in angles:
        dtheta_values.append(d_theta(angle_v, paramv['k'],paramv['R'],
                                         paramv['alpha'], An))
    
    beamshape = [20*mpmath.log10(abs(each/dzero_value)) for each in dtheta_values]
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(angles, beamshape, '-*',label='calculated')
    #plt.plot(angles, beamshape_nonpll, label='serial')
    plt.ylim(-40,0);plt.yticks(np.arange(-40,10,10))
    plt.xticks(np.arange(0,2*np.pi,np.pi/6))
    # load digitised textbook data
    plt.plot(angles, ka5['relonaxis_db'], '*', label='actual')
    plt.savefig(f'ka{ka_val}_pistoninasphere.png')
    # Also compare the error between prediction and textbook values
    plt.figure()
    plt.plot(angles, ka5['relonaxis_db'],'-',label='ground truth') # textbook
    plt.plot(angles, beamshape,'-*',label='calculated') # calculated
    plt.plot(angles, beamshape-ka5['relonaxis_db'],'-*',label='error') # relative error
    plt.yticks(np.arange(-36,4,2))
    plt.grid();plt.legend();plt.title('gauss-legendre')
# plt.savefig(f'ka{ka_val}_pistoninasphere_error.png')

    error = beamshape-ka5['relonaxis_db']
    median_error = np.median(np.abs(error))
    avg_error = np.mean(np.abs(error))
    rms_error = np.sqrt(np.mean(np.square(error)))
    print(median_error, avg_error, rms_error)
