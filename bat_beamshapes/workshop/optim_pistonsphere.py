"""
This example shows how to parallelise the Imn term by 
converting mpmath objects to strings --> allow them to 
be pickled - run the calculation, and then output the strings 
back -- which can then in term be converted to mpmath objects.

"""

#%%
import bat_beamshapes
import bat_beamshapes.piston_in_sphere as pins
import copy
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import mpmath
mpmath.mp.dps = 50
import numpy as np 
import tqdm

v_sound = mpmath.mpf(330.0)
frequency = mpmath.mpf(51000.0)
k = mpmath.fdiv(2*mpmath.pi, (v_sound/frequency)) # AKA wavenumber


a = mpmath.mpf(0.0065) # meters
ka = mpmath.fmul(k,a)
alpha = mpmath.pi/8.0 # report of gape angle of 90 degrees --> 45degrees half angle
R = mpmath.fdiv(a, mpmath.sin(alpha))

params = {}
params['R'] = R
params['alpha'] = alpha
params['k'] = k
#params['ka'] = ka
params['a'] = a


# %% 
#%%time 

# params['m'] = 10
# params['n'] = 10
# Imn_value = pins.Imn_func(params['m'], params['n'],params['k'], params['R'],params['alpha'])
# #%%
# %%time 
# Kmn_value = pins.Kmn_func(params['m'],params['n'],params['alpha'])
# #%%
# %%time 
# numerator_hankels = pins.mmn_hankels_func(j,params['k'],params['R'])

#%%
import sympy
from sympy import symbols, legendre, sin, cos, tan, summation, I, diff, pi, sqrt
from sympy import Matrix, besselj, bessely, Piecewise
from sympy import  lambdify, integrate, expand,Integral
import tqdm
from joblib import Parallel, delayed
import time

#%% 
from sympy import symbols, sin, lambdify
# x = symbols('x')

# hjh = lambdify([x], sin(x)*x**2, 'mpmath')

# def imn_joblib(**args):
#     mpmath_args = {key: mpmath.mpmathify(value) for key,value in args.items()}
#     output = pins.Imn_func(**mpmath_args)
#     backto_str = str(output)
#     return backto_str

# imn_vals = {'m':mpmath.mpf(params['m']),
#             'n':params['n'],
#             'k': params['k'],
#             'R': params['R'],
#             'alpha': params['alpha']}

def args_to_str(**args):
    return {key: str(value) for key,value in args.items()}

def args_to_mpmath(**args):
    return {key: mpmath.mpmathify(value) for key,value in args.items()}


# imn_strvals = args_to_str(**imn_vals)
# %% 
# non-parallelised

# %%time 
# for each in range(4):
#     pins.Imn_func(**imn_vals)
    
#%% 

# paramset = [imn_strvals]*4

# t_start = time.time()

# uty = Parallel(n_jobs=4, backend='multiprocessing')(delayed(imn_joblib)(**each) for each in paramset)
# print("With pickle from stdlib and wrapper: {:.3f}s"
#       .format(time.time() - t_start))

#%%

# %%time 
def calc_one_Mmn_term(**params):
    '''
    '''
    Imn_value = pins.Imn_func(params['m'], params['n'],params['k'],
                         params['R'],params['alpha'])
    Kmn_value = pins.Kmn_func(int(params['m']),int(params['n']),params['alpha'])
    numerator_hankels = pins.mmn_hankels_func(params['n'],params['k'],params['R'])
    numerator = Imn_value+ numerator_hankels*Kmn_value
    denom = 2*params['n']+1
    return numerator/denom

params_str = args_to_str(**params)

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
    Nv = 12 + int(2*params['ka']/sin(params['alpha']))
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
    M_mn_out = Parallel(n_jobs=-1, backend='multiprocessing')(delayed(parallel_calc_one_Mmn_term)(**inputs) for inputs in tqdm.tqdm(multi_paramset_str))
    M_matrix = format_Mmn_to_matrix(M_mn_out)
    return M_matrix
    

#%% 
params = {}

params['alpha'] = alpha
params['k'] = k
#params['ka'] = ka
params['a'] = a
params['m'] = 10
params['n'] = 10
paramst_str = args_to_str(**params)
paramst_mpm  = args_to_mpmath(**paramst_str)

#pins.Kmn_func(int(paramst_mpm['m']),int(paramst_mpm['n']),paramst_mpm['alpha'])
#pins.mmn_hankels_func(paramst_mpm['n'],paramst_mpm['k'],paramst_mpm['R'])
parallel_calc_one_Mmn_term(**paramst_str)

#%%
# This is what the function call should look like 



# Nv = 20
# multi_paramsets = []
# for i in range(Nv):
#     for j in range(Nv):
#         this_paramset = copy.deepcopy(params)
#         this_paramset['m'] = i
#         this_paramset['n'] = j
#         multi_paramsets.append(this_paramset)
#%%
%%time 

# multi_paramset_str = [args_to_str(**each) for each in multi_paramsets]
# M_mn_out = Parallel(n_jobs=-1, backend='multiprocessing')(delayed(parallel_calc_one_Mmn_term)(**inputs) for inputs in tqdm.tqdm(multi_paramset_str))
# M_mn_matrix = format_Mmn_to_matrix(M_mn_out)

compute_Mmn_parallel(params)





