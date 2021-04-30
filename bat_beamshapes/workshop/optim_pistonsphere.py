"""
This example shows how to parallelise the Imn term by 
converting mpmath objects to strings --> allow them to 
be pickled - run the calculation, and then output the strings 
back -- which can then in term be converted to mpmath objects.

"""

#%%
import bat_beamshapes
import bat_beamshapes.piston_in_sphere as pins
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import mpmath
mpmath.mp.dps = 80
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
params['ka'] = ka
params['a'] = a

#onlyparams = {[]}

aucalc = params['R']*mpmath.sin(params['alpha'])



# angles = mpmath.linspace(0,mpmath.pi,100)
# beamshape = bat_beamshapes.piston_in_sphere_directionality(angles, params)
# plt.figure()
# a0 = plt.subplot(111, projection='polar')
# plt.plot(angles, beamshape)
# plt.ylim(-40,0);plt.yticks(np.arange(-40,10,10))
# plt.xticks(np.arange(0,2*np.pi,np.pi/6))



# i,j = 3,3
# params['m'],params['n'] = i,j
# # pins.compute_Mmn(params)
# %% 
%%time 

params['m'] = 10
params['n'] = 10
Imn_value = pins.Imn_func(params['m'], params['n'],params['k'], params['R'],params['alpha'])
# #%%
# %%time 
# Kmn_value = pins.Kmn_func(params['m'],params['n'],params['alpha'])
# #%%
# %%time 
# numerator_hankels = pins.mmn_hankels_func(j,params['k'],params['R'])

# #%%
# numerator = Imn_value+ numerator_hankels*Kmn_value
# denom = 2*params['n']+1



# #%% 
# # replicting the joblib example 
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects
import time

# @delayed
# @wrap_non_picklable_objects
# def func_async_wrapped(i,*args): 
#     return 2 * i

# large_list = list(range(1000000))

# t_start = time.time()
# Parallel(n_jobs=2)(func_async_wrapped(21, large_list) for _ in range(1))
# print("With pickle from stdlib and wrapper: {:.3f}s"
#       .format(time.time() - t_start))

#%% 
from sympy import symbols, sin, lambdify
x = symbols('x')

hjh = lambdify([x], sin(x)*x**2, 'mpmath')

def imn_joblib(**args):
    mpmath_args = {key: mpmath.mpmathify(value) for key,value in args.items()}
    output = pins.Imn_func(**mpmath_args)
    backto_str = str(output)
    return backto_str
    


# def joblib_hjh_wrapper(xas_str):
#     inval_as_mpf = mpmath.mpf(xas_str)
#     mpmathout = hjh(inval_as_mpf)
#     return str(mpmathout)

# t_start = time.time()
# Parallel(n_jobs=2, backend='multiprocessing')(joblib_hjh_wrapper(each) for each in invals_str)
# print("With pickle from stdlib and wrapper: {:.3f}s"
#       .format(time.time() - t_start))


imn_vals = {'m':mpmath.mpf(params['m']),
            'n':params['n'],
            'k': params['k'],
            'R': params['R'],
            'alpha': params['alpha']}

def args_to_str(**args):
    return {key: str(value) for key,value in args.items()}

imn_strvals = args_to_str(**imn_vals)
# %% 
# non-parallelised

%%time 
for each in range(4):
    pins.Imn_func(**imn_vals)
    
#%% 

paramset = [imn_strvals]*4

t_start = time.time()

uty = Parallel(n_jobs=4, backend='multiprocessing')(delayed(imn_joblib)(**each) for each in paramset)
print("With pickle from stdlib and wrapper: {:.3f}s"
      .format(time.time() - t_start))


# #%%

#%%