# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# End-to-end tests to check the piston in a sphere output. 

# Here I check that the symbolically calculated output from Sympy is the same
# as that expected from substitutions in Beranek & Mellow 2012
# """


import os 
try:
    os.chdir('beamshapes/tests')
except:
    pass

import copy
import unittest
from joblib import Parallel, delayed
import tqdm
import mpmath
import numpy as np 
import pandas as pd
from beamshapes.piston_in_sphere import piston_in_sphere_directivity, index
from beamshapes.piston_in_sphere import relative_directivity_db
from sympy import symbols, lambdify, legendre, sin, cos, diff, Sum, Piecewise
# #%% 
# alpha, m, n, theta, z = symbols('alpha m n theta z')

# #%%
# # equation 12.107
# Kmn_expr = legendre(n, cos(theta))*legendre(m, cos(theta))*sin(theta) # the integrl of this expression 
# # has a solution given in Appendix II, eqn 70
# legendre_1stderiv = diff(legendre(n,z),z).subs(z, cos(alpha))

# # when m != n -- The 'typo' version
# num_legendre_term1 = legendre(m,cos(alpha))*legendre_1stderiv
# num_legendre_term2 = legendre(n,cos(alpha))*legendre_1stderiv.subs({'n':m})
# eqn70_mnoteqn = sin(alpha)*(num_legendre_term1-num_legendre_term2)/(m*(m+1)-n*(n+1))

# # when m==n
# # substituting 'j' for 'index' 
# # because j is also used for sqrt(-1) all through the book!!
# summn_funcn = legendre(index,cos(alpha))*(legendre(index,cos(alpha))*cos(alpha)-legendre(index+1,cos(alpha)))

# meqn_sumterm = 2*Sum(summn_funcn, (index,1,m-1))
# eqn70_meqn = (1+ cos(alpha)*legendre(m,cos(alpha))**2 + meqn_sumterm)/(2*m+1)
# Kmn = Piecewise((eqn70_mnoteqn,m>n),
#                 (eqn70_mnoteqn,m<n),
#                 (eqn70_meqn,True), )
# Kmn_func3 = lambdify([m,n,alpha],Kmn,'mpmath')



# ### parallel zone 
# def calc_one_Mmn_term3(**params):
#     '''The 'typo' version
#     '''
#     Imn_value = Imn_func(params['m'], params['n'],params['k'],
#                          params['R'],params['alpha'])
#     Kmn_value = Kmn_func3(int(params['m']),int(params['n']),params['alpha'])
#     numerator_hankels = mmn_hankels_func(params['n'],params['k'],params['R'])
#     numerator = Imn_value+ numerator_hankels*Kmn_value
#     denom = 2*params['n']+1
#     return numerator/denom

# def parallel_calc_one_Mmn_term3(**args):
#     mpmath_args = args_to_mpmath(**args)
#     output = calc_one_Mmn_term3(**mpmath_args)
#     backto_str = str(output)
#     return backto_str

# def compute_Mmn_parallel3(params):
#     '''
#     The 'typo' version
#     '''
#     params['ka'] = mpmath.fmul(params['k'], params['a'])
#     Nv = calc_N(params)#12 + int(2*params['ka']/sin(params['alpha']))
#     M_matrix = mpmath.matrix(Nv,Nv)

    
#     # create multiple paramsets with changing m,n
#     multi_paramsets = []
#     for i in range(Nv):
#         for j in range(Nv):
#             this_paramset = copy.deepcopy(params)
#             this_paramset['m'] = i
#             this_paramset['n'] = j
#             multi_paramsets.append(this_paramset)
            
#     multi_paramset_str = [args_to_str(**each) for each in multi_paramsets]
#     num_cores = int(params.get('num_cores',-1))
#     M_mn_out = Parallel(n_jobs=num_cores, backend='multiprocessing')(delayed(parallel_calc_one_Mmn_term3)(**inputs) for inputs in tqdm.tqdm(multi_paramset_str))
#     M_matrix = format_Mmn_to_matrix(M_mn_out)
#     return M_matrix#%%


# def piston_in_sphere_directivity3(angles, params, **kwargs):
#     '''
#    The 'typo' version
    
#     '''
#     A_n = kwargs.get('A_n', None)
#     if A_n is None:
#         Mmatrix = compute_Mmn_parallel3(params)
#         bmatrix = compute_b(params)
#         A_n = compute_a(Mmatrix, bmatrix)
#     directivity = []
#     for angle_v in angles:
#         rel_dirnlty = relative_directivity_db(angle_v,
#                                             params['k'],
#                                             params['R'], 
#                                             params['alpha'],
#                                             A_n)
#         directivity.append(rel_dirnlty)
    
#     directivity = np.array(directivity,'float32')
#     return A_n, directivity



#%%
class PistonInSphereT2(unittest.TestCase):
    '''
    The test data assumes that alpha = pi/3
    '''
    
    def test_substitution_and_automatedversion(self):
        kavals = [1,3]
        
        frequency = mpmath.mpf(50*10**3) # kHz
        vsound = mpmath.mpf(330) # m/s
        wavelength = vsound/frequency
        alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
        k_value = 2*mpmath.pi/(wavelength)
        
        
        paramv = {}
        paramv['alpha'] = alpha_value
        paramv['k'] = k_value
        
        angles = mpmath.linspace(0, mpmath.pi, 20)
        
        
        for kaval in kavals:
            ka = mpmath.mpf(kaval)
            a_value = ka/paramv['k']
            R_value = a_value/mpmath.sin(paramv['alpha'])  # m
            paramv['R'] = R_value        
            paramv['a'] = a_value
            paramv['a'] = kaval/paramv['k']        
            _, output_dirnlty_actual = piston_in_sphere_directivity(angles,
                                                                      paramv)
            
            self.assertTrue(np.unique(output_dirnlty - output_dirnlty_actual)==0)

if __name__=='__main__':
    unittest.main()