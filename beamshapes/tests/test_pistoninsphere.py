#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end tests to check the piston in a sphere output. 

Here the test isn't as straightforward as the substitution of Kmn's solution 
(eqn. 12.107) has a typo in the original Mathematica code. The typo is the 
result of a wrong substitution of P'n(z), where the sin(theta) term should be
-sin^2(theta).

Aside from the typo, the rest of the code seems to be fine. 
The only way to 'test' that this code is working correctly is to 
replicate the 'bug'/typo and check that the implementation matches the expected.

"""


import os 
try:
    os.chdir('beamshapes/tests')
except:
    pass


import unittest
import mpmath
import numpy as np 
import pandas as pd
from beamshapes.piston_in_sphere import *
del Kmn_func, calc_one_Mmn_term

#%% 
    

#%%
# equation 12.107
Kmn_expr = legendre(n, cos(theta))*legendre(m, cos(theta))*sin(theta) # the integrl of this expression 
# has a solution given in Appendix II, eqn 70
legendre_1stderiv = diff(legendre(n,z),z)

legendre_cosalpha = (n*(n+1)/((sin(alpha))*(2*n+1)))*(legendre(n+1,cos(alpha))-legendre(n-1, cos(alpha)))

# when m != n -- The 'typo' version
num_legendre_term1 = legendre(m,cos(alpha))*legendre_cosalpha
num_legendre_term2 = legendre(n,cos(alpha))*legendre_cosalpha.subs({'n':m})
eqn70_mnoteqn = sin(alpha)*(num_legendre_term2-num_legendre_term1)/(m*(m+1)-n*(n+1))

# when m==n
# substituting 'j' for 'index' 
# because j is also used for sqrt(-1) all through the book!!
summn_funcn = legendre(index,cos(alpha))*(legendre(index,cos(alpha))*cos(alpha)-legendre(index+1,cos(alpha)))

meqn_sumterm = 2*Sum(summn_funcn, (index,1,m-1))
eqn70_meqn = (1+ cos(alpha)*legendre(m,cos(alpha))**2 + meqn_sumterm)/(2*m+1)
Kmn = Piecewise((eqn70_mnoteqn,m>n),
                (eqn70_mnoteqn,m<n),
                (eqn70_meqn,True), )
Kmn_func2 = lambdify([m,n,alpha],Kmn,'mpmath')



### parallel zone 
def calc_one_Mmn_term2(**params):
    '''The 'typo' version
    '''
    Imn_value = Imn_func(params['m'], params['n'],params['k'],
                         params['R'],params['alpha'])
    Kmn_value = Kmn_func2(int(params['m']),int(params['n']),params['alpha'])
    numerator_hankels = mmn_hankels_func(params['n'],params['k'],params['R'])
    numerator = Imn_value+ numerator_hankels*Kmn_value
    denom = 2*params['n']+1
    return numerator/denom

def parallel_calc_one_Mmn_term2(**args):
    mpmath_args = args_to_mpmath(**args)
    output = calc_one_Mmn_term2(**mpmath_args)
    backto_str = str(output)
    return backto_str

def compute_Mmn_parallel2(params):
    '''
    The 'typo' version
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
    M_mn_out = Parallel(n_jobs=num_cores, backend='multiprocessing')(delayed(parallel_calc_one_Mmn_term2)(**inputs) for inputs in tqdm.tqdm(multi_paramset_str))
    M_matrix = format_Mmn_to_matrix(M_mn_out)
    return M_matrix#%%


def piston_in_sphere_directivity2(angles, params, **kwargs):
    '''
   The 'typo' version
    
    '''
    A_n = kwargs.get('A_n', None)
    if A_n is None:
        Mmatrix = compute_Mmn_parallel2(params)
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



#%%
class PistonInSphere(unittest.TestCase):
    '''
    The test data assumes that alpha = pi/3
    '''
    
    def setUp(self):
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/pistoninsphere.csv')        
        self.by_ka = plotdata.groupby('ka')
        frequency = mpmath.mpf(50*10**3) # kHz
        vsound = mpmath.mpf(330) # m/s
        wavelength = vsound/frequency
        alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
        k_value = 2*mpmath.pi/(wavelength)
        self.paramv = {}
        self.paramv['alpha'] = alpha_value
        self.paramv['k'] = k_value
    
        
    def perform_ka_match(self, kaval):
        #print(f'Starting piston in sphere for ka={ka_val}')
        ka = mpmath.mpf(kaval)
        a_value = ka/self.paramv['k']
        R_value = a_value/mpmath.sin(self.paramv['alpha'])  # m
        self.paramv['R'] = R_value        
        self.paramv['a'] = a_value
        self.paramv['a'] = kaval/self.paramv['k']
        angles = np.radians(self.by_ka.get_group(kaval)['theta_deg'])
        actual_dirnlty = self.by_ka.get_group(kaval)['relonaxis_db'].to_numpy()
        
        _, output_dirnlty = piston_in_sphere_directivity2(angles,
                                                         self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        kavals = [1,3,5,10]
        max_error_allowed = [1,1,1,3.5]
        max_abs_errors = np.zeros(len(kavals)) # incorporates digitisation error
        for i, (each, allowed_error)  in enumerate(zip(kavals,max_error_allowed)):
            max_abs_errors[i] = self.perform_ka_match(each)
            print(max_abs_errors[i])
            self.assertTrue(max_abs_errors[i]<=allowed_error)

if __name__=='__main__':
    unittest.main()
    # #%% 
    # import numpy as np 
    # import matplotlib.pyplot as plt 
    # #%% 
    # plotdata = pd.read_csv('plots_data/pistoninsphere.csv')        
    # by_ka = plotdata.groupby('ka')
    # frequency = mpmath.mpf(50*10**3) # kHz
    # vsound = mpmath.mpf(330) # m/s
    # wavelength = vsound/frequency
    # alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
    # k_value = 2*mpmath.pi/(wavelength)

    # kaval = 10
    # paramv = {}
    # paramv['alpha'] = alpha_value
    # paramv['k'] = k_value
    # ka = mpmath.mpf(kaval)
    # a_value = ka/paramv['k']
    # R_value = a_value/mpmath.sin(paramv['alpha'])  # m
    # paramv['R'] = R_value        
    # paramv['a'] = a_value
    # paramv['a'] = kaval/paramv['k']

    # angles = np.radians(by_ka.get_group(kaval)['theta_deg'])
    # actual_dirnlty = by_ka.get_group(kaval)['relonaxis_db'].to_numpy()
    
    # _, output_dirnlty = piston_in_sphere_directivity2(angles,
    #                                                               paramv)
    # error = np.abs(output_dirnlty-actual_dirnlty)
    # error_deg = pd.DataFrame(data={'theta':angles, 
    #                                 'r':actual_dirnlty,
    #                                 'rout':output_dirnlty})
    # error_deg = error_deg.sort_values(by='theta')
    # error_deg['diff'] = np.abs(error_deg['r']-error_deg['rout'])
    # error_deg['thetadeg'] = by_ka.get_group(kaval)['theta_deg']
    # print(np.max(error))
    # #%%   
    # plt.figure()
    # a0 = plt.subplot(111, projection='polar')
    # plt.plot(angles, output_dirnlty, '*',label='mine')
    # plt.plot(angles, actual_dirnlty, '*',label='textbook')
    # plt.ylim(-40,10);plt.yticks([-20,-10,0])
    # plt.xticks(np.arange(0,2*np.pi,np.radians(15)))
    # plt.legend()