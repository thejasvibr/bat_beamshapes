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
from beamshapes.piston_in_sphere import piston_in_sphere_directivity

#%% 
    


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
        
        _, output_dirnlty = piston_in_sphere_directivity(angles,
                                                         self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        #kavals = [1,3,5,10]
        # max_error_allowed = [1,1,1,3.5]
        kavals = [1,3]
        max_error_allowed = [1,1]
        max_abs_errors = np.zeros(len(kavals)) # incorporates digitisation error
        for i, (each, allowed_error)  in enumerate(zip(kavals,max_error_allowed)):
            max_abs_errors[i] = self.perform_ka_match(each)
            print(max_abs_errors[i])
            self.assertTrue(max_abs_errors[i]<=allowed_error)

if __name__=='__main__':
    #unittest.main()
    #%% 
    import numpy as np 
    import matplotlib.pyplot as plt 
    #%% 
    plotdata = pd.read_csv('plots_data/pistoninsphere.csv')        
    by_ka = plotdata.groupby('ka')
    frequency = mpmath.mpf(50*10**3) # kHz
    vsound = mpmath.mpf(330) # m/s
    wavelength = vsound/frequency
    alpha_value = mpmath.pi/3 # 60 degrees --> pi/3
    k_value = 2*mpmath.pi/(wavelength)

    kaval = 3
    paramv = {}
    paramv['alpha'] = alpha_value
    paramv['k'] = k_value
    ka = mpmath.mpf(kaval)
    a_value = ka/paramv['k']
    R_value = a_value/mpmath.sin(paramv['alpha'])  # m
    paramv['R'] = R_value        
    paramv['a'] = a_value
    paramv['a'] = kaval/paramv['k']

    angles = np.radians(by_ka.get_group(kaval)['theta_deg'])
    actual_dirnlty = by_ka.get_group(kaval)['relonaxis_db'].to_numpy()
    
    _, output_dirnlty = piston_in_sphere_directivity(angles,
                                                                  paramv)
    error = np.abs(output_dirnlty-actual_dirnlty)
    error_deg = pd.DataFrame(data={'theta':angles, 
                                    'r':actual_dirnlty,
                                    'rout':output_dirnlty})
    error_deg = error_deg.sort_values(by='theta')
    error_deg['diff'] = np.abs(error_deg['r']-error_deg['rout'])
    error_deg['thetadeg'] = by_ka.get_group(kaval)['theta_deg']
    print(np.max(error))
    #%%   
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(angles, output_dirnlty, '*',label='beamshapes implementation')
    plt.plot(angles, actual_dirnlty, '*',label='textbook')
    plt.ylim(-40,10);plt.yticks([-20,-10,0]);plt.title(f'ka = {kaval}')
    plt.xticks(np.arange(0,2*np.pi,np.radians(15)))
    plt.legend()