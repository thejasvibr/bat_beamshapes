#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end tests to check the piston in a sphere output. 

"""


import os 
try:
    os.chdir('beamshapes/tests')
except:
    pass


import unittest
# installing FLINT with Travis-CI is tricky, and so avoid doing it for now. 
try:
    import flint
    import numpy as np 
    import pandas as pd
    from beamshapes.piston_in_sphere_flint import piston_in_sphere_directivity
    
    cos = flint.acb.cos
    sin = flint.acb.sin
    tan = flint.acb.tan
    pi = flint.acb.pi()
    arb = flint.arb
    acb = flint.acb
    
    
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
            frequency = acb(50*10**3) # kHz
            vsound = acb(330) # m/s
            wavelength = vsound/frequency
            alpha_value = pi/3 # 60 degrees --> pi/3
            k_value = 2*pi/(wavelength)
            self.paramv = {}
            self.paramv['alpha'] = alpha_value
            self.paramv['k'] = k_value
        
            
        def perform_ka_match(self, kaval):
            print(f'Starting piston in sphere for ka={kaval}')
            ka = acb(kaval)
            a_value = ka/self.paramv['k']
            R_value = a_value/sin(self.paramv['alpha'])  # m
            self.paramv['R'] = R_value        
            self.paramv['a'] = a_value
            self.paramv['a'] = ka/self.paramv['k']
            angles = np.radians(self.by_ka.get_group(kaval)['theta_deg'])
            angles_acb = [acb(each) for each in angles]
            actual_dirnlty = self.by_ka.get_group(kaval)['relonaxis_db'].to_numpy()
            
            _, output_dirnlty = piston_in_sphere_directivity(angles_acb,
                                                             self.paramv)
            error = np.abs(output_dirnlty-actual_dirnlty)
            return np.max(error)
    
        def test_kamatches(self):
            kavals = [3,10]
            max_error_allowed = [1,1.3]
            max_abs_errors = np.zeros(len(kavals)) # incorporates digitisation error
            for i, (each, allowed_error)  in enumerate(zip(kavals,max_error_allowed)):
                max_abs_errors[i] = self.perform_ka_match(each)
                print(max_abs_errors[i])
                self.assertTrue(max_abs_errors[i]<=allowed_error)
except:
    pass

if __name__=='__main__':
    unittest.main()