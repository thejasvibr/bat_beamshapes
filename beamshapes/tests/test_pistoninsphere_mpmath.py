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
import mpmath
import numpy as np 
import pandas as pd
from beamshapes.piston_in_sphere import piston_in_sphere_directivity

#%%
class PistonInSphereMpmath(unittest.TestCase):
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

        ka = kaval
        a_value = ka/self.paramv['k']
        R_value = a_value/mpmath.sin(self.paramv['alpha'])  # m
        self.paramv['R'] = R_value        
        self.paramv['a'] = a_value
        self.paramv['a'] = ka/self.paramv['k']
        angles = np.radians(self.by_ka.get_group(kaval)['theta_deg'])
        actual_dirnlty = self.by_ka.get_group(kaval)['relonaxis_db'].to_numpy()
        
        _, output_dirnlty = piston_in_sphere_directivity(angles,
                                                         self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        #kavals = [1,3,5,10]
        kavals = [1,3,5]
        max_error_allowed = [1,1,1]
        max_abs_errors = np.zeros(len(kavals)) # incorporates digitisation error
        for i, (each, allowed_error)  in enumerate(zip(kavals,max_error_allowed)):
            max_abs_errors[i] = self.perform_ka_match(each)
            print(max_abs_errors[i])
            self.assertTrue(max_abs_errors[i]<=allowed_error) 

if __name__=='__main__':
    unittest.main()