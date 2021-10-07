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
from beamshapes.piston_in_sphere import compute_Mmn_parallel, compute_a, compute_b
from beamshapes.piston_in_sphere import d_zero

#%%
class PistonInSphereMpmath(unittest.TestCase):
    '''
    The test data assumes that alpha = pi/3
    '''

    
    def setUp(self):
        # load the ground-truth data and parameters used to generate them
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
        # runs and outputs results for a parameter set across the azimuth
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
        
        kavals = [1,3,5]
        max_error_allowed = [1,1,1]
        max_abs_errors = np.zeros(len(kavals)) # incorporates digitisation error
        for i, (each, allowed_error)  in enumerate(zip(kavals,max_error_allowed)):
            max_abs_errors[i] = self.perform_ka_match(each)
            print(max_abs_errors[i])
            self.assertTrue(max_abs_errors[i]<=allowed_error) 
#%%
class PistonInSphereOnAxisResponse(unittest.TestCase):
    '''
    Run tests only with alpha=pi/3 because everything takes so long to run
    '''
    def setUp(self):
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/onaxis_respons_pistoninsphere.csv')        
        self.by_alpha = plotdata.groupby('alpha_deg')
        frequency = mpmath.mpf(50*10**3) # Hz
        vsound = mpmath.mpf(330) # m/s
        wavelength = vsound/frequency
        k_value = 2*mpmath.pi/(wavelength)
        self.alpha_values = plotdata['alpha_deg'].unique()
        self.paramv['k'] = k_value
    
    def dB(self, X):
        return 20*np.log10(np.abs())
    
    def compute_An(self, params):
        Mmatrix = compute_Mmn_parallel(params)
        bmatrix = compute_b(params)
        A_n = compute_a(Mmatrix, bmatrix)
        return A_n

    def compute_onaxis_level(self, params):
        An = self.compute_An(params)
        return d_zero(params['k'], params['R'], params['alpha'], An)
    
    def test_perform_onaxiscomparison(self):
        '''
        Check the onaxis response curve match to alpha=60deg
        as it has the 0 dB referene in its line!
        '''
        for alpha in self.alpha_values:
            subdf = self.by_alpha.get_group(alpha)
            ground_truth = subdf['normalised_onaxis_db']
            obtained = []
            for ka in subdf['ka']:
                print(ka, alpha)
                self.paramv['a'] = ka/self.paramv['k']
                self.paramv['alpha'] = mpmath.mpf(np.radians(alpha))
                self.paramv['R'] = self.paramv['a']/mpmath.sin(self.paramv['alpha'])  
                onaxis_level = self.compute_onaxis_level(self.paramv)
                obtained.append(float(onaxis_level))
        obtained_dB = 20*np.log10(np.array(obtained)/np.min(obtained))
        difference = obtained_dB - ground_truth
        self.assertTrue(np.max(np.abs(difference))<0.5) # set a 0.5 dB threshold 

if __name__=='__main__':
    unittest.main()