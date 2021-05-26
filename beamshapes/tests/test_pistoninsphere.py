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
from beamshapes import piston_in_sphere_directivity


class PistonInSphere(unittest.TestCase):
    '''
    The test data assumes that alpha = pi/3
    '''
    
    def setUp(self):
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/piston_in_sphere_fig12-23.csv')        
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
        angles = np.radians(self.by_ka.get_group(kaval)['angle_deg'])
        actual_dirnlty = self.by_ka.get_group(kaval)['relonaxis_db'].to_numpy()
        
        _, output_dirnlty = piston_in_sphere_directivity(angles,
                                                         self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        kavals = [1, 3, 5]
        max_abs_errors = np.zeros(len(kavals))
        for i, each in enumerate(kavals):
            max_abs_errors[i] = self.perform_ka_match(each)
        print(max_abs_errors)
        #self.assertTrue(np.all(max_abs_errors<=1))

if __name__=='__main__':
    unittest.main()
    # #%% 
    # import numpy as np 
    # paramv = {}
    # plotdata = pd.read_csv('plots_data/piston_in_sphere_fig12-23.csv')               
    # by_ka = plotdata.groupby('ka')
    # kaval = 10
    # paramv['k'] = 10
    # paramv['a'] = kRval/paramv['k']
    # angles = np.radians(by_ka.get_group(kRval)['theta'])
    # actual_dirnlty = by_ka.get_group(kRval)['r'].to_numpy()
    
    # _, output_dirnlty = point_source_on_a_sphere_directivity(angles,
    #                                                               paramv)
    # error = np.abs(output_dirnlty-actual_dirnlty)
    # error_deg = pd.DataFrame(data={'theta':angles, 
    #                                 'r':actual_dirnlty,
    #                                 'rout':output_dirnlty})
    # error_deg = error_deg.sort_values(by='theta')
    # error_deg['diff'] = np.abs(error_deg['r']-error_deg['rout'])
    # error_deg['thetadeg'] = by_ka.get_group(kRval)['theta']
    # print(np.max(error))
    # #%%
    
    # import matplotlib.pyplot as plt 
    # plt.figure()
    # a0 = plt.subplot(111, projection='polar')
    # plt.plot(angles, output_dirnlty, '*')
    # plt.plot(angles, actual_dirnlty, '*')
    # plt.ylim(-50,10)
    
        