#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end tests to check if various beamshape models are generating the
predicted outputs. 

"""

import os 
try:
    os.chdir('bat_beamshapes/tests')
except:
    pass

import unittest
import numpy as np 
import pandas as pd
from bat_beamshapes import piston_in_infinite_baffle_directivity


class PistonInInfBaffle(unittest.TestCase):
    
    def setUp(self):
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/piston_in_infbaffle_fig13-5_2019edn.csv')        
        self.by_ka = plotdata.groupby('ka')
    
    def perform_ka_match(self, kaval):
        self.paramv['k'] = 10
        self.paramv['a'] = kaval/self.paramv['k']
        angles = np.radians(self.by_ka.get_group(kaval)['theta'])
        actual_dirnlty = self.by_ka.get_group(kaval)['r'].to_numpy()
        
        _, output_dirnlty = piston_in_infinite_baffle_directivity(angles,
                                                                     self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        kavals = [1, 3, 5, 10]
        max_abs_errors = np.zeros(len(kavals))
        for i, each in enumerate(kavals):
            max_abs_errors[i] = self.perform_ka_match(each)
        print(max_abs_errors)
        self.assertTrue(np.all(max_abs_errors<=1))

if __name__=='__main__':
    unittest.main()
    # #%% 
    # paramv = {}
    # plotdata = pd.read_csv('plots_data/piston_in_infbaffle_fig13-5_2019edn.csv')        
    # by_ka = plotdata.groupby('ka')
    # kaval = 5
    # paramv['k'] = 10
    # paramv['a'] = kaval/paramv['k']
    # angles = np.radians(by_ka.get_group(kaval)['theta'])
    # actual_dirnlty = by_ka.get_group(kaval)['r'].to_numpy()
    
    # _, output_dirnlty = piston_in_infinite_baffle_directivity(angles,
    #                                                               paramv)
    # error = np.abs(output_dirnlty-actual_dirnlty)
    # error_deg = pd.DataFrame(data={'theta':angles, 
    #                                 'r':actual_dirnlty,
    #                                 'rout':output_dirnlty})
    # error_deg = error_deg.sort_values(by='theta')
    # error_deg['diff'] = np.abs(error_deg['r']-error_deg['rout'])
    # error_deg['thetadeg'] = by_ka.get_group(kaval)['theta']
    # print(np.max(error))
    # #%%
    
    # import matplotlib.pyplot as plt 
    # plt.figure()
    # a0 = plt.subplot(111, projection='polar')
    # plt.plot(angles, output_dirnlty, '*')
    # plt.plot(angles, actual_dirnlty, '*')
    # plt.ylim(-50,10)
    
        