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
from bat_beamshapes import point_source_on_a_sphere_directionality


class PointSourceSphere(unittest.TestCase):
    
    def setUp(self):
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/pointonsphere.csv')        
        self.by_ka = plotdata.groupby('kR')
    
    def perform_kR_match(self, kRval):
        self.paramv['k'] = 10
        self.paramv['R'] = kRval/self.paramv['k']
        angles = np.radians(self.by_ka.get_group(kRval)['theta'])
        actual_dirnlty = self.by_ka.get_group(kRval)['r'].to_numpy()
        
        _, output_dirnlty = point_source_on_a_sphere_directionality(angles,
                                                                     self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        kRvals = [1, 3, 5, 10]
        max_abs_errors = np.zeros(len(kRvals))
        for i, each in enumerate(kRvals):
            max_abs_errors[i] = self.perform_kR_match(each)
        print(max_abs_errors)
        self.assertTrue(np.all(max_abs_errors<=1))

if __name__=='__main__':
    unittest.main()
    # #%% 
    # import numpy as np 
    # paramv = {}
    # plotdata = pd.read_csv('plots_data/pointonsphere.csv')               
    # by_ka = plotdata.groupby('kR')
    # kRval = 10
    # paramv['k'] = 10
    # paramv['R'] = kRval/paramv['k']
    # angles = np.radians(by_ka.get_group(kRval)['theta'])
    # actual_dirnlty = by_ka.get_group(kRval)['r'].to_numpy()
    
    # _, output_dirnlty = point_source_on_a_sphere_directionality(angles,
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
    
        