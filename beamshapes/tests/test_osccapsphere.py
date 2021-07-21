#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end tests to check if various beamshape models are generating the
predicted outputs. 

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
from beamshapes import cap_in_sphere_directivity
from beamshapes.cap_in_sphere import d_zero

#%%

# class CapOfSphere(unittest.TestCase):
    
#     def setUp(self):
#         self.paramv = {}
#         plotdata = pd.read_csv('plots_data/cap-sphere_fig12-27.csv')        
#         self.by_ka = plotdata.groupby('ka')
    
#     def perform_ka_match(self, kaval):
#         self.paramv['k'] = 10
#         self.paramv['R'] = 0.01
#         self.paramv['alpha'] = mpmath.pi/3
#         self.paramv['a'] = self.paramv['R']*mpmath.sin(self.paramv['alpha'])
#         self.paramv['k'] = kaval/self.paramv['a']
#         self.paramv['a'] = kaval/self.paramv['k']
#         angles = np.radians(self.by_ka.get_group(kaval)['theta'])
#         actual_dirnlty = self.by_ka.get_group(kaval)['r'].to_numpy()
        
#         _, output_dirnlty = cap_in_sphere_directivity(angles,
#                                                                      self.paramv)
#         error = np.abs(output_dirnlty-actual_dirnlty)
#         return np.max(error)
    
#     def test_kamatches(self):
#         kavals = [1, 3, 10, 30]
#         max_abs_errors = np.zeros(len(kavals))
#         for i, each in enumerate(kavals):
#             max_abs_errors[i] = self.perform_ka_match(each)
#         print(max_abs_errors)
#         self.assertTrue(np.all(max_abs_errors<=1.05))


class OnAxisCapOfSphere(unittest.TestCase):
    '''
    Check the on-axis levels over different cap apertures and ka values
    '''
    
    def setUp(self):
        self.paramv = {}
        self.plotdata = pd.read_csv('plots_data/onaxis_response_capofsphere.csv')        
        self.alpha_values = self.plotdata['alpha_deg'].unique()
        self.by_alpha = self.plotdata.groupby('alpha_deg')
        
    def test_onaxis_match(self):
        self.paramv['R'] = 0.1
        groundtruth_levels = self.plotdata['normalised_onaxis_db']
        all_obtained = []
        
        for alpha in self.alpha_values:
            subdf = self.by_alpha.get_group(alpha)
            this_alpha_set = []
            for kaval in subdf['ka']:
                self.paramv['alpha'] = np.radians(alpha)
                self.paramv['k'] = kaval/self.paramv['R']
                onaxis_value = d_zero(self.paramv['k'], self.paramv['R'], 
                                      self.paramv['alpha'])
                this_alpha_set.append(float(np.abs(onaxis_value)))

            this_alpha_normalised = np.array(this_alpha_set)
            this_alpha_normalised *= 1/np.min(this_alpha_normalised)
            all_obtained_raw.append(this_alpha_normalised)

        all_obtained_raw =  np.concatenate(all_obtained_raw)
        
        all_obtained = all_obtained_raw.copy()
        all_obtained_db = 20*np.log10(all_obtained)
        difference = all_obtained_db - groundtruth_levels

        print(difference.tolist())
        print('\n', all_obtained_db)
        print('\n', groundtruth_levels)
        
        df = pd.DataFrame(data={'groundtruth_db':groundtruth_levels,
                                'obtained_raw': all_obtained_raw,
                                'obtained_db': all_obtained_db,
                                'alpha': self.plotdata['alpha_deg'],
                                'ka': self.plotdata['ka']})
        df.to_csv('miaow_cap.csv')        
        self.assertTrue(np.max(np.abs(difference))<1)

if __name__=='__main__':
    unittest.main()
    #%% 
    # paramv = {}
    # plotdata = pd.read_csv('plots_data/cap-sphere_fig12-27.csv')        
    # by_ka = plotdata.groupby('ka')
    # kaval = 1
    # paramv['R'] = 0.01
    # paramv['alpha'] = np.pi/3
    # paramv['a'] = paramv['R']*np.sin(paramv['alpha'])
    # paramv['k'] = kaval/paramv['a']
    
    # angles = np.radians(by_ka.get_group(kaval)['theta'])
    # actual_dirnlty = by_ka.get_group(kaval)['r'].to_numpy()
    
    # _, output_dirnlty = cap_in_sphere_directivity(angles,
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
    
        