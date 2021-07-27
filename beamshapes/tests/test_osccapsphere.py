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
import tqdm
from beamshapes import cap_in_sphere_directivity
from beamshapes.cap_in_sphere import d_zero

class CapOfSphere(unittest.TestCase):
    
    def setUp(self):
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/cap-sphere_fig12-27.csv')        
        self.by_ka = plotdata.groupby('ka')
    
    def perform_ka_match(self, kaval):
        self.paramv['k'] = 10
        self.paramv['R'] = 0.01
        self.paramv['alpha'] = mpmath.pi/3
        self.paramv['a'] = self.paramv['R']*mpmath.sin(self.paramv['alpha'])
        self.paramv['k'] = kaval/self.paramv['a']
        self.paramv['a'] = kaval/self.paramv['k']
        angles = np.radians(self.by_ka.get_group(kaval)['theta'])
        actual_dirnlty = self.by_ka.get_group(kaval)['r'].to_numpy()
        
        _, output_dirnlty = cap_in_sphere_directivity(angles,
                                                                      self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)
    
    def test_kamatches(self):
        kavals = [1, 3, 10, 30]
        max_abs_errors = np.zeros(len(kavals))
        for i, each in enumerate(kavals):
            max_abs_errors[i] = self.perform_ka_match(each)
        print(max_abs_errors)
        self.assertTrue(np.all(max_abs_errors<=1.05))


class OnAxisCapOfSphere(unittest.TestCase):
    '''
    Check the on-axis levels over different cap apertures and ka values
    '''
    
    def setUp(self):
        self.paramv = {}
        self.plotdata = pd.read_csv('plots_data/onaxis_capofsphere_v2.csv')        
        self.alpha_values = self.plotdata['alpha_deg'].unique()
        self.by_alpha = self.plotdata.groupby('alpha_deg')
        
    def test_onaxis_match(self):
        self.paramv['R'] = 0.1
        groundtruth_levels = self.plotdata['onaxis_db'].tolist()
        all_obtained_raw = []
        
        for alpha in self.alpha_values:
            subdf = self.by_alpha.get_group(alpha)
            this_alpha_set = []
            for kaval in tqdm.tqdm(subdf['ka']):
                self.paramv['alpha'] = np.radians(alpha)
                self.paramv['a'] = self.paramv['R']*mpmath.sin(self.paramv['alpha'])
                self.paramv['k'] = kaval/self.paramv['a']
                onaxis_value = d_zero(self.paramv['k'], self.paramv['R'], 
                                      self.paramv['alpha'])
                this_alpha_set.append(float(np.abs(onaxis_value)))

            this_alpha_normalised = np.array(this_alpha_set)
            #this_alpha_normalised *= 1/np.min(this_alpha_normalised)
            all_obtained_raw.append(this_alpha_normalised)
        
        all_obtained_raw =  np.concatenate(all_obtained_raw)
        all_obtained_db = 20*np.log10(all_obtained_raw)
        difference = all_obtained_db - groundtruth_levels

        self.assertTrue(np.max(np.abs(difference))<0.5)

if __name__=='__main__':
    unittest.main()

    # #%% 
    # import matplotlib.pyplot as plt
    
    # paramv = {}
    # paramv['R'] = 0.1
    # plotdata = pd.read_csv('plots_data/onaxis_capofsphere_v2.csv')
    # groundtruth_levels = plotdata['onaxis_db']
    # all_obtained = []
    # all_obtained_raw = {}
    # alpha_values = plotdata['alpha_deg'].unique()
    # by_alpha = plotdata.groupby('alpha_deg')
    # all_kas = {}
    # for alpha in alpha_values:
    #     this_alpha_set = []
    #     ka_values = by_alpha.get_group(alpha)['ka']
    #     all_kas[alpha] = ka_values
    #     for kaval in tqdm.tqdm(ka_values):
    #         paramv['alpha'] = np.radians(alpha)
    #         paramv['a'] = paramv['R']*mpmath.sin(paramv['alpha'])
    #         paramv['k'] = kaval/paramv['a']
    #         onaxis_value = d_zero(paramv['k'], paramv['R'], 
    #                               paramv['alpha'])
    #         this_alpha_set.append(float(np.abs(onaxis_value)))
    #     all_obtained_raw[alpha] = this_alpha_set

    # #%%
    # db = lambda X: 20*np.log10(X)
    
    # plt.figure()
    # for i, each in enumerate(alpha_values):

    #     subdf = by_alpha.get_group(each)
    #     plt.plot(subdf['ka'], subdf['onaxis_db'],'-*', label=f'{each}')
    #     plt.plot(all_kas[each], db(all_obtained_raw[each]), label=f'{each} -calculated')
    #     plt.xscale('log');
    #     plt.grid(which='both',axis='both');plt.title('ground truth')
    #     plt.xlim(0.1,100);plt.ylim(-30,7)
    #     plt.legend()

    
    
    
    
    
    
    
    
    
    