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


class OnAxisCapOfSphere(unittest.TestCase):
    '''
    Check the on-axis levels over different cap apertures and ka values
    '''
    
    def setUp(self):
        # load the ground-truth data and parameters used for the results in the plots
        self.paramv = {}
        self.plotdata = pd.read_csv('plots_data/onaxis_capofsphere_v2.csv')        
        self.alpha_values = self.plotdata['alpha_deg'].unique()
        self.by_alpha = self.plotdata.groupby('alpha_deg')
        
    def test_onaxis_match(self):
        # set fixed radius
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
            all_obtained_raw.append(this_alpha_normalised)
        
        all_obtained_raw =  np.concatenate(all_obtained_raw)
        all_obtained_db = 20*np.log10(all_obtained_raw)
        difference = all_obtained_db - groundtruth_levels
        print(np.max(np.abs(difference)))
        self.assertTrue(np.max(np.abs(difference))<1)

if __name__=='__main__':
    unittest.main()
