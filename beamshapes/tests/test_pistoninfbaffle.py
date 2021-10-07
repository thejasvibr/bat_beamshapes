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
import numpy as np 
import pandas as pd
from beamshapes import piston_in_infinite_baffle_directivity


class PistonInInfBaffle(unittest.TestCase):

    def setUp(self):
        # load the ground-truth data and the parameters used to generate them
        self.paramv = {}
        plotdata = pd.read_csv('plots_data/piston_in_infbaffle_fig13-5_2019edn.csv')        
        self.by_ka = plotdata.groupby('ka')
    
    def perform_ka_match(self, kaval):
        '''
        Function to actually generate the radiation pattern given a ka value
        '''
        self.paramv['k'] = 10
        self.paramv['a'] = kaval/self.paramv['k']
        angles = np.radians(self.by_ka.get_group(kaval)['theta'])
        actual_dirnlty = self.by_ka.get_group(kaval)['r'].to_numpy()
        
        _, output_dirnlty = piston_in_infinite_baffle_directivity(angles,
                                                                     self.paramv)
        error = np.abs(output_dirnlty-actual_dirnlty)
        return np.max(error)

    def test_kamatches(self):
        # check the match across different ka
        kavals = [1, 3, 5, 10]
        max_abs_errors = np.zeros(len(kavals))
        for i, each in enumerate(kavals):
            max_abs_errors[i] = self.perform_ka_match(each)
        print(max_abs_errors)
        self.assertTrue(np.all(max_abs_errors<=1))

if __name__=='__main__':
    unittest.main()