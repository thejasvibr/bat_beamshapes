#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests to check point-source on a sphere

"""


import os 
try:
    os.chdir('beamshapes/tests')
except:
    pass


import unittest
import numpy as np 
import pandas as pd
from beamshapes import point_source_on_a_sphere_directivity
from beamshapes.point_source_on_a_sphere import d_zero_func 


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
        
        _, output_dirnlty = point_source_on_a_sphere_directivity(angles,
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

class OnAxisPointSourceOnSphere(unittest.TestCase):
    
    def setUp(self):
        self.R = 0.1
        self.plotdata = pd.read_csv('plots_data/onaxis_response_pointsourceonsphere.csv')
        self.kr_values = self.plotdata['kR'].tolist()
        self.tolerance = 0.1
    
    def test_onaxis(self):
        dzero_outputs = []
        for each_kr in self.kr_values:
            k_val = each_kr/self.R
            NN = int(10 + 2*each_kr)
            dzero_outputs.append(np.abs(d_zero_func(k_val, self.R, NN)))
        db_dzero = 20*np.log10(dzero_outputs)
        # check all values are < deviation
        deviation = db_dzero-self.plotdata['onaxis_db']
        less_than_threshold = np.all(np.abs(deviation)<self.tolerance)
        self.assertTrue(less_than_threshold)


if __name__=='__main__':
    unittest.main()        
    
    
    
        