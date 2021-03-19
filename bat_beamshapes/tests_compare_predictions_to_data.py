from beamshape_predictions import piston_in_infinite_baffle
from compare_predictions_to_data import sum_absolute_error, dbeam_by_dtheta_error

import unittest
import pandas as pd
import numpy as np

class TestSumAbsoluteError(unittest.TestCase):

    def setUp(self):
        angles = np.linspace(0+0.01, 2*np.pi, 100)
        self.beamshape = pd.DataFrame(columns=['theta','pred_relative_pressure'],
                                 index=range(angles.size))
        k = 2000
        a = 0.01
        for row, theta in enumerate(angles):
            self.beamshape['pred_relative_pressure'][row] = piston_in_infinite_baffle(theta,k,a=a)
            self.beamshape['theta'][row] = theta

        self.observed_beamshape = pd.DataFrame(columns=['theta', 'obs_relative_pressure'],
                                                           index=range(angles.size))

    def test_basic_sumabsolute_error(self):
        '''
        Check if a simple case with the same number of emission angles between
        predictions and real data matches up.

        In this test, the predictions and the observed are the *same* data.
        '''
        self.observed_beamshape['theta'] = self.beamshape['theta']
        self.observed_beamshape['obs_relative_pressure'] = self.beamshape['pred_relative_pressure']

        sum_abs_error  = sum_absolute_error(self.observed_beamshape, self.beamshape)
        self.assertEqual(sum_abs_error, 0.0)

    def test_1_sumabsolute_error(self):
        '''
        Observed data is +0.5  than the prediction all throughout
        '''
        self.observed_beamshape['theta'] = self.beamshape['theta']
        self.observed_beamshape['obs_relative_pressure'] = self.beamshape['pred_relative_pressure']
        self.observed_beamshape['obs_relative_pressure'] += 0.5
        predicted_sum_error = 0.5*len(self.beamshape['theta'])

        sum_abs_error  = sum_absolute_error(self.observed_beamshape, self.beamshape)
        self.assertEqual(sum_abs_error, predicted_sum_error)

if __name__ == '__main__':
    unittest.main()
