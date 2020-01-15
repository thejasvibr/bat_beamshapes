"""
Test beamshapes predictions given model data.

Author: Thejasvi Beleyur, Acoustic and Functional Ecology,
        Max Planck Institute for Ornithology, Seewiesen

License : This code is released under an MIT License.
 Copyright 2020, Thejasvi Beleyur

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
  persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

def calculate_model_and_data_error(real_data, predictions, **kwargs):
    """ This function compares and calculates the error between model prediction
    and data.

    Parameters
    ----------
    real_data : pd.DataFrame with 2 columns
                theta: float. Angle in radians.
                obs_relative_pressure: float>0. Observed relative pressure in comparison to on-axis.

    predictions : pd.DataFrame with 2 columns.
                  theta: float. Angle in radians.
                  pred_relative_pressure: float>0. Predicted relative pressure in comparison to on-axis.

    Keyword Arguments
    -----------------
    error_function : function that calculates the mis/match between data and predictions.
                     Defaults to the sum of absolute error observed.

    Returns
    --------
    prediction_error : output format depends on the exact error function used.
    """

    prediction_error = kwargs.get('error_function', sum_absolute)(real_data, predictions)
    return(prediction_error)


def sum_absolute_error(real_data, predictions):
    """
    Calculates the absolute difference between the predictions and the observed data
    and outputs the sum.

    sum_absolute_error = Sum(real_data - prediction)


    """

    pass

def dbeam_by_dtheta(real_data, predictions):
    '''
    Calculates the first order derivative of the beamshape with reference to
    the angle of emission.

    If the overall error in dbeam/dtheta is low it means that the real data and
    predictions match well in their shape, but not so much in their exact values.

    Parameters
    ----------
    real_data
    predictions

    Returns
    -------
    error_dbeam_by_dtheta


    '''





