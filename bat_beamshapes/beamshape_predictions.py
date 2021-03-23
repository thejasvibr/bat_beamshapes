"""
Generate beamshape predictions given a model and a set of input parameters.

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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import jv as bessel_firstkind
from sympy import besselj, symbols, hankel2, legendre, sin, cos, summation, I, oo, diff, pi
from sympy import Abs, lambdify


def simple_cardiod(theta,**kwargs):
    '''
    From Giuggioli et al. 2015 : (A+2)(cos(theta)-1)

    Parameters
    ----------
    theta : np.array
        Angles in radians
    A : float>0, optional
        Asymmetry parameter. Defaults to 7.3 as in Giuggioli et al. 2015.
        
    Returns
    -------
    rel_onaxis

    References
    ----------
    Giuggioli, L., McKetterick, T. J., & Holderied, M. (2015).
    Delayed response and biosonar perception explain movement coordination
    in trawling bats. PLoS Comput Biol, 11(3), e1004089.
    '''
    rel_onaxis = (kwargs.get('A',7.3))*(np.cos(theta)-1)
    return rel_onaxis

def omnidirn_model(theta,**kwargs):
    '''

    Parameters
    ----------
    theta : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    return np.tile(0, theta.size)

def piston_in_infinite_baffle(theta, k, **kwargs):
    '''
    Generate predictions for the classic piston vibrating in an infinite baffle.

    Parameters
    ----------
    theta: 0<=float.<=2*pi
           The angle in azimuth/elevation at which the relative sound pressure is calculated.
           The relative soundpressure is calculated in comparison to the sound pressure level
           on-axis (0degrees/radians)
    k : float>0.
        The wavenumber of the sound (2*pi/wavelength)

    Keyword Arguments
    ------------------
    a : float>0.
        The piston radius.

    Returns
    --------
    rp_theta : float >=0.
               The relative sound pressure level at the given theta
    '''
    ka_sintheta = k*kwargs['a']*np.sin(theta)
    numerator = 2*besselj(1, ka_sintheta).evalf()#2*bessel_firstkind(1, ka_sintheta)

    try:
        rp_theta = np.abs(numerator/ka_sintheta)
    except:
        rp_theta = 0
    return(rp_theta)


##### Vibrating cap on sphere
Wn, Pn, h2n, h2n_prime, theta, w0, n, theta0, k, R, H_theta, kR = symbols(
        'Wn Pn h2n h2n_prime theta w0 n theta0 k R H_theta kR')
# take equation 6 to define Wn now
Wn = 0.5 * w0 * (legendre(n - 1, cos(theta0)) - legendre(n + 1, cos(theta0)))

h2n = hankel2(n, kR)
# h2n_prime_term = diff(hankel2(n, k*R), k*R)
h2n_prime_term = diff(h2n, kR)

# equation 23
common_term = ((I ** n) * Wn) / h2n_prime_term
denominator = summation(common_term, (n, 0, 20))

numerator_term = common_term * legendre(n, cos(theta))
numerator = summation(numerator_term, (n, 0, 20))

H_theta = Abs(numerator / denominator)

calc_H_theta = lambdify((theta, k, R, theta0, kR,w0), H_theta)

def vibrating_cap_of_sphere(theta, k, **kwargs):
    '''
    Parameters
    ----------
    theta: 0<=float<=2*pi
           The angle in azimuth/elevation at which the relative sound pressure is calculated.
           The relative soundpressure is calculated in comparison to the sound pressure level
           on-axis (0degrees/radians)

    k : float>0
        Wavenumber ( 2pi/wavelength)

    Keyword Arguments
    ------------------
    R : float>0
        Radius of the sphere

    theta_0 : float>0
              Aperture size of the vibrating cap in radians.

    w0 : float>0
        Defaults to 1 if unspecified.

    Returns
    --------
    H_theta : float>0
              The relative sound pressure in comparison to the on-axis.
    '''
    H_theta = calc_H_theta(theta, k, kwargs['R'], kwargs['theta_0'],
                           k*kwargs['R'], kwargs.get('w0',1))

    return(H_theta)


if __name__ == '__main__':
    thetas = np.arange(0.0, 2 * np.pi, 0.01)
    q = np.array([vibrating_cap_of_sphere(theta, 2000, R=0.01, theta_0=np.pi/8) for theta in thetas], dtype='float64')

    plt.figure()
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_zero_location('N')
    ax.plot(thetas, 20 * np.log10(q))


