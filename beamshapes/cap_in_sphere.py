"""
Oscillating cap of a sphere
============================

References
----------
Beranek, L. L., & Mellow, T. (2012). Acoustics: sound fields and transducers.
Academic Press.

Notes
-----
Experience has shown it's best to leave the dps un-specified while using this
module. I guess this is because of the summation terms.

To Do
-----
* Implement the special case of alpha = pi/2, described by equations 12.60
"""
import copy
from joblib import Parallel, delayed
import numpy as np
import mpmath
from sympy import symbols,  I, cos, sin, legendre
from sympy import lambdify
import tqdm
import warnings
from beamshapes.special_functions import sph_hankel2, legendre_mvz
from beamshapes.utilities import args_to_mpmath, args_to_str

n, z, k, R, alpha, theta = symbols('n z k R alpha theta')

def sphhankel2(nv,zv):
    return sph_hankel2.subs({'n': nv, 'z': zv})

# equation 12.59
# split the big parenthesis into 3 parts (j/sph_hankel, cos(theta) term and the summation)
d_theta_term1 = I/(2*sphhankel2(1,k*R))

d_theta_term2_num = 3*(1-cos(alpha)**3)*cos(theta)
d_theta_term2_denom = (sin(alpha)**2)*(sphhankel2(0,k*R)-2*sphhankel2(2,k*R))
d_theta_term2 = d_theta_term2_num/d_theta_term2_denom
P_1ncosalpha = legendre_mvz.subs({'m': 1, 'v': n, 'z': cos(alpha)}).doit()

# the summation term 
dtheta_t3_num = (I**(n+1))*((2*n+1)**2)*(sin(alpha)*legendre(n, cos(alpha)) + cos(alpha)*P_1ncosalpha)
dtheta_t3_denom = (n-1)*(n+2)*sin(alpha)*(n*sphhankel2( n-1,k*R)-(n+1)*sphhankel2(n+1,k*R))
dtheta_t3_oneterm = (dtheta_t3_num/dtheta_t3_denom)*legendre(n, cos(theta))

d_theta_t1_func = lambdify([k, R, alpha, theta], d_theta_term1, 'mpmath')
d_theta_t2_func = lambdify([k, R, alpha, theta], d_theta_term2, 'mpmath')

#  the problem is here that the summation doesn't happen as
# 'n' is somehow not treated as an internal integer variable by mpmath
# the alternative is to then to create a function out of it and then
# perform summation through mpmath.


def dtheta_t3_func(k_v, R_v, alpha_v, theta_v):
    '''
    '''
    version_with_freen = dtheta_t3_oneterm.subs({'k': k_v,
                                                 'R': R_v,
                                                 'alpha': alpha_v,
                                                 'theta': theta_v})
    freen_func = lambdify([n], version_with_freen, 'mpmath')
    return mpmath.nsum(freen_func, [2, mpmath.inf])
    #return mpmath.nsum(freen_func, [2, int(5+2*k_v*R_v)])


#  In eqn. 12.61, only term 2 differs by the absence of a cos(theta)
d_0_term2_num = 3*(1-cos(alpha)**3)
d_0_term2 = d_0_term2_num/d_theta_term2_denom
d_0_t2_func = lambdify([k, R, alpha, theta], d_0_term2, 'mpmath')


#  term 3 of eqn. 12.61 doesn't have a Pn(cos(theta)) at the end.
d_0_t3_oneterm = dtheta_t3_num/dtheta_t3_denom


def d_0_t3_func(k_v, R_v, alpha_v, theta_v):
    
    version_with_freen = d_0_t3_oneterm.subs({'k': k_v, 'R': R_v,
                                              'alpha': alpha_v,
                                              'theta': theta_v})
    return mpmath.nsum(lambdify([n],version_with_freen, 'mpmath'), [2,mpmath.inf])
    #return mpmath.nsum(lambdify([n],version_with_freen, 'mpmath'), [2, int(5+2*k_v*R_v)])

def d_theta(**param):
    '''Calculates off-axis levels for cap in a sphere. 
    
    Parameters
    ----------
    param : dictionary
        With the following keys and entries:
        k : mpmath.mpf      
            Wavenumber
        R : mpmath.mpf
            Sphere radius
        alpha: mpmath.mpf
            Cap half-aperture in radians
        theta : mpmath.mpf  
            Azimuth/elevation at which level is to be calculated
            in radians.

    Returns
    -------
    final_d_theta: mpmath.mpf
        Calculated level at angle :math:`\theta`
    '''
    brackets_term1 = d_theta_t1_func(param['k'], param['R'],
                                     param['alpha'], param['theta'])
    brackets_term2 = d_theta_t2_func(param['k'], param['R'],
                                     param['alpha'], param['theta'])
    brackets_term3 = dtheta_t3_func(param['k'], param['R'],
                                         param['alpha'], param['theta'])
    kR_square = (param['k']*param['R'])**2
    final_d_theta = -(2/kR_square)*(brackets_term1+brackets_term2+brackets_term3)
    return final_d_theta


def d_theta_pll(**param):
    mpmath_args = args_to_mpmath(**param)
    output = d_theta(**mpmath_args)
    str_output = str(output)
    return str_output


def d_zero(kv, Rv, alphav, thetav=0):
    '''Calculates the on-axis level :math:`D_{0}` for cap in sphere.
    
    Parameters
    ----------
    kv : mpmath.mpf 
        Wavenumber
    Rv : mpmath.mpf
        Sphere radius
    alphav : mpmath.mpf
        Cap aperture.

    Returns
    -------
    final_d_0 : mpmath.mpf  
        On-axis level :math:`D_{0}`
    '''
    brackets_term1 = d_theta_t1_func(kv, Rv, alphav, thetav)
    brackets_term2 = d_0_t2_func(kv, Rv, alphav, thetav)
    brackets_term3 = d_0_t3_func(kv, Rv, alphav, thetav)
    kR_square = (kv*Rv)**2
    final_d_0 = -(2/kR_square)*(brackets_term1+brackets_term2+brackets_term3)
    return final_d_0

#
def relative_directivity_db(angle, k_v, R_v, alpha_v):
    '''
    Calculates the :math:`D_{\\theta}/D_{0}` in dB for all angles
    in the cap in a sphere model.
    '''
    off_axis = d_theta(k_v, R_v, alpha_v, angle)
    on_axis = d_zero(k_v, R_v, alpha_v)
    rel_level = 20*mpmath.log10(abs(off_axis/on_axis))
    return rel_level


def cap_in_sphere_directivity(angles, params):
    '''
    Calculates relative directivity dB (D(theta)/D(0))
    of an oscillating cap in a rigid sphere.

    Parameters
    ----------
    angles : array-like
        Angles at which the directivity is to be calculated in radians.
    params : dictionary
        Dictionary with at least the following keys:
            k : mpmath.mpf>0
                Wavenumber.
            R : mpmath.mpf>0
                Radius of sphere
            alpha: 0<mpmath.mpf<pi/2
                Half-angular aperture of cap.

    
    Returns
    -------
    _ : None
    directivity : np.array
        Array with relative directionalities in (20log10) dB.
        The number of items is equal to the number of angles.


    Notes
    -----
    While the dictionary entries can be normal floats, it is best
    to use mpmath.mpf float objects to maintain digit precision.
    eg. when specifying an angle of 60degrees for the
    aperture alpha, use ```mpmath.pi/3``` rather than ```np.pi/3```
    for instance.

    Warning
    -------
    The validity of the results at alpha>=pi/2
    haven't been checked. Proceed after verifying  it yourself. The function
    will throw a warning and give an output nonetheless. On pg. 506, the authors
    write 'When alpha = 1/2pi, the second term simplifies to that for an oscillating
    sphere, ...so that Eq. (12.59) simplifies to' and proceed to give Equation 12.6.

    It seems like 12.59 still holds, and at alpha = pi/2 it changes, but this
    hasn't been checked this numerically. It seems like it should work though.

    '''
    if params['alpha'] >= mpmath.pi/2.0:
        warnings.warn('The validity of this function for alpha>= pi/2 has not \
                      been checked. Proceed with caution.')

    params['a'] = params['R']*mpmath.sin(params['alpha'])  # effective piston radius
    ka = params['k']*params['a']
    if ka >= 5:
        num_cores = -1
    else:
        num_cores = 1

    dzero_value =  d_zero(params['k'],
                          params['R'],
                          params['alpha'])
    paramset = []
    for each in angles:
        one_paramset = copy.deepcopy(params)
        one_paramset['theta'] = each
        paramset.append(args_to_str(**one_paramset))

    dtheta_out = Parallel(n_jobs=num_cores, backend='multiprocessing')(delayed(d_theta_pll)(**params_str) for params_str in tqdm.tqdm(paramset))
    dtheta_values = [mpmath.mpmathify(each) for each in dtheta_out]

    directivity = np.array([20*mpmath.log10(abs(each/dzero_value)) for each in dtheta_values],
                            'float32')

    return None, directivity
