'''
Piston in a sphere (FLINT implementation)
=========================================
Very fast implementation of piston in a sphere directivity using python-Flint. 

Installation/requirements
-------------------------
This implementation requires the installation of the 'python-flint' package.
See `here <https://github.com/fredrik-johansson/python-flint/>`_ or the `docs <https://fredrikj.net/python-flint/>`_ .

While the installation is a bit involved due to the packages dependencies (GMP,MPFR, FLINT etc),
the speed-up afterwards is worth it. This implementation currently only works on 
Linux machines (as of python-flint 0.3.0 - but perhaps later releases may be OS-agnostic!).

References
----------
Beranek, L. L., & Mellow, T. (2012). Acoustics: sound fields and transducers.
Academic Press.

See Also 
--------
beamshapes.piston_in_sphere

'''


from joblib import Parallel, delayed
import flint 
from flint import ctx
import numpy as np 
import tqdm
from beamshapes.flint_parallelisation import conv_acb_to_str, conv_str_acb, interchange_params_and_str
ctx.dps = 50
cos = flint.acb.cos
sin = flint.acb.sin
tan = flint.acb.tan
pi = flint.acb.pi()
arb = flint.arb
acb = flint.acb
acb_mat = flint.acb_mat
good = flint.good
integral = flint.acb.integral

 

def conv_acb_to_int(X):
    try:
        int_out = int(float(X.abs_upper()))
    except:
        int_out = int(float(str(X.abs_upper())))
    return int_out


## !! ATTENTION -- the FLINT legendre has diff conventions than the 
# mpmath one !! 
legendre_p = lambda n,z : flint.acb.legendre_p(z,n) 

# spherical bessel and hankel functions
besselj = lambda n,z : acb.bessel_j(z, n) # opposite order of mpmath!
bessely = lambda n,z : acb.bessel_y(z, n)
sph_besselj = lambda n,z : acb.sqrt(pi/(2*z))*besselj(n+acb(0.5),z)
sph_bessely = lambda n,z : acb.sqrt(pi/(2*z))*bessely(n+acb(0.5),z)
sph_hankel2 = lambda n,z : sph_besselj(n, z) - 1j*sph_bessely(n,z)


# P-prime-cos(alpha) in Appendix II eqn.70 and NOT in eqn. 12.98 in 2012edn.
def pprime_cosalpha(n, alpha):
    legendre_term = legendre_p(n+1,cos(alpha))-legendre_p(n-1,cos(alpha))
    numerator = (n*(n+1))*legendre_term
    denominator = (1+2*n)*(-1+cos(alpha)**2)
    final_term = -sin(alpha)*(numerator/denominator)
    return final_term

 Calculating N -  using the formula from T. Mellow's code.
def calc_defaultNN(param):
    ka = param['k']*param['a']
    return conv_acb_to_int(12+2*ka/sin(param['alpha']))



r1 = lambda theta,alpha,R : R*cos(alpha)/cos(theta) 

# eqn. 12.108 Lm 
# r1/R --> cos(alpha)/cos(theta)
Lm_term = lambda theta,m,alpha : legendre_p(m, cos(theta))*((cos(alpha)/cos(theta))**2)*tan(theta)

def Lm_func(mv, alphav):
    '''
    '''
    only_theta_lm = lambda theta, a : legendre_p(mv, cos(theta))*((cos(alphav)/cos(theta))**2)*tan(theta)
    return integral(only_theta_lm, acb(0.0), alphav)

# eqn 12.107 solution for the Kmn integral - See eqn. 70 in Appendix II
# the 2012 edition has a typo for eqn.70 App II, and so be sure to check
# the Errata. The order of denominator term's 'n' and 'm' are different.
def m_noteq_n(alpha,m,n):
    term1 = legendre_p(m,cos(alpha))*pprime_cosalpha(n, alpha)
    term2 = legendre_p(n, cos(alpha))*pprime_cosalpha(m, alpha)
    value = (sin(alpha)*(term1-term2))/(n*(n+1)-m*(m+1))
    return value

def meqn_legendre_summn_term(jj,alpha):
    first_term = legendre_p(jj,cos(alpha))
    second_term = (legendre_p(jj,cos(alpha))*cos(alpha)-legendre_p(jj+1,cos(alpha)))
    return first_term*second_term

def meqn_summterm(m,alpha):
    m = int(float(str(m)))
    terms = acb_mat([[meqn_legendre_summn_term(jj, alpha) for jj in range(1,m)]])
    return 2*sum(terms)
     
def m_eq_n(alpha,m,n):
    terms_wo_summation = 1+cos(alpha)*legendre_p(m, cos(alpha))**2
    summation_term = meqn_summterm(m,alpha)
    numerator = terms_wo_summation + summation_term
    denominator = 2*m+1
    return numerator/denominator

def kmn_func(alpha,m,n):
    if m!=n:
        output = m_noteq_n(alpha, m, n)
    else:
        output = m_eq_n(alpha, m, n)
    return output

 eqn. 12.106
alternate_hankels = lambda n,kr1 : n*sph_hankel2(n-1, kr1) - (n+1)*sph_hankel2(n+1, kr1)


def imn_term(theta,m,n,k,R,alpha):
    kr1 = k*r1(theta,alpha,R)
    r1byR_squared = (cos(alpha)/cos(theta))**2
    subterm_pt1 = alternate_hankels(n,kr1)*legendre_p(n, cos(theta))*cos(theta)
    subterm_pt2 = n*(n+1)*sph_hankel2(n, kr1)*(legendre_p(n-1, cos(theta)) - legendre_p(n+1, cos(theta)))/kr1 
    final_term = (subterm_pt1 + subterm_pt2)*legendre_p(m,cos(theta))*r1byR_squared*tan(theta)
    return final_term

def Imn_func(mv,nv,kv,Rv,alphav):
    imn_wrapper = lambda theta,_: imn_term(theta, mv,nv,kv,Rv,alphav)
    return integral(imn_wrapper, acb(0.0), alphav)

 
def calc_one_Mmn_term(param):
    m,n = param['m'], param['n']
    k,R,alpha = [param[each] for each in ['k','R','alpha']]
    kR = k*R
    
    Imn_value = Imn_func(m,n,k,R,alpha)
    hankels_value = n*sph_hankel2(n-1, kR) - (n+1)*sph_hankel2(n+1, kR)
    Kmn_value = kmn_func(alpha, m, n)
    numerator = Imn_value + hankels_value*Kmn_value
    denominator = 2*n + 1 
    return numerator/denominator

def calc_one_Mmn_term_pll(param):
    '''
    '''
    # convert from string to acb
    param_acb = interchange_params_and_str(param, to_str=False)
    # obtain output
    one_Mmn_term = calc_one_Mmn_term(param_acb)
    # convert back to string
    one_Mmn_term_str = conv_acb_to_str(one_Mmn_term)
    return one_Mmn_term_str



def make_Mmn(param):
    NN = param.get('NN', calc_defaultNN(param))
    Mmn_matrix = acb_mat(NN,NN)
    for mm in tqdm.trange(NN):
        for nn in range(NN):
            param['m'] = mm
            param['n'] = nn
            Mmn_matrix[mm,nn] = calc_one_Mmn_term(param)
    return Mmn_matrix

def make_Mmn_pll(param):
    '''
    Wrapper to calculate one Mmn term in a parallel loop. 
    This function interchanges `acb` and `str` objects 
    in both directions. For proper parallelisation, only
    string-converted acb objects are used, converted back to
    acb objects for calculations, and then sent back as
    string again. 
    
    Parameters
    ----------
    param : dictionary
        Dictionary with the following entries
        
    Returns 
    -------
    Mmn_matrix : flint.acb_mat
        A matrix with all coefficients that satisfy model 
        conditions. 
    
    See Also
    --------
    calc_one_Mmn_term_pll, calc_one_Mmn_term
    '''
    NN = param.get('NN', calc_defaultNN(param))
    Mmn_matrix = acb_mat(NN,NN)
    n_cores = param.get('n_cores', -1)
    
    Mmn_as_str = []
    params_as_str = {}
    
    for mm in range(NN):
        for nn in range(NN):
            param['m'] = acb(mm)
            param['n'] = acb(nn)
            params_as_str[(mm,nn)] = interchange_params_and_str(param,
                                                                to_str=True)

    # run parallel calc
    Mmn_as_str = Parallel(n_jobs=n_cores, backend='multiprocessing')(delayed(calc_one_Mmn_term_pll)(paramset)   for position, paramset in tqdm.tqdm(params_as_str.items()))
    Mmn_acb = [conv_str_acb(each) for each in Mmn_as_str]

    entry_num = 0
    for rowpos in range(NN):
        for colpos in range(NN):
            Mmn_matrix[rowpos,colpos] = Mmn_acb[entry_num]
            entry_num += 1 
    return Mmn_matrix


def make_bm(param):
    NN = param.get('NN', calc_defaultNN(param))
    bm_mat = acb_mat(NN,1)
    
    alpha = param['alpha']
    for mm in range(NN):
        bm_mat[mm,0] = -1j*Lm_func(acb(mm),alpha)
    return bm_mat

def dtheta(theta, An, param):
    k,R,alpha = (param[each] for each in ['k','R','alpha'])
    NN = max([An.ncols(),An.nrows()])
    I_n_mat = acb_mat(1,NN)
    legendre_mat = acb_mat(1,NN)
    for nn in range(NN):
        I_n_mat[0,nn] = acb(1j)**nn
        legendre_mat[0,nn] = legendre_p(acb(nn), cos(theta))
    
    an_in_legendre = acb_mat(1,NN)
    for nn in range(NN):
        an_in_legendre[0,nn] = An[nn,0]*I_n_mat[0,nn]*legendre_mat[0,nn]
    sum_term = sum(an_in_legendre)
    
    pt1 = -acb(4)/(k**2*R**2*sin(alpha)**2)
    
    return pt1*sum_term

def dzero(An, param):
    k,R,alpha = (param[each] for each in ['k','R','alpha'])
    NN = max([An.ncols(),An.nrows()])
    I_n_mat = acb_mat(1,NN)
    for nn in range(NN):
        I_n_mat[0,nn] = acb(1j)**nn
    an_in = acb_mat([[An[nn,0]*I_n_mat[0,nn] for nn in range(NN)]])
    sum_term = sum(an_in)
    pt1 = -acb(4)/(k**2*R**2*sin(alpha)**2)
    
    return pt1*sum_term

def directivity(thetas, An, param):
    d_zero_val = dzero(An, param)
    ratios = []
    for each in thetas:
        dtheta_val = dtheta(each, An, param)
        ratio = abs(dtheta_val)/abs(d_zero_val)
        ratios.append(ratio)
    return ratios

def piston_in_sphere_directivity(thetas, param,**kwargs):
    '''Calculates directivity for piston in a sphere using the
    fast python-flint package. 

    Parameters
    ----------
    thetas : list with acb entries
        List with angles in radians.
    param: dictionary
    Dictionary with at least the following keys:
            k : mpmath.mpf>0
                Wavenumber. 
            R : mpmath.mpf>0
                Radius of sphere
            alpha: 0<mpmath.mpf<pi
                Half-angular aperture of piston. 
            num_cores : int, optional
                The number of cores to be used for the Mmn matrix computation. 
                Defaults to using all cores.
        

    Keyword Arguments
    -----------------
    An : optional, acb_mat
        A complex matrix from a previous calculation. 

    Returns
    -------
    A_n : acb_mat
        The 'An' term required for calculating directionalities. Using the pre-calculated `An` 
        saves time for repeated use if it was not already available. 
    dB_directivity : np.array
        Array with 20log10(on-axis/off-axis) values. 
    '''
    A_n = kwargs.get('A_n', None)
    if A_n is None:
        mmn_mat  = make_Mmn_pll(param)
        b_mat = make_bm(param)
        A_n = mmn_mat.solve(b_mat)
    
    ratios = directivity(thetas, A_n, param)
    dB_directivity = 20*np.log10(np.float32(ratios))
    return A_n, dB_directivity