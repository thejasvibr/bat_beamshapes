'''I have had problems running the piston in a finite baffle model by Mello and Kaerkaeinnen 2005.


'''


from decimal import Decimal
from mpmath import besseljzero
from scipy import linalg
from sympy import symbols, Matrix
from sympy import I, summation, gamma, besselj, factorial, sqrt, pi
from sympy.solvers import linsolve
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

Phi_q, tau_m, b, m_S_q, m_B_q, = symbols('Phi_q, tau_m, b, m_S_q  m_B_q')
# The velocity distribution function Phi_q
# calculating the values directly and numerically because you can't use the besseljzero function with a symbolic variable.
# Gaurav suggested the hack to compute the values first and then proceed to solve the equations.


def calc_Phi_q_values(a, b, q, N=40):
    '''Calculate Phi_q, the values for equation 50.

    Parameters
    ----------
    a : float>0. Piston radius

    b : float>a. Baffle radius. b must be > a

    q : integer. External summation index.

    N : integer. Upper limit for summations over the variable n.
        Defaults to 40 as used in Mellow & Karkainnen 2005.

    Returns
    --------
    Phi_q : float.
            The value of the Phi_q function.

    '''

    if b < a:
        raise ValueError('The value of b cannot be less than a')

    j0_n = lambda n: besseljzero(0, n)
    n_limits = range(1, N + 1)
    j0_n_values = [j0_n(nth_zero) for nth_zero in n_limits]  # delte it ?

    all_sum_terms = []
    for n in n_limits:
        j0_n_value = j0_n(n)

        multip_factor = (j0_n_value / 2) ** (2*q - 1)
        bessels_numerator = besselj(1, (j0_n_value * a) / b).evalf()
        numerator = ((-1) ** q) * bessels_numerator

        bessels_denom = (besselj(1, j0_n_value).evalf())**2
        denominator = (factorial(q) ** 2) * bessels_denom
        this_sum_term = (numerator / denominator) * multip_factor

        all_sum_terms.append(this_sum_term)

    Phi_q = (a / b) * np.sum(all_sum_terms)
    # this much has been checked by two pairs of eyes -- the algebra has been coded correctly TWICE. T
    # the resulting values of phi_Q seem to get very neative ...and so was wondering if this is normal/expected.
    return (Phi_q)


# The Bouwkamp function mBq(kb) (eq. 48)
def calc_mBQ_for_given_m(k, b, m, q, M=200):
    '''Calculate m_B_q values for the Bouwkamp function for a given
    set of values. The Bouwkamp funciton is defined in equation 48 of
    Mellow & Karkainnen 2005, JASA


    Parameters
    ----------
    k : 0> float. The wavenumber

    b : 0> float. The baffle radius

    m : 0> integer.An external index.

    q : an external index.

    M : 0> integer. The upper limit over which to perform the summation.
        Defaults to 200 as set by Mellow and Karkainnen 2005.

    Returns
    -------
    final_mBq_value : float.
                    The number resulting from summing up the term over index variable r. The upper limit
                    for r is set by M (see Parameters).
    this function's algebra has been checked by GOGO + TBR on
     11/12/2019, Galeria Kaufhaus Cafe
     AND
     17/01/2020 via google hangouts
    '''
    summation_loop = []
    for r in range(M + 1):
        # the Decimal object is used because the exponent-ing can lead to very large numbers and this
        # raises an "OverflowError: (34, 'Result too large')" message
        multiplication_factor = Decimal((k * b / 2)) ** (Decimal(2 *(q + r) + 3))  # on the rightside of eq 48
        numerator = ((-1)**(q + r)) * gamma(m + 2.5) * gamma(q + r + 1)
        denominator = factorial(r) * (factorial(q) ** 2) * gamma(r + m + 2.5) * gamma(q + r + 2.5)

        this_term = (numerator/denominator) * multiplication_factor
        summation_loop.append(this_term)
    final_mBq_value = np.sqrt(np.pi) * np.sum(summation_loop)
    return (final_mBq_value)

def calc_mSQ_for_given_m(k, b, m, q, M=200):
    '''Calculate m_S_q values of the Streng function for a given
    set of values. The Streng funciton is defined in equation 49 of
    Mellow & Karkainnen 2005, JASA

    Parameters
    ----------
    k : 0> float. The wavenumber

    b : 0> float. The baffle radius

    m : 0> integer.An external index.

    q : an external index.

    M : 0> integer. The upper limit over which to perform the summation.
        Defaults to 200 as set by Mellow and Karkainnen 2005.

    Returns
    -------
    final_mSq_value : float.
                    The number resulting from summing up the term over index variable r. The upper limit
                    for r is set by M (see Parameters).
   This function's algebra has been checked by GOGO + TBR on 11/12/2019, Galeria Kaufhaus Cafe
   AND
   over Google hangouts on 17/01/2020

    '''
    summation_loop = []
    for r in range(M + 1):
        multiplication_factor = Decimal((k * b / 2)) ** (Decimal(2 * (q + r - m)))  # on the rightside of eq 49
        numerator = ((-1)**(q + r + m)) * gamma(m + 2.5) * gamma(q + r - m - 0.5)
        denominator = factorial(r) * (factorial(q)**2) * gamma(r - m - 0.5) * gamma(q + r - m + 1)

        this_term = (numerator / denominator) * multiplication_factor
        summation_loop.append(this_term)
    final_mSq_value = np.sqrt(np.pi) * np.sum(summation_loop)
    return (final_mSq_value)

################# - now trying to estimate tau_m ###################

a_value = 1
b_value = 8*a_value
k_value = np.pi/2
M = 50
Q = M


Phi_q = np.zeros((Q+1,1))
q_index_range = range(0,Q+1)
print('Calculating Phi_q')
for q in tqdm(q_index_range):
    Phi_q[q] = calc_Phi_q_values(a_value, b_value, q)


m_index_range = range(0,M+1)
m_B_q = np.zeros((len(q_index_range),len(m_index_range)))
m_S_q = np.zeros(m_B_q.shape)

print('Calculating m_B_q and m_S_q, progress bar over q')
for row,q in enumerate(tqdm(q_index_range)):
    for col,m in enumerate(m_index_range):
        m_B_q[row, col] = calc_mBQ_for_given_m(k_value, b_value, m, q, M)
        m_S_q[row, col] = calc_mSQ_for_given_m(k_value, b_value, m, q, M)

m_BS_q = m_B_q - (1j*m_S_q)

tau_m = linalg.solve(m_BS_q, -Phi_q)
Phi1_verifying = np.dot(m_BS_q, tau_m)
error = Phi1_verifying/-Phi_q

print('tau_m solved...')
compare_results_and_Phiq = pd.DataFrame(data={'RHS':-Phi_q.flatten(),
                                              'prediction':Phi1_verifying.flatten()})

######## trying out SVD w code based on https://sukhbinder.wordpress.com/2013/03/26/solving-axb-by-svd/

def svdsolve(a,b):
    u,s,v = np.linalg.svd(a)
    c = np.dot(u.T,b)
    w = np.linalg.solve(np.diag(s),c)
    x = np.dot(v.T,w)
    return x

tau_m_w_svd = svdsolve(m_BS_q, -Phi_q)
svd_Phi1_verifying = np.dot(m_BS_q, tau_m_w_svd)
svd_error = svd_Phi1_verifying/-Phi_q

svd_compare_results_and_Phiq = pd.DataFrame(data={'RHS':-Phi_q.flatten(),
                                                  'linalg_prediction':Phi1_verifying.flatten(),
                                              'svd_prediction':svd_Phi1_verifying.flatten()})

# ################## - Gogo's calculations to estimate tau manually.
# new_Phi = np.concatenate((-Phi_q.flatten(), np.zeros(M+1)))
# #new_matrix = np.array(([m_B_q, m_S_q],[-m_S_q, m_B_q]))
# new_matrix_row1 = np.column_stack((m_B_q, m_S_q))
# new_matrix_row2 = np.column_stack((-m_S_q, m_B_q))
# new_matrix = np.row_stack((new_matrix_row1, new_matrix_row2))

# solve this new_matrix for new_Phi
# new_tau = linalg.solve(new_matrix, new_Phi)
# tau_real, tau_imag = new_tau[:M+1], new_tau[M+1:]
# tau = tau_real+1j*tau_imag
#
# output_Phiq = np.dot(m_BS_q, tau)
#
# error = -Phi_q - output_Phiq
# print(-Phi_q)
# print(output_Phiq)

############ - - the error in solving for tau_m probably comes from the floating point error.
# a 64 bit number can at the most handle 2^-16 and our numbers are way below that or close to it.
# What if we use a Decimal representation of the numbers.
from mpmath import matrix, lu_solve, mpf, mp, fdot, svd, mpmathify,diag, residual,eye
mp.dps = 200

mp_tau_m = lu_solve(m_BS_q, -Phi_q)
test = m_BS_q*mp_tau_m


def svdsolve_mpmath(a,b):
    u,s,v = svd(a)
    c = u.transpose()*b
    w = lu_solve(diag(s), c)
    x = v.transpose()*w
    return x

def convert_numpy_array_to_mp_matrix(np_array):
    '''

    Parameters
    ----------
    np_array

    Returns
    -------

    '''
    rows, columns = np_array.shape
    mp_matrix  = matrix(rows, columns)
    for row in range(rows):
        for column in range(columns):
            mp_matrix[row,column] = np_array[row, column]
    return mp_matrix

print('....now calculating the high precision lu_solve for tau_m')
m_BS_q_MP = convert_numpy_array_to_mp_matrix(m_BS_q)
Phi_q_MP = convert_numpy_array_to_mp_matrix(Phi_q)

tau_m_MP = lu_solve(m_BS_q_MP, -Phi_q_MP)
residual(m_BS_q_MP, tau_m_MP, -Phi_q_MP)

tau_m_mp_np = np.array(tau_m_MP.tolist())

## open piston in finite baffle:

def disk_in_finite_baffle(theta, k, **kwargs):
    """ function to calculate the directivity of a closed piston in a finite baffle as per
    Mellow and Kaerkkkaeinen, JASA 2005 equantion 67

    Parameters
    -----------
    theta: float in radians.
            Angle of emission, where 0 radians is on-axis, and increasing values are
            increasingly off-axis.
    k: float>0.
        wavenumber.

    Keyword Arguments
    -----------------
    a : float>0.
        Piston radius
    b : b>a
        Baffle radius.
    tau_m : ?? dimensions
        coefficients obtained by solving equation 47.

    Note to self:
    the multiplication term is *INSIDE* the summation!!!

    The algebra for this function has been checked by TBR and GD on 17/1/2020
    over google hangouts
    """
    # TERM 1 with the kasintheta fraction
    a = kwargs['a']
    b = kwargs['b']

    # TERM 2 with the summation
    m = np.arange(0, tau_m.size)
    summation_over_m_values = np.zeros(m.size, dtype=np.complex)
    kb_sintheta = k*b*np.sin(theta)
    two_by_kbsintheta = mpf(2/kb_sintheta)
    # calculate values over each m
    for each_m in m:
        part1 = tau_m[each_m]*gamma(each_m+2.5)*two_by_kbsintheta**(each_m+1.5)
        part2 = besselj(each_m+1.5, kb_sintheta).evalf()
        summation_over_m_values[each_m] = part1[0]*part2

    term2 = k*b*(b**2/a**2)*np.cos(theta)*np.sum(summation_over_m_values)

    D_theta = -term2 # there is a typo in the paper -- it should be 2*J1(kasin(theta))
    D_zero = k*b*(b**2/a**2)*np.sum(tau_m)

    return D_theta, D_zero




#######
def closed_piston_in_finite_baffle(theta, k, **kwargs):
    """ function to calculate the directivity of a closed piston in a finite baffle as per
    Mellow and Kaerkkkaeinen, JASA 2005 equantion 92.

    Parameters
    -----------
    theta: float in radians.
            Angle of emission, where 0 radians is on-axis, and increasing values are
            increasingly off-axis.
    k: float>0.
        wavenumber.

    Keyword Arguments
    -----------------
    a : float>0.
        Piston radius
    b : b>a
        Baffle radius.
    tau_m : ?? dimensions
        coefficients obtained by solving equation 47.

    Note to self:
    the multiplication term is *INSIDE* the summation!!!

    The algebra for this function has been checked by TBR and GD on 17/1/2020
    over google hangouts
    """
    # TERM 1 with the kasintheta fraction
    a = kwargs['a']
    b = kwargs['b']
    ka_sintheta = k*a*np.sin(theta)
    term1 = 2*besselj(1, ka_sintheta).evalf()/ka_sintheta

    # TERM 2 with the summation
    m = np.arange(0, tau_m.size)
    summation_over_m_values = np.zeros(m.size, dtype=np.complex)
    kb_sintheta = k*b*np.sin(theta)
    two_by_kbsintheta = mpf(2/kb_sintheta)
    # calculate values over each m
    for each_m in m:
        part1 = tau_m[each_m]*gamma(each_m+2.5)*two_by_kbsintheta**(each_m+1.5)
        part2 = besselj(each_m+1.5, kb_sintheta).evalf()
        summation_over_m_values[each_m] = part1[0]*part2

    term2 = k*b*(b**2/a**2)*np.cos(theta)*np.sum(summation_over_m_values)

    D_theta = 0.5*(term1 - term2) # there is a typo in the paper -- it should be 2*J1(kasin(theta))
    D_zero = 0.5*(1 - k*b*(b**2/a**2)*np.sum(tau_m))

    return D_theta, D_zero

######### --- calculating the values and plotting the directivity function
all_thetas = np.arange(0.001, np.pi, 0.01)
d_theta_values_closed_finite = np.zeros(all_thetas.size, dtype=np.complex)
d_theta_values_open_finite = np.zeros(all_thetas.size, dtype=np.complex)

print('calculating the beamshape now....')
for i, theta in enumerate(tqdm(all_thetas)):
    d_theta_values_closed_finite[i], d_zero_closed = closed_piston_in_finite_baffle(theta,
                                                       k=k_value, a=a_value,b=b_value,
                                                       tau_m=tau_m_mp_np)

    d_theta_values_open_finite[i], d_zero_open = disk_in_finite_baffle(theta,
                                                                       k=k_value, a=a_value, b=b_value,
                                                                       tau_m=tau_m_mp_np)

directivity_closed = np.float64((np.abs(d_theta_values_closed_finite)/np.abs(d_zero_closed)))
directivity_open = np.float64((np.abs(d_theta_values_open_finite)/np.abs(d_zero_open)))

dB = lambda X : 20*np.log10(X)

plt.figure()
ax = plt.subplot(111, projection='polar')
ax.set_theta_zero_location('E')
plt.plot(all_thetas, dB(directivity_closed))
plt.plot(all_thetas, dB(directivity_open))

ax.set_rgrids([-30,-20,-10,0,10,20], angle=60)
ax.set_ylim(-40,20)
#ax.set_xticks(np.linspace(0,2*np.pi,24))
ax.set_title('k:'+str(np.round(k_value,3))+' '+'a:'+str(a_value)+' '+'b'+str(b_value))
# hangouts with Gogo on 16/1/2020
