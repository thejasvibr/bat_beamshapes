'''I have had problems running the piston in a finite baffle model by Mello and Kaerkaeinnen 2005.


'''


from decimal import Decimal
from mpmath import besseljzero
from scipy import linalg
from sympy import symbols, Matrix
from sympy import I, summation, gamma, besselj, factorial, sqrt, pi
from sympy.solvers import linsolve
import numpy as np

Phi_q, tau_m, b, k, m_S_q, m_B_q, n, q, a, b = symbols('Phi_q, tau_m, b, k, m_S_q  m_B_q n q a b ')
# The velocity distribution function Phi_q
# calculating the values directly and numerically because you can't use the besseljzero function with a symbolic variable.
# Gaurav suggested the hack to compute the values first and then proceed to solve the equations.

a = 0.1
b_factor = 2.0

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

        multip_factor = (j0_n_value / 2) ** (2 * q - 1)
        bessels_numerator = besselj(1, j0_n_value * a / b).evalf()
        numerator = ((-1) ** q) * bessels_numerator

        bessels_denom = (besselj(1, j0_n_value).evalf()) ** 2
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
    set of values. The Bouwkamp funciton is defined in equation 49 of
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
    '''
    summation_values = []
    for r in range(M + 1):
        multiplication_factor = Decimal((k * b / 2)) ** (Decimal(2 * (q + r) + 3))  # on the rightside of eq 48

        numerator = ((-1) ** (q + r)) * gamma(m + 5 / 2) * gamma(q + r + 1)
        denominator = factorial(r) * (factorial(q) ** 2) * gamma(r + m + 5 / 2) * gamma(q + r + 5 / 2)

        this_term = (numerator / denominator) * multiplication_factor
        summation_values.append(this_term)
    final_mBq_value = np.sqrt(np.pi) * np.sum(summation_values)
    return (final_mBq_value)
## this function's algebra has been checked by GOGO + TBR on 11/12/2019, Galeria Kaufhaus Cafe

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
    '''
    summation_values = []
    for r in range(M + 1):
        multiplication_factor = Decimal((k * b / 2)) ** (Decimal(2 * (q + r - m)))  # on the rightside of eq 49
        numerator = ((-1) ** (q + r + m)) * gamma(m + 5 / 2) * gamma(q + r - m - 1 / 2)
        denominator = factorial(r) * (factorial(q) ** 2) * gamma(r - m - 1 / 2) * gamma(q + r - m + 1)

        this_term = (numerator / denominator) * multiplication_factor
        summation_values.append(this_term)
    final_mSq_value = np.sqrt(np.pi) * np.sum(summation_values)

    return (final_mSq_value)
## this function's algebra has been checked by GOGO + TBR on 11/12/2019, Galeria Kaufhaus Cafe

################# - now trying to estimate tau_m ###################



freq = 35000.0
v_sound = 330.0
wavelength = v_sound/freq
k = 2*np.pi/wavelength
a = 0.005
b = 4*a
M = 3
Q = M


Phi_q = np.zeros((Q+1,1))
q_index_range = range(0,Phi_q.shape[0])
for i,q in enumerate(q_index_range):
    Phi_q[i] = calc_Phi_q_values(a,b,q)
Phi_q = np.concatenate(Phi_q)


m_index_range = range(0,M+1)
m_B_q = np.zeros((len(q_index_range),len(m_index_range)))
m_S_q = np.zeros(m_B_q.shape)
for row,q in enumerate(q_index_range):
    for col,m in enumerate(m_index_range):
        m_B_q[row, col] = calc_mBQ_for_given_m(k,b,m,q,M=M)
        m_S_q[row, col] = calc_mSQ_for_given_m(k,b,m,q,M=M)

m_BS_q = m_B_q - (1j*m_S_q)


tau_m = linalg.solve(m_BS_q, -Phi_q)
Phi1_verifying = np.dot(m_BS_q, tau_m)
error = Phi1_verifying/-Phi_q

new_Phi = np.concatenate((-Phi_q, np.zeros(M+1)))
#new_matrix = np.array(([m_B_q, m_S_q],[-m_S_q, m_B_q]))
new_matrix_row1 = np.column_stack((m_B_q, m_S_q))
new_matrix_row2 = np.column_stack((-m_S_q, m_B_q))
new_matrix = np.row_stack((new_matrix_row1, new_matrix_row2))

# solve this new_matrix for new_Phi
new_tau = linalg.solve(new_matrix, new_Phi)
tau_real, tau_imag = new_tau[:M+1], new_tau[M+1:]
tau = tau_real+1j*tau_imag

output_Phiq = np.dot(m_BS_q, tau)

error = -Phi_q - output_Phiq

############## Equation 92: directivity function of a closed piston in a finite baffle

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

    Trial 1 : multiplicative term inside of summation.
    """
    # TERM 1 with the kasintheta fraction
    ka_sintheta = k*a*np.sin(theta)
    term1 = besselj(1, ka_sintheta).evalf()/ka_sintheta

    # TERM 2 with the summation
    m = np.arange(0, tau_m.size)
    summation_over_m_values = np.zeros(m.size,dtype=np.complex_)
    kb_sintheta = k*b*np.sin(theta)
    two_by_kbsintheta = 2/kb_sintheta
    # calculate values over each m
    for i,each_m in enumerate(m):
        part1 = tau_m[i]*gamma(each_m+5/2)*two_by_kbsintheta**(each_m+3/2)
        part2 = besselj(each_m+3/2, kb_sintheta).evalf()
        summation_over_m_values[i] = part1*part2

    term2 = k*b*(b**2/a**2)*np.cos(theta)*np.sum(summation_over_m_values)

    D_theta = term1 - term2
    return D_theta

#all_thetas = np.arange(0.001, np.pi, 0.1)
#d_theta_values = np.zeros(all_thetas.size, dtype=np.complex)

#for i, theta in enumerate(all_thetas):
#    d_theta_values[i] = closed_piston_in_finite_baffle(theta, k=500, a=0.01,b=0.2,tau_m=tau_m)


# hangouts with Gogo on 16/1/2020