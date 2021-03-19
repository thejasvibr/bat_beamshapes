from time import time
import numpy as np
import matplotlib.pyplot as plt
from mpmath import lu_solve, qr_solve, residual, log, mp
from sympy import I, summation, Sum, gamma, besselj, factorial, sqrt, pi, lambdify
from sympy import KroneckerDelta, conjugate, symbols, Matrix, cos, sin, N
from tqdm import tqdm
from joblib import Parallel, delayed

k, beta, NN, P, Q, n, m, p, q, r = symbols('k, beta, NN, P, Q, n, m, p, q, r')
mp.dps = 300
a, b = symbols('a, b')
dB = lambda X : 20*log(abs(X),10)
precision = 300

k_value = 1
b_value = 2
a_value = 1

beta = b / a
NN = int(10 + 2*k_value*b_value)  #
P = 2*NN
Q = 2*NN

# equation 13.218 - checked by TB and GD 23/1/2020
kronecker = N(KroneckerDelta(r, 1), maxn=precision, strict=True)
numerator_Streng = gamma(r/2 - 1/2 + kronecker)
denominator_Streng = gamma(r/2+1)*gamma(r/2-m-1/2)*gamma(r/2+n-m+1)
Streng = (numerator_Streng/denominator_Streng).evalf(precision)
calc_Streng = lambdify(('n','m','r'), Streng, modules='sympy')

# diff between lambdified and evalf with subs -- a dramatic 10-20 time decrease!
# start = time(); [calc_Streng(0,0,0) for i in range(100)]; print(time()-start)
# start = time(); [Streng.evalf(subs={'r':0, 'm':0, 'n':0}) for i in range(100)]; print(time()-start)

# equation 13.217 - checked by Tb and Gogo 23/1/2020
ka = k*a
term1 = -I*sqrt(pi)*gamma(n+5/2)*(1/factorial(m)**2)
function_to_be_summed = Sum(calc_Streng(n,m,r)*((-I*ka/2)**r), (r, 0, NN)).doit()
n_B_m = (term1*function_to_be_summed).evalf(precision)

calc_Bouwkamp = lambdify(('k','a','n','m'), n_B_m, modules='sympy')

# comparing the drop in calculation times after lambdifying
# start = time(); [calc_Bouwkamp(1,1,0,0) for i in range(100)]; print(time()-start)
# start = time(); [n_B_m.evalf(subs={'k':1,'a':1, 'm':0, 'n':0}) for i in range(100)]; print(time()-start)

print('13.226....')
# equation 13.226 -- check by TB and Gd 23/1/2020
# calculates the M term for one value in the matrix
def calc_M_m_n(args, precision=300):
    '''

    Parameters
    ----------
    args :  tuple/list with
            k_value, b_value, m_value, n_value, P, Q
    precision : int>0
                Defaults to 300

    Returns
    -------
    M_m_n : multiple precision sympy number
    '''
    k_value, b_value, m_value, n_value, P, Q = args
    numerator_Mmn = conjugate(calc_Bouwkamp(k_value, b_value,m_value,p))*calc_Bouwkamp(k_value, b_value,n_value,q)
    denominator_Mmn = p + q + 1
    sumover_q_M_m_n = summation(numerator_Mmn/denominator_Mmn, (q, 0, Q)).evalf(precision)
    M_m_n = summation(sumover_q_M_m_n, (p,0,P)).evalf(precision)
    return M_m_n

b_m_left = conjugate(calc_Bouwkamp(k,b,m,p))/(p+1)
b_m_right = (a/b)**(2*p)
b_m = lambdify(('k','b','a','m','P'), -summation(b_m_left*b_m_right, (p,0,P)))
#               modules='sympy')


def calc_D_theta(theta, k_value, a_value, b_value, An, precision=300):
    '''Implementing Equation 13.252

    Parameters
    ----------
    theta
    k_value
    a_value
    b_value
    An

    Returns
    -------

    '''

    D_theta_inf_baffle = (besselj(1, k_value*a_value*sin(theta))/(k_value*a_value*sin(theta))).evalf(precision)

    summation_terms = Matrix.zeros(len(An)+1,1)
    for n_index in range(0, len(An)):
        summation_term_left1 = An[n_index]*gamma(n_index+5/2)
        summation_term_left2 = (2/(k_value*b_value*sin(theta)))**(n_index+3/2)
        summation_term_left = summation_term_left1*summation_term_left2
        summation_term_right = besselj(n_index+3/2, k_value*b_value*sin(theta))
        summation_terms[n_index] = summation_term_left*summation_term_right
    D_theta_open_finite_baffle = (k_value*b_value*0.5*cos(theta)*sum(summation_terms)).evalf(precision)

    D_theta = D_theta_inf_baffle + D_theta_open_finite_baffle

    return D_theta

def calc_D_0(k_value, b_value,An):
    '''

    Parameters
    ----------
    k_value
    b_value
    An

    Returns
    -------

    '''

    D_0 = (0.5*(1 + k_value*b_value*sum(An)))
    return D_0


if __name__ == '__main__':

    M_matrix = Matrix.zeros(NN + 1, NN + 1)
    params = []
    for m_v in range(0, NN + 1):
        print('\n n loop...')
        for n_v in range(0, NN + 1):
            params.append((k_value, b_value, m_v, n_v, P, Q))

    print('now calculating M_m_n matrix..')
    M_values = Parallel(n_jobs=4, verbose=1)(delayed(calc_M_m_n)(param) for param in tqdm(params))

    i = 0
    for m_v in range(0, NN + 1):
        for n_v in range(0, NN + 1):
            M_matrix[m_v, n_v] = M_values[i]
            i += 1

    b_matrix = Matrix.zeros(NN + 1, 1)
    print('now calculate  b_m matrix')
    for m_v in tqdm(range(0, NN + 1)):
        b_matrix[m_v] = b_m(k_value, b_value, a_value, m_v, P).evalf(precision)

    print('... now solving for An matrix')
    An, error = qr_solve(M_matrix, b_matrix)
    #Anq = qr_solve(M_matrix, b_matrix)

    print('...now calculating beam shapes')
    angles = np.linspace(0.01, 2 * np.pi, 1000)
    directivity = np.zeros(angles.size)
    for i, theta in enumerate(tqdm(angles)):
        d_theta = calc_D_theta(theta, k_value, a_value, b_value, An)
        d_zero = calc_D_0(k_value, b_value, An)
        directivity[i] = abs(d_theta / d_zero)

    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    a0.plot(angles, 20 * np.log10(directivity))

