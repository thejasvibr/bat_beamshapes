'''This code snippet is heavily based on Mathematica code sent by Tim Mellow.

'''
import matplotlib.pyplot as plt
import numpy as np
from sympy import I, summation, gamma, besselj, factorial, sqrt, pi, Sum
from sympy import KroneckerDelta, conjugate, symbols, Matrix, cos, sin, log
from tqdm import tqdm

log10 = lambda X: log(X,10)
dB = lambda X : 20*log10(X).evalf()
k = 1
b = 2
a = 1


def Streng(n,m,r,precision=300):
    '''Equation 13.218 in Chapter 13 of Mellow & Beranek 2012

    Parameters
    ----------
    n_value
    m_value
    r_value

    Returns
    -------
    n_S_m

    '''
    #r, n, m = symbols('r,n,m')
    numerator = gamma(r/2- 1/2 + KroneckerDelta(r,1))
    denominator = (gamma(r/2 + 1)*gamma(r/2-m-1/2)*gamma(r/2+n-m+1))
    n_S_m = (numerator/denominator).evalf(precision)
    return n_S_m

def Bouwkamp(ka,n,m,b,a,precision=300):
    '''Chapter 13, Equation 13.217 of Mellow & Beranek 2012

    Parameters
    ----------
    ka : the product of k*a
    n
    m

    Returns
    -------
    n_B_m

    '''

    beta = b / a
    NN = int(10 + 2 * ka * beta)
    r = symbols('r')
    term1 = -1j*sqrt(pi)*gamma(n+5/2)*(1/factorial(m)**2)
    summation_function = Streng(n,m,r)*(-I*ka/2)**r
    summation_term = summation(summation_function, (r, 0, 2*NN))
    n_B_m = (term1*summation_term).evalf(precision)
    return n_B_m

def calculate_M_b_An(k,b,a,precision=300):
    '''Equations 13.226 and 13.227

    Parameters
    ----------
    N

    Returns
    -------
    M

    '''

    beta = b / a
    NN = int(20 + 2*k*a*beta) # just changed to 20
    P = 2 * NN
    Q = 2 * NN

    p,q = symbols('p,q')
    M_matrix = Matrix.zeros(NN+1,NN+1)
    b_matrix = Matrix.zeros(NN+1,1)

    for m in tqdm(range(0, NN+1)):
        eq13227_rightside = (conjugate(Bouwkamp(k*b, m, p, b, a))/(p + 1)) * (a/b)**(2*p)
        b_matrix[m] = summation(eq13227_rightside,(p,0,P)).evalf(precision)
        for n in range(0,NN+1):
            rightside_term_numerator = conjugate(Bouwkamp(k*b,n,m,b,a)) * Bouwkamp(k*b, n, m,b,a)
            rightside_term_denominator = p + q + 1
            rightside_term = rightside_term_numerator / rightside_term_denominator

            summation_over_q = summation(rightside_term,(q,0,Q))
            M_matrix[m,n] = summation(summation_over_q, (p,0,P)).evalf(precision)
    An = M_matrix.LUsolve(b_matrix)
    return M_matrix, b_matrix, An

def D_disc_in_finite_open_baffle(theta,k,a,b,An):
    '''Equation 13.235, Chapter 13, Beranek & Mellow 2012
    Parameters
    ----------
    theta
    k
    a
    b
    An

    Returns
    -------

    '''
    D_theta_terms = []
    for n, An_n in enumerate(An):
        left_term = k*b*cos(theta)
        term1 = An_n*gamma(n+5/2)
        term2 = (2/(k*b*sin(theta)))**(n+3/2)
        term3 = besselj(n+3/2, k*b*sin(theta))
        eq13235_rightside = term1*term2*term3
        D_theta_terms.append(eq13235_rightside)
    D_theta = sum(D_theta_terms)

    D_0 = k*b*sum(An) # equation 13.236

    return D_theta, D_0

if __name__ == '__main__':
    k =1
    b = 2
    a = 1

    M_matrix, b_matrix, An = calculate_M_b_An(k,b,a)
    angles = np.linspace(0,np.pi/2,100)
    directivity = np.zeros(angles.size)
    for i,theta in enumerate(angles):
        d_theta, d_0 = D_disc_in_finite_open_baffle(theta, k, a, b, An)
        directivity[i] = dB(abs(d_theta)/abs(d_0)) +6

    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    a0.set_theta_zero_location('E')
    plt.plot(angles, directivity)
    a0.set_rgrids([-30,-20, -10, 0], angle=60)

