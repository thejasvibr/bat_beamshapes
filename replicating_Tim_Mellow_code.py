from time import time
import matplotlib.pyplot as plt
import numpy as np
from sympy import I, summation, gamma, besselj, factorial, sqrt, pi, Sum, lambdify
from sympy import KroneckerDelta, conjugate, symbols, Matrix, cos, sin, log, N
from tqdm import tqdm

k, beta, NN, P, Q, n, m, p, q, r = symbols('k, beta, NN, P, Q, n, m, p, q, r')
a, b = symbols('a, b')

beta = (b/a).evalf()
NN = (10 + 2*k*a*beta).evalf()
P = 2*NN.evalf()
Q = 2*NN.evalf()

# equation 13.218
kronecker = N(KroneckerDelta(r, 1))
numerator_Streng = gamma(r/2 - 1/2 +kronecker)
denominator_Streng = gamma(r/2+1)*gamma(r/2-m-1/2)*gamma(r/2+n-m+1)
Streng = N(numerator_Streng/denominator_Streng, 300)
calc_Streng = lambdify(('n','m','r'), Streng, modules='sympy')

# diff between lambdified and evalf with subs -- a dramatic 10-20 time decrease!
# start = time(); [calc_Streng(0,0,0) for i in range(100)]; print(time()-start)
# start = time(); [Streng.evalf(subs={'r':0, 'm':0, 'n':0}) for i in range(100)]; print(time()-start)

# equation 13.217
ka = k*a
term1 = -I*sqrt(pi)*gamma(n+5/2)*(1/factorial(m)**2)
summation_function = calc_Streng(n,m,r)*(-I*ka/2)**r
summation_term = summation(summation_function, (r, 0, 2*NN))
n_B_m = N(term1*summation_term, 300)

calc_Bouwkamp = lambdify(('k','a','n','m','r'), n_B_m, modules='sympy')

# comparing the drop in calculation times after lambdifying
# start = time(); [calc_Bouwkamp(1,1,0,0,0) for i in range(100)]; print(time()-start)
# start = time(); [n_B_m.evalf(subs={'k':1,'a':1,'r':0, 'm':0, 'n':0}) for i in range(100)]; print(time()-start)

# equation 13.226

# calculates the M term for one value in the matrix
numerator_Mmn = conjugate(calc_Bouwkamp(k,b,p,m,r))*calc_Bouwkamp(k,b,n,q,r)
denomnator_Mmn = p + q + 1
sumover_p_M_m_n = summation(numerator_Mmn/denomnator_Mmn, (p, 0, P.evalf()))
M_m_n = summation(sumover_p_M_m_n, (q,0,Q.evalf()))
calc_M_m_n = lambdify( ('k','b','m','n','r'), M_m_n, modules='sympy')



