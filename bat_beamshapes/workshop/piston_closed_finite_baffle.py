#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A piston in a closed finite baffle. 
Model described in Chp.13 'Radiation from a rigid circular piston in a finite
circular closed baffle'.

Contributors:  Thejasvi Beleyur, Gaurav Dhariwal

References
==========
Chp. 13, Beranek, L., & Mellow, T. (2019). Acoustics: Sound Fields, Transducers and 
Vibration. Academic Press.
"""

from gmpy2 import *
from symengine import * 
import mpmath
#mpmath.mp.dps = 15
from sympy import expand,symbols, Sum,summation, I, cos, sin, legendre, besselj
from sympy import conjugate,sqrt,  lambdify,pi, factorial, gamma, KroneckerDelta
from sympy import oo
a, b, B, m, n,NN, p, P, Q, r, z,k,R,alpha,theta = symbols('a b B m n NN p P Q r z k R alpha theta')

ka = k*a
kb =  k*b
beta = b/a
NN = 10 + 2*ka*beta
P = 2*NN
Q = 2*NN

# The 'streng' portion of the Bouwkamp-Streng term
numerator_n_S_mr = gamma(r/2 - 1/2 + KroneckerDelta(r,1)) #eqn 13.218
denominator_n_S_mr  = gamma(r/2 +1)*gamma(r/2-m-1/2)*gamma(r/2+n-m+1)
n_S_mr = numerator_n_S_mr/denominator_n_S_mr


S_sumterm = Sum( n_S_mr*(-I*ka/2)**r  ,(r, 0, oo))
# here the term is treated as nBm(z)
n_B_mka = -I*sqrt(pi)*gamma(n+(5/2))*(1/(factorial(m)**2))*S_sumterm
                                      
# eqn 13.231 for B term
m_B_p_kb = n_B_mka.subs({'n':m, 'm':p,'ka':kb})
b_matrix = Sum(conjugate(m_B_p_kb)*(a/b)**(2*p),  (p,0,P))


# eqn. 13.226
M_mn = Sum()



eq13_252_pt1 = besselj(1, ka*sin(theta))/( ka*sin(theta))


# d_0 
# eqn13_253 = (1/2)*(1+kb*Sum(A_n, (n,0,N)))