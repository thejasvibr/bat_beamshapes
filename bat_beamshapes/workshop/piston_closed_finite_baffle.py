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


Troubleshooting
===============
2021-05-13 : 
    > Unable to run the summation on eqn. 13.217 because the Gamma term goes to 0 every 
    now and then, which creates poles. The gamma function has poles at 0 and negative
    integers. The summation term needs to skip these points. Is this what's happening in the 
    Mathematica code? (somewhat related : https://mathematica.stackexchange.com/questions/148203/how-to-elegantly-tell-mathematica-where-are-the-poles-of-an-integrand)
    > 
"""

from gmpy2 import *
from symengine import * 
import mpmath
#mpmath.mp.dps = 15
from sympy import expand,symbols, Sum,summation, I, cos, sin, legendre, besselj
from sympy import conjugate,sqrt,  lambdify,pi, factorial, gamma, KroneckerDelta
from sympy import oo, Piecewise

a, b, B, m, n, p, P, q, Q, z,k,R,alpha,theta = symbols('a b B m n p P q Q z k R alpha theta')
A_n, N = symbols('A_n, N')
r = symbols('r')
NN = symbols('NN', integer=True)
ka = k*a
kb =  k*b
beta = b/a



#%% The 'streng' portion of the Bouwkamp-Streng term
#eqn 13.218
numerator_n_S_mr = gamma(r*0.5 - 0.5 + KroneckerDelta(r,1)) 
denominator_n_S_mr  = gamma(r*0.5 +1)*gamma(0.5*r-m-0.5)*gamma(r*0.5 + n - m +1)
n_S_mr = Piecewise((numerator_n_S_mr/denominator_n_S_mr, 

nsmr_func = lambdify([m,n,r], n_S_mr,'mpmath')

#%% 
# here the term is treated as nBm(z)
# eqn. 13.217, 13.218
S_sumterm = Sum( n_S_mr*(-I*ka/2)**r  ,(r, 0, 2*NN))
ssum_f = lambdify([m,n,k,a], S_sumterm)
ssum_f(1,1,30.0,0.1)
#%%
n_B_mka = -I*sqrt(pi)*gamma(n+ 5/2)*(1/(factorial(m)**2))*S_sumterm

nbmka_f = lambdify([m,n,k,a], n_B_mka, 'mpmath')
                                    
# eqn 13.231 for B term
m_B_p_kb = n_B_mka.subs({'n':m, 'm':p,'ka':kb})


#%% 
# eqn. 13.226
n_B_q_kb = n_B_mka.subs({'m':n, 'n':q, 'ka':kb})
term_13_226 = (conjugate(m_B_p_kb)*n_B_q_kb)/(p+q+1)
M_mn = Sum(term_13_226, (q,0,Q),(p,0,P))

#%% 
# eqn. 13.227

b_term = - Sum((conjugate(m_B_p_kb)/(p+1))*(a/b)**(2*p),  (p,0,20))
bterm_func = lambdify([a,b,k,m,n], b_term)
#%% 
# piston in an infinite baffle portion
eq13_252_pt1 = besselj(1, ka*sin(theta))/( ka*sin(theta))
# piston in a finite open baffle portion
pt2_term = A_n*gamma(n+5/2)*((2/kb*sin(theta))**(n+3/2))*besselj(n+3/2, kb*sin(theta))
eq13_252_pt2 = (kb/2)*cos(theta)*Sum(pt2_term , (n,0,N))

# eqn. 13.252
d_theta = eq13_252_pt1 + eq13_252_pt2

# eqn 13.253 
d_0 = (1/2)*(1+kb*Sum(A_n, (n,0,N)))


if __name__ == '__main__':
    
    k_v = 2*mpmath.pi/(50000.0/330.0)
    ka = 1
    a_v = ka/k_v
    b_v = 2*a_v
    def M_mn_func(mv,nv,kv,av,bv):
        NN_val = NN.subs({'k':kv,'b':bv,})
        P = Integer(2*NN_val)
        Q = Integer(2*NN_val)
        M_mn_lambda  = lambdify([m,n,k,a,b], n_B_q_kb, 'mpmath')
        return M_mn_lambda(mv,nv,kv,av,bv)
        
        
    
    M_mn_func(0,0,k_v,a_v,b_v)






