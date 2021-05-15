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
    
    > But it works when using SymPy, as somehow a 0 is output.
    >
"""

import gmpy2
import tqdm
from symengine import * 
import mpmath
#mpmath.mp.dps = 15
from sympy import expand,symbols, Sum,summation, I, cos, sin, legendre, besselj
from sympy import conjugate,sqrt,  lambdify,pi, factorial, gamma, KroneckerDelta
from sympy import  Piecewise, IndexedBase,limit
a, b, B, m, n, p, P, q, Q, z,k,R,alpha,theta = symbols('a b B m n p P q Q z k R alpha theta')
A_n, N, r = symbols('A_n N r')
r = symbols('r')
NN = symbols('NN', integer=True)
P = 2*NN
Q = 2*NN
ka = k*a
kb =  k*b
beta = b/a
# %%
def not_0_or_negint(X):
    is_positive = X>0
    if is_positive:
        return True
    else:
        is_not_negint = X % 1 !=0
        if is_not_negint:
            return True
        else: 
            return False

#%% The 'streng' portion of the Bouwkamp-Streng term
#eqn 13.218
term1, denom_term1, denom_term2, denom_term3 = symbols('term1 denom_term1 denom_term2 denom_term3')

term1 = r*0.5 - 0.5 + KroneckerDelta(r,1)
denom_term1 = r/2 + 1
denom_term2 = r/2 - m - 1/2
denom_term3 = r/2 + n - m + 1

# The summation has been altered to result in a 0 in case of a pole. 
numerator_n_S_mr = gamma(term1) 
denominator_n_S_mr  = gamma(denom_term1)*gamma(denom_term2)*gamma(denom_term3)

nSmr = numerator_n_S_mr/denominator_n_S_mr

onnsmr = lambdify([m,n,r], nSmr,'scipy')



# eg. m=2, n=1, r=0 --> use solveset(denomterm1,m) to get the pole points
# at r = 2m +1 
# at r = 2m - 2n -2 

# n_S_mr = Piecewise( (one_n_S_mr, (denomterm1>0) & (denomterm2>0)),
#                     (one_n_S_mr, (denomterm1%-1.0!=0) & (denomterm2%-1!=0)),
#                     (0,True))

#inputs = {'m':5,'n':2,'r':5}
#n_S_mr = one_n_S_mr

#%% 
# here the term is treated as nBm(z)
# eqn. 13.217, 13.218
S_sumterm = Sum(nSmr*(-I*ka/2)**r  ,(r, 0, 2*NN))
ssum_f = lambdify([m,n,k,a,NN], S_sumterm,'scipy')

#ssum_f(5,2,10,0.01,20)

#%%
n_B_mka = -I*sqrt(pi)*gamma(n+ 5/2)*(1/(factorial(m)**2))*S_sumterm

nbmka_f = lambdify([m,n,k,a,NN], n_B_mka, 'scipy')
                                    
# eqn 13.231 for B term
m_B_p_kb = n_B_mka.subs({'n':m, 'm':p,'ka':kb})


#%% 
# eqn. 13.226
n_B_q_kb = n_B_mka.subs({'m':n, 'n':q, 'ka':kb})
term_13_226 = (conjugate(m_B_p_kb)*n_B_q_kb)/(p+q+1)
M_mn = Sum(term_13_226, (q,0,Q),(p,0,P))

one_Mmn_term = lambdify([a,b,k,m,n,NN], M_mn,'scipy')

#%% 
# eqn. 13.227

b_term = - Sum((conjugate(m_B_p_kb)/(p+1))*(a/b)**(2*p),  (p,0,P))
bterm_func = lambdify([a,b,k,m,NN], b_term,'scipy')


#%% 
# piston in an infinite baffle portion
eq13_252_pt1 = besselj(1, ka*sin(theta))/( ka*sin(theta))
# piston in a finite open baffle portion
A_n = IndexedBase('A_n')
jj = symbols('jj',integer=True)
pt2_term = A_n[jj]*gamma(jj+5/2)*((2/kb*sin(theta))**(jj+3/2))*besselj(jj+3/2, kb*sin(theta))
eq13_252_pt2 = (kb/2)*cos(theta)*Sum(pt2_term , (jj,0,NN-1))

# eqn. 13.252
d_theta = eq13_252_pt1 + eq13_252_pt2
d_thetaf = lambdify([k,a,b,theta,A_n,NN], d_theta,'sympy')

# eqn 13.253 
def d_0f(kv,bv,Anmat):
    return (1/2)*(1+kv*bv*sum(Anmat))


if __name__ == '__main__':
    from sympy import Matrix, zeros, pi
    #%%
    k_v = 2*pi/(50000.0/330.0)
    ka = 1
    a_v = ka/k_v
    b_v = 2*a_v
    beta = b_v/a_v
    
    N_max = 7#int(10 + 2*ka*beta )
    #%%
    Mmn_mat = zeros(N_max, N_max)
    for mm in tqdm.trange(N_max):
        for nn in range(N_max):
            Mmn_mat[mm,nn] = one_Mmn_term(a_v,b_v,k_v,mm,nn,N_max)
    #%% 
    b_mat = zeros(N_max,1)
    for mm in tqdm.trange(N_max):
        b_mat[mm] = bterm_func(a_v,b_v,k_v,mm, N_max)

    #%% 
    Ann = Mmn_mat.QRsolve(b_mat)
    
    #%% 
    import numpy as np 
    angles = np.linspace(0,np.pi,50)
    angles_rad = [each for each in angles]
    
    at_theta = []    
    for angle in angles:
        at_theta.append( d_thetaf(k_v,a_v,b_v,angle,Ann,N_max).evalf())
    zeroval = d_0f(k_v, b_v, Ann)
    db_theta = 20*np.log10(np.array(np.abs(at_theta)/np.abs(zeroval),'float32'))
    
    import matplotlib.pyplot as plt
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(angles_rad, db_theta)
    



