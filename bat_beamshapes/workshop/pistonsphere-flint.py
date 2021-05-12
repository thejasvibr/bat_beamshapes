'''
Re-working of piston in a sphere using python-Flint

It seemed like there might be an issue with the numerical routines in 
mpmath - and so this is trying to see if there's a chance this package
can overcome the problems

The model is implemented from Beranek & Mellow 2012 Chp. 12
All equations refer to those there.

'''
#%%
import flint 
from flint import ctx
ctx.dps = 500
cos = flint.acb.cos
sin = flint.acb.sin
tan = flint.acb.tan
pi = flint.acb.pi()
arb = flint.arb
acb = flint.acb
acb_mat = flint.acb_mat
good = flint.good
## !! ATTENTION -- the FLINT legendre has diff conventions than the 
# mpmath one !! 
legendre_p = lambda n,z : flint.acb.legendre_p(z,n) 


# P-prime-cos(alpha) in Appendix II eqn.70 and NOT in eqn. 12.98 in 2012edn.

def pprime_cosalpha(n, alpha):
    legendre_term = legendre_p(n+1,cos(alpha))-legendre_p(n+1,cos(alpha))
    numerator = n*(n+1)*legendre_term
    denominator = (1+2*n)*sin(alpha)
    return -numerator/denominator

integral = flint.acb.integral

#%%
r1 = lambda theta,alpha,R : R*cos(alpha)/cos(theta) 

# eqn. 12.018 Lm 
# r1/R --> cos(alpha)/cos(theta)
Lm_term = lambda theta,m,alpha : legendre_p(m, cos(theta))*((cos(alpha)/cos(theta))**2)*tan(theta)

def Lm_func(mv, alphav):
    '''
    '''
    only_theta_lm = lambda theta, a : legendre_p(mv, cos(theta))*((cos(alphav)/cos(theta))**2)*tan(theta)
    return integral(only_theta_lm, acb(0.0), alphav)

#Lm_func(acb(10.0), pi/3.0)

# %% eqn 12.107 solution for the Kmn integral 
def m_noteq_n(alpha,m,n):
    term1 = legendre_p(m,cos(alpha))*pprime_cosalpha(n, alpha)
    term2 = legendre_p(n, cos(alpha))*pprime_cosalpha(m, alpha)
    value = sin(alpha)*(term1-term2)/(m*(m+1)-n*(n+1))
    return value

def meqn_legendre_summn_term(jj,alpha):
    return legendre_p(jj,cos(alpha))*(legendre_p(jj,cos(alpha))*cos(alpha)-legendre_p(jj+1,cos(alpha)))

def meqn_summterm(m,alpha):
    terms = acb_mat([[meqn_legendre_summn_term(jj, alpha) for jj in range(1,m)]])
    return 2*sum(terms)
     
def m_eq_n(alpha,m,n):
    terms_wo_summation = 1+cos(alpha)*legendre_p(m, cos(alpha))**2
    summation_term = meqn_summterm(m,alpha)
    numerator = terms_wo_summation + summation_term
    denominator = 2*m+1
    return numerator/denominator

    






