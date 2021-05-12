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
import tqdm
ctx.dps = 150
cos = flint.acb.cos
sin = flint.acb.sin
tan = flint.acb.tan
pi = flint.acb.pi()
arb = flint.arb
acb = flint.acb
acb_mat = flint.acb_mat
good = flint.good

#%%
## !! ATTENTION -- the FLINT legendre has diff conventions than the 
# mpmath one !! 
legendre_p = lambda n,z : flint.acb.legendre_p(z,n) 

# spherical bessel and hankel functions
besselj = lambda n,z : acb.bessel_j(z, n) # opposite order of mpmath!
bessely = lambda n,z : acb.bessel_y(z, n)
sph_besselj = lambda n,z : acb.sqrt(pi/(2*z))*besselj(n+acb(0.5),z)
sph_bessely = lambda n,z : acb.sqrt(pi/(2*z))*bessely(n+acb(0.5),z)

sph_hankel2 = lambda n,z : sph_besselj(n, z) - 1j*sph_bessely(n,z)

#%%
# P-prime-cos(alpha) in Appendix II eqn.70 and NOT in eqn. 12.98 in 2012edn.

def pprime_cosalpha(n, alpha):
    legendre_term = legendre_p(n-1,cos(alpha))-legendre_p(n+1,cos(alpha))
    numerator = n*(n+1)*legendre_term
    denominator = (1+2*n)*(-sin(alpha)**2)
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

# %% eqn 12.107 solution for the Kmn integral 
def m_noteq_n(alpha,m,n):
    term1 = legendre_p(m,cos(alpha))*pprime_cosalpha(n, alpha)
    term2 = legendre_p(n, cos(alpha))*pprime_cosalpha(m, alpha)
    value = sin(alpha)*(term1-term2)/(m*(m+1)-n*(n+1))
    return value

def meqn_legendre_summn_term(jj,alpha):
    return legendre_p(jj,cos(alpha))*(legendre_p(jj,cos(alpha))*cos(alpha)-legendre_p(jj+1,cos(alpha)))

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

#%% eqn. 12.06
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

#%% 
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

#%%
def conv_acb_to_int(X):
    try:
        int_out = int(float(X.abs_upper()))
    except:
        int_out = int(float(str(X.abs_upper())))
    return int_out

def make_Mmn(param):
    ka = params['k']*params['a']
    NN = param.get('NN', conv_acb_to_int(12+2*ka/sin(param['alpha'])))
    Mmn_matrix = acb_mat(NN,NN)
    for mm in tqdm.trange(NN):
        for nn in range(NN):
            param['m'] = mm
            param['n'] = nn
            Mmn_matrix[mm,nn] = calc_one_Mmn_term(param)
    return Mmn_matrix
#%%

def make_bm(param):
    ka = params['k']*params['a']
    NN = param.get('NN', conv_acb_to_int(12+2*ka/sin(param['alpha'])))
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

def directionality(thetas, An, param):
    d_zero_val = dzero(An, param)
    ratios = []
    for each in thetas:
        dtheta_val = dtheta(each, An, param)
        ratio = abs(dtheta_val)/abs(d_zero_val)
        ratios.append(ratio)
    return ratios


#%%


if __name__ == '__main__':
    import numpy as np 
    ka = acb(10)
    kv = 2*pi/(330.0/50000.0)
    av = ka/kv
    alphav = pi/3
    Rv = av/sin(alphav)
    params = {'k':kv,
              'a':av,
              'alpha': alphav,
              'R': Rv}
    print(f'The DPS is: {ctx.dps}')
    mmn_mat = make_Mmn(params)
    b_mat = make_bm(params)
    an = mmn_mat.solve(b_mat)
    #%%
    ratios = directionality([pi*0, pi/6,pi/4, pi/2,pi], an, params)
    dbratios = 20*np.log10(np.float32(ratios))
