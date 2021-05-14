'''
Rectangular cap of a sphere. From Chapter 12, Beranek & Mellow 2012

'''
# alpha is the elevation half-angle, beta is the azimuthal half-angle
# this formula is valid only when alpha,beta  are <pi/4!!
S = 4*R**2*sin(alpha)*sin(beta)


#%% There are a set of 6 integral terms in here - naming them accordingly  as Imn_int1 for Integral 1,2,3...etc.
Imn_int1 = Integral(cos(2*m*phi), (phi,0,atan(tan(beta)/tan(alpha))))
Imn_int2 = Integral(legendre(n,cos(theta))**2*m*sin(theta), (theta,0,atan(tan(alpha)/cos(phi))),())

first_integral_pair = *
Imn = 


Amn_num =  (2*n+1)**2 * factorial(n-2*m)*Imn
Amn_denom = I*2*pi*factorial(n+2*m)*(n*sph_hankel2(n-1,k*R)-(n+1)*sph_hankel2(n+1, k*R))



Amat = MatrixSymbol('Amat',N,N)
dtheta_summ_term = Sum(Amat[m,n]*I**(2*m)*(legendre(n,cos(theta))**(2*m))*cos(2*m*phi), (m,0,n/2),(n,0,N-1) )
d_theta =   -(4*pi/(k**2*S))*dtheta_summ_term
dzero_summ_term = Sum(Amat[0,n]*I**n, (n,0,N-1))
d_zero =    -(4*pi/(k**2*S))*dtheta_summ_term    



