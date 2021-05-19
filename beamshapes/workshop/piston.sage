
x, alpha, index = var('x alpha index')
k, m,n,p, r1, R, theta, y, z = var('k m n p r1 R theta y z')


r1 = (R*cos(alpha))/cos(theta)

lm_term = legendre_P(m, cos(theta))*(r1/R)^2
assume(alpha>0)

Lm = integral(lm_term, (theta,0,alpha), algorithm='sympy')
