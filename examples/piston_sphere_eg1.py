"""
Example to show the calculation of a piston in a sphere.
"""
import bat_beamshapes
import matplotlib.pyplot as plt
import mpmath
mpmath.mp.dps = 50
import numpy as np 


#%%
# Simple Example: Piston in a rigid sphere
# ========================================

#%% 
# Parametrising piston diameter with 'gape height'
# ------------------------------------------------
# Echolocating bats are pretty tiny animals. Jakobsen, Ratcliffe, Surlykke (2012)
# performed a series of skull and video measurements to estimate gape angle
# and gape height. Here, let's take gape height as the piston diameter. The 
# range is between 4-9 mm.
#
# Here, we'll use parameters for the well-studied *Myotis daubentonii* bat. 
# We focus on it's peak frequency at 51 kHz and gape heigh of 6.5mm.
# Looking at the product *ka*, this represents a value of about ~6. 

v_sound = mpmath.mpf(330.0)
frequency = mpmath.mpf(51000.0)
k = mpmath.fdiv(2*mpmath.pi, (v_sound/frequency)) # AKA wavenumber


a = mpmath.mpf(0.0065) # meters
ka = mpmath.fmul(k,a)
alpha = mpmath.pi/8.0 # report of gape angle of 90 degrees --> 45degrees half angle
R = mpmath.fdiv(a, mpmath.sin(alpha))

paramv = {}
paramv['R'] = R
paramv['alpha'] = alpha
paramv['k'] = k
paramv['ka'] = ka
paramv['a'] = a

aucalc = paramv['R']*mpmath.sin(paramv['alpha'])
# %%
# Let's now run the model, and plot the results. Here we will use
# the default 50 dps for all models. 

angles = mpmath.linspace(0,mpmath.pi,100)
beamshape = bat_beamshapes.piston_in_sphere_directionality(angles, paramv)
plt.figure()
a0 = plt.subplot(111, projection='polar')
plt.plot(angles, beamshape)
plt.ylim(-40,0);plt.yticks(np.arange(-40,10,10))
plt.xticks(np.arange(0,2*np.pi,np.pi/6))
plt.savefig(f'ka-{np.round(float(ka),2)}_dps={mpmath.mp.dps}.png')



# %%
# References
# ----------
# Jakobsen, L., Ratcliffe, J. M., & Surlykke, A. (2013). Convergent acoustic 
# field of view in echolocating bats. `Nature`, 493(7430), 93-96.