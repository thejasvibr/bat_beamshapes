"""
Piston in a rigid sphere: the bat version
=========================================

"""
import beamshapes
import matplotlib.pyplot as plt
import mpmath
import numpy as np 



#%% 
# Parametrising piston diameter with 'gape height'
# ------------------------------------------------
# Echolocating bats are pretty tiny animals. Jakobsen, Ratcliffe, Surlykke (2012)
# performed a series of skull and video measurements to estimate gape angle
# and gape height. Here, let's take gape height as the piston diameter. The 
# range is between 4-9 mm.
#
# Here, we'll use parameters for the well-studied *Myotis daubentonii* bat. 
# We focus on it's peak frequency at 51 kHz and gape height of ~`4mm. This is equivalent to a piston radius of 3.25mm.
# Looking at the product *ka*, this represents a value of about ~2. We'll then go onto 
# see how the beamshape looks like for the next harmonic at 102 kHz, representing a `ka` of
# about 4.  

v_sound = mpmath.mpf(330.0)
frequency = mpmath.mpf(51000.0)
k = mpmath.fdiv(2*mpmath.pi, (v_sound/frequency)) # AKA wavenumber

alpha = mpmath.pi/8.0 # report of gape angle of 90 degrees --> 45degrees half angle
R = 6e-3  
#%%
# Here I'm assuming a total skull 'diameter' of 1.2 cm, and thus a radius of 6mm
# The value is admittedly a 'data-free' choice, and perhaps the actual effective 
# 'sphere' is really the inside of the mouth cavity!!


paramv = {}
paramv['R'] = R # This is admittedly a 'data-free' choice
paramv['alpha'] = alpha
paramv['k'] = k

a_calc = paramv['R']*mpmath.sin(paramv['alpha'])
paramv['a'] = a_calc # about 2 mm, which broadly matches the 5mm diameter in Jakobsen, Ratcliffee, Surlykke

# %%
# Let's now run the model, and plot the results.

angles = mpmath.linspace(0,2*mpmath.pi,100)
An_51khz, beamshape_51khz = beamshapes.piston_in_sphere_directivity(angles, paramv)

#%% 
# Let's also run the same for the next harmonic t 102 kHz - by doubling the `k` value. 

k_2ndharmonic = k*2
paramv['k'] = k_2ndharmonic
An_102khz, beamshape_102khz = beamshapes.piston_in_sphere_directivity(angles, paramv)


#%%
# Bat-inspired piston in a sphere beamshape
# -----------------------------------------


plt.figure()
a0 = plt.subplot(111, projection='polar')
plt.plot(angles, beamshape_51khz ,'b', label='51 kHz')
plt.plot(angles, beamshape_102khz,'k', label='102 kHz')

plt.ylim(-36,0);plt.yticks(np.arange(-24,0,12))

plt.xticks(np.arange(0,2*np.pi,np.pi/6))
a0.set_theta_zero_location('N')
plt.legend()
plt.savefig('../docs/source/_static/peak_and_higher_harmonic.png')

#%%
#.. image:: ../_static/peak_and_higher_harmonic.png
#  :width: 400
#  :alt: bat beamshape

# %%
# References
# ----------
# Jakobsen, L., Ratcliffe, J. M., & Surlykke, A. (2013). Convergent acoustic 
# field of view in echolocating bats. `Nature`, 493(7430), 93-96.