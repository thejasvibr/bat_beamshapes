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
# We focus on it's peak frequency at 51 kHz and gape height of 6.5mm. This is equivalent to a piston radius of 3.25mm.
# Looking at the product *ka*, this represents a value of about ~3. 

v_sound = mpmath.mpf(330.0)
frequency = mpmath.mpf(51000.0)
k = mpmath.fdiv(2*mpmath.pi, (v_sound/frequency)) # AKA wavenumber

alpha = mpmath.pi/8.0 # report of gape angle of 90 degrees --> 45degrees half angle
R = 6e-3  # assume a total skull 'diameter' of 1.2 cm, and thus a radius of 6mm


paramv = {}
paramv['R'] = R
paramv['alpha'] = alpha
paramv['k'] = k

aucalc = paramv['R']*mpmath.sin(paramv['alpha'])
# %%
# Let's now run the model, and plot the results.

if __name__ == '__main__':
    # the name == main part is required on Win10 for some reason - may 
    # not be required on other OS's.
    angles = mpmath.linspace(0,mpmath.pi,50)
    An_51khz, beamshape = beamshapes.piston_in_sphere_directivity(angles, paramv)
    
#%% Let's also calculate the pattern for the next harmonic at 102 kHz
# We implement this by doubling the k value. 
    k_2ndharmonic = k*2
    paramv['k'] = k_2ndharmonic
    #%%
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    plt.plot(angles, beamshape,'b')
    plt.plot(-np.array(angles), beamshape,'b')
    plt.ylim(-40,0);plt.yticks(np.arange(-48,0,12))
    
    plt.xticks(np.arange(0,2*np.pi,np.pi/6))
    a0.set_theta_zero_location('N')
    #plt.savefig(f'ka-{np.round(float(ka),2)}_dps={mpmath.mp.dps}.png')

#%%
#.. image:: ../_static/ka-6.31_dps=100.png
#  :width: 400
#  :alt: bat beamshape

# %%
# References
# ----------
# Jakobsen, L., Ratcliffe, J. M., & Surlykke, A. (2013). Convergent acoustic 
# field of view in echolocating bats. `Nature`, 493(7430), 93-96.