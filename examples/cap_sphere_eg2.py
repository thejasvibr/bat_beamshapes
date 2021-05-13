"""
Oscillating cap of a sphere
===========================

"""

import mpmath
mpmath.mp.dps = 50
from mpmath import sin 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import bat_beamshapes
from bat_beamshapes import cap_in_sphere_directionality

# %% 
# This model assumes a curved portion of a sphere (the 'cap') oscillates
# to produce sound. One cool thing about this model is that it produces somewhat
# uniform-ish beams that are frontally biased.
# For some `ka` values, the intensity is actually higher a bit off-axis (at higher `ka`'s). 
#
# Below, let's reproduce an example to demonstrate these lobes that peak off-axis.

wavelength = (mpmath.mpf(330.0)/mpmath.mpf(50000))
ka = 10
k_v = 2*mpmath.pi/wavelength
a_v = ka/k_v
alpha_v = mpmath.pi/3
R_v = a_v/sin(alpha_v) #mpmath.mpf(0.01)

angles = mpmath.linspace(0,mpmath.pi,50)
input_params = {'angle': angles, 'k':k_v, 'R':R_v, 'alpha':alpha_v}
rato = cap_in_sphere_directionality(angles, input_params)
print('kooo')
plt.figure()
a0 = plt.subplot(111, projection='polar')
plt.plot(np.array(angles), rato, label='calculated')
plt.ylim(-40,10);plt.yticks(np.arange(-40,20,10))

a0.set_xticks(np.arange(0,2*np.pi,np.pi/6))
df = pd.read_csv('fig12-17.csv')
plt.plot(np.radians(df['deg']), df['rel_db_0deg'],'*',
         label='Beranek & Mellow 2012') 
# Data digitised from figure 12.17 using WebPlotDigitizer (Ankit Rohatgi)
plt.legend()

plt.savefig(f'capsphere_ka-{np.round(float(ka),2)}_dps={mpmath.mp.dps}.png')


#%%
#.. image:: ../_static/capsphere_ka-10.0_dps=50.png
#  :width: 400
#  :alt: cap of a sphere
