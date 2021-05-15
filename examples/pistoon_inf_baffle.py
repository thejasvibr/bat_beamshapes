#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Piston in an infinite baffle example
====================================
"""

#%% The piston in an infinite baffle is a classic model used to describe
# the beamshapes of mouth-emitting bat species [1,2] (among other animals eg. [3]). 
# It's easily to compute, but however can only describe the beamshape upto :math:`\pm90^{\circ}`
# off-axis. Let's see an example of what the piston model outputs:


from bat_beamshapes import piston_in_infinite_baffle_directionality
import matplotlib.pyplot as plt
import numpy as np 

# %% 
# The piston in an infinite baffle has two parameters `k` and `a` 
# (for more on these)parameters see :doc:`here<../general_intro>`.

dB = lambda X: 20*np.log10(abs(X))

angles = np.linspace(0,np.pi/2,100)
k = 10.0
ka_values = np.array([1,3,5,10])
a_values = ka_values/k

parameters = {'k':k}

plt.figure()
a0 = plt.subplot(111, projection='polar')

for ka_v, a_v in zip(ka_values, a_values):
    parameters['a'] = a_v
    _, dirn = piston_in_infinite_baffle_directionality(angles, parameters)
    plt.plot(angles, dB(np.array(dirn)), label=str(ka_v))
    angles *= -1 # switch between L & R of 
plt.legend(title='ka')
plt.savefig('PIB.png')


#%%
#.. image:: ../_static/PIB.png
#  :width: 400
#  :alt: piston in infinite baffle





#%% References
#   ~~~~~~~~~~
#   
#  1. Mogensen, F., & Møhl, B. (1979). Sound radiation patterns in the frequency domain of cries from a Vespertilionid bat. Journal of comparative physiology, 134(2), 165-171.
#  2. Jakobsen, L., Ratcliffe, J. M., & Surlykke, A. (2013). Convergent acoustic field of view in echolocating bats. Nature, 493(7430), 93-96. 
#  3. Macaulay, J. D., Malinka, C. E., Gillespie, D., & Madsen, P. T. (2020). High resolution three-dimensional beam radiation pattern of harbour porpoise clicks with implications for passive acoustic monitoring. The Journal of the Acoustical Society of America, 147(6), 4175-4188.