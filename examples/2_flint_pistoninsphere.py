# -*- coding: utf-8 -*-
"""
FLINT piston in sphere
======================
This example shows how to run the FLINT piston in sphere implementation. 
The ```python-flint``` package is the python port for the FLINT package in C, 
and runs extremely fast in comparison to ```mpmath```. However, it does
have its own conventions. This run-through will highlight. 

Important
~~~~~~~~~
The Python ```python-flint``` package currently only works on Linux. See [here](https://github.com/fredrik-johansson/python-flint/)
and [here too](https://fredrikj.net/python-flint/).
"""

import matplotlib.pyplot as plt
import numpy as np 
from beamshapes import piston_in_sphere_flint as piston_sphere # to speed things up!!
pistonsphere_directivity = piston_sphere.piston_in_sphere_directivity
from flint import acb # import the python-flint package

#%% Let's calculate the directivity for 3 equivalent `ka` values.
ka_values  = [1, 5, 10]
a = 0.01 # piston of 1cm radius 
angles = np.linspace(0,-np.pi,100) # angles from 0-180 degrees

#%% The ````flint``` implementation of piston in a sphere only accepts an ```acb``` 
# complex number (see `here <https://fredrikj.net/python-flint/acb.html?highlight=acb>`),
# and not an inbuilt float or np.float object.  Also note how even functions
# like `sin` and :math:`\pi` are called from the ```acb``` module (```acb.pi```, ```acb.sin```). 

angles_as_acb = [acb(each) for each in angles] 
R = acb(0.1) # radius of a sphere
alpha = acb.pi()/acb(6)
a_value = R*acb.sin(alpha)
piston_in_sphere_D = np.zeros(angles.size)

for j,each_ka in enumerate(ka_values):
    k_value = acb(each_ka)/a_value
    parameters = {'k':k_value, 'a':a_value, 'R': R, 'alpha':alpha}
    A_n, piston_in_sphere_D[:,j] = pistonsphere_directivity(angles_as_acb,
                                                                    parameters)

line_types = ['dotted', 'dashed', 'solid']
plt.figure()
a1 = plt.subplot(111, projection='polar')
for j,ka in enumerate(ka_values):
    plt.plot(angles, piston_in_sphere_D[:,j], color='k', linestyle=line_types[j])
    plt.plot(-angles, piston_in_sphere_D[:,j], color='k', label=f'$ka$={ka}', linestyle=line_types[j]) 
a1.set_theta_zero_location('N')
plt.yticks([-12,-24,-36,-48], fontsize=6)
a1.set_rlabel_position(0)  # Move radial labels away from plotted line
tick_angles = [0,45,90,135, 180,225,270,315]
tick_labels = []
for each in tick_angles:
    tick_labels.append(f'{each}'+'$^{\circ}$')


plt.yticks([-12,-24,-36, -48],['',-24,'', -48], fontsize=6)        

plt.xticks(np.radians(tick_angles) ,tick_labels , fontsize=6)
#a1.tick_params(axis='both', pad=0.5)
a1.xaxis.set_tick_params(pad=-3.5)
plt.legend(loc=(-0.27,-0.23), frameon=False, fontsize=7)
