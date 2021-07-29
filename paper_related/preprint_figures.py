#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes a comparison of the piston in a sphere and piston in an infinite baffle.
The plot is meant to go in the preprint associated with the paper
"""
import numpy as np 
import matplotlib.pyplot as plt
from beamshapes import piston_in_infinite_baffle as piston_infbaf
from beamshapes import piston_in_sphere_flint as piston_sphere
from flint import acb

#%% set up the infinite baffle directivity on the left and in sphere on the right
ka_values  = [1, 5, 10]
a = 0.01 # piston of 1cm radius 


angles = np.linspace(0,-np.pi,100)
piston_infbaf_D = np.zeros((angles.size, len(ka_values)))



for j,each_ka in enumerate(ka_values):
    k_value = each_ka/a
    parameters = {'k':k_value, 'a':a}
    _, piston_infbaf_D[:,j] = piston_infbaf.piston_in_infinite_baffle_directivity(angles,
                                                                               parameters)
    
#%% 
# Now calculate the piston in sphere 
pistonsphere_directivity = piston_sphere.piston_in_sphere_directivity
angles_as_acb = [acb(each) for each in angles]

R = acb(0.1) # radius of a sphere
alpha = acb.pi()/acb(6)
a_value = R*acb.sin(alpha)
piston_in_sphere_D = np.zeros(piston_infbaf_D.shape)

for j,each_ka in enumerate(ka_values):
    k_value = acb(each_ka)/a_value
    parameters = {'k':k_value, 'a':a_value, 'R': R, 'alpha':alpha}
    A_n, piston_in_sphere_D[:,j] = pistonsphere_directivity(angles_as_acb,
                                                                    parameters)




#%%

line_types = ['dotted', 'dashed', 'solid']

plt.figure(figsize=(4,3.5))
a0 = plt.subplot(111, projection='polar')
for j,ka in enumerate(ka_values):
    plt.plot(angles, piston_infbaf_D[:,j], color='k', linestyle=line_types[j])
    plt.plot(-angles, piston_in_sphere_D[:,j], color='k', label=f'$ka$={ka}', linestyle=line_types[j]) 
a0.set_theta_zero_location('N')
plt.yticks([-12,-24,-36,-48], fontsize=6)
a0.set_rlabel_position(0)  # Move radial labels away from plotted line
tick_angles = [0,45,90,135, 180,225,270,315]
tick_labels = []
for each in tick_angles:
    if each==135:
        tick_labels.append('')
    else:
        tick_labels.append(f'{each}'+'$^{\circ}$')
plt.xticks(np.radians(tick_angles) ,tick_labels , fontsize=6)
#plt.legend()
plt.legend(loc=(-0.175,-0.07), frameon=False)
plt.text(0.1, 1.02, 'piston in a sphere', transform=a0.transAxes, fontsize=8)
plt.text(0.55, 1.02, 'piston in an infinite baffle', transform=a0.transAxes,
                         multialignment='right', fontsize=8)
plt.tight_layout()
plt.savefig('piston_sphere_baffle.png')