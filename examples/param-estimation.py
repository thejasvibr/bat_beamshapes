'''
Source parameter estimation
===========================
Sounds recorded in nature are often directional (bird songs, bat calls, etc.). This
directionality means a set of observations (reduction in call level) may be only due
to the animal turning one way or the other - and not necessarily because of an actual 
change in sound production. Fitting a `source model <../general_intro.rst>` helps in accounting for 
sound directionality and actually being able to estimate many other aspects of sound production. 

Let's simulate an example to show the power of actually using a source model to study 
animal sound production. Of course, any type of sound can be studied this way!! 

'''
import matplotlib.pyplot as plt
import numpy as np 


#%% 
# Let's simulate a bat that emits five calls as it flies through a room with 12 microphones. 
# We'll pretend our simulated bat is a flying piston in a sphere for now. In general, here we'll
# simulate the following progressive trends 1) reduction in call level 2) 'widening' of the beam. 
# These two trends mimic what a bat coming in to capture an insect will do. 

# Define mic positions  in 2D at flight height of the bat
mic_posns = np.zeros((12,2))
# define x-positions
mic_posns[:4,0] = 0 # first three mics on left wall
mic_posns[4:8,0] = np.linspace(0.5,2,4)
mic_posns[-4:,0] = 2.5

# define y-positions 
mic_posns[:4,1] = np.linspace(0,2,4)
mic_posns[4:8,1] = 2.5
mic_posns[-4:,1] = np.linspace(0,2,4)[::-1]

# Define on-axis level for each of the 5 calls (dB SPL re 20 muPa at 1m)
# These levels are semi-realistic!
onaxis_level = [100, 96, 90, 84, 84]

# Define the directivity function broadly by altering ka for each call
overall_ka = [10, 5, 4, 3, 2]

# Assign 
