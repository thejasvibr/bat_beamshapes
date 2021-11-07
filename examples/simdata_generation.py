'''
Source parameter estimation - generating the data
=================================================
Sounds recorded in nature are often directional (bird songs, bat calls, etc.). This
directionality means a set of observations (reduction in call level) may be only due
to the animal turning one way or the other - and not necessarily because of an actual 
change in sound production. Fitting a `source model <../general_intro.rst>` helps in accounting for 
sound directionality and actually being able to estimate many other aspects of sound production. 

Let's make simulated data to show the power of actually using a source model to study 
animal sound production. Of course, any type of sound can be studied this way!! 

This example will go through the following steps:

    1. Define mic and bat positions 
    2. Define call direction, on-axis level and parameters of the source model
    3. Calculate received levels at mics assuming an omnidirectional source
    4. Calculate mic to call-direction relative angles 
    5. Generate beamshape from source model 
    6. Add/subtract the relative change for off-axis angles. 
    

'''
# sphinx_gallery_thumbnail_path = '_static/simulated_data.png'

#%% the if __name__ == '__main__' is required when running on Windows systems
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    import numpy as np 
    import pandas as pd
    import beamshapes.piston_in_sphere
    import beamshapes.sim_reclevels as simlevels
    import mpmath

#%% 
# Define mic and bat positions 
# ----------------------------
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

#%% 
# Define call direction, on-axis level and parameters of the source model
# -----------------------------------------------------------------------
    # Define on-axis level for each of the 5 calls (dB SPL re 20 muPa at 1m)
    # These levels are semi-realistic!
    onaxis_levels = [100, 96, 90, 84, 84]
    
    # Define the directivity function broadly by altering ka for each call
    # Set 'k' for a sound of 80 kHz frequency
    overall_ka = [5, 4, 3, 2, 1]
    
#%%
# Where was the bat at each of the 5 call emissinos?
    flight_path_x =  np.array([0.45, 1.0, 1.5, 2.0, 2.4])
    flight_path_y = np.array([0.3, 0.8, 1.2, 0.6, 0.4])
    flight_path = np.column_stack((flight_path_x, flight_path_y))
    
    # Where was the bat aiming it's call? 
    call_directions = np.pi/2-np.deg2rad(np.array([15, 60, 140, 160, 200]))
    
    
    plt.figure(figsize=(5,4))
    plt.plot(flight_path[:,0], flight_path[:,1], '*')
    plt.plot(mic_posns[:,0], mic_posns[:,1],'r*')
    
    for i, direction in enumerate(call_directions):
        arrow_dx, arrow_dy = np.cos(direction), np.sin(direction)
        arrow_dx *= 0.5
        arrow_dy *= 0.5
        plt.arrow(flight_path[i,0], flight_path[i,1],arrow_dx, arrow_dy, head_width=0.05)
        plt.text(flight_path[i,0]-0.5, flight_path[i,1]+0.1,f'direction: {np.round(np.degrees(call_directions[i]),2)}')
        plt.text(flight_path[i,0]-0.1, flight_path[i,1]+0.2,f'SL: {onaxis_levels[i]}')
        plt.text(flight_path[i,0]-0.1, flight_path[i,1]+0.3,f'ka: {overall_ka[i]}')
    plt.savefig('../docs/source/_static/simulated_data.png')

#%%
#.. image:: ../_static/simulated_data.png
#  :width: 400
#  :alt: simulated data

#%% 
#
# `An overview of the call parameters with respect to mic positions (red stars). Call direction is shown by an arrow, while
# other parameters are shown in text`

#%%
# Calculate received levels at mics assuming an omnidirectional source
# --------------------------------------------------------------------
    
    received_omnidirn_level = np.zeros((mic_posns.shape[0], 5))
    for i, call_level in enumerate(onaxis_levels):
        received_omnidirn_level[:,i] = simlevels.calc_mic_level_nobeamshape(call_level,
                                                           flight_path[i,:].reshape(-1,2),
                                                           mic_posns)
#%%
# Calculate mic to call-direction relative angles 
# -----------------------------------------------
# Let's calculate the relative mic to bat angles to get the relative addition or reduction
# of call level. 
    comparitive_angle = np.zeros((mic_posns.shape[0], 5))
    for i, call_direction in enumerate(call_directions):
        comparitive_angle[:,i] = simlevels.angle_bet_bat_and_mic(flight_path[i,:],
                                                                 mic_posns,
                                                                     call_directions[i])


#%% 
# Generate beamshape from source model 
# ------------------------------------
# Having calculated the received levels assuming an omnidirectional call, let's now
# include the beamshapes arising from a piston in a sphere. Let's start by calculating the 
# relative bat to mic 
    alpha = mpmath.pi/8 # half-aperture of 22.
    R = 6e-3
    paramv = {}
    paramv['R'] = R
    paramv['alpha'] = alpha
    
    a_calc = paramv['R']*mpmath.sin(paramv['alpha'])
    k = overall_ka[0]/a_calc
    
    
    
    all_beamshapes = []
    for i, ka in enumerate(overall_ka):
        k = overall_ka[i]/a_calc
        paramv['k'] = k
        paramv['a'] = a_calc
        _, relative_levels = beamshapes.piston_in_sphere_directivity(comparitive_angle[:,i], paramv)
        all_beamshapes.append(relative_levels)

#%%
    relative_off_axis = np.zeros((mic_posns.shape[0],5))
    for i in range(5):
        relative_off_axis[:,i] += np.array(all_beamshapes[i]).flatten()
#%%
# Add/subtract the relative change for off-axis angles 
# ----------------------------------------------------
# save the relative levels
    relative_off_axis_csv = pd.DataFrame(relative_off_axis)
    relative_off_axis_csv.to_csv('../docs/source/_static/call_beamshapes_dthetadzero-db.csv')
    
    
    # relative_off_axis = pd.read_csv('../docs/source/_static/call_beamshapes_dthetadzero-db.csv')
    # relative_off_axis = relative_off_axis.to_numpy()
    
    # include the directivity into the omnidirectional model received levels 
    
    received_with_beamshapes = np.zeros((received_omnidirn_level.shape))
    for i in range(5):
        received_with_beamshapes[:,i] = received_omnidirn_level[:,i] + relative_off_axis[:,i]

#%% 
# Let's now save the parameters used to generate the data along with the final 
# simulated data. 

# Received levels
pd.DataFrame(received_with_beamshapes).to_csv('../docs/source/_static/bat_as_pistoninsphere.csv')
# Original parameters: ka, source level and call position

call_data = pd.DataFrame(data={'ka':overall_ka, 'source_level':onaxis_levels,
                                   'call_directions':call_directions})

position_data = pd.DataFrame(data={'position_x':flight_path[:,0],
                                   'position_y':flight_path[:,1]})

position_data.to_csv('../docs/source/_static/simdata_positions.csv')
call_data.to_csv('../docs/source/_static/simdata_calls.csv')
