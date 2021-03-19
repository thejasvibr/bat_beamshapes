#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module which incorporates beam-shapes into call received levels 
when bat position and microphone positions are known. 

The default direction of the bat is 0 degrees, pointing North/upwards
of the graph. To simulate call-direction alterations, add a few degrees
here and there. 

@author: Thejasvi Beleyur
"""

import numpy as np 
import scipy.spatial as spatial
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from bat_beamshapes.beamshape_predictions import simple_cardiod

     
def calc_relangle(a,b,c):
    '''
     Here 'b' is the 'vertex', and c is the point +1 y unit above the 
     bat. 'a' is the microphone position.

    Parameters
    ----------
    a,b,c : TYPE
        DESCRIPTION.
    
    Returns
    -------
    relangle
    
    Notes
    -----
    thanks to @6502 https://stackoverflow.com/a/31334882

    '''
    theta = np.arctan2(c[1]-b[1],c[0]-b[0])-np.arctan2(a[1]-b[1],a[0]-b[0])
    return theta 

def angle_bet_bat_and_mic(batpos, micpos, call_direction=0):
    '''Outputs the 2D angle between the xy-axis, the bat
    and the microphone positions. 
    
    The angle calculated is the internal angle created
    by 3 points: a point in line with the call direction, 
    the bat itself, and the microphone position. 
    
    Parameters
    ----------
    batpos : 1x2 np.array
        XY position
    call_direction : float, optinal 
        Defaults to 0 radians (North)
    micpos : Nx2 np.array
        XY positions of 1-N mics
    Returns
    -------
    rel_angle : Nx1 np.array
        The relative angle in radians
    
    '''
    call_vector = np.concatenate((batpos,[0])).flatten()
    call_vector[1] += 1
    # now rotate the call_vector in the call direction
    rotmat = R.from_euler('z',call_direction)
    call_vector = rotmat.apply(call_vector)[:-1]
    
    rel_angle = np.apply_along_axis(calc_relangle,1,micpos,batpos,call_vector)
    return rel_angle

def calc_mic_level_nobeamshape(sl,batposition,micpositions,refdist=1):
    '''
    Calculates the received levels at different microphones, assuming the 
    bat is emitting an omnidirectional call.  

    Parameters
    ----------
    sl : float
        source level in dB SPL
    micpositions : Nx2 np.array
        Microphone postiions
    refdist : float, optional
        The reference distance for source level. The default is 1m. 

    Returns
    -------
    miclevels : Nx2 np.array
        The received levels at each microphone

    '''
    bat_and_micpos = np.row_stack((batposition, micpositions))
    r_mic = spatial.distance_matrix(bat_and_micpos,bat_and_micpos)[1:,0]
    
    miclevels = sl - 20*np.log10(r_mic/refdist)
    return miclevels

beamshapename_and_model = {'cardiod':simple_cardiod}


def sim_miclevel_beamshape(sl, call_direction,
                               batposition,micpositions,refdist=1,**kwargs):
    '''
    Calculate microphone received levels with beamshapes, in 2D

    Parameters
    ----------
    sl : float
        Source level in dB
    call_direction: 1xM array-like
        The call direction with reference to North (0 degrees).         
    batposition : 1x2 array-like
        XY positions
    micpositions : Nx2 np.array
        XY positions
    refdist : float, optional
        The reference distance for the source level. The default is 1.
    **kwargs : dictionary
        The keywords provide information on the beamshape pattern 
        and necessary parameters to use. 
        beamshape_model : str.
            Right now, only 'cardiod' is supported. More to come.
        model_params : dict.
            The parameters required for the beamshape model in question
            within a dictionary. 
    Returns
    -------
    miclevel_beamshape

    '''
    # first generate the received levels assuming omnidirectionality
    mic_levels_omni = calc_mic_level_nobeamshape(sl, 
                                                 batposition,
                                                 micpositions,refdist)
    # now calculate the relative angles bet bat and mic a
    rel_mic_angles = angle_bet_bat_and_mic(batposition, micpositions, 
                                           call_direction)
  
    beamshape_model = beamshapename_and_model[kwargs['beamshape_model']]
    rel_onaxis = beamshape_model(rel_mic_angles, **kwargs['model_params'])
    
    miclevel_beamshape = mic_levels_omni + rel_onaxis
    return miclevel_beamshape

if __name__=='__main__':
    print('yay')
    mic_positions = np.array([[-1,1],
                              [0,1],
                              [1,1],
                              [1.5,1.5],
                              [2,0.75],
                              [-1,-1]])
    bat_position = np.array([0,0.0])
    angles = angle_bet_bat_and_mic(bat_position, mic_positions)
    #%%
    plt.figure()
    plt.plot(bat_position[0],bat_position[1],'*')
    plt.plot(mic_positions[:,0],mic_positions[:,1],'*')
    plt.arrow(bat_position[0],bat_position[1],0,0.5,width=0.01)
    
    #%% Showing the effect of the asymmetry parameter
    thetas = np.linspace(0,2*np.pi,100)
    plt.figure()
    a0 = plt.subplot(111, projection='polar')
    a0.set_theta_zero_location('N')
    a0.set_theta_direction(1)
    plt.plot(thetas, simple_cardiod(thetas,A=5),label='A=5')
    plt.plot(thetas, simple_cardiod(thetas,A=10),label='A=10')
    plt.plot(thetas, simple_cardiod(thetas,A=20),label='A=20')
    plt.legend()

    # %% 
    rel = simple_cardiod(thetas,A=20)
    rel_x = np.cos(thetas)
    rel_y = np.sin(rel)
    plt.figure()
    plt.plot(thetas, rel)
    #%% with cardiod pattern
    calldirn = 0
    op = sim_miclevel_beamshape(100, calldirn, bat_position, 
                           mic_positions,beamshape_model='cardiod',
                           model_params={'A':7.3})
    print(op)