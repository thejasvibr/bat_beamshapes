#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions 
=================
Functions that helps with formatting shape and converting between object types.

"""
import mpmath
import numpy as np 


dB = lambda X : 20*np.log10(np.abs(X))

def args_to_str(**args):
    '''
    Converts dictionary entries into string.
    '''
    return {key: str(value) for key,value in args.items()}

def args_to_mpmath(**args):
    '''
    Converts dictionary entries into mpmath objects (mpmath.mpf or mpmath.mpc)
    '''
    return {key: mpmath.mpmathify(value) for key,value in args.items()}

