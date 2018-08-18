#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 10:14:38 2018

@author: sfielder
"""

class Error(Exception):
   """Base class for other exceptions"""
   pass

class YTErrorReshape(Error):
    """
    Coordinate of the velocity array cannot be reshaped into the 256x256
    structure of the original data simulation file
    """
    pass

class YTErrorValue(Error):
    """
    Signifies that velocity_array_reducer function detected a length of 0
    on the v_positions variable, thus computations will break.
    """
    pass

class YTRuntimeError(Error):
    """
    Signifies that the clump found has no dimensional length along one of its
    edges. This is made clear in the FITS file itself.
    """
    pass
    
class YTPassThrough(Error):
    """
    Signifies we need to pass through code to not kick any more errors.
    """
    pass

class YTNoClump(Error):
    '''
    Happens when YT finds a clump, but that clump has no data associated 
    with it, thus program crashes out instantly.
    '''
    pass