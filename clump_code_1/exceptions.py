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