#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 06:54:23 2018
Copied from github: 
    https://github.com/low-sky/py-low-sky/blob/master/plane_fit.py
@author: Erik Rosolowski
"""

from scipy.optimize import least_squares as lsq
import numpy as np

def myplane(p, x, y, z):
    return(p[0] + p[1] * x + p[2] * y - z)


def plane_fit(x, y, z, robust=False):
    """
    Fits a plane to data without any given uncertainties or weighting.
    Arguments:
    ----------
    x, y: float
       x and y coordinates of the data
    z: float
       z-coordinates to which the values are fit
    Returns:
    --------
    coefficients:
       3-element vector with components [ z0 (constant offset) , grad_x, grad_y]
    
    """

    x0, y0 = np.median(x), np.median(y)
    dataz = np.c_[np.ones(x.size), 
                  x-x0, 
                  y-y0]

    lsqcoeffs, residuals, _, _ = np.linalg.lstsq(dataz,z)
    if robust:
        outputs = lsq(myplane, np.r_[lsqcoeffs],
                      args=([x-x0,
                             y-y0, z]),
                      loss = 'soft_l1')
        lsqcoeffs = outputs.x
    return(lsqcoeffs,residuals)