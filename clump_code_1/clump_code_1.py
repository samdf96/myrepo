#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 19:53:52 2018

@author: sfielder
"""

#Importing packages here
import yt
import numpy as np
import astropy
from yt.analysis_modules.level_sets.api import *
import math as m
import sys


#filename = sys.argv[1]

#Defining Global Variables Here
l=10	#Length of Original Data Set
cmin = 5e-21	#Minimum Density Threshold for Clump Finding
cmax = 5e-22	#Maximum Density Threshold for Clump Finding
step = 1000	#Step-size multiplier for Clump Finding


def octant_split(data_object):
    """
    Splits cubic data object into a sub-set of 8 data objects, i.e. octants.
    	
    Input: Data_object (usually representing whole simulation space)
    Output: List of 8 data-objects.
    Given Global Parameters: Length of data_object along the side, in pc.
    """
    dbox_array = [] #Creating Empty List to Store Data Objects
    #Creating Boundaries for Octants Here
    x1 = np.array((-l/2,0,-l/2,0,-l/2,0,-l/2,0))
    x2 = np.array((0,l/2,0,l/2,0,l/2,0,l/2))
    y1 = np.array((-l/2,-l/2,0,0,-l/2,-l/2,0,0))
    y2 = np.array((0,0,l/2,l/2,0,0,l/2,l/2))
    z1 = np.array((-l/2,-l/2,-l/2,-l/2,0,0,0,0))
    z2 = np.array((0,0,0,0,l/2,l/2,l/2,l/2))
    for i in range(0,8):
        dbox_array.append(ds.r[(x1[i],'pc'):(x2[i],'pc'), (y1[i],'pc'):(y2[i],'pc'), (z1[i],'pc'):(z2[i],'pc')])	
    return(dbox_array)

def master_clump_maker(data_object):
    '''
    Creates master clump for imputted data object.
    Input: data_object
    Output: master_clump object
    '''
    #Defining Field for Contouring
    field = ("gas", "density")
    master_clump = Clump(data_object, field)
    return(master_clump)

def clump_finder(master_clump,clump_sizing):
    '''
    Finds clumps that match certain criteria.
    Input: master_clump object, clump sizing parameter
    Output: Lowest Tree Values (leaves) of the clump data.
    Given Global Variables: c_min, c_max, step.
    '''
    #Setting Parameters for Clump Finding
    master_clump.add_validator("min_cells", clump_sizing)
    master_clump.add_info_item("center_of_mass") #Adds Center of Mass info for Clumps
    #Finding Clumps Here
    find_clumps(master_clump, cmin, cmax, step)
    lc = get_lowest_clumps(master_clump)
    return(lc)

def center_of_mass(lowest_clump):
    com = lowest_clump.quantities.center_of_mass()
    return(com)

ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")
ad = ds.all_data()
octants = octant_split(ad)
master_clump_1 = master_clump_maker(octants[0])
clump_1 = clump_finder(master_clump_1,150)
#center_of_mass_1 = center_of_mass(clump_1)


