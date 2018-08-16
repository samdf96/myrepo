#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:12:57 2018

@author: sfielder
"""

import clump_code_1 as cc
import header_printer as hp
import glob
import os
import yt
from definitions import OctantSplit
from definitions import MasterClumpMaker
from definitions import ClumpFinder
from definitions import ClumpBox
from definitions import BoundingBox
from yt.utilities.lib.misc_utilities import gravitational_binding_energy

from yt.utilities.physical_constants import \
    gravitational_constant_cgs as G

filename = '/mnt/bigdata/erosolow/Orion2/Fiducial00/data.0070.3d.hdf5'

## Imported Values
l = 10
clump_sizing = 30
cmin = 5.0e-21
step = 100

ds = yt.load(filename)

ad = ds.all_data()

octant = OctantSplit(ad,ds,l)

clumps = [] #Defining Empty List for loop

for i in range(0,len(octant)):

    master_clump_main = MasterClumpMaker(octant[i])
    cmax = octant[i]["gas", "density"].max()
    lc = ClumpFinder(master_clump_main,clump_sizing,cmin,cmax,step)
    for j in range(0,len(lc)):
        clumps.append(lc[j])

#%%


clump = clumps[8] # This will pull the 9th clump which gives error.

import pdb; pdb.set_trace()




