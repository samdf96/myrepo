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



flist = glob.glob('/Users/sfielder/Documents/Astro_Data/**/*.hdf5', recursive=True)

filename = flist[0]

## Imported Values
l = 10
clump_sizing = 30
cmin = 5.0e-21
step = 100

ds = yt.load(filename)

ad = ds.all_data()

octant = OctantSplit(ad,ds,l)

clumps = [] #Defining Empty List for loop

for i in range(2,3):

    master_clump_main = MasterClumpMaker(octant[i])
    cmax = octant[i]["gas", "density"].max()
    lc = ClumpFinder(master_clump_main,clump_sizing,cmin,cmax,step)
    for j in range(0,len(lc)):
        clumps.append(lc[j])
#%%


clump = clumps[0]

bregion = BoundingBox(clump)
data_object_clump = ClumpBox(ds,bregion)


bulk_velocity = data_object_clump.quantities.bulk_velocity(use_gas=True)
x_arr = data_object_clump['x']
y_arr = data_object_clump['y']
z_arr = data_object_clump['z']

kinetic = 0.5 * (data_object_clump['cell_mass'] *
        ((bulk_velocity[0] - data_object_clump["gas", "velocity_x"])**2 +
         (bulk_velocity[1] - data_object_clump["gas", "velocity_y"])**2 +
          (bulk_velocity[2] - data_object_clump["gas", "velocity_z"])**2)).sum()


# Units won't make sense because gravitational_binding_energy returns as just 
# value with no associated unit value to make it a quantity
grav_bind = G * gravitational_binding_energy(data_object_clump['cell_mass'],
                                             x_arr,
                                             y_arr,
                                             z_arr,
                                             True,
                                             (kinetic/G).in_cgs())
if grav_bind.value > kinetic.value:
    gravitationally_bound = True
else:
    gravitationally_bound = False
