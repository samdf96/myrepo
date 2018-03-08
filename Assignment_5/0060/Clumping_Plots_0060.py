#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 19:35:29 2018

@author: sfielder
"""

# Coding Clumping and Fill-Null Areas

import yt
from yt.analysis_modules.level_sets.api import *

#Loading Data Set into Script
ds = yt.load("~/bigdata/Fiducial00/data.0060.3d.hdf5")

#All Data Object for Clumping Variable
ad = ds.all_data()

#Creating the first Clump
master_clump = Clump(ad, ("gas", "density"))

#Setting the Gravitationally Bound Indicator to the Clump
master_clump.add_validator("gravitationally_bound", use_particles=False)

#Clump Finding
c_min = ad["gas", "density"].min()
c_max = ad["gas", "density"].max()
step = 2.0
find_clumps(master_clump, c_min, c_max, step)

#Saving as a yt object
fn = master_clump.save_as_dataset(fields=["density", "particle_mass"])

#Extracting the Leaf Clumps
leaf_clumps = get_lowest_clumps(master_clump)

prj = yt.ProjectionPlot(ds, 2, 'density', center='c',width=(10,'pc'))



