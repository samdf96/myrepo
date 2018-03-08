#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 03:45:14 2018

@author: sfielder
"""

#Advanced Plotting Script for Experimentation

import yt as yt
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
ds = yt.load("~/bigdata/Fiducial00/data.0060.3d.hdf5")
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''
# Make a density projection.
p_annotate = yt.ProjectionPlot(ds, "z", "density")
p_annotate.annotate_grids()
p_annotate.save('z_density_annotated_grids.png')