#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 05:01:24 2018

@author: sfielder
"""

import yt as yt

#Loading Data from file.
ds = yt.load("~/bigdata/data.0453.3d.hdf5")

#Creating Slice Plots
'''
yt.SlicePlot(ds, 'x', "density").save("x_slice_plot.png")
yt.SlicePlot(ds, 'y', "density").save("y_slice_plot.png")
yt.SlicePlot(ds, 'z', "density").save("z_slice_plot.png")
'''

'''
#Creating Projection Plots in x-direction
#print(sorted(ds.field_list))
print(sorted(ds.derived_field_list))
yt.ProjectionPlot(ds, "x", "density").save("x_mass_density_plot.png")
yt.ProjectionPlot(ds, "x", "X-momentum").save("x_momentum_density_plot.png")
'''


#Creating 1D Probability Density functions

