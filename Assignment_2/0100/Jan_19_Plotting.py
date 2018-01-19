#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 19:53:35 2018

@author: sfielder
"""

# Import Packages here.
import yt as yt

# Call different data sets here, possibly create loop.
ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")

#Loading all data into one object.
ad = ds.all_data()

#Printing Derived Variables to be able to use (comment in if needed)
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''


# Basic Projection Plots

#Density Projection along z
yt.ProjectionPlot(ds, "z", "density").save("z_mass_density_plot.png")

#Angular Momentum (x,y) Projected along z
z_projected_angular_momentum_y = yt.ProjectionPlot(ds, "z", "angular_momentum_y",data_source=ad)
z_projected_angular_momentum_y.set_cmap(field="angular_momentum_y", cmap='bwr')
z_projected_angular_momentum_y.save("z_projected_angular_momentum_y.png")
z_projected_angular_momentum_x = yt.ProjectionPlot(ds, "z", "angular_momentum_x",data_source=ad)
z_projected_angular_momentum_x.set_cmap(field="angular_momentum_x", cmap='bwr')
z_projected_angular_momentum_x.save("z_projected_angular_momentum_x.png")

#Velocity along z direction of z-velocity
z_projected_z_velocity = yt.ProjectionPlot(ds, "z", "velocity_z")
z_projected_z_velocity.set_cmap(field="velocity_z", cmap='bwr')
z_projected_z_velocity.save("z_velocity_density_plot.png")

#Basic Probability Density Plots

plot_1D_PDF_Density = yt.ProfilePlot(ad, "density", "ones", weight_field=None)
plot_1D_PDF_Density.save("1D_PDF_Density.png")