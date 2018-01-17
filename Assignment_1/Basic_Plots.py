#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 05:01:24 2018

@author: sfielder
"""

import yt as yt
#Loading Data from file.
ds = yt.load("~/bigdata/data.0453.3d.hdf5")
ds2 = yt.load("~/bigdata/chk.0453.3d.hdf5")


#Printing Derived Variables to be able to use
for i in sorted(ds.derived_field_list):
    print(i)



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
# Create a data object that represents the whole box.
ad = ds.all_data()
#Creating my Speed Derived Variable for Computations and Plotting

def rms_speed(field, data):
    return ((data['gas', 'velocity_x']**2) + (data['gas', 'velocity_y']**2) + (data['gas', 'velocity_z']**2))**(1/2)

ds.add_field(("gas", "rms_speed"), units="cm/s", function=rms_speed)

#Re-defining all data to include added rms_speed

ad = ds.all_data()

#Plotting 2D PDF
'''
plot_2D_PDF = yt.PhasePlot(ad, "density", "rms_speed", "cell_mass",
                    weight_field=None, fractional=True)

plot_2D_PDF.save("2D_PDF_Density_RMS_Speed.png")
'''

#Plotting 1D PDF
plot_1D_energy = yt.ProfilePlot(ad,"kinetic_energy","density")
plot_1D_energy.save("1D_PDF_Energy.png")
plot_1D_mass = yt.ProfilePlot(ad,"cell_mass","density")
plot_1D_mass.save("1D_PDF_Mass.png")

#Plotting Mean Velocity Projection Weighted
'''
plot_weighted_projection_z = yt.ProjectionPlot(ds,"z","rms_speed",weight_field="cell_mass")
plot_weighted_projection_z.save("2D_Projection_Weighted_RMS_Speed_z.png")
plot_weighted_projection_x = yt.ProjectionPlot(ds,"x","rms_speed",weight_field="cell_mass")
plot_weighted_projection_x.save("2D_Projection_Weighted_RMS_Speed_x.png")
plot_weighted_projection_y = yt.ProjectionPlot(ds,"y","rms_speed",weight_field="cell_mass")
plot_weighted_projection_y.save("2D_Projection_Weighted_RMS_Speed_y.png")
'''
#Computations and Data Analysis

#Computing Total Mass in Domain
#print("Total Mass of Gas in Domain = ", ad.quantities.total_quantity([("gas", "cell_mass")]))
#print("Total Kinetic Energy in Domain = ", ad.quantities.total_quantity([("gas", "kinetic_energy")]))
#print("Total Momentum in x-direction in Domain = ", ad.quantities.total_quantity([("gas", "momentum_x")]))



