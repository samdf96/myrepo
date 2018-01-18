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
#Slice Plot with Velocity Vectors overlayed
'''
plot_vector = yt.SlicePlot(ds, 'x', "density")
plot_vector.annotate_velocity(factor = 16)
plot_vector.save("x_slice_plot_vector_field_overlay.png")
'''

#Slice plot with density contours overlayed
'''
plot_contour = yt.SlicePlot(ds, 'x', "density")
plot_contour.annotate_contour("density")
plot_contour.save("x_slice_plot_density_contours.png")
'''

#Projection Plots x-direction.
'''
#Creating Projection Plots in x-direction
#print(sorted(ds.field_list))
print(sorted(ds.derived_field_list))
yt.ProjectionPlot(ds, "x", "density").save("x_mass_density_plot.png")
'''
'''
plot_x_momentum_density = yt.ProjectionPlot(ds, "x", "X-momentum")
plot_x_momentum_density.set_cmap(field="X-momentum", cmap='bwr')
plot_x_momentum_density.save("x_momentum_density_plot.png")
'''


#Creating 1D Probability Density functions
# Create a data object that represents the whole box.
ad = ds.all_data()
#Creating my Speed Derived Variable for Computations and Plotting

def rms_speed(field, data):
    return ((data['gas', 'velocity_x']**2) + (data['gas', 'velocity_y']**2) + (data['gas', 'velocity_z']**2))**(1/2)

ds.add_field(("gas", "rms_speed"), units="cm/s", function=rms_speed)

#Re-defining all data to include added rms_speed
#Not sure if I need to do this.

ad = ds.all_data()

#Plotting 2D PDF
'''
plot_2D_PDF = yt.PhasePlot(ad, "density", "rms_speed", "cell_mass",
                    weight_field=None, fractional=True)

plot_2D_PDF.save("2D_PDF_Density_RMS_Speed.png")
'''

#Plotting 1D PDF
'''
plot_1D_PDF_Density = yt.ProfilePlot(ad, "density", "ones", weight_field=None)
plot_1D_PDF_Density.save("1D_PDF_Density.png")
plot_1D_PDF_Energy = yt.ProfilePlot(ad, "kinetic_energy", "ones", weight_field=None)
plot_1D_PDF_Energy.save("1D_PDF_Energy.png")
'''
#Plotting Mean Velocity Projection Weighted
'''
plot_weighted_projection_z = yt.ProjectionPlot(ds,"z","velocity_z",weight_field="cell_mass")
plot_weighted_projection_z.set_cmap(field="velocity_z", cmap='bwr')
plot_weighted_projection_z.save("2D_Projection_Weighted_RMS_Velocity_z.png")
'''

#In x and y directions

'''
plot_weighted_projection_x = yt.ProjectionPlot(ds,"x","velocity_x",weight_field="cell_mass")
plot_weighted_projection_x.set_cmap(field="velocity_x", cmap='bwr')
plot_weighted_projection_x.save("2D_Projection_Weighted_RMS_Velocity_x.png")
'''

#Computations and Data Analysis

#Computing Total Mass in Domain
#print("Total Mass of Gas in Domain = ", ad.quantities.total_quantity([("gas", "cell_mass")]))
#print("Total Kinetic Energy in Domain = ", ad.quantities.total_quantity([("gas", "kinetic_energy")]))
#print("Total Momentum in x-direction in Domain = ", ad.quantities.total_quantity([("gas", "momentum_x")]))

# 3D Rendering Plot just for fun

#im, sc = yt.volume_render(ds, field=('gas', 'density'))
'''
plot_density_dusk = yt.ProjectionPlot(ds, "z", "density")
plot_density_dusk.set_cmap(field="density", cmap='dusk')
plot_density_dusk.save("plot_density_dusk.png")
'''




#Projection Plots of Angular Momentum
dd = ds.sphere([0,0,0],(4,'pc'))
dbox = ds.box([0,0,0],[5,5,5])


z_projected_angular_momentum_x_box = yt.ProjectionPlot(ds,"z","angular_momentum_y",data_source=dbox)
z_projected_angular_momentum_x_box.save("z_projected_angular_momentum_x_box.png")
'''
z_projected_angular_momentum_y = yt.ProjectionPlot(ds, "z", "angular_momentum_y",data_source=dd)
z_projected_angular_momentum_y.set_cmap(field="angular_momentum_y", cmap='bwr')
z_projected_angular_momentum_y.save("z_projected_angular_momentum_y.png")
z_projected_angular_momentum_x = yt.ProjectionPlot(ds, "z", "angular_momentum_x",data_source=dd)
z_projected_angular_momentum_x.set_cmap(field="angular_momentum_x", cmap='bwr')
z_projected_angular_momentum_x.save("z_projected_angular_momentum_x.png")
'''
