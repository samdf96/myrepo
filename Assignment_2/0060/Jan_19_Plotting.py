#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 19:53:35 2018

@author: sfielder
"""

# Import Packages here.
import yt as yt

# Call different data sets here, possibly create loop.
ds = yt.load("~/bigdata/Fiducial00/data.0060.3d.hdf5")

#Loading all data into one object.
ad = ds.all_data()

#Printing Derived Variables to be able to use (comment in if needed)
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''

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
'''



import matplotlib.pyplot as plt

# This creates a projection weighting by density
vz = ad.integrate('velocity_z',weight='density',axis='z')

# The projection has fields 'px' and 'py' for the position
# and whatever quantity you average over like velocity_z.

plt.scatter(vz['px'],vz['py'],c=vz['velocity_z'])
plt.colorbar()
plt.show()


arr = vz['velocity_z']
arr.shape = (256,256)  # Since the underlying data are 256^3
import numpy as np
arr = np.array(arr)  # Make this into a plain numpy array if needed.
plt.imshow(arr)      
plt.show()

#Generate Subsections of Arrays (16 regions of 64x64).
def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))


#This creates a list of arrays, specifically a list of 16 subregions.
    # This makes the index number go left to right and top to bottom.
    
array_split_size = 64 #Value for blocks 
arr_reshape = blockshaped(arr,array_split_size,array_split_size)
#Callable Sizing Value
array_size = np.size(arr_reshape[0])
#Plane Fitting Process
    
    #Creating x,y coordintes for blocks of data
x_coordinates = np.arange(0,array_size)
y_coordinates = np.arange(0,array_size)
z_velocity_mapping = np.empty(4096)


z_velocity_mapping_sec_1 = np.ndarray.flatten(arr_reshape[0])
z_velocity_mapping_sec_2 = np.ndarray.flatten(arr_reshape[1])
z_velocity_mapping_sec_3 = np.ndarray.flatten(arr_reshape[2])
z_velocity_mapping_sec_4 = np.ndarray.flatten(arr_reshape[3])



#Importing Plane Fit Script here.
from plane_fit import plane_fit

#Making z of zeros for plane fit script

first_plane_sec_1 = plane_fit(x_coordinates,y_coordinates,z_velocity_mapping_sec_1)
first_plane_sec_2 = plane_fit(x_coordinates,y_coordinates,z_velocity_mapping_sec_2)
first_plane_sec_3 = plane_fit(x_coordinates,y_coordinates,z_velocity_mapping_sec_3)
first_plane_sec_4 = plane_fit(x_coordinates,y_coordinates,z_velocity_mapping_sec_4)

