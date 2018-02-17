#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 19:53:35 2018

@author: sfielder
"""
# This script will create the following:
#   Projection Plots:
#                   Mass Density along z-axis
#                   Angular momentum (x and y) along z-axis
#                   Velocity (z direction) along z-axis
#   PDF's:
#                   Mass Density
#
# This script will also map the velocity field onto a plane to extract
# the x and y gradients, and their magnitude is extracted and saved
# in a csv which can be ported into Latex directly into a table
# which follows exactly the subregions that the original data was split to

'''
Will need to set a few values in the script for differing data sets:
    Filename for importing data
    Data Length (Number of pixels along a given axis for original data set)
    Array_split_size - number of pixels per subregion along a given axis
    Array_split_size_sub - number of subregions in total
'''


# Import Packages here.
import yt as yt
import matplotlib.pyplot as plt
import numpy as np

# Call different data sets here, possibly create loop.
ds = yt.load("~/bigdata/Fiducial00/data.0060.3d.hdf5")

#Loading all data into one object.
ad = ds.all_data()

#Printing Derived Variables to be able to use (comment in if needed).

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
'''
##############################################################################

#Analysis of Data Section


# This creates a projection weighting by density
vz = ad.integrate('velocity_z',weight='density',axis='z')

# The projection has fields 'px' and 'py' for the position.

plt.scatter(vz['px'],vz['py'],c=vz['velocity_z'])
plt.colorbar()
plt.show()

#Defining Value for underlying data set length
data_length = 256
arr = vz['velocity_z']
arr.shape = (data_length,data_length)  # Since the underlying data are 256^3

arr = np.array(arr)  # Make this into a plain numpy array if needed.
plt.imshow(arr)      
plt.show()

#Generate Subsections of Arrays.
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

#This creates a list of arrays, specifically a list of 'x' subregions.
# This makes the index number go left to right and top to bottom.
 
array_split_size = 64 #Value for blocks
array_split_size_sub = 16
array_split_size_sub_tot = array_split_size // array_split_size_sub
arr_reshape = blockshaped(arr,array_split_size,array_split_size)
#Callable Sizing Value
array_size = np.size(arr_reshape[0])
#Plane Fitting Process

#Creating x,y coordintes for blocks of data
xvalues = np.arange(0,data_length);
yvalues = np.arange(0,data_length);

xx, yy = np.meshgrid(xvalues, yvalues)

xx = blockshaped(xx,array_split_size,array_split_size)
yy = blockshaped(yy,array_split_size,array_split_size)

#Importing Plane Fit Script here.
from plane_fit import plane_fit

#Making z of zeros for plane fit script
first_plane_sec = np.zeros((array_split_size_sub,1,3))


for i in range(0,array_split_size_sub):
    first_plane_sec[i] = plane_fit(np.ndarray.flatten(xx[i]),np.ndarray.flatten(yy[i]),np.ndarray.flatten(arr_reshape[i]))

#Need Conversion for values in km/s/pc instead of cm/s/pixel
conv_length = 1/(10000) #cm to km
conv_distance = 1/0.039 #pixel to pc

#Grabs matrix / row / column.


for i in range(0,array_split_size_sub):
    for j in range(1,3):
        first_plane_sec[i,0,j] = first_plane_sec[i,0,j]*conv_length*conv_distance
    
#Now computing magnitudes of gradient into empty superset array as defined below

gradients = np.zeros((16,1))

for i in range(0,array_split_size_sub):
    gradients[i,0] = ((first_plane_sec[i,0,1]**2) + (first_plane_sec[i,0,2]**2))**(1/2)
    
#Reshaping matrix for correct position to original plot scheme
gradients_table = np.reshape(gradients,(array_split_size_sub_tot,array_split_size_sub_tot))

np.savetxt("gradient_magnitudes_0060.csv",gradients_table,delimiter=' & ', fmt='%.4g', newline=' \\\\\n')

#New Steps for Feb 2nd Assignment

#Finding mass contribution for all subregions

cell_mass = ad.sum('cell_mass',axis='z') #Sums mass plane over z
arr_2 = cell_mass['cell_mass']          #defining variable
arr_2.shape = (data_length,data_length)
arr_2 = np.array(arr_2)                 #Turning into Numpy Array
#Resizing into superset of arrays
arr_2 = blockshaped(arr_2,array_split_size,array_split_size)

#Defining Empty array to store mass values.
mass_sub_region = np.zeros((array_split_size_sub,1,1))

for i in range(0,array_split_size_sub):
    mass_sub_region[i] = np.sum(arr_2[i])

#Now must compute Implied angular momentum for each subregion

#Defining Constants
beta = 1
#Subbox size in pc
length_pc_km = (1.25**2)*(3.086e13) #Gives L^2 as km*pc to cancel the gradient term
km_to_m_squared = 1000**2

#Computing implied angular momentum for each subbox

angular_momentum_implied = np.zeros((array_split_size_sub,1,1))

for i in range(0,array_split_size_sub):
    angular_momentum_implied[i] = beta*((length_pc_km)**2) * mass_sub_region[i] * gradients[i] * km_to_m_squared


#Total Actual Angular Momentum Magnitudes in the regions

angular_momentum_magnitude = ad.integrate('angular_momentum_magnitude',weight=None,axis='z')
angular_momentum_magnitude = angular_momentum_magnitude['angular_momentum_magnitude']
angular_momentum_magnitude.shape = (data_length,data_length)  # Since the underlying data are 256^3

angular_momentum_magnitude = np.array(angular_momentum_magnitude)  #Making Numpy array

angular_momentum_magnitude = blockshaped(angular_momentum_magnitude,array_split_size,array_split_size)

#Defining Empty Array to store values
angular_momentum_magnitude_sub_regions = np.zeros((array_split_size_sub,1,1))

for i in range(0,array_split_size_sub):
    angular_momentum_magnitude_sub_regions[i] = np.sum(angular_momentum_magnitude[i])
    
#Plotting Values Obtaining

#define x axis quadrant values
x = np.linspace(1,16,16)
y1 = np.zeros((array_split_size_sub,1))
y2 = np.zeros((array_split_size_sub,1))
for i in range(0,array_split_size_sub):
    y1[i] = angular_momentum_magnitude_sub_regions[i]
    y2[i] = angular_momentum_implied[i]
    
np.savetxt("Actual_J.csv",y1,delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
np.savetxt("Implied_J.csv",y2,delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
    
plt.plot(x,y1,'b.',label='Angular Momentum')
plt.plot(x,y2,'r.',label='Implied Angular Momentum')
plt.legend(loc=1)
plt.xlabel('Qaudrant Numbers')
plt.ylabel(r'Angular Momentum $(m^2 g \ s^{-2})$')
plt.title('Actual Angular Momentum vs Implied Angular Momentum')
plt.savefig('Ang_Mom_Comparison.pdf')