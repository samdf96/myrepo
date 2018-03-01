#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:29:55 2018

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
    Data Width - Length of one side of simulation space in pc (original data set)
    Beta - Constant for Computing Implied Angular Momentum
    Subregion_edge_length (Implied Section) - Subregion edge length in pc.
'''


# Import Packages here.
import yt as yt
import matplotlib.pyplot as plt
import numpy as np

# Call different data sets here, possibly create loop.
ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")

#Loading all data into one object.
ad = ds.all_data()

#Printing Derived Variables to be able to use (comment in if needed).
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''

# Basic Projection Plots
'''
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
plt.show() #Rid of comment to show plot.

#Defining Value for underlying data set length
data_length = 256
data_length_pc = 10


arr = vz.to_frb((data_length_pc,'pc'),[data_length,data_length]) #Forces 256x256 mesh
arr = np.array(arr['velocity_z'])   #Imports velocity data into np array

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
 


for n in range(0,5):

    
    #Need to make it so that it runs over a certain array of sizes:   
    #Length of the edge of a box
    length_variable = np.array([10,5,2.5,1.25,0.625])
    #Total Number of Sub regions
    array_split_size_sub = np.array([1,4,16,64,256], dtype=np.int64)
    array_split_size = np.array([256,128,64,32,16], dtype=np.int64)
    gradients_split_size = np.array([1,2,4,8,16], dtype=np.int64)

    array_split_size = array_split_size[n]
    array_split_size_sub = array_split_size_sub[n]
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
    
    gradients = np.zeros((array_split_size_sub,1))
    
    for i in range(0,array_split_size_sub):
        gradients[i,0] = ((first_plane_sec[i,0,1]**2) + (first_plane_sec[i,0,2]**2))**(1/2)
        
    #Reshaping matrix for correct position to original plot scheme
    gradients_table = np.reshape(gradients,(gradients_split_size[n],gradients_split_size[n]))
    
    np.savetxt("gradient_magnitudes_0060_{n}.csv".format(n=n),gradients_table,delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
    
    
    # Computing Mass for each region
    
    # This next line integrates over a line of sight to create a Line-of-sight projection weighted by density
    mass = ad.integrate('cell_mass',axis='z',weight=None)
    # Now we have to make this into a gridded set of pixel data.
    data_width = data_length_pc
    mass_cell = mass.to_frb((data_width,'pc'),[data_length,data_length])
    mass_cell = np.array(mass_cell['cell_mass'])
    
    mass_cell = blockshaped(mass_cell,array_split_size,array_split_size)
    mass_cell_total = np.zeros((array_split_size_sub,1,1))
    for i in range(0,array_split_size_sub):
        mass_cell_total[i] = np.sum(mass_cell[i])
        
    #Defining Constants
    beta = 1
    #Subbox size in pc
    subregion_edge_length = length_variable[n] #Given in pc
    #unit_conv_2 = (1/(1.988e33)) old unit conversion for full J
    
    #Computing implied angular momentum for each subbox
    
    angular_momentum_implied_specific = np.zeros((array_split_size_sub,1,1))
    
    for i in range(0,array_split_size_sub):
        angular_momentum_implied_specific[i] = beta * gradients[i] * subregion_edge_length**2
    
        
    # Computing Actual Angular Momentum Weighted by Density?
        
    angular_momentum_actual = ad.integrate('angular_momentum_magnitude',axis='z',weight=None)
    # Now we have to make this into a gridded set of pixel data.
    angular_momentum_actual = angular_momentum_actual.to_frb((data_width,'pc'),[data_length,data_length])
    angular_momenutm_actual = np.array(angular_momentum_actual['angular_momentum_magnitude'])
    angular_momentum_actual = blockshaped(angular_momenutm_actual,array_split_size,array_split_size)
    
    angular_momentum_actual_specific_sum = np.zeros((array_split_size_sub,1,1))
    unit_conv = (1/(3.086e18)) * (1e-5)
    for i in range(0,array_split_size_sub):
        angular_momentum_actual_specific_sum[i] = np.sum(angular_momentum_actual[i]) * unit_conv / mass_cell_total[i]
    
    #Plotting Values Obtaining
    
    #define x axis quadrant values
    x_plotting = np.linspace(2,11,1000)
    y1 = np.zeros((array_split_size_sub,1))
    y2 = np.zeros((array_split_size_sub,1))
    for i in range(0,array_split_size_sub):
        y1[i] = angular_momentum_actual_specific_sum[i]
        y2[i] = angular_momentum_implied_specific[i]
    
    np.savetxt("Actual_J_specific_{n}.csv".format(n=n),y1,delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
    np.savetxt("Implied_J_specific_{n}.csv".format(n=n),y2,delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
        
    #Plotting - Loglog (residual based with the unity line)
    plt.loglog(y1,y2,'r.',label='Specific Angular Momentum')
    plt.loglog(x_plotting, x_plotting, 'k-', alpha=0.75, zorder=0,label='Line of Unity')
    plt.legend(bbox_to_anchor=(1, 0.5))
    #plt.grid(True, which='both')
    plt.xlabel(r'Actual Specific Angular Momentum $(pc^2 \ Myr^{-1})$')
    plt.ylabel(r'Implied Specific Angular Momentum $(pc^2 \ Myr^{-1})$')
    plt.title('Actual Specific Angular Momentum vs Implied Specific Angular Momentum', y=1.08)
    plt.savefig("Ang_Mom_Specific_Comparison_{n}.pdf".format(n=n), bbox_inches='tight')
    plt.gcf().clear()