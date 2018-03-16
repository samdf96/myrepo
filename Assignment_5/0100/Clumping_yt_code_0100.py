#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 19:52:29 2018

@author: sfielder
"""

import numpy as np

import yt
from yt.analysis_modules.level_sets.api import *
#from yt.utilities.physical_constants import kboltz, mh
#from yt.units import Kelvin
import matplotlib.pyplot as plt
import math as m
from plane_fit import plane_fit

#Loading Data Set into Script
ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")

#All Data Object for Clumping Variable
ad = ds.all_data()
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''
#Comment In to add in a Thermal Energy Density field for clumping
'''
def _therm(field, data):
   return(1.5 * (kboltz) * 10 * Kelvin / (2.32 * mh))
yt.add_field(("gas","thermal_energy"), function=_therm, units="erg/g", force_override=True)
'''


#Used for naming scheme of dbox_x areas
namespace = globals()
#Splitting Arrays for Octant Beginning and Ending Values for L length space
l = 10 #Data Length of One edge of original ds data set

x1 = np.array((-l/2,0,-l/2,0,-l/2,0,-l/2,0))
x2 = np.array((0,l/2,0,l/2,0,l/2,0,l/2))
y1 = np.array((-l/2,-l/2,0,0,-l/2,-l/2,0,0))
y2 = np.array((0,0,l/2,l/2,0,0,l/2,l/2))
z1 = np.array((-l/2,-l/2,-l/2,-l/2,0,0,0,0))
z2 = np.array((0,0,0,0,l/2,l/2,l/2,l/2))
octant = 8

#Creating 4 Loop for 8 Octants of data space to reduce computation overload
dbox_array = [] #Creating Empty List

for i in range(0,octant):
    dbox_array.append(ds.r[(x1[i],'pc'):(x2[i],'pc'), (y1[i],'pc'):(y2[i],'pc'), (z1[i],'pc'):(z2[i],'pc')])

    # the field to be used for contouring
    field = ("gas", "density")

    master_clump = Clump(dbox_array[i], ("gas", "density")) #Makes the first big clump
    clump_sizing = 50  #Any Clump size smaller than this value get eliminated
    master_clump.add_validator("min_cells", clump_sizing)
    '''
    master_clump.add_validator("gravitationally_bound",
                           use_particles=False,
                           use_thermal_energy=False)
    '''
    
    master_clump.add_info_item("center_of_mass") #Adds Center of Mass info for Clumps
    #Setting Limits for Clump Finding 
    c_min = 5e-21
    c_max = dbox_array[i]["gas", "density"].max()
    step = 1000
    find_clumps(master_clump, 5e-21, c_max, step) #Finds Clumps

    prj = yt.ProjectionPlot(ds,"z",field, data_source=dbox_array[i])    #Projcetion Plot
    prj.set_zlim('density', 1e-3, 3e-2) #Setting limits for all pictures to be the same
    prj.annotate_clumps(get_lowest_clumps(master_clump))    #Draws Contours on top of projection
    prj.save('Clump_0100_ClumpSizing_' + str(clump_sizing) + '_Octant_' + str(i) + '.png')  #Saves Fig

    #Grabs Lowest Tree Values (Leaves)
    lc = get_lowest_clumps(master_clump)
    
    center_of_mass = np.zeros((len(lc),3))
    box_region = np.zeros((len(lc),3,2))
    max_x_distance = np.zeros((len(lc),1))
    master_dist_array = np.array(ad.quantities.extrema([('gas','x'),('gas','y'),('gas','z')]))
    master_dist = master_dist_array[0,1] * 2
    master_dist_data = 256    #Original Data Set underlying data structure
    master_dist_pc = 10         #Original Data set length of edge in pc
    for j in range(0,len(lc)):
        com = np.array(lc[j].quantities.center_of_mass())   #Writes current com into array
        center_of_mass[j] = com     #Writes to center_of_mass array
        #print(lc[j].quantities.angular_momentum_vector())
        bregions = np.array(lc[j].quantities.extrema([('gas','x'),('gas','y'),('gas','z')]))
        box_region[j] = bregions # Creates Boundaries for clump boxes
        
    for k in range(0,len(lc)):
        dbox_clump = ds.r[(box_region[k,0,0],'cm'):(box_region[k,0,1],'cm'), (box_region[k,1,0],'cm'):(box_region[k,1,1],'cm'), (box_region[k,2,0],'cm'):(box_region[k,2,1],'cm')]
        #test_proj = yt.ProjectionPlot(ds,"z",field, data_source=dbox_clump)
        #test_proj.save('test_proj_' + str(k) + '.png')
        
        # This creates a projection weighting by density for clump
        vz = dbox_clump.integrate('velocity_z',weight='density',axis='z')
        
        test_image = plt.scatter(vz['px'],vz['py'],c=vz['velocity_z'])
        plt.colorbar()
        plt.show()
        plt.clf()
        
        # For x-axis
        dist_x = box_region[k,0,1] - box_region[k,0,0]  #Finds x-axis distance
        scale_x = dist_x / master_dist  #Gets Scale factor for x
        data_length_x_pc = scale_x * master_dist_pc
        data_length_x_non_rounded = scale_x * master_dist_data  #Data_length for forcing mesh
        data_length_x = m.floor(data_length_x_non_rounded)
        
        # For y-axis
        dist_y = box_region[k,1,1] - box_region[k,1,0]  #Finds y-axis distance
        scale_y = dist_y / master_dist  #Gets Scale factor for y
        data_length_y_pc = scale_y * master_dist_pc
        data_length_y_non_rounded = scale_y * master_dist_data  #Data_length for forcing mesh
        data_length_y = m.floor(data_length_y_non_rounded)
        
        #Data of z_velocity into np.arrray
        arr = vz.to_frb((data_length_x_pc,'pc'),(data_length_x,data_length_y), height=(data_length_y_pc,'pc'))
        arr = np.array(arr['velocity_z'])
        
        #Plane Fitting Process Start
        #Creating x,y coordintes for blocks of data
        xvalues = np.arange(0,data_length_x);
        yvalues = np.arange(0,data_length_y);
    
        xx, yy = np.meshgrid(xvalues, yvalues)  #Creates x,y for planefitting
        
        first_plane_sec = np.zeros((1,3))       #Creates empty array
        
        #Runs the Plane fitting to get a gradient value
        first_plane_sec = plane_fit(np.ndarray.flatten(xx),np.ndarray.flatten(yy),np.ndarray.flatten(arr))

        #Need Conversion for values in km/s/pc instead of cm/s/pixel
        conv_length = 1/(10000) #cm to km
        conv_distance = 1/0.039 #pixel to pc
        #Correct For units
        first_plane_sec = first_plane_sec * conv_length * conv_distance
        
        #Creating Empty Array for Graient term
        gradient = np.zeros((1,1))
        #Computing Gradient Term
        gradient = ((first_plane_sec[1]**2) + (first_plane_sec[2]**2))**(1/2)
        
        #Computing Mass for the region:
        mass_cell_total = np.array(dbox_clump.quantities.total_mass())[0]
        
        
        #Making empty array for Observational Estimator
        angular_momentum_implied_specific = np.zeros((1,1))
        
        #Defining Constants
        beta = 1
    
        #Computing implied angular momentum for each subbox
        angular_momentum_implied_specific = beta * gradient * data_length_x_pc * data_length_y_pc
        
        #Computing Actual Angular Momentum
        
        #Integration on line of sight
        angular_momentum_actual = dbox_clump.integrate('angular_momentum_magnitude',axis='z',weight=None)
        angular_momentum_actual = angular_momentum_actual.to_frb((data_length_x_pc,'pc'),[data_length_x,data_length_y],height=(data_length_y_pc,'pc'))
        angular_momentum_actual = np.array(angular_momentum_actual['angular_momentum_magnitude'])
        
        #Defining Empty Array
        angular_momentum_actual_specific_sum = np.zeros((1,1))
        unit_conv = (1/(3.086e18)) * (1e-5)
        angular_momentum_actual_specific_sum = np.sum(angular_momentum_actual) * unit_conv / mass_cell_total
        