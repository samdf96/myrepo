#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 19:52:29 2018

@author: sfielder
"""

import numpy as np
np.set_printoptions(threshold=np.inf)

import yt
from yt.analysis_modules.level_sets.api import *
#from yt.utilities.physical_constants import kboltz, mh
#from yt.units import Kelvin
import matplotlib.pyplot as plt
import math as m
from plane_fit import plane_fit

################# INPUTS TO FILE GO HERE #############################
'''
l = int(input("Enter length of the data set in pc (leave blank to enter default of 10: " or 10))
print("Defined Edge Length of Data Set as: ",str(l)," pc", sep='')

step = input("Enter Step Size for clump finding (leave blank to enter default of 1000):" or 1000)

'''
######################################################################
l=10
step=1000

#Loading Data Set into Script
ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")

#All Data Object for Clumping Variable
ad = ds.all_data()

#Used for naming scheme of dbox_x areas
namespace = globals()

#Splitting Arrays for Octant Beginning and Ending Values for l length space
x1 = np.array((-l/2,0,-l/2,0,-l/2,0,-l/2,0))
x2 = np.array((0,l/2,0,l/2,0,l/2,0,l/2))
y1 = np.array((-l/2,-l/2,0,0,-l/2,-l/2,0,0))
y2 = np.array((0,0,l/2,l/2,0,0,l/2,l/2))
z1 = np.array((-l/2,-l/2,-l/2,-l/2,0,0,0,0))
z2 = np.array((0,0,0,0,l/2,l/2,l/2,l/2))
octant = 8

#Creating Loop for 8 Octants of data space to reduce computation overload

#Creating Empty Lists for data input later on
dbox_array = []
angular_momentum_list_x_projection = []
angular_momentum_list_y_projection = []
angular_momentum_list_z_projection = []

for i in range(0,octant):   #Master Loop for Octants
    #Creating Boxed Regions for Octants of Data Simulation Set Here
    dbox_array.append(ds.r[(x1[i],'pc'):(x2[i],'pc'), (y1[i],'pc'):(y2[i],'pc'), (z1[i],'pc'):(z2[i],'pc')])

    #Defining Field for Contouring
    field = ("gas", "density")

    #Creating Master Clump, and Clump Parameters
    master_clump = Clump(dbox_array[i], ("gas", "density")) #Makes the first big clump
    clump_sizing = 50  #Any Clump size smaller than this value get eliminated
    
    #Validator is the minimum clump size:
    master_clump.add_validator("min_cells", clump_sizing)
    
    master_clump.add_info_item("center_of_mass") #Adds Center of Mass info for Clumps
    
    #Setting Limits for Clump Finding
    c_min = 5e-21
    c_max = dbox_array[i]["gas", "density"].max()
    step = 1000
    find_clumps(master_clump, c_min, c_max, step) #Finds Clumps

#Projection plots for Clumping Overlay
    prj = yt.ProjectionPlot(ds,"z",field, data_source=dbox_array[i])    #Projcetion Plot
    prj.set_zlim('density', 1e-3, 3e-2) #Setting limits for all pictures to be the same
    prj.annotate_clumps(get_lowest_clumps(master_clump))    #Draws Contours on top of projection
    prj.save('Clump_0100_ClumpSizing_' + str(clump_sizing) + '_Octant_' + str(i) + '.png')  #Saves Figures

#Grabs Lowest Tree Values (Leaves)
    lc = get_lowest_clumps(master_clump)
    
# Sets up parameters needed for further calculations    
    center_of_mass = np.zeros((len(lc),3))
    box_region = np.zeros((len(lc),3,2))
    max_x_distance = np.zeros((len(lc),1))
    master_dist_array = np.array(ad.quantities.extrema([('gas','x'),('gas','y'),('gas','z')]))
    master_dist = master_dist_array[0,1] * 2
    master_dist_data = int(ds.domain_dimensions[0]) #Original Data Set underlying data structure
    
#This creates new center_of_mass and new box_regions for the octant
    for j in range(0,len(lc)):
        com = np.array(lc[j].quantities.center_of_mass())   #Writes current com into array
        center_of_mass[j] = com     #Writes to center_of_mass array
        #print(lc[j].quantities.angular_momentum_vector())
        bregions = np.array(lc[j].quantities.extrema([('gas','x'),('gas','y'),('gas','z')]))
        box_region[j] = bregions # Creates Boundaries for clump boxes
        

#Main Loop for All the clump regions, to determine parameters
    for k in range(0,len(lc)):
        dbox_clump = ds.r[(box_region[k,0,0],'cm'):(box_region[k,0,1],'cm'), (box_region[k,1,0],'cm'):(box_region[k,1,1],'cm'), (box_region[k,2,0],'cm'):(box_region[k,2,1],'cm')]

######### Section for finding length of each box in pc ###############

#Also finds scale factors needed for angular momentum computation
        # i.e. Length of box in pc of each axis.

        # For x-axis
        dist_x = box_region[k,0,1] - box_region[k,0,0]  #Finds x-axis distance
        dist_x_scale_1 = np.abs(com[0]/master_dist)
        scale_x = dist_x / master_dist  #Gets Scale factor for x
        data_length_x_pc = scale_x * l
        data_length_x_non_rounded = scale_x * master_dist_data  #Data_length for forcing mesh
        data_length_x = m.floor(data_length_x_non_rounded)
        
        # For y-axis
        dist_y = box_region[k,1,1] - box_region[k,1,0]  #Finds y-axis distance
        scale_y = dist_y / master_dist  #Gets Scale factor for y
        data_length_y_pc = scale_y * l
        data_length_y_non_rounded = scale_y * master_dist_data  #Data_length for forcing mesh
        data_length_y = m.floor(data_length_y_non_rounded)
        
        # For y-axis
        dist_z = box_region[k,2,1] - box_region[k,2,0]  #Finds z-axis distance
        scale_z = dist_z / master_dist  #Gets Scale factor for z
        data_length_z_pc = scale_z * l
        data_length_z_non_rounded = scale_z * master_dist_data  #Data_length for forcing mesh
        data_length_z = m.floor(data_length_z_non_rounded)
        
######################################################################  
        
        # This creates projections weighting by density for clump
        vx = dbox_clump.integrate('velocity_x', weight='density', axis='x')
        vy = dbox_clump.integrate('velocity_y', weight='density', axis='y')
        vz = dbox_clump.integrate('velocity_z', weight='density', axis='z')
        
        #Data of x_velocity into arrray
        arr_x = vx.to_frb((l,'pc'),(master_dist_data,master_dist_data))
        arr_x = np.array(arr_x['velocity_x'])
        arr_x = np.nan_to_num(arr_x) #Gets rid of nan entries
        
        #Data of y_velocity into arrray
        arr_y = vy.to_frb((l,'pc'),(master_dist_data,master_dist_data))
        arr_y = np.array(arr_y['velocity_y'])
        arr_y = np.nan_to_num(arr_y) #Gets rid of nan entries
        
        #Data of z_velocity into arrray
        arr_z = vz.to_frb((l,'pc'),(master_dist_data,master_dist_data))
        arr_z = np.array(arr_z['velocity_z'])
        arr_z = np.nan_to_num(arr_z) #Gets rid of nan entries

        #Plane Fitting Process Start
        #Creating i,j coordintes for blocks of data
        ivalues = np.arange(0,master_dist_data);
        jvalues = np.arange(0,master_dist_data);
    
        ii, jj = np.meshgrid(ivalues, jvalues)  #Creates i,j for planefitting
        
        #Creating Empty Arrays for Gradient Values
        first_plane_sec_x = np.zeros((1,3))
        first_plane_sec_y = np.zeros((1,3))
        first_plane_sec_z = np.zeros((1,3))
        
        #Runs the Plane fitting to get a gradient values
        first_plane_sec_x = plane_fit(np.ndarray.flatten(ii),np.ndarray.flatten(jj),np.ndarray.flatten(arr_x))
        first_plane_sec_y = plane_fit(np.ndarray.flatten(ii),np.ndarray.flatten(jj),np.ndarray.flatten(arr_y))
        first_plane_sec_z = plane_fit(np.ndarray.flatten(ii),np.ndarray.flatten(jj),np.ndarray.flatten(arr_z))

        #Need Conversion for values in km/s/pc instead of cm/s/pixel
        conv_length = 1/(10000) #cm to km
        conv_distance = master_dist_data/l #pixel to pc
        #Correct For units
        first_plane_sec_x = first_plane_sec_x * conv_length * conv_distance
        first_plane_sec_y = first_plane_sec_y * conv_length * conv_distance
        first_plane_sec_z = first_plane_sec_z * conv_length * conv_distance

        #Creating Empty Arrays for gradient value
        gradient_x = np.zeros((1,1))
        gradient_y = np.zeros((1,1))
        gradient_z = np.zeros((1,1))
        #Computing Gradient Term
        gradient_kms_x = ((first_plane_sec_x[1]**2) + (first_plane_sec_x[2]**2))**(1/2)
        gradient_kms_y = ((first_plane_sec_y[1]**2) + (first_plane_sec_y[2]**2))**(1/2)
        gradient_kms_z = ((first_plane_sec_z[1]**2) + (first_plane_sec_z[2]**2))**(1/2)
        
        #Correcting for units (km/pc to unity, i.e. per second)
        gradient_x = gradient_kms_x*(3.24078e-14)
        gradient_y = gradient_kms_y*(3.24078e-14)
        gradient_z = gradient_kms_z*(3.24078e-14)
        
        #Computing Mass for the region:
        mass_cell_total = np.array(dbox_clump.quantities.total_mass())[0]
           
        #Making empty array for Observational Estimator
        angular_momentum_implied_specific_x_projection = np.zeros((1,1))
        angular_momentum_implied_specific_y_projection = np.zeros((1,1))
        angular_momentum_implied_specific_z_projection = np.zeros((1,1))
        
        #Defining Constants
        beta = 1
        #Converts pc^2 to cm^2 in numerator to keep cgs units
        convert_factor = (1/(3.24078e-19))**2
        #Computing implied angular momentum for each subbox
        angular_momentum_implied_specific_x_projection = beta * gradient_x * data_length_y_pc * data_length_z_pc * convert_factor
        angular_momentum_implied_specific_y_projection = beta * gradient_y * data_length_x_pc * data_length_z_pc * convert_factor
        angular_momentum_implied_specific_z_projection = beta * gradient_z * data_length_x_pc * data_length_y_pc * convert_factor
        
        #Computing Actual Angular Momentum
        
        #Integration on line of sight
        angular_momentum_actual_x = dbox_clump.sum('angular_momentum_x')
        angular_momentum_actual_y = dbox_clump.sum('angular_momentum_y') 
        angular_momentum_actual_z = dbox_clump.sum('angular_momentum_z')
        angular_momentum_actual_xy = ((angular_momentum_actual_x**2)+(angular_momentum_actual_y**2))**(1/2)
        angular_momentum_actual_xz = ((angular_momentum_actual_x**2)+(angular_momentum_actual_z**2))**(1/2)
        angular_momentum_actual_yz = ((angular_momentum_actual_y**2)+(angular_momentum_actual_z**2))**(1/2)
        angular_momentum_actual_total = ((angular_momentum_actual_x**2)+(angular_momentum_actual_y**2)+(angular_momentum_actual_z**2))**(1/2)
        
        #Defining Empty Array
        angular_momentum_actual_total_specific = angular_momentum_actual_total / mass_cell_total
        angular_momentum_actual_xy_specific = angular_momentum_actual_xy / mass_cell_total
        angular_momentum_actual_xz_specific = angular_momentum_actual_xz / mass_cell_total
        angular_momentum_actual_yz_specific = angular_momentum_actual_yz / mass_cell_total
    
        angular_momentum_list_x_projection.append([gradient_kms_x,
                                      gradient_x,
                                      mass_cell_total,
                                      angular_momentum_implied_specific_x_projection,
                                      angular_momentum_actual_total_specific,
                                      angular_momentum_actual_yz_specific
                                      ])
        angular_momentum_list_y_projection.append([gradient_kms_y,
                                      gradient_y,
                                      mass_cell_total,
                                      angular_momentum_implied_specific_y_projection,
                                      angular_momentum_actual_total_specific,
                                      angular_momentum_actual_xz_specific
                                      ])
        angular_momentum_list_z_projection.append([gradient_kms_z,
                                      gradient_z,
                                      mass_cell_total,
                                      angular_momentum_implied_specific_z_projection,
                                      angular_momentum_actual_total_specific,
                                      angular_momentum_actual_xy_specific
                                      ])

# Creating Array of Values for this list, for plotting below.
angular_momentum_list_x_projection = np.array(angular_momentum_list_x_projection)
angular_momentum_list_y_projection = np.array(angular_momentum_list_y_projection)
angular_momentum_list_z_projection = np.array(angular_momentum_list_z_projection)

#Saving Data to files.
np.savetxt("Angular_Momentum_0100_x_projection.csv",
           angular_momentum_list_x_projection,
           delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
np.savetxt("Angular_Momentum_0100_y_projection.csv",
           angular_momentum_list_y_projection,
           delimiter=' & ', fmt='%.4g', newline=' \\\\\n')
np.savetxt("Angular_Momentum_0100_z_projection.csv",
           angular_momentum_list_z_projection,
           delimiter=' & ', fmt='%.4g', newline=' \\\\\n')

#Plotting Comparison Graphs

# =============================================================================
# Unit Conversion for cgs [cm^2 / s] to Astro units [pc^2 / Myr]
cgs_to_astro = ((3.24078e-19)**2) * (3.154e13) 
# =============================================================================
    
################### For X-Projection #################################

qx = len(angular_momentum_list_x_projection)
#Redfining Empty Arrays
y1 = np.zeros((qx,1))
y2 = np.zeros((qx,1))
y3 = np.zeros((qx,1))
for i in range(0,qx):
    y1[i] = angular_momentum_list_x_projection[i,4]
    y2[i] = angular_momentum_list_x_projection[i,5]
    y3[i] = angular_momentum_list_x_projection[i,3]
    
#Switching to Astro Units
y1 = np.multiply(y1,cgs_to_astro)
y2 = np.multiply(y2,cgs_to_astro)
y3 = np.multiply(y3,cgs_to_astro)

y1_min = np.around(np.min(y1),2)
y2_min = np.around(np.min(y2),2)
y1_max = np.around(np.max(y1),1)
y2_max = np.around(np.max(y2),1)

x_min = np.min([y1_min,y2_min])
x_max = np.max([y1_max,y2_max])

#Plotting - Loglog (residual based with the unity line)
plt.figure(1)
plt.loglog(y1,y3,'r.',label='Implied vs Total')
plt.loglog(y2,y3,'b*', label='Implied vs Partial')
x_plotting = np.linspace(x_min,x_max,1000)
plt.loglog(x_plotting, x_plotting, 'k-', alpha=0.75, zorder=0,label='Line of Unity')
plt.legend(bbox_to_anchor=(1, 0.5))
#plt.grid(True, which='both')
#This is not working correctly.
#plt.axes().set_aspect('equal', 'datalim')
plt.xlabel(r'Actual Specific Angular Momentum $(pc^2 / Myr)$')
plt.ylabel(r'Implied Specific Angular Momentum $(pc^2 / Myr)$')
plt.title('X-Axis Line of Sight', y=1.08)
plt.savefig("Ang_Mom_Specific_Comparison_0100_X-LOS.pdf", bbox_inches='tight')
plt.savefig("Ang_Mom_Specific_Comparison_0100_X-LOS.png", bbox_inches='tight')
plt.show()


################### For Y-Projection #################################

qy = len(angular_momentum_list_y_projection)
#Redfining Empty Arrays
y1 = np.zeros((qy,1))
y2 = np.zeros((qy,1))
y3 = np.zeros((qy,1))
for i in range(0,qy):
    y1[i] = angular_momentum_list_y_projection[i,4]
    y2[i] = angular_momentum_list_y_projection[i,5]
    y3[i] = angular_momentum_list_y_projection[i,3]

#Switching to Astro Units
y1 = np.multiply(y1,cgs_to_astro)
y2 = np.multiply(y2,cgs_to_astro)
y3 = np.multiply(y3,cgs_to_astro)

y1_min = np.around(np.min(y1),2)
y2_min = np.around(np.min(y2),2)
y1_max = np.around(np.max(y1),1)
y2_max = np.around(np.max(y2),1)

x_min = np.min([y1_min,y2_min])
x_max = np.max([y1_max,y2_max])
                              
#Plotting - Loglog (residual based with the unity line)
plt.figure(2)
plt.loglog(y1,y3,'r.',label='Implied vs Total')
plt.loglog(y2,y3,'b*', label='Implied vs Partial')
x_plotting = np.linspace(x_min,x_max,1000)
plt.loglog(x_plotting, x_plotting, 'k-', alpha=0.75, zorder=0,label='Line of Unity')
plt.legend(bbox_to_anchor=(1, 0.5))
#plt.grid(True, which='both')
#This is not working correctly.
#plt.axes().set_aspect('equal', 'datalim')
plt.xlabel(r'Actual Specific Angular Momentum $(pc^2 / Myr)$')
plt.ylabel(r'Implied Specific Angular Momentum $(pc^2 / Myr)$')
plt.title('Y-Axis Line of Sight', y=1.08)
plt.savefig("Ang_Mom_Specific_Comparison_0100_Y-LOS.pdf", bbox_inches='tight')
plt.savefig("Ang_Mom_Specific_Comparison_0100_Y-LOS.png", bbox_inches='tight')
plt.show()
    
################### For Z-Projection #################################

qz = len(angular_momentum_list_z_projection)
#Redfining Empty Arrays
y1 = np.zeros((qz,1))
y2 = np.zeros((qz,1))
y3 = np.zeros((qz,1))
for i in range(0,qz):    
    y1[i] = angular_momentum_list_z_projection[i,4]
    y2[i] = angular_momentum_list_z_projection[i,5]
    y3[i] = angular_momentum_list_z_projection[i,3]

#Switching to Astro Units
y1 = np.multiply(y1,cgs_to_astro)
y2 = np.multiply(y2,cgs_to_astro)
y3 = np.multiply(y3,cgs_to_astro)

y1_min = np.around(np.min(y1),2)
y2_min = np.around(np.min(y2),2)
y1_max = np.around(np.max(y1),1)
y2_max = np.around(np.max(y2),1)

x_min = np.min([y1_min,y2_min])
x_max = np.max([y1_max,y2_max])
                              
#Plotting - Loglog (residual based with the unity line)
plt.figure(3)
plt.loglog(y1,y3,'r.',label='Implied vs Total')
plt.loglog(y2,y3,'b*', label='Implied vs Partial')
x_plotting = np.linspace(x_min,x_max,1000)
plt.loglog(x_plotting, x_plotting, 'k-', alpha=0.75, zorder=0,label='Line of Unity')
plt.legend(bbox_to_anchor=(1, 0.5))
#plt.grid(True, which='both')
#This is not working correctly.
#plt.axes().set_aspect('equal', 'datalim')
plt.xlabel(r'Actual Specific Angular Momentum $(pc^2 / Myr)$')
plt.ylabel(r'Implied Specific Angular Momentum $(pc^2 / Myr)$')
plt.title('Z-Axis Line of Sight', y=1.08)
plt.savefig("Ang_Mom_Specific_Comparison_0100_Z-LOS.pdf", bbox_inches='tight')
plt.savefig("Ang_Mom_Specific_Comparison_0100_Z-LOS.png", bbox_inches='tight')
plt.show() 
       