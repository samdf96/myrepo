
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 19:53:52 2018

@author: sfielder
"""


#Importing packages here
import yt
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# To be used with YAML file which we will create later on
#import sys
#filename = sys.argv[1]


#Defining Global Variables Here - Move to File After
l=10	#Length of Original Data Set
cmin = 5e-21	#Minimum Density Threshold for Clump Finding
#cmax = 5e-20	#Maximum Density Threshold for Clump Finding
step = 100	 #Step-size multiplier for Clump Finding
beta = 1        # For Implied Angular Momentum Calculation.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Section for importing all the definitions from defitions.py file
from definitions import octant_split
from definitions import master_clump_maker
from definitions import clump_finder
from definitions import center_of_mass
from definitions import bounding_box
from definitions import clump_box
from definitions import velocity_array
from definitions import velocity_array_reducer
from definitions import array_flattener
#Add in if needed
#from definitions import myplane
#from definitions import plane_fit_visualization
from definitions import plane_fit
from definitions import gradient
from definitions import angular_momentum_implied
from definitions import angular_momentum_actual
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#Loads data into File
ds = yt.load("~/Documents/Astro_AngMom/Astro-Clump/data.0060.3d.hdf5")
master_dist_data = int(ds.domain_dimensions[0])

#Creates a Data Object containing all the Simulation Data
ad = ds.all_data()

#Splits the Data into Octants
octant = octant_split(ad,ds,l)

#Grabs each of the octants and runs them through the Clump Finding algorithm
clumps = [] #Defining Empty List for loop
for i in range(4,len(octant)):
    master_clump_main = master_clump_maker(octant[i])
    cmax = octant[i]["gas", "density"].max()
    print("Now Finding clumps for Octant:",i)
    lc = clump_finder(master_clump_main,30,cmin,cmax,step)
    for j in range(0,len(lc)):
        clumps.append(lc[j])

# =============================================================================
# Making empty arrays to store the data, variable to allow for different
# size of clumps list. Actual array (1,3) and (3,2) will not change as
# that is how they are written with no difference between clump objects.
# =============================================================================

#Creating Arrays for later export (see bottom of script for this step)
com = np.zeros((len(clumps),1,3))
bregion = np.zeros((len(clumps),3,2))

# Creates lists for center_of_mass of clump list -> stores in com
# Creates list for bounding boxes of clump list -> stores in bregion
for i in range(0,len(clumps)):
    com[i] = center_of_mass(clumps[i])
    bregion[i] = bounding_box(clumps[i])



# =============================================================================
# CHECK IN POINT: What we have so far here in the script
#
#     com: arrays that contain x,y,z center of mass values for clumps
#     bregion: x,y,z min/max value to build the boxes for clumps
# =============================================================================


data_object_clump = [] #Setting Empty list for loop
mass_clump_g = []       #Setting Empty list for loop
dist_span_cm = np.zeros((len(bregion),3,1))
for i in range(0,len(bregion)):
    data_object_clump.append(clump_box(ds,bregion[i]))
    
    mass_clump_g.append(np.array(data_object_clump[i].quantities.total_mass())[0])
    
    dist_x = bregion[i,0,1] - bregion[i,0,0]
    dist_y = bregion[i,1,1] - bregion[i,1,0]
    dist_z = bregion[i,2,1] - bregion[i,2,0]
    
    dist_tot = np.array([[dist_x],[dist_y],[dist_z]])
    
    # Already in array format for export later
    dist_span_cm[i] = dist_tot

#Turn into Array for later export (see bottom of script for this step)
mass_clump_g = np.array(mass_clump_g)


# =============================================================================
# CHECK IN POINT: What we have so far here in the script
#
# list known as data_object_clump which has each of the regions found by the 
# clumping algorithm.
# =============================================================================

#Defining Empty Lists to store data into, will convert into arrays later on
grad_x = []
grad_y = []
grad_z = []
am_implied_x = []
am_implied_y = []
am_implied_z = []
am_actual_total = []
am_actual_partial_xy = []
am_actual_partial_xz = []
am_actual_partial_yz = []


for i in range(0,len(bregion)):
    clump = data_object_clump[i]

    # =============================================================================
    # Computing Integrated Velocity Arrays
    arr_z, vz = velocity_array(clump,'velocity_z','z',master_dist_data,l)
    arr_x, vx = velocity_array(clump,'velocity_x','x',master_dist_data,l)
    arr_y, vy = velocity_array(clump,'velocity_y','y',master_dist_data,l)
    # =============================================================================
    
    # =============================================================================
    # Computing Reduced Velocity Arrays
    arr_x_red, vx_py, vx_pz = velocity_array_reducer(arr_x,vx,'x',master_dist_data)
    arr_y_red, vy_px, vy_pz = velocity_array_reducer(arr_y,vy,'y',master_dist_data)
    arr_z_red, vz_px, vz_py = velocity_array_reducer(arr_z,vz,'z',master_dist_data)
    # =============================================================================
    
    # =============================================================================
    # Flattening Arrays for Plane Fitting Process
    arr_x_red_flat = array_flattener(arr_x_red)
    arr_y_red_flat = array_flattener(arr_y_red)
    arr_z_red_flat = array_flattener(arr_z_red)
    
    vx_py_flat = array_flattener(vx_py)
    vx_pz_flat = array_flattener(vx_pz)
    vy_px_flat = array_flattener(vy_px)
    vy_pz_flat = array_flattener(vy_pz)
    vz_px_flat = array_flattener(vz_px)
    vz_py_flat = array_flattener(vz_py)
    # =============================================================================
    
    #Plane Fitting Process
    result_x = plane_fit(vx_py_flat,vx_pz_flat,arr_x_red_flat)
    result_y = plane_fit(vy_px_flat,vy_pz_flat,arr_y_red_flat)
    result_z = plane_fit(vz_px_flat,vz_py_flat,arr_z_red_flat)
    
    #Computing Gradients from Results
    gradient_x = gradient(result_x)
    gradient_y = gradient(result_y)
    gradient_z = gradient(result_z)

    #Computing Specific Angular Momentum
    # BETA VALUE INPUT TO DEFINITION WILL DEFAULT TO ONE
    ang_mom_implied_x = angular_momentum_implied(gradient_x,
                                                  dist_span_cm[i,1,0],
                                                  dist_span_cm[i,2,0])
    ang_mom_implied_y = angular_momentum_implied(gradient_y,
                                                  dist_span_cm[i,0,0],
                                                  dist_span_cm[i,2,0])
    ang_mom_implied_z = angular_momentum_implied(gradient_z,
                                                  dist_span_cm[i,0,0],
                                                  dist_span_cm[i,1,0])
    
    ang_mom_actual_total, ang_mom_actual_xy, ang_mom_actual_xz, ang_mom_actual_yz = angular_momentum_actual(clump,mass_clump_g[i])
    
    
    # =========================================================================
    #     STORAGE OF VALUES FOR LATER USE
    grad_x.append(gradient_x)
    grad_y.append(gradient_y)
    grad_z.append(gradient_z)
    am_implied_x.append(ang_mom_implied_x)
    am_implied_y.append(ang_mom_implied_y)
    am_implied_z.append(ang_mom_implied_z)
    am_actual_total.append(ang_mom_actual_total)
    am_actual_partial_xy.append(ang_mom_actual_xy)
    am_actual_partial_xz.append(ang_mom_actual_xz)
    am_actual_partial_yz.append(ang_mom_actual_yz)
    # =========================================================================


# Turning all the data into numpy arrays to export
clump_number = np.arange(1,len(bregion)+1)
am_actual_total = np.array(am_actual_total)
am_actual_partial_xy = np.array(am_actual_partial_xy)
am_actual_partial_xz = np.array(am_actual_partial_xz)
am_actual_partial_yz = np.array(am_actual_partial_yz)
am_implied_x = np.array(am_implied_x)
am_implied_y = np.array(am_implied_y)
am_implied_z = np.array(am_implied_z)
grad_x = np.array(grad_x)
grad_y = np.array(grad_y)
grad_z = np.array(grad_z)

#Exporting into archive with the name in the string for first input
np.savez('clump_data',
         clump_number=clump_number,
         center_of_mass=com,
         mass=mass_clump_g,
         spanning_distance=dist_span_cm,
         gradient_x_los=grad_x,
         gradient_y_los=grad_y,
         gradient_z_los=grad_z,
         am_imp_x_los=am_implied_x,
         am_imp_y_los=am_implied_y,
         am_imp_z_los=am_implied_z,
         am_actual_total=am_actual_total,
         am_actual_xy=am_actual_partial_xy,
         am_actual_xz=am_actual_partial_xz,
         am_actual_yz=am_actual_partial_yz)



