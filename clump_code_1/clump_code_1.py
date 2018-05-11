
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 19:53:52 2018

@author: sfielder
"""

#%%
#Importing packages here
import yt
import numpy as np
import astropy
from yt.analysis_modules.level_sets.api import *
import math as m
import sys
import matplotlib.pyplot as plt


#filename = sys.argv[1]

#Defining Global Variables Here - Move to File After
l=10	#Length of Original Data Set
cmin = 5e-21	#Minimum Density Threshold for Clump Finding
#cmax = 5e-20	#Maximum Density Threshold for Clump Finding
step = 100	#Step-size multiplier for Clump Finding



def octant_split(data_object):
    """
    Splits cubic data object into a sub-set of 8 data objects, i.e. octants.
    	
    Input: Data_object (usually representing whole simulation space)
    Output: List of 8 data-objects.
    Given Global Parameters: Length of data_object along the side, in pc.
    """
    dbox_array = [] #Creating Empty List to Store Data Objects
    #Creating Boundaries for Octants Here
    x1 = np.array((-l/2,0,-l/2,0,-l/2,0,-l/2,0))
    x2 = np.array((0,l/2,0,l/2,0,l/2,0,l/2))
    y1 = np.array((-l/2,-l/2,0,0,-l/2,-l/2,0,0))
    y2 = np.array((0,0,l/2,l/2,0,0,l/2,l/2))
    z1 = np.array((-l/2,-l/2,-l/2,-l/2,0,0,0,0))
    z2 = np.array((0,0,0,0,l/2,l/2,l/2,l/2))
    for i in range(0,8):
        dbox_array.append(ds.r[(x1[i],'pc'):(x2[i],'pc'), (y1[i],'pc'):(y2[i],'pc'), (z1[i],'pc'):(z2[i],'pc')])	
    return(dbox_array)

def master_clump_maker(data_object):
    '''
    Creates master clump for imputted data object.
    Input: data_object
    Output: master_clump object
    '''
    #Defining Field for Contouring
    field = ("gas", "density")
    master_clump = Clump(data_object, field)
    return(master_clump)

def clump_finder(master_clump,clump_sizing):
    '''
    Finds clumps that match certain criteria.
    Input: master_clump object, clump sizing parameter
    Output: Lowest Tree Values (leaves) of the clump data.
    Given Global Variables: c_min, c_max, step.
    '''
    #Setting Parameters for Clump Finding
    master_clump.add_validator("min_cells", clump_sizing)
    master_clump.add_info_item("center_of_mass") #Adds Center of Mass info for Clumps
    #Finding Clumps Here
    find_clumps(master_clump, cmin, cmax, step)
    leaf_clumps = get_lowest_clumps(master_clump)
    return(leaf_clumps)

def center_of_mass(lc):
    '''
    Computes the center of mass of a imputted data object
    Input: Lowest Leaf Clump Data Object
    Output: 3-list of x,y,z coordinates of center of mass
    '''
    com = np.array(lc.quantities.center_of_mass())
    return(com)

def bounding_box(lc):
    '''
    Writes the lower and upper bounds for each axis of a box, given a data
    clump object
    
    Input: Lowest Lead Clump Data Object
    Output: (3,2)-list of x,y,z min/max values.
    '''
    bounding_box = np.array(lc.quantities.extrema([('gas','x'),('gas','y'),('gas','z')]))
    return(bounding_box)
    
def clump_box(data_object,br):
    '''
    Takes original data set, and constructs a data object according to
    the bounding box given.
    Input:  Original Data in File (ds)
            br: bounding region of the clump
    Output: Data object with is constructed according to bounding box,
            represents a data object around a clump
    '''
    dbox_clump = data_object.r[(br[0,0],'cm'):(br[0,1],'cm'), (br[1,0],'cm'):(br[1,1],'cm'), (br[2,0],'cm'):(br[2,1],'cm')]
    return(dbox_clump)
    
def velocity_array(data_object,velocity,axis):
    '''
    Takes Data object and integrates velocity along chosen line of sight
    and returns an array that matches data structure of original data set
    with only the clump having values in the array, all others are nan
    
    Input:  Data object (clump or otherwise)
        Both next imputs have to be in strings in a specific way
        Velocity must be: 'velocity_x', 'velocity_y', 'velocity_z'
        Axis must be: 'x', 'y', 'z'
            Velocity: the component of velocity to be integrated
            Axis: Line of sight along which velocity component is integrated
    '''
    v_object = data_object.integrate(velocity, weight='density', axis=axis)
    v_reform = v_object.to_frb((l,'pc'),(master_dist_data,master_dist_data))
    v_arr = np.array(v_reform[velocity])
    v_arr = np.reshape(v_arr,(master_dist_data,master_dist_data))
    return(v_arr,v_object)
    
def velocity_array_reducer(velocity_int_array,velocity_int_data_object,axis):
    '''
    Takes velocity integrated array that is in the original dimension size
    with imbedded data values surrounded by nan values
    and reduces that array into just the imbedded data,
    keeping the generally shape of the data intact.
    Input: Integrated Velocity Array
    Output:
        velocity_integrated_array_reduced: just the significant data, no nan values
        v_pi: First perpendicular coordinate values to the integrated array
        v_pj: Second perpendicular coordinate values to the integrated array
        
        Values for coordinate system are always defined such that the relationship
        between the integrated array and the perpendicular coordinates
        follow a cyclical pattern
        e.g. Calling Z-LOS array, gives first coordinate x, and second coordinate y
    '''
    #Find where the data is non-nan valued
    v_positions = np.argwhere(~np.isnan(velocity_int_array))
    
    v_positions_ij = []    #Creating list for array slicing
    v_positions_ij.append(v_positions[0,0]) #First Row Value
    v_positions_ij.append(v_positions[-1,0]) #Last Row Value
    v_positions_ij.append(v_positions[0,1])    #First Column Value
    v_positions_ij.append(v_positions[-1,1])   #Last Column Value
    
    
    if axis == 'z':
        # Finds appropriate px coordinates and py coordinates
        vx_coordinates = np.array(velocity_int_data_object['px'])
        vx_coordinates = np.reshape(vx_coordinates,(master_dist_data,master_dist_data))
        vz_px = vx_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        vy_coordinates = np.array(velocity_int_data_object['py'])
        vy_coordinates = np.reshape(vy_coordinates,(master_dist_data,master_dist_data))
        vz_py = vy_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                        v_positions_ij[2]:v_positions_ij[3]+1]
        
        return(velocity_int_array_reduced,vz_px,vz_py)
    
    if axis == 'y':
        # Finds appropriate px coordinates and pz coordinates
        vx_coordinates = np.array(velocity_int_data_object['px'])
        vx_coordinates = np.reshape(vx_coordinates,(master_dist_data,master_dist_data))
        vy_px = vx_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        vz_coordinates = np.array(velocity_int_data_object['pz'])
        vz_coordinates = np.reshape(vz_coordinates,(master_dist_data,master_dist_data))
        vy_pz = vz_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                        v_positions_ij[2]:v_positions_ij[3]+1]
        return(velocity_int_array_reduced,vy_px,vy_pz)
        
    
    
    

#Loads data into File
ds = yt.load("~/Documents/Astro_AngMom/Astro-Clump/data.0060.3d.hdf5")
master_dist_data = int(ds.domain_dimensions[0])

#Creates a Data Object containing all the Simulation Data
ad = ds.all_data()

#Splits the Data into Octants
octant = octant_split(ad)

#Grabs each of the octants and runs them through the Clump Finding algorithm
clumps = [] #Defining Empty List for loop
for i in range(4,len(octant)):
    master_clump_main = master_clump_maker(octant[i])
    cmax = octant[i]["gas", "density"].max()
    print("Now Finding clumps for Octant:",i)
    lc = clump_finder(master_clump_main,30)
    for j in range(0,len(lc)):
        clumps.append(lc[j])

# =============================================================================
# Making empty arrays to store the data, variable to allow for different
# size of clumps list. Actual array (1,3) and (3,2) will not change as
# that is how they are written with no difference between clump objects.
# =============================================================================

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
#Creates list of data objects for all the clumps found through algorithm
for i in range(0,len(bregion)):
    data_object_clump.append(clump_box(ds,bregion[0]))


# =============================================================================
# CHECK IN POINT: What we have so far here in the script
#
# list known as data_object_clump which has each of the regions found by the 
# clumping algorithm.
# =============================================================================
    


#%%

#Testing: Calling the velocity_array definition
arr_z, vz = velocity_array(data_object_clump[0],'velocity_z','z')
arr_x, vx = velocity_array(data_object_clump[0],'velocity_x','x')
arr_y, vy = velocity_array(data_object_clump[0],'velocity_y','y')


arr_z_reduced, vz_px, vz_py = velocity_array_reducer(arr_z,vz,'z')

x_position_vx = np.array(vx['px'])
y_position_vy = np.array(vy['py'])
z_position_vz = np.array(vz['pz'])







