
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
from scipy.optimize import least_squares as lsq

#For Visualization of Data (plane_fit_visualization definition)
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick

#filename = sys.argv[1]

#Defining Global Variables Here - Move to File After
l=10	#Length of Original Data Set
cmin = 5e-21	#Minimum Density Threshold for Clump Finding
#cmax = 5e-20	#Maximum Density Threshold for Clump Finding
step = 100	#Step-size multiplier for Clump Finding
beta = 1        # For Implied Angular Momentum Calculation.



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
        
    Relationship between line-of-sight and 'px' and 'py' values were found by
    comparting 'px' and 'py' to 'x','y','z' values for an arbitrary data
    object. The relationship is hard coded in this definition, as the data
    extracted from simulations will not be different in this sense.
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
        vx_coordinates = np.array(velocity_int_data_object['py'])
        vx_coordinates = np.reshape(vx_coordinates,(master_dist_data,master_dist_data))
        vz_px = vx_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        vy_coordinates = np.array(velocity_int_data_object['px'])
        vy_coordinates = np.reshape(vy_coordinates,(master_dist_data,master_dist_data))
        vz_py = vy_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                        v_positions_ij[2]:v_positions_ij[3]+1]
        
        return(velocity_int_array_reduced,vz_px,vz_py)
    
    if axis == 'y':
        # Finds appropriate px coordinates and py coordinates
        vx_coordinates = np.array(velocity_int_data_object['px'])
        vx_coordinates = np.reshape(vx_coordinates,(master_dist_data,master_dist_data))
        vy_px = vx_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        vz_coordinates = np.array(velocity_int_data_object['py'])
        vz_coordinates = np.reshape(vz_coordinates,(master_dist_data,master_dist_data))
        vy_pz = vz_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                        v_positions_ij[2]:v_positions_ij[3]+1]
        return(velocity_int_array_reduced,vy_px,vy_pz)
        
    if axis == 'x':
        # Finds appropriate px coordinates and py coordinates
        vy_coordinates = np.array(velocity_int_data_object['py'])
        vy_coordinates = np.reshape(vy_coordinates,(master_dist_data,master_dist_data))
        vx_py = vy_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        vz_coordinates = np.array(velocity_int_data_object['px'])
        vz_coordinates = np.reshape(vz_coordinates,(master_dist_data,master_dist_data))
        vx_pz = vz_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,v_positions_ij[2]:v_positions_ij[3]+1]
        velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                        v_positions_ij[2]:v_positions_ij[3]+1]
        return(velocity_int_array_reduced,vx_py,vx_pz)

def array_flattener(x):
    '''
    Flattens arrays using .ndarray.flatten command
    '''
    arr_flat = np.ndarray.flatten(x)
    return(arr_flat)
        
def myplane(p, x, y, z):
    '''
    Imported from Erik Rosolowsky
    '''
    return(p[0] + p[1] * x + p[2] * y - z)

def plane_fit_visualization(x,y,z,k):
    '''
    Takes flattened data, plane fits in, and visualizes it in 3d.
    NOTE: This method of fitting the data is different than the plane_fit
    definition used for computation of gradients. The resulting coefficients
    are different from those in plane_fit, thus this should only be used as 
    a rough indication of how the plane fitting process works, only use for 
    visualization of the data if it is needed at all.
    
    Input:
            x,y - mesh coordinates (flattened)
            z - coordinate data points to fit plane to
    Outputs:
            PDF and PNG of matplotlib figure
    '''
    data = np.c_[x,y,z]
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
    
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

    # evaluate it on grid
    Z = C[0]*X + C[1]*Y + C[2]
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=1)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.zaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax.set_zlabel('Z')
    ax.axis('equal')
    ax.axis('tight')
    plt.savefig("Plane_Fitting_Clump_{k}.pdf".format(k=k), bbox_inches='tight')
    plt.savefig("Plane_Fitting_Clump_{k}.png".format(k=k), bbox_inches='tight')        
    plt.show()
    return()


def plane_fit(x, y, z, robust=False):
    """
    Imported from Erik Rosolowsky
    Fits a plane to data without any given uncertainties or weighting.
    Arguments:
    ----------
    x, y: float
       x and y coordinates of the data
    z: float
       z-coordinates to which the values are fit
    Returns:
    --------
    coefficients:
       3-element vector with components [ z0 (constant offset) , grad_x, grad_y]
    
    """

    x0, y0 = np.median(x), np.median(y)
    dataz = np.c_[np.ones(x.size), 
                  x-x0, 
                  y-y0]

    lsqcoeffs, _, _, _ = np.linalg.lstsq(dataz,z)
    if robust:
        outputs = lsq(myplane, np.r_[lsqcoeffs],
                      args=([x-x0,
                             y-y0, z]),
                      loss = 'soft_l1')
        lsqcoeffs = outputs.x
    return(lsqcoeffs)
    
def gradient(results):
    gradient = ((results[1]**2) + (results[2]**2))**(1/2)
    return(gradient)


def angular_momentum_implied(gradient,dist_perp_1,dist_perp_2):
    '''
    Computes the specific angular momentum given input parameters
    '''
    ang_mom_implied = beta * gradient * dist_perp_1 * dist_perp_2
    return(ang_mom_implied)

def angular_momentum_actual(data_object):
    '''
    Takes input data object and computes actual angular momentum values
    NON SPECIFIC
    '''
    angular_momentum_x = data_object.sum('angular_momentum_x')
    angular_momentum_y = data_object.sum('angular_momentum_y') 
    angular_momentum_z = data_object.sum('angular_momentum_z')
    angular_momentum_xy = ((angular_momentum_x**2)+
                           (angular_momentum_y**2))**(1/2)
    angular_momentum_xz = ((angular_momentum_x**2)+
                           (angular_momentum_z**2))**(1/2)
    angular_momentum_yz = ((angular_momentum_y**2)+
                           (angular_momentum_z**2))**(1/2)
    angular_momentum_total = ((angular_momentum_x**2)+
                              (angular_momentum_y**2)+
                              (angular_momentum_z**2))**(1/2)
    return(angular_momentum_total,
           angular_momentum_xy,
           angular_momentum_xz,
           angular_momentum_yz)


#%%
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
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


#%%
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
    
    dist_span_cm[i] = dist_tot

mass_clump_g = np.array(mass_clump_g)
cm_to_pc = 3.0857e18
g_to_sol = 1.989e33
mass_clump_solar = np.divide(mass_clump_g,g_to_sol)
dist_span_pc = np.divide(dist_span_cm,cm_to_pc)
# =============================================================================
# CHECK IN POINT: What we have so far here in the script
#
# list known as data_object_clump which has each of the regions found by the 
# clumping algorithm.
# =============================================================================



#%%
#Defining Empty Lists to store data into
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
    arr_z, vz = velocity_array(clump,'velocity_z','z')
    arr_x, vx = velocity_array(clump,'velocity_x','x')
    arr_y, vy = velocity_array(clump,'velocity_y','y')
    # =============================================================================
    
    # =============================================================================
    # Computing Reduced Velocity Arrays
    arr_x_red, vx_py, vx_pz = velocity_array_reducer(arr_x,vx,'x')
    arr_y_red, vy_px, vy_pz = velocity_array_reducer(arr_y,vy,'y')
    arr_z_red, vz_px, vz_py = velocity_array_reducer(arr_z,vz,'z')
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
    ang_mom_specific_x = angular_momentum_implied(gradient_x,
                                                  dist_span_cm[i,1,0],
                                                  dist_span_cm[i,2,0])
    ang_mom_specific_y = angular_momentum_implied(gradient_y,
                                                  dist_span_cm[i,0,0],
                                                  dist_span_cm[i,2,0])
    ang_mom_specific_z = angular_momentum_implied(gradient_z,
                                                  dist_span_cm[i,0,0],
                                                  dist_span_cm[i,1,0])
    
    ang_mom_actual_total, ang_mom_actual_xy, ang_mom_actual_xz, ang_mom_actual_yz = angular_momentum_actual(clump)
    
    
    # =========================================================================
    #     STORAGE OF VALUES FOR LATER USE
    grad_x.append(gradient_x)
    grad_y.append(gradient_y)
    grad_z.append(gradient_z)
    am_implied_x.append(ang_mom_specific_x)
    am_implied_y.append(ang_mom_specific_y)
    am_implied_z.append(ang_mom_specific_z)
    am_actual_total.append(ang_mom_actual_total)
    am_actual_partial_xy.append(ang_mom_actual_xy)
    am_actual_partial_xz.append(ang_mom_actual_xz)
    am_actual_partial_yz.append(ang_mom_actual_yz)
    # =========================================================================


# Turning all the data into numpy arrays to use later
am_actual_total = np.array(am_actual_total)
am_actual_partial_xy = np.array(am_actual_partial_xy)
am_actual_partial_xz = np.array(am_actual_partial_xz)
am_actual_partial_yz = np.array(am_actual_partial_yz)
am_implied_x = np.array(am_implied_x)
am_implied_y = np.array(am_implied_y)
am_implied_z = np.array(am_implied_z)


