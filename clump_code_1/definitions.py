#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 14:43:26 2018

@author: sfielder

This file contains all the definitions used for the clump finding coding.
Default coding structure for comments inside modules:
"""    
"""

Arguments:
----------

Global: for global parameter inputs

Returns:
--------


"""

#Importing Modules used in the definitions
#Importing packages here
import numpy as np
import yt as yt
from yt.analysis_modules.level_sets.api import Clump
from yt.analysis_modules.level_sets.api import find_clumps
from yt.analysis_modules.level_sets.api import get_lowest_clumps



import matplotlib.pyplot as plt
from scipy.optimize import least_squares as lsq

#For Visualization of Data (plane_fit_visualization definition)
import scipy.linalg
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick





def octant_split(data_object,ds,l):
    """
    Splits cubic data object into a sub-set of 8 data objects, i.e. octants.
    	
    Arguments:
    ---------- 
    data_object: YTRegion
    (usually representing whole simulation space)
    ds: yt object
        original data set given from importing file
    l: float
        length of data_object along the side, in units of pc.
        
    Returns:
    -------- 
    dbox_array: list of YTRegion 
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
        dbox_array.append(ds.r[(x1[i],'pc'):(x2[i],'pc'),
                               (y1[i],'pc'):(y2[i],'pc'),
                               (z1[i],'pc'):(z2[i],'pc')])	
    return(dbox_array)

def master_clump_maker(data_object):
    '''
    Creates master clump for imputted data object.
    
    Arguments:
    ---------- 
    data_object: YTRegion
    
    Returns:
    -------- 
    master_clump: object in yt
    '''
    #Defining Field for Contouring
    field = ("gas", "density")
    master_clump = Clump(data_object, field)
    return(master_clump)

def clump_finder(master_clump,clump_sizing,cmin,cmax,step):
    '''
    Finds clumps that match certain criteria.
    
    Arguments:
    ---------- 
    master_clump: yt object
        from master clump maker
    clump sizing parameter: float
        minumum cell count to make a clump
    c_min: float
        minimum density threshold
    c_max: float
        maximum density threshold
    step: float
        increment to work through clump algorithm
    
    Returns:
    -------- 
    leaf_clumps: yt Clump object
        Lowest Tree Values (leaves) of the clump data
    '''
    #Setting Parameters for Clump Finding
    master_clump.add_validator("min_cells", clump_sizing)
    #Adds Center of Mass info for Clumps
    master_clump.add_info_item("center_of_mass") 
    #Finding Clumps Here
    find_clumps(master_clump, cmin, cmax, step)
    leaf_clumps = get_lowest_clumps(master_clump)
    return(leaf_clumps)

def center_of_mass(lc):
    '''
    Computes the center of mass of a imputted data object
    
    Arguments:
    ----------
    lc: yt Clump object
        Lowest Leaf Clump Data Object
    
    Returns:
    -------- 
    com: float
        3-list of x,y,z coordinates of center of mass
    '''
    com = np.array(lc.quantities.center_of_mass())
    return(com)

def bounding_box(lc):
    '''
    Writes the lower and upper bounds for each axis of a box, given a data
    clump object
    
    Arguments:
    ----------
    lc: yt Clump object
        Lowest Lead Clump Data Object
        
    Returns:
    -------- 
    bounding_box: float
        (3,2)-list of x,y,z min/max values
    '''
    bounding_box = np.array(lc.quantities.extrema([('gas','x'),
                                                   ('gas','y'),
                                                   ('gas','z')]))
    return(bounding_box)
    
def clump_box(data_object_original,br):
    '''
    Takes original data set, and constructs a data object according to
    the bounding box given.
    Arguments:
    ----------
    data_object_original: YTRegion object
        Original Data in File (ds, usually called)
    br: float
        bounding region of the clump
        
    Returns:
    -------- 
    dbox_clump: YTRegion
        Data object which is constructed according to bounding box,
        represents a data object around a clump region.
    '''
    dbox_clump = data_object_original.r[(br[0,0],'cm'):(br[0,1],'cm'),
                                        (br[1,0],'cm'):(br[1,1],'cm'),
                                        (br[2,0],'cm'):(br[2,1],'cm')]
    return(dbox_clump)
    
def velocity_array(data_object,velocity,axis,master_dist_data,l):
    '''
    Takes Data object and integrates velocity along chosen line of sight
    and returns an array that matches data structure of original data set
    with only the clump having values in the array, all others are nan
    
    Arguments:
    ----------
    data_object: YTRegion object
        Data object (clump or otherwise)
    velocity: string
        the component of velocity to be integrated
            options include: 
                - 'velocity_x'
                - 'velocity_y'
                - 'velocity_z'
    axis: string
        line of sight along which velocity component is integrated
            options include: (usually same component as velocity)
                - 'x'
                - 'y'
                - 'z'
    master_dist_data: float
        dimensional number along one axis for the data set
    l: float
        distance the original data set along one axis in pc
    
    Returns:
    -------- 
    v_arr: array
        velocity integrated array, in the dimensions of the original dataset
    v_object: YTQuadTree
        velocity integrated data object
    '''
    v_object = data_object.integrate(velocity, weight='density', axis=axis)
    v_reform = v_object.to_frb((l,'pc'),(master_dist_data,master_dist_data))
    v_arr = np.array(v_reform[velocity])
    v_arr = np.reshape(v_arr,(master_dist_data,master_dist_data))
    return(v_arr,v_object)
    
def velocity_array_reducer(velocity_int_array,
                           velocity_int_data_object,
                           axis,
                           master_dist_data):
    '''
    Takes velocity integrated array that is in the original dimension size
    with imbedded data values surrounded by nan values
    and reduces that array into just the imbedded data,
    keeping the general shape of the data intact.
    Arguments:
    ----------
    velocity_int_array: array
        this is the velocity integrated array from velocity_array module
    velocity_int_data: YTQuadTree
        associated data object with the aformentioned array
    axis: string
        should be in the style 'i', where axis matches the axis chosen for 
        the two previously mentioned inputs.
    master_dist_data: float
        dimensional number along one axis for the data set
    
    Returns:
    -------- 
    velocity_int_array_reduced: array
        reduced array with all non nan values ridden
    vi_pj: array
        first perpendicular component for coordinate system
    vu_pw: array
        second perpendicular component for coordinate system
    broken: int
        signifies broken code:
            value == 0: not broken - passes normally
            value == 1: broken due to v_positions having length 0 (empty)
            value == 2: broken due to failure to reshape

    Extra Notes:
    
    Relationship between line-of-sight and 'px' and 'py' values were found by
    comparting 'px' and 'py' to 'x','y','z' values for an arbitrary data
    object. The relationship is hard coded in this definition, as the data
    extracted from simulations will not be different in this case.
    '''
    #Setting transferable quantity to exit loop if error occurs
    broken=0
    
    #Find where the data is non-nan valued
    v_positions = np.argwhere(~np.isnan(velocity_int_array))
    
    #Exit loop with broken=1 if v_positions is empty
    if len(v_positions) == 0:
        broken=1
        arr = np.empty((1,1))
        perp_coord_1 = np.empty((1,1))
        perp_coord_2 = np.empty((1,1))
        return(arr,perp_coord_1,perp_coord_2,broken)
        
    v_positions_ij = []    #Creating list for array slicing
    v_positions_ij.append(v_positions[0,0]) #First Row Value
    v_positions_ij.append(v_positions[-1,0]) #Last Row Value
    v_positions_ij.append(v_positions[0,1])    #First Column Value
    v_positions_ij.append(v_positions[-1,1])   #Last Column Value
    
    while True:
        try:
            if axis == 'z':
                # Finds appropriate px coordinates and py coordinates
                vx_coordinates = np.array(velocity_int_data_object['py'])
                vx_coordinates = np.reshape(vx_coordinates,
                                            (master_dist_data,master_dist_data))
                vz_px = vx_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,
                                       v_positions_ij[2]:v_positions_ij[3]+1]
                vy_coordinates = np.array(velocity_int_data_object['px'])
                vy_coordinates = np.reshape(vy_coordinates,
                                            (master_dist_data,master_dist_data))
                vz_py = vy_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,
                                       v_positions_ij[2]:v_positions_ij[3]+1]
                velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                                v_positions_ij[2]:v_positions_ij[3]+1]
               
                return(velocity_int_array_reduced,vz_px,vz_py,broken)
            
            if axis == 'y':
                # Finds appropriate px coordinates and py coordinates
                vx_coordinates = np.array(velocity_int_data_object['px'])
                vx_coordinates = np.reshape(vx_coordinates,
                                            (master_dist_data,master_dist_data))
                vy_px = vx_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,
                                       v_positions_ij[2]:v_positions_ij[3]+1]
                vz_coordinates = np.array(velocity_int_data_object['py'])
                vz_coordinates = np.reshape(vz_coordinates,
                                            (master_dist_data,master_dist_data))
                vy_pz = vz_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,
                                       v_positions_ij[2]:v_positions_ij[3]+1]
                velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                                v_positions_ij[2]:v_positions_ij[3]+1]
                return(velocity_int_array_reduced,vy_px,vy_pz,broken)
                
            if axis == 'x':
                # Finds appropriate px coordinates and py coordinates
                vy_coordinates = np.array(velocity_int_data_object['py'])
                vy_coordinates = np.reshape(vy_coordinates,
                                            (master_dist_data,master_dist_data))
                vx_py = vy_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,
                                       v_positions_ij[2]:v_positions_ij[3]+1]
                vz_coordinates = np.array(velocity_int_data_object['px'])
                vz_coordinates = np.reshape(vz_coordinates,
                                            (master_dist_data,master_dist_data))
                vx_pz = vz_coordinates[v_positions_ij[0]:v_positions_ij[1]+1,
                                       v_positions_ij[2]:v_positions_ij[3]+1]
                velocity_int_array_reduced = velocity_int_array[v_positions_ij[0]:v_positions_ij[1]+1,
                                                                v_positions_ij[2]:v_positions_ij[3]+1]
            return(velocity_int_array_reduced,vx_py,vx_pz,broken)
        
        except ValueError:
            broken=2
            arr = np.empty((1,1))
            perp_coord_1 = np.empty((1,1))
            perp_coord_2 = np.empty((1,1))
            return(arr,perp_coord_1,perp_coord_2,broken)
        

def array_flattener(x):
    '''
    Flattens arrays using .ndarray.flatten command
    
    Arguments:
    ----------
    x: array
    
    Returns:
    -------- 
    arr_flat: array
        flattened array
    '''
    arr_flat = np.ndarray.flatten(x)
    return(arr_flat)
        
def myplane(p, x, y, z):
    '''
    Imported from Erik Rosolowsky. Creates data that encompases a plane
    fitted from data
    '''
    return(p[0] + p[1] * x + p[2] * y - z)

def plane_fit_visualization(x,y,z,k):
    '''
    Takes flattened data, plane fits in, and visualizes it in 3d.
    
    Arguments:
    ----------
    x,y: array 
        mesh coordinates (flattened)
    z: array
        coordinate data points to fit plane to
    k: integer
        used for naming the output files, should match clump number being
        worked on
    Returns:
    -------- 
    Special: PDF and PNG of matplotlib figure
    
    Notes:
        
    This method of fitting the data is different than the plane_fit
    definition used for computation of gradients. The resulting coefficients
    are different from those in plane_fit, thus this should only be used as 
    a rough indication of how the plane fitting process works, only use for 
    visualization of the data if it is needed at all.
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
    lsqcoeffs: list
       3-element vector with components[ z0 (constant offset) , grad_x, grad_y]
    
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
    '''
    Takes output from lsqcoeffs and creates a gradient value from it.
    
    Arguments:
    ----------
    result: array
        must be a array of size 3 with [z0,grad_x,grad_y] format
        same as the output from lsqcoeffs
        
    Returns:
    --------
    gradient: float
    '''
    gradient = ((results[1]**2) + (results[2]**2))**(1/2)
    return(gradient)


def angular_momentum_implied(gradient,dist_perp_1,dist_perp_2,beta=1):
    '''
    Computes the specific angular momentum given input parameters
    
    Arguments:
    ----------
    gradient: float
        from gradient module
    dist_perp_1: float
        distance of the length perpendicular to line of sight (1st value)
    dist_perp_2: float
        2nd Value (see above)
    beta: float
        global input for coefficient to implied angular momentum
            - Will default to one if no argument given
    
    Returns:
    --------
    ang_mom_implied: float
        implied specific angular momentum
    '''
    ang_mom_implied = beta * gradient * dist_perp_1 * dist_perp_2
    return(ang_mom_implied)

def angular_momentum_actual(data_object,mass):
    '''
    Takes input data object and computes specific actual angular 
    momentum values.
    
    Arguments:
    ----------
    data_object: YTRregion
    
    Returns:
    --------
    angular_momentum_total: float
        total specific angular momentum of data region
    angular_momentum_xy: float
        total specific angular momemtum (only two components, (x,y))
    angular_momentum_xz: float
        total specific angular momemtum (only two components, (x,z))
    angular_momentum_yz: float
        total specific angular momemtum (only two components, (y,z))
    '''
    angular_momentum_x = data_object.sum('angular_momentum_x')
    angular_momentum_y = data_object.sum('angular_momentum_y') 
    angular_momentum_z = data_object.sum('angular_momentum_z')
    angular_momentum_xy = (((angular_momentum_x**2)+
                           (angular_momentum_y**2))**(1/2))/mass
    angular_momentum_xz = (((angular_momentum_x**2)+
                           (angular_momentum_z**2))**(1/2))/mass
    angular_momentum_yz = (((angular_momentum_y**2)+
                           (angular_momentum_z**2))**(1/2))/mass
    angular_momentum_total = (((angular_momentum_x**2)+
                              (angular_momentum_y**2)+
                              (angular_momentum_z**2))**(1/2))/mass
    return(angular_momentum_total,
           angular_momentum_xy,
           angular_momentum_xz,
           angular_momentum_yz)

def proj_clump_com(ds,axis,com):
    prj = yt.ProjectionPlot(ds,axis,("gas","density"),
                            center='c', width = (10,'pc'))
    return(prj)
