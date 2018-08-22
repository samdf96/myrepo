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
#Definitions for Gravitational Energy Calculations
from yt.utilities.lib.misc_utilities import gravitational_binding_energy

from yt.utilities.physical_constants import \
    gravitational_constant_cgs as G

import matplotlib
matplotlib.use('agg') # This is to default matplitlib to draw off screen and not throw errors when running in screen
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as lsq

#For Visualization of Data (plane_fit_visualization definition)
import scipy.linalg
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick
import itertools

#For Plotting Definition
from astropy import units as u
from astropy.io import fits

import logging
 
module_logger = logging.getLogger("initialize.definitions")




def OctantSplit(data_object,ds,l):
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
    logger = logging.getLogger("initialize.definitions.OctantSplit")
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
    logger.debug("OctantSplit has been run successfully.")
    return(dbox_array)

def MasterClumpMaker(data_object):
    '''
    Creates master clump for imputted data object.
    
    Arguments:
    ---------- 
    data_object: YTRegion
    
    Returns:
    -------- 
    master_clump: object in yt
    '''
    logger = logging.getLogger("initialize.definitions.MasterClumpMaker")
    #Defining Field for Contouring
    field = ("gas", "density")
    master_clump = Clump(data_object, field)
    logger.debug("MasterClumpMaker has been run successfully.")
    return(master_clump)

def ClumpFinder(master_clump,clump_sizing,cmin,cmax,step):
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
    logger = logging.getLogger("initialize.definitions.ClumpFinder")
    #Setting Parameters for Clump Finding
    master_clump.add_validator("min_cells", clump_sizing)
    #Adds Center of Mass info for Clumps
    master_clump.add_info_item("center_of_mass") 
    #Finding Clumps Here
    while True:
        try:
            find_clumps(master_clump, cmin, cmax, step)
            leaf_clumps = get_lowest_clumps(master_clump)
            logger.debug("ClumpFinder has been run successfully.")
            print('ClumpFinder has been run successfully.')
        except TypeError:
            leaf_clumps = []
            logger.debug('ClumpFinder has not been run successfully.')
            print('ClumpFinder has not been run successfully.')
            break
        break

    return(leaf_clumps)

def CenterOfMass(lc):
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
    logger = logging.getLogger("initialize.definitions.CenterOfMass")
    com = np.array(lc.quantities.center_of_mass())
    logger.debug("Center Of Mass found to be: %s", str(com))
    logger.debug("CenterOfMass has run successfully.")
    return(com)

def BoundingBox(lc):
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
    logger = logging.getLogger("initialize.definitions.BoundingBox")
    bounding_box = np.array(lc.quantities.extrema([('gas','x'),
                                                   ('gas','y'),
                                                   ('gas','z')]))
    logger.debug("Bounding Box found to be: %s", str(bounding_box))
    logger.debug("BoundingBox has run successfully.")
    return(bounding_box)
    
def ClumpBox(data_object_original,br):
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
    logger = logging.getLogger("initialize.definitions.ClumpBox")
    dbox_clump = data_object_original.r[(br[0,0],'cm'):(br[0,1],'cm'),
                                        (br[1,0],'cm'):(br[1,1],'cm'),
                                        (br[2,0],'cm'):(br[2,1],'cm')]
    logger.debug("ClumpBox has been run successfully.")
    return(dbox_clump)
    
def VelocityArray(data_object,velocity,axis,master_dist_data,l):
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
    broken: float
        used for error detection
    '''
    logger = logging.getLogger("initialize.definitions.VelocityArray")
    broken=0    #Setting default value for error detection
    logger.debug("Broken Export Value has initialized to: %s", str(broken))
    while True:
        try:
            v_object = data_object.integrate(velocity, weight='density', axis=axis)
            v_reform = v_object.to_frb((l,'pc'),(master_dist_data,master_dist_data))
            v_arr = np.array(v_reform[velocity])
            v_arr = np.reshape(v_arr,(master_dist_data,master_dist_data))
            logger.debug("VelocityArray has been run successfully.")
            return(v_arr,v_object,broken)
        except RuntimeError:
            logging.exception("RunTimeError has been caught.")
            broken=1
            logging.warning("Broken Export Value has been overwritten to: %s",
                            broken)
            logging.debug("Setting v_object to False, and v_arr to zeros.")
            v_object = False
            v_arr = np.zeros((master_dist_data,master_dist_data))
            logger.debug("VelocityArray has not been run successfully.")
            return(v_arr,v_object,broken)

    
def VelocityArrayReducer(velocity_int_array,
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
    logger = logging.getLogger("initialize.definitions.VelocityArrayReducer")
    #Setting transferable quantity to exit loop if error occurs
    broken=0
    logger.debug("Broken Export Value has initialized to: %s", str(broken))
    
    #Find where the data is non-nan valued
    v_positions = np.argwhere(~np.isnan(velocity_int_array))
    
    #Exit loop with broken=1 if v_positions is empty
    if len(v_positions) == 0:
        logger.debug("v_positions check has returned that it's length is zero.")
        broken=1
        logger.warning("Broken Export Value is being overwritten to: %s", broken)
        logger.debug("Setting all export values of this function to empty.")
        arr = np.empty((1,1))
        perp_coord_1 = np.empty((1,1))
        perp_coord_2 = np.empty((1,1))
        logger.debug("VelocityArrayReducer has not been run successfully.")
        return(arr,perp_coord_1,perp_coord_2,broken)
    
    logger.debug("Setting v_positions.")
    v_positions_ij = []                         #Creating list for array slicing
    v_positions_ij.append(v_positions[0,0])     #First Row Value
    v_positions_ij.append(v_positions[-1,0])    #Last Row Value
    v_positions_ij.append(v_positions[0,1])     #First Column Value
    v_positions_ij.append(v_positions[-1,1])    #Last Column Value
    logger.debug("v_positions set to: %s", str(v_positions))
    while True:
        try:
            if axis == 'z':
                logger.debug("axis string input detected as: %s", axis)
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
                logger.debug("VelocityArrayReducer has been run successfully on axis: %s",
                            str(axis))
                return(velocity_int_array_reduced,vz_px,vz_py,broken)
            
            if axis == 'y':
                logger.debug("axis string input detected as: %s", axis)
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
                logger.debug("VelocityArrayReducer has been run successfully on axis: %s",
                            str(axis))
                return(velocity_int_array_reduced,vy_px,vy_pz,broken)
                
            if axis == 'x':
                logger.debug("axis string input detected as: %s", str(axis))
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
                logger.debug("VelocityArrayReducer has been run successfully on axis: %s",
                            str(axis))
            return(velocity_int_array_reduced,vx_py,vx_pz,broken)
        
        except ValueError:
            logger.exception("ValueError has been caught.")
            broken=2
            logging.warning("Broken Export Value has been overwritten to: %s",
                            broken)
            logging.debug("Setting arr, perp_coord1a and perp_coord2 to zeros.")
            arr = np.empty((1,1))
            perp_coord_1 = np.empty((1,1))
            perp_coord_2 = np.empty((1,1))
            logging.debug("VelocityArrayReducer has not been run successfully.")
            return(arr,perp_coord_1,perp_coord_2,broken)
        

def ArrayFlattener(x):
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
    logger = logging.getLogger("initialize.definitions.ArrayFlattener")
    arr_flat = np.ndarray.flatten(x)
    logger.debug("ArrayFlattener has been run successfully.")
    return(arr_flat)
        
def MyPlane(p, x, y, z):
    '''
    Imported from Erik Rosolowsky. Creates data that encompases a plane
    fitted from data
    '''
    logger = logging.getLogger("initialize.definitions.MyPlane")
    logger.debug("Returning Flattened plane array.")
    return(p[0] + p[1] * x + p[2] * y - z)

def PlaneFitVisualization(x,y,z,k):
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
    logger = logging.getLogger("initialize.definitions.PlaneFitVisualization")
    data = np.c_[x,y,z]
    logger.debug("data has been set to: %s", str(data))
    mn = np.min(data, axis=0)
    logger.debug("Minimum Value for Data found to be: %s", str(mn))
    mx = np.max(data, axis=0)
    logger.debug("Maximum Value for Data found to be: %s", str(mn))
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
    logger.debug("Meshgrid Created.")
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
    logger.debug("Results Coefficients from linalg.lstsq: %s", str(C))
    
    # evaluate it on grid
    Z = C[0]*X + C[1]*Y + C[2]
    logger.debug("Z values grid set as: %s", str(Z))
    
    logger.debug("Creating Figure Object.")
    fig = plt.figure()
    logger.debug("Setting ax object to 3D axes.")
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=1)
    logger.debug("Setting Labels")
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.zaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax.set_zlabel('Z')
    logger.debug("Setting Axes Properties.")
    ax.axis('equal')
    ax.axis('tight')
    logger.debug("Saving Figure Object to PDF and PNG.")
    plt.savefig("Plane_Fitting_Clump_{k}.pdf".format(k=k), bbox_inches='tight')
    plt.savefig("Plane_Fitting_Clump_{k}.png".format(k=k), bbox_inches='tight')        
    plt.show()
    logger.idebug("PlaneFitVisualization has been run successfully.")
    return()


def PlaneFit(x, y, z, robust=False):
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
    logger = logging.getLogger("initialization.definitions.PlaneFit")
    logger.debug("Setting Median Values.")
    x0, y0 = np.median(x), np.median(y)
    logger.debug("x median value found to be: %s", str(x0))
    logger.debug("y median value found to be: %s", str(y0))
    dataz = np.c_[np.ones(x.size), 
                  x-x0, 
                  y-y0]
    logger.debug("dataz to fit set as: %s", str(dataz))
    lsqcoeffs, _, _, _ = np.linalg.lstsq(dataz,z)
    logger.debug("Coefficients from fit found to be: %s", str(lsqcoeffs))
    if robust:
        outputs = lsq(MyPlane, np.r_[lsqcoeffs],
                      args=([x-x0,
                             y-y0, z]),
                      loss = 'soft_l1')
        lsqcoeffs = outputs.x
    logger.debug("PlaneFit has been run successfully.")
    return(lsqcoeffs)
    
def Gradient(results):
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
    logger = logging.getLogger("initialize.definitions.Gradient")
    gradient = ((results[1]**2) + (results[2]**2))**(1/2)
    logger.debug("Gradient has been run successfully.")
    return(gradient)


def AngularMomentumImplied(gradient,dist_perp_1,dist_perp_2,beta=1):
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
    logger = logging.getLogger("initialize.definitions.AngularMomentumImplied")
    ang_mom_implied = beta * gradient * dist_perp_1 * dist_perp_2
    logger.debug("ang_mom_implied found to be: %s", str(ang_mom_implied))
    logger.debug("AngularMomentumImplied has been run successfully.")
    return(ang_mom_implied)

def AngularMomentumActual(data_object,mass):
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
    logger = logging.getLogger("initialize.definitions.AngularMomentumActual")
    logger.debug("Computing the Sum of all components of angular momentum.")
    angular_momentum_x = data_object.sum('angular_momentum_x')
    angular_momentum_y = data_object.sum('angular_momentum_y')
    angular_momentum_z = data_object.sum('angular_momentum_z')
    logger.debug("X component found to be: %s", str(angular_momentum_x))
    logger.debug("Y component found to be: %s", str(angular_momentum_y))
    logger.debug("Z component found to be: %s", str(angular_momentum_z))
    logger.debug("Combining Sums on: xy, xz, yz, and total.")
    angular_momentum_xy = (((angular_momentum_x**2)+
                          (angular_momentum_y**2))**(1/2))/mass
    angular_momentum_xz = (((angular_momentum_x**2)+
                           (angular_momentum_z**2))**(1/2))/mass
    angular_momentum_yz = (((angular_momentum_y**2)+
                           (angular_momentum_z**2))**(1/2))/mass
    angular_momentum_total = (((angular_momentum_x**2)+
                              (angular_momentum_y**2)+
                              (angular_momentum_z**2))**(1/2))/mass
    logger.debug("xy sum found to be: %s", str(angular_momentum_xy))
    logger.debug("xz sum found to be: %s", str(angular_momentum_xz))
    logger.debug("yz sum found to be: %s", str(angular_momentum_yz))
    logger.debug("total found to be: %s", str(angular_momentum_total))
    logger.debug("AngularMomentumActual has been run successfully.")
    return(angular_momentum_total,
           angular_momentum_xy,
           angular_momentum_xz,
           angular_momentum_yz)

def KineticEnergy(data_object,
                  bulk_velocity):
    """
    Calculates the Kinetic Energy of a YTRegion
    
    Inputs:
        - data_object: YTRegion
        - bulk_velocity: array (quantity)
            - 3 components (x,y,z) in cm/s
            
    Outputs:
        - kinetic: quantity
            - Value + unit of energy in cgs
    """
    logger = logging.getLogger("initialize.definitions.KineticEnergy")
    kinetic = 0.5 * (data_object['cell_mass'] *
        ((bulk_velocity[0] - data_object["gas", "velocity_x"].value)**2 +
         (bulk_velocity[1] - data_object["gas", "velocity_y"].value)**2 +
          (bulk_velocity[2] - data_object["gas", "velocity_z"].value)**2)).sum()
    logger.debug("KineticEnergy has been run successfully.")
    return(kinetic)
    
def GravitationalEnergy(data_object, kinetic_energy):
    """
    Calculates the Gravitational Energy of a YTRegion
    Inputs:
        - data_object: YTRegion
        - kinetic_energy: quantity
    Outputs:
        - grav_energy: quantity
    
    This will return grav_energy with the wrong units, because gravitational_binding_energy does not output with units (or as quantity).
    """
    logger = logging.getLogger("initialize.definitions.GravitationalEnergy")
    grav_energy = G * gravitational_binding_energy(data_object['cell_mass'],
                                                   data_object['x'],
                                                   data_object['y'],
                                                   data_object['z'],
                                                   True,
                                                   (kinetic_energy/G))
    logger.debug("GravitationalEnergy has been run successfully.")
    return(grav_energy)


def ProjCreator(ds,
                 data_object,
                 com,
                 com_x,
                 com_y,
                 com_z,
                 save_directory,
                 call_str,
                 time_stamp):
    """
    Takes the center of mass coordinates, and overlays markers on the simulation
    projection plots to where those clumps are located.
    
    Arguments:
    ----------
    ds: data set
        This is the original data set file object
    data_object: YTDataContainer
        This is the all_data() object created for the simulation file
    com_x: array
        This has the center of mass data for x LOS (y,z coodinates)
    com_y: array
        This has the center of mass data for y LOS (x,z coodinates)
    com_z: array
        This has the center of mass data for z LOS (x,y coodinates)
    save_directory: string
        Has the specific location which plots are to be saved to
    call_str: string
         Number for simulation number
    time_stamp: string
        Time Stamp for Data simulation
    
    Returns:
    --------
    For Testing, can return the prj object, but this definition will auto save
    the plot into the correct directory
    """
    logger = logging.getLogger("initialize.definitions.ProjCreator")
    
    #Defining centimetre unit for quantity
    logger.debug("Creating Quantities out of com_x/y/z values: set to u.cm")
    com_x = com_x * u.cm
    com_y = com_y * u.cm
    com_z = com_z * u.cm
    logger.debug("com_x set to: %s", str(com_x))
    logger.debug("com_y set to: %s", str(com_y))
    logger.debug("com_z set to: %s", str(com_z))
    
    #Converting to Parsec for plot overlay
    logger.debug("Converting Quantities to parsec units.")
    com_x_pc = com_x.to(u.parsec)
    com_y_pc = com_y.to(u.parsec)
    com_z_pc = com_z.to(u.parsec)
    logger.debug("com_x reset to: %s", str(com_x_pc))
    logger.debug("com_y reset to: %s", str(com_y_pc))
    logger.debug("com_z reset to: %s", str(com_z_pc))
    
    ## For x LOS
    logger.debug("Working on X LOS.")
    prj_x = yt.ProjectionPlot(ds,
                            'x',
                            ("gas","density"),
                            center='c',
                            width = (10,'pc'),
                            data_source=data_object)
    logger.debug("Projection Plot created.")
    for i in range(0,len(com_x)):
        prj_x.annotate_marker(com_x_pc[i],
                            coord_system='plot',
                            plot_args={'color':'red','s':500})
    logger.debug("COM markers have been overlayed.")
    prj_x.save(save_directory+call_str+"_"+time_stamp+"_x_los_clump_marker.pdf")
    logger.debug("Projection Plot has been saved.")
    
    ## For y LOS
    logger.debug("Working on Y LOS.")
    prj_y = yt.ProjectionPlot(ds,
                            'y',
                            ("gas","density"),
                            center='c',
                            width = (10,'pc'),
                            data_source=data_object)
    logger.debug("Projection Plot created.")
    for i in range(0,len(com_y)):
        prj_y.annotate_marker(com_y_pc[i],
                            coord_system='plot',
                            plot_args={'color':'red','s':500})
    logger.debug("COM markers have been overlayed.")
    prj_y.save(save_directory+call_str+"_"+time_stamp+"_y_los_clump_marker.pdf")
    logger.debug("Projection Plot has been saved.")
    
    ## For z LOS
    logger.debug("Working on Z LOS.")
    prj_z = yt.ProjectionPlot(ds,
                            'z',
                            ("gas","density"),
                            center='c',
                            width = (10,'pc'),
                            data_source=data_object)
    logger.debug("Projection Plot created.")
    for i in range(0,len(com_z)):
        prj_z.annotate_marker(com_z_pc[i],
                            coord_system='plot',
                            plot_args={'color':'red','s':500})
    logger.debug("COM markers have been overlayed.")
    prj_z.save(save_directory+call_str+"_"+time_stamp+"_z_los_clump_marker.pdf")
    logger.debug("Projection Plot has been saved.")
    
    logger.debug("All Figures have been saved.")
    logger.debug("ProjCreator has been run successfully.")
    return()

def DataGrabber(data):
    """
    Function that takes in the hdu_table.data object, and makes a dictionary
    with all the data pertaining to the plotting that will be done
    
    Inputs:
        data: hdu_table.data object
    
    Outputs:
        dict: dictionary
            Has all the relevant quantities inside
    
    Notes: Can expand this to include all the data inside the object if needed
    """
    logger = logging.getLogger("initialize.definitions.DataGrabber")
    
    
    clump_number = data['Clump Number']
    mass = data['Mass'] 
    
    actual_tot = data['Actual Total Angular Momentum']
    actual_par_xy = data['Actual Partial Angular Momentum (x+y components)']
    actual_par_xz = data['Actual Partial Angular Momentum (x+z components)']
    actual_par_yz = data['Actual Partial Angular Momentum (y+z components)']
    
    grad_x = data['Plane Fitting Gradient (x LOS)']
    grad_y = data['Plane Fitting Gradient (y LOS)']
    grad_z = data['Plane Fitting Gradient (z LOS)']
    
    imp_x_los = data['Implied Total Angular Momentum (x LOS)']
    imp_y_los = data['Implied Total Angular Momentum (y LOS)']
    imp_z_los = data['Implied Total Angular Momentum (z LOS)']
    
    volume_pix = data['Volume (pixels)']
    
    logger.debug("All the Data extracted.")
    logger.debug("Making Dictionary with all key/value pairs.")
    dict = {'clump_number': clump_number,
            'mass': mass,
            'actual_tot': actual_tot,
            'actual_par_xy': actual_par_xy,
            'actual_par_xz': actual_par_xz,
            'actual_par_yz': actual_par_yz,
            'imp_x_los': imp_x_los,
            'imp_y_los': imp_y_los,
            'imp_z_los': imp_z_los,
            'volume_pix': volume_pix,
            'grad_x': grad_x,
            'grad_y': grad_y,
            'grad_z': grad_z}
    
    logger.debug("DataGrabber has been run successfully.")
    return(dict)
    
def jCompPlotter(x, y1, y2, axis_str, tit_str):
    """
    Function that plots j/j for a specific in a time step in a specific Simulation
    Run.
    
    Inputs:
        - x: array
            specific angular momentum corresponding to axis given
        - y1: array
            Full component data set for actual specific angular momentum
        - y2: array
            Partial component data set for actual specific angular momentum
        - axis_str: string
            line of sight argument: 'X' , 'Y' , 'Z' .
        - tit_str: string
            Either 'Regular' or 'Colormapped' - sets title, and savename
    Outputs:
        - fig: matplotlib object
            This is directed back out into previously tree file, and saved
            approprately to a save location pre-defined in the code.
    """
    logger = logging.getLogger("initialize.definitions.jCompPlotter")
    #Insert Actual and Partial Data Here in LOGLOG style
    
    logger.debug("Setting Axis/Title Strings.")
    x_str = 'Implied Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    y_str = 'Actual Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    title = 'j Comparison - '+axis_str+ ' LOS'
    
    logger.debug("Creating Figure Environment.")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    logger.debug("Adding Main Data.")
    ax.loglog(x, y1, 'b.', label='Full')
    ax.loglog(x, y2, 'r.', label='Partial')
    
    #Grabs Min and Max for x-axis (searches all imputted data however)
    x_min = min(x)
    x_max = max(x)
    length = len(x)
    logger.debug("x_min found to be: %s", str(x_min))
    logger.debug("x_max found to be: %s", str(x_max))
    
    logger.debug("Creating Line of Constant Slope.")
    unity_x = np.linspace(x_min,x_max,1000) #Makes line of unity values
    ax.loglog(unity_x,unity_x, 'k--', label='Line of Unity') #Adds to plot here
    
    if length > 10:
        #Best Fit Calculations here
        logger.debug("Starting Line of Best Fit Calculations.")
        log_x = np.log10(x)
        log_y1 = np.log10(y1)
        log_y2 = np.log10(y2)
        slope_1, intercept_1 = np.polyfit(log_x, log_y1, 1)
        slope_2, intercept_2 = np.polyfit(log_x, log_y2, 1)
        
        logger.debug("Creating Data for Lines.")
        fit_1 = 1e1**(slope_1 * log_x + intercept_1)
        fit_2 = 1e1**(slope_2 * log_x + intercept_2)
        
        #Insert on Plot Here
        logger.debug("Inserting on Plot.")
        ax.loglog(x,
                fit_1,
                'b-',
                alpha=0.5,
                label='Line of Best Fit - Full')
        ax.loglog(x,
                fit_2,
                'r-',
                alpha=0.5,
                label='Line of Best Fit - Partial')
    
    
    logger.debug("Tweaking plot settings.")
    ax.set_aspect('equal')   
    ax.set_facecolor('#f2f2f2')
    ax.grid()
    
    #Labelling and Legend
    logger.debug("Setting plot labels.")
    ax.legend(bbox_to_anchor=(1.02, 1)) 
    ax.set_xlabel(x_str)
    ax.set_ylabel(y_str)
    ax.set_title(title)

    logger.debug("jCompPlotter has been run successfully.")
    return(fig)

def jCompPlotterColormap(x, y1, y2, mass, axis_str, tit_str):
    """
    Function that plots j/j for a specific in a time step in a specific Timestep
    Run - colormapped by mass argument given.
    
    Inputs:
        - x: array
            specific angular momentum corresponding to axis given
        - y1: array
            Full component data set for actual specific angular momentum
        - y2: array
            Partial component data set for actual specific angular momentum
        - mass: array
            Corresponding mass terms for the inputted data set.
        - axis_str: string
            line of sight argument: 'X' , 'Y' , 'Z' .
        - tit_str: string
            Either 'Regular' or 'Colormapped' - sets title, and savename
    Outputs:
        - fig: matplotlib object
            This is directed back out into previously tree file, and saved
            approprately to a save location pre-defined in the code.
    """
    logger = logging.getLogger("initialize.definitions.jCompPlotterColormap")
    
    logger.debug("Setting Axis/Title Strings.")
    x_str = 'Implied Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    y_str = 'Actual Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    title = 'j Comparison - '+axis_str+ ' LOS - Colormapped'
    
    logger.debug("Creating Figure Environment.")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    logger.debug("Adding Main Data.")
    data1 = ax.scatter(x,
                       y1,
                       c=np.log10(mass),
                       s=100,
                       marker='.',
                       cmap='viridis',
                       linewidth=0,
                       alpha=1,
                       label='Full')
    ax.scatter(x,
               y2,
               c=np.log10(mass),
               s=40,
               marker='v',
               cmap='viridis',
               linewidth=0,
               alpha=1,
               label='Partial')
    
    #Grabs Min and Max for x-axis
    x_min = min(x)
    x_max = max(x)
    length = len(x)
    logger.debug("x_min found to be: %s", str(x_min))
    logger.debug("x_max found to be: %s", str(x_max))
    logger.debug("Creating Line of Constant Slope.")
    unity_x = np.linspace(x_min,x_max,1000) #Makes line of unity values
    ax.plot(unity_x,
            unity_x,
            linestyle='--',
            linewidth=1,
            alpha=0.75,
            c='k',
            label='Line of Unity') #Adds to plot here
    
    if length > 10:
        #Line Fitting Linear in Loglog Scaling
        logger.debug("Starting Line of Best Fit Calculations.")
        log_x = np.log10(x)
        log_y1 = np.log10(y1)
        log_y2 = np.log10(y2)
        slope_1, intercept_1 = np.polyfit(log_x, log_y1, 1)
        slope_2, intercept_2 = np.polyfit(log_x, log_y2, 1)
    
        fit_1 = 1e1**(slope_1 * log_x + intercept_1)
        fit_2 = 1e1**(slope_2 * log_x + intercept_2)
        
        #Insert on Plot Here
        logger.debug("Inserting on Plot.")
        ax.plot(x,
                fit_1,
                linestyle='-',
                linewidth=1,
                alpha=0.75,
                c='b',
                label='Line of Best Fit - Full')
        ax.plot(x,
                fit_2,
                linestyle='-',
                linewidth=1,
                alpha=0.75,
                c='r',
                label='Line of Best Fit - Partial')
    
    #Setting LOG LOG Scale for Scatter Plots
    logger.debug("Tweaking plot settings.")
    ax.set_aspect('equal')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_facecolor('#f2f2f2') # Grey Background
    ax.grid()
    
    #Labelling and Legend
    logger.debug("Setting plot labels.")
    ax.legend(bbox_to_anchor=(1.02, 1),
               bbox_transform=plt.gcf().transFigure) 
    ax.set_xlabel(x_str)
    ax.set_ylabel(y_str)
    ax.set_title(title)
    logger.debug("Setting Colorbar.")
    fig.colorbar(data1,label=r'$\log_{10}$(mass) [g]')
    
    logger.debug("jCompPlotterColormap has been run succesfully.")
    return(fig)

def DictionarySifter(d):
    """
    Takes in a dictionary of values corresponding to clump analysis data
    and sifts through the act_tot data set, and finds nan values. Once found,
    all other data in the dictionary is filtered through that, and nan rows
    are taken out to ensure proper plotting. This does not affect the fits file
    and the fits file remains with all comments and nan values for safekeeping.
    
    Inputs:
        - d: dictionary
            This is the dictionary with all the values in the fits file
    
    Outputs:
        - sifted_dict: dictionary
            See above, non-nan valued dictionary.
    """
    logger = logging.getLogger("initialize.definitions.DictionarySifter")
    
    logger.debug("Masking Data using act_tot nan values as filter.")
    
    mask_array_actual_tot = np.isnan(d['actual_tot'])
    
    mask_array_grad_x = np.equal(d['grad_x'],0)
    mask_array_grad_y = np.equal(d['grad_y'],0)
    mask_array_grad_z = np.equal(d['grad_z'],0)
    
    mask_master = np.zeros((len(mask_array_actual_tot)), dtype=bool)
    
    for i in range(0,len(mask_array_actual_tot)):
        if (mask_array_actual_tot[i]==True) or (
                mask_array_grad_x[i]==True
                ) or (
                        mask_array_grad_y[i]==True
                ) or (
                        mask_array_grad_z[i]==True):
            mask_master[i] = True
    
    imp_x_los_mask = np.ma.masked_where(mask_master == True,
                                        d['imp_x_los'])
    imp_y_los_mask = np.ma.masked_where(mask_master == True,
                                        d['imp_y_los'])
    imp_z_los_mask = np.ma.masked_where(mask_master == True,
                                        d['imp_z_los'])
    act_tot_mask = np.ma.masked_where(mask_master == True,
                                      d['actual_tot'])
    act_par_xy_mask = np.ma.masked_where(mask_master == True,
                                         d['actual_par_xy'])
    act_par_xz_mask = np.ma.masked_where(mask_master == True,
                                         d['actual_par_xz'])
    act_par_yz_mask = np.ma.masked_where(mask_master == True,
                                         d['actual_par_yz'])
    mass_mask = np.ma.masked_where(mask_master == True,
                                   d['mass'])
    
    #Compress all the masked array to array with correct data
    logger.debug("Compressing Masked array to sift out masked values.")
    imp_x_los_mask_compressed = np.ma.compressed(imp_x_los_mask)
    imp_y_los_mask_compressed = np.ma.compressed(imp_y_los_mask)
    imp_z_los_mask_compressed = np.ma.compressed(imp_z_los_mask)
    act_tot_mask_compressed = np.ma.compressed(act_tot_mask)
    act_par_xy_mask_compressed = np.ma.compressed(act_par_xy_mask)
    act_par_xz_mask_compressed = np.ma.compressed(act_par_xz_mask)
    act_par_yz_mask_compressed = np.ma.compressed(act_par_yz_mask)
    mass_mask_compressed = np.ma.compressed(mass_mask)
    
    logger.debug("Creating new dictionary with sifted data.")
    sifted_dict = {}
    sifted_dict['imp_x_los'] = imp_x_los_mask_compressed
    sifted_dict['imp_y_los'] = imp_y_los_mask_compressed
    sifted_dict['imp_z_los'] = imp_z_los_mask_compressed
    sifted_dict['act_tot'] = act_tot_mask_compressed
    sifted_dict['act_par_xy'] = act_par_xy_mask_compressed
    sifted_dict['act_par_xz'] = act_par_xz_mask_compressed
    sifted_dict['act_par_yz'] = act_par_yz_mask_compressed
    sifted_dict['mass'] = mass_mask_compressed
    
    logger.debug("DicationarySifter has been run successfully.")
    return(sifted_dict)

def jComparisonPlotter(current_file):
    """
    Master Function for the Comparison Plotter by specific timestep
    
    Inputs:
        - current_file: string
            Corresponds to the .fits file being plotted
    
    Output:
        - No Outputs, self contained.
    """
    logger = logging.getLogger("initialize.definitions.jComparisonPlotter")
    
    logger.debug("Opening inputted current file: %s", str(current_file))
    hdu = fits.open(current_file)
    logger.debug("hdu set as: %s", str(hdu))
    #Print Statements for File being worked on:
    filename_printing = current_file.split("/")[-1]
    logger.debug("Current File being worked on: %s", str(filename_printing))
    
    #Grabbing BinHDUTable object here, should always be second object
    hdu_table = hdu[1]
    logger.debug("hdu_table extracted as: %s", str(hdu_table))
    logger.debug("Extracting data out of the HDU Table Object.")
    data = hdu_table.data #Grabs data stored in the table -> FITS REC   
    logger.debug("Invoking DataGrabber function.")
    d = DataGrabber(data)
    logger.debug("Closing hdu object.")
    hdu.close()
    
    #Making the Output Directory for current file:
    save_dir_list = current_file.split("/")[:-1]
    save_dir = '/'.join(save_dir_list)
    logger.debug("Save_directory (General) set to: %s", save_dir)
    
    #Make Masked arrays for data that is nan valued
    """
    Will use the act_tot_mask array for master mask condition
    Below creates masked data from the dictionary entries using the condition 
    stated above to parse through the data and extract non-nan values
    """
    #Call the Data Sifter Function from definitions. See Note Above
    logger.debug("Invoking DictionarySifter function.")
    data_sifted = DictionarySifter(d)
    
    logger.debug("Creating title string and axes string to loop over.")
    tit_str = ['Regular','Colormapped']    #Loop over these two values
    axis_str = ['X','Y','Z']        #Loop over these LOS axes
    logger.debug("Looping over Title String Now.")
    for i in range(0,len(tit_str)): #Loop over full and partial data sets for actual sim
        logger.debug("Currently Working on tit_str: %s", str(tit_str[i]))
        logger.debug("Looping over Axis String Now.")
        for j in range(0,len(axis_str)): #Loop over all the LOS axes
            logger.debug("Currently Working on axis_str: %s", str(axis_str[j]))
            #Looping to set the correct data for the los axis
            if axis_str[j] == 'X':
                x_string = 'imp_x_los'
                y_string1 = 'act_tot'
                y_string2 = 'act_par_yz'
            if axis_str[j] == 'Y':
                x_string = 'imp_y_los'
                y_string1 = 'act_tot'
                y_string2 = 'act_par_xz'
            if axis_str[j] == 'Z':
                x_string = 'imp_z_los'
                y_string1 = 'act_tot'
                y_string2 = 'act_par_xy'
            
            logger.debug("x_string set as: %s", str(x_string))
            logger.debug("y_string1 set as: %s", str(y_string1))
            logger.debug("y_string2 set as: %s", str(y_string2))
            # Setting the proper title string and axis los string
            title_string = tit_str[i]
            axis_string = axis_str[j]

            #Calling Function Here for Plotting - detects which function to call
            if tit_str[i] == 'Regular':
                logger.debug("tit_str detected as: %s", str(tit_str[i]))
                logger.debug("Invoking jCompPlotter function.")
                fig = jCompPlotter(data_sifted[x_string],
                                     data_sifted[y_string1],
                                     data_sifted[y_string2],
                                     axis_string,
                                     title_string)
                plt.tight_layout()
                save_fig_dir = save_dir + '/' + axis_string + '_LOS_' + title_string + '.pdf'
                logger.debug("Save Directory set to: %s", str(save_fig_dir))
                fig.savefig(save_fig_dir, bbox_inches='tight')
                plt.close(fig)
            else:
                logger.debug("tit_str detected as: %s", str(tit_str[i]))
                logger.debug("Invoking jCompPlotterColormap function.")
                fig = jCompPlotterColormap(data_sifted[x_string],
                                              data_sifted[y_string1],
                                              data_sifted[y_string2],
                                              data_sifted['mass'], #CALLS MASS HERE
                                              axis_string,
                                              title_string)
                plt.tight_layout()
                save_fig_dir = save_dir + '/' + axis_string + '_LOS_' + title_string + '.pdf'
                logger.debug("Save Directory set to: %s", str(save_fig_dir))
                fig.savefig(save_fig_dir, bbox_inches='tight')
                plt.close(fig)
    
    logger.debug("jComparisonPlotter has been run successfully.")
    return()
    

def jTimestepPlotter(dict_list,
                     simulation_list,
                     save_dir,
                     x_data_key,
                     y_data_key,
                     axis_str,
                     tit_str,
                     fitting_overlay): 
    """
    Plots the Timestep Comparison Plots for each timestep in all of the data
    acquired by the analyzer function.
    
    Inputs:
        - dict_list: list of dictionaries
            Has all the data sorted by simulation run for each timestep (one at
            a time however)
        - simulation_list: list of strings
            Simulation strings corresponding directly to the dict_list items
        - save_dir: string
            Save directory for plot
        - x_data_key: string
            Used to call the proper dataset in the dict_list
        - y_data_key: string
            Used to call the proper dataset in the dict_list
        - axis_str: string
            Used to set the LOS value: Options: 'X', 'Y', 'Z'.
        - tit_str: Either 'Full' or 'Partial' - Sets the title and datasets used
        - fitting_overlay: boolean
            Run the fits on top of the plot or not.
    
    Outputs:
        - fig: matplotlib object
            - Used to step back in the directory tree and save later.
    """
    logger = logging.getLogger("initialize.definitions.jTimestepPlotter")
    
    logger.debug("Setting Axis Strings.")
    x_str = 'Implied Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    y_str = 'Actual Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    
    if tit_str == 'Full':
        logger.debug("tit_str has been detected as: %s", str(tit_str))
        title = 'j Comparison (x+y+z components) '+axis_str+ ' LOS'
    else:
        logger.debug("tit_str has been detected as: %s", str(tit_str))
        if axis_str == 'X':
            title = 'j Comparison (y+z components) '+axis_str+ ' LOS'
        if axis_str == 'Y':
            title = 'j Comparison (x+z components) '+axis_str+ ' LOS'
        if axis_str == 'Z':
            title = 'j Comparison (x+y components) '+axis_str+ ' LOS'
    
    logger.debug("title has been set to: %s", str(title))
    #Starting Plotting Here
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #Setting Empty Lists for unity line measurements
    x_min = []
    x_max = []
    y_min = []
    y_max = []
    
    #Setting Combination Tuples for marker and line style calling
    logger.debug("Setting tuples for markers.")
    tuples_max = 12
    tuples_max_list = np.arange(0,tuples_max,1)
    marker_color = itertools.cycle(['r','k','b'])
    marker_style = itertools.cycle(['>','.','^','*','>'])
    marker_tuples = zip(tuples_max_list,marker_style,marker_color)
    marker_tuples_list = list(marker_tuples)
    logger.debug("Marker List Tuple has been set to: %s", str(marker_tuples_list))

    logger.debug("Setting tuples for best fit lines.")
    line_color = itertools.cycle(['r','k','b'])
    line_style = itertools.cycle(['-','--','-.',':'])
    line_tuples = zip(tuples_max_list,line_style,line_color)
    line_tuples_list = list(line_tuples)
    logger.debug("Best Fit Lines List Tuple has been set to: %s", str(line_tuples_list))

    #Here will start the Plotting of the data
# =============================================================================       
    for i in range(0,len(dict_list)):
        logger.debug("Working on Entry: %s in dict_list", str(i))
        x = dict_list[i][x_data_key]
        y = dict_list[i][y_data_key]
        label_str = simulation_list[i]
        ax.loglog(x,y,
                  marker = marker_tuples_list[i][1],
                  markerfacecolor = marker_tuples_list[i][2],
                  markeredgewidth = 0,
                  linestyle='',
                  label=label_str)
        x_min.append(min(x))
        x_max.append(max(x))
        y_min.append(min(y))
        y_max.append(max(y))
        
        if fitting_overlay == True:
            #Fitting Here
            log_x = np.log10(x)
            log_y = np.log10(y)
            slope, intercept = np.polyfit(log_x, log_y, 1)
            fit_y = 1e1**(slope * log_x + intercept) #Creating Fit Line
            label_fit = 'Line of Best Fit: ' + simulation_list[i]
            ax.plot(x,fit_y,
                    linestyle = line_tuples_list[i][1],
                    color = line_tuples_list[i][2],
                    alpha=0.5,
                    label = label_fit)
        logger.debug("Data for entry: %s has been plotted.", str(i))

# =============================================================================
    logger.debug("All Data inputted on graph.")
    #Doing Unity Line after to avoid overwrite by function above
    #Grabs Min and Max for x-axis
    x_min = min(x_min)
    x_max = max(x_max)
    unity_x = np.linspace(x_min,x_max,1000) #Makes line of unity values
    logger.debug("Setting Line of Constant Slope.")
    ax.loglog(unity_x,unity_x, 'k--', label='Line of Unity')
        
    #Axes Specifications
    logger.debug("Setting Plot Specifications.")
    ax.set_aspect('equal')
    #ax.set_facecolor('#f2f2f2')
    ax.grid()
    
    #Labelling and Legend
    logger.debug("Setting Labels.")
    ax.legend(bbox_to_anchor=(1.02, 1)) #Sets Legend off to the side of the plot
    ax.set_xlabel(x_str)
    ax.set_ylabel(y_str)
    ax.set_title(title)
    
    logger.debug("jTimestepPlotter has been run successfully.")
    return(fig)

def TimestepPlotter(flist, save_dir, timestamp,fitting_overlay):
    """
    Master Function for starting the Timestep Plotting
    
    Input:
        - flist: list of strings
            Corresponds to a list of fits files that share the same timestamp
        - save_dir: string
            Top level save directory for this step
        - timestamp: string
            Used to keep track of which timestamp we are working on, used for 
            filename for saving also.
    
    Output:
        - No Outputs, self contained.
    """
    logger = logging.getLogger("initialize.definitions.TimestepPlotter")

    #Set empty lists for append statements
    dict_list = []
    simulation_list = []
    logger.debug("Starting Loop over all files inputted.")
    for i in range(0,len(flist)):
        current_file = flist[i]
        logger.debug("Currently working on file: %s", str(current_file))
        hdu = fits.open(current_file)
        timestep_stamp = current_file.split("/")[-3]
        #Grabbing BinHDUTable object here, should always be second object
        hdu_table = hdu[1]
        data = hdu_table.data #Grabs data stored in the table -> FITS REC   
        logger.debug("Invoking DataGrabber function.")
        d = DataGrabber(data)
        hdu.close()
        #Append the sifted dictionaries into the dict_list for computation.
        logger.debug("Invoking DictionarySifter function.")
        dict_list.append(DictionarySifter(d))
        simulation_list.append(timestep_stamp)
    #Call Function Here
    logger.debug("Creating tit_str and axis_str to loop over.")
    tit_str = ['Full','Partial']    #Loop over these two values
    axis_str = ['X','Y','Z']        #Loop over these LOS axes
    for i in range(0,len(tit_str)): #Loop over full and partial data sets for actual sim
        logger.debug("Currently working on tit_str: %s", str(tit_str[i]))
        logger.debug("Looping over Axis String now.")
        for j in range(0,len(axis_str)): #Loop over all the LOS axes
            if tit_str == 'Full':
                #Setting actual data to be the full amount
                y_string = 'act_tot'
            else:
                #Looping for partial data sets for actual data
                if axis_str[j] == 'X':
                    y_string = 'act_par_yz'
                if axis_str[j] == 'Y':
                    y_string = 'act_par_xz'
                if axis_str[j] == 'Z':
                    y_string = 'act_par_xy'
            logger.debug("y_string has been set to: %s", str(y_string))
            #Looping to set the correct data for the los axis
            if axis_str[j] == 'X':
                x_string = 'imp_x_los'
            if axis_str[j] == 'Y':
                x_string = 'imp_y_los'
            if axis_str[j] == 'Z':
                x_string = 'imp_z_los'
            
            logger.debug("x_string has been set to: %s", str(x_string))
            # Setting the proper title string and axis los string
            title_string = tit_str[i]
            axis_string = axis_str[j]
            #Call Main Function Here with all included inputs.
            logger.debug("Invoking jTimestepPlotter function.")
            jTimestepPlotter(dict_list,
                             simulation_list,
                             save_dir,
                             x_string,
                             y_string,
                             axis_string,
                             title_string,
                             fitting_overlay)
            #Make the filename for the pdf saved
            save_string = save_dir + timestamp + '_'+title_string+'_j_comp_'+axis_string+'_LOS.pdf'
            logger.debug("Save Directory set to: %s", str(save_string))
            plt.savefig(save_string, bbox_inches='tight')
            logger.debug("Plot Saved.")
            plt.close()
    logger.debug("TimestepPlotter has been run successfully.")
    return()
