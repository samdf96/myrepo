#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 15:40:39 2018

@author: sfielder
"""

"""
Fiducial Plotter Function

Input:
    - data_dir: string
        - Will be the otuput directory of the data.
    - data_check_list: list of strings
        - Will be the timestamps that we are looking for, are going to be the 
        same ones as the data is sorted through.
"""

#Setting Local inputs for now:
data_dir = '/Users/sfielder/Documents/Astro_Data/'
data_check_list = ['0060','0070','0080','0090','0100']

import glob
from astropy.io import fits
import matplotlib.pylab as plt

#From Definitions Here:
from definitions import data_grabber


#Find all data sets with one time stamp here - will loop
for i in range(0,1): #len(data_check_list)
    flist = glob.glob(data_dir+'**/*'+data_check_list[i]+'.fits', recursive=True)
    flist.sort()
    
    #For X-LOS
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for j in range(0,len(flist)):
        hdu = fits.open(flist[j])
        hdu_table = hdu[1]
        data = hdu_table.data
        d = data_grabber(data)
        
        ax.scatter(x,
                   y1,
                   c=np.log10(mass),
                   s=100,
                   marker='.',
                   cmap='viridis',
                   linewidth=0,
                   alpha=1,
                   label='Full'))
        
        




