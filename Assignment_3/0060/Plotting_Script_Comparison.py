#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 23:04:27 2018

@author: sfielder
"""
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt


#Importing Data For Actual Angular Momentum
actual_j_0 = genfromtxt('Actual_J_specific_0.csv', delimiter='')
actual_j_0 = actual_j_0[np.logical_not(np.isnan(actual_j_0))]
actual_j_1 = genfromtxt('Actual_J_specific_1.csv', delimiter='')
actual_j_1 = actual_j_1[np.logical_not(np.isnan(actual_j_1))]
actual_j_2 = genfromtxt('Actual_J_specific_2.csv', delimiter='')
actual_j_2 = actual_j_2[np.logical_not(np.isnan(actual_j_2))]
actual_j_3 = genfromtxt('Actual_J_specific_3.csv', delimiter='')
actual_j_3 = actual_j_3[np.logical_not(np.isnan(actual_j_3))]
actual_j_4 = genfromtxt('Actual_J_specific_4.csv', delimiter='')
actual_j_4 = actual_j_4[np.logical_not(np.isnan(actual_j_4))]
#Importing Data for Implied Angular Momentum
implied_j_0 = genfromtxt('Implied_J_specific_0.csv', delimiter='')
implied_j_0 = implied_j_0[np.logical_not(np.isnan(implied_j_0))]
implied_j_1 = genfromtxt('Implied_J_specific_1.csv', delimiter='')
implied_j_1 = implied_j_1[np.logical_not(np.isnan(implied_j_1))]
implied_j_2 = genfromtxt('Implied_J_specific_2.csv', delimiter='')
implied_j_2 = implied_j_2[np.logical_not(np.isnan(implied_j_2))]
implied_j_3 = genfromtxt('Implied_J_specific_3.csv', delimiter='')
implied_j_3 = implied_j_3[np.logical_not(np.isnan(implied_j_3))]
implied_j_4 = genfromtxt('Implied_J_specific_4.csv', delimiter='')
implied_j_4 = implied_j_4[np.logical_not(np.isnan(implied_j_4))]


#Plotting Scale vs j graphs

scale = [1,2,3,4,5]

plt.figure(1)
plt.plot(np.full((np.size(actual_j_0),1),scale[4]),actual_j_0,'b.',label='0.625pc')
plt.plot(np.full((np.size(actual_j_1),1),scale[3]),actual_j_1,'g.',label='1.25pc')
plt.plot(np.full((np.size(actual_j_2),1),scale[2]),actual_j_2,'m.',label='2.5pc')
plt.plot(np.full((np.size(actual_j_3),1),scale[1]),actual_j_3,'b.',label='5pc')
plt.plot(np.full((np.size(actual_j_4),1),scale[0]),actual_j_4,'c.',label='10pc')
plt.yscale('log')
plt.xticks(np.array([1,2,3,4,5]), ('0.625', '1.25', '2.5', '5', '10') )
plt.xlabel('pc (Scale Factor, length of edge of subregion)')
plt.ylabel(r'Actual Specific Angular Momentum $(pc^2 \ Myr^{-1})$')
plt.title('Specific Actual Angular Momentum vs Scale Factor')
plt.legend(bbox_to_anchor=(1.25, 1.0))
plt.savefig("j_vs_scale_factor_actual.pdf", bbox_inches='tight')


plt.figure(2)
plt.plot(np.full((np.size(implied_j_0),1),scale[4]),implied_j_0,'b.',label='0.625pc')
plt.plot(np.full((np.size(implied_j_1),1),scale[3]),implied_j_1,'g.',label='1.25pc')
plt.plot(np.full((np.size(implied_j_2),1),scale[2]),implied_j_2,'m.',label='2.5pc')
plt.plot(np.full((np.size(implied_j_3),1),scale[1]),implied_j_3,'b.',label='5pc')
plt.plot(np.full((np.size(implied_j_4),1),scale[0]),implied_j_4,'c.',label='10pc')
plt.yscale('log')
plt.xticks(np.array([1,2,3,4,5]), ('0.625', '1.25', '2.5', '5', '10') )
plt.xlabel('pc (Scale Factor, length of edge of subregion)')
plt.ylabel(r'Implied Specific Angular Momentum $(pc^2 \ Myr^{-1})$')
plt.title('Specific Implied Angular Momentum vs Scale Factor')
plt.legend(bbox_to_anchor=(1.25, 1.0))
plt.savefig("j_vs_scale_factor_implied.pdf", bbox_inches='tight')

#Plotting Comparison Graphs with line of unity for all data



