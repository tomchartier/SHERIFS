#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import numpy as np

def calculation_of_percentile(values,weights):
    weighted_data = []
    for i in range(len(values)):
        weighted_data.append(values[i]*weights[i]/np.sum(weights))
    mean_value = np.sum(weighted_data)
    
    #calcultion of the weighted percentiles
    sorted_inds = np.array(values).argsort()
    
    data_weight_sorted = []
    data_weight_sorted_sum = []
    data_value_sorted = []
    
    per_16 = 0
    per_50 = 0
    per_84 = 0
    
    i_16 = 0
    i_50 = 0
    i_84 = 0        
    
    for i in sorted_inds :
        data_weight_sorted.append(weights[i]/np.sum(weights))
        data_value_sorted.append(values[i])
        
        if (np.sum(data_weight_sorted) > 0.16) & (i_16 == 0):
            i_16 = 1
            per_16 = values[i]
            
        if (np.sum(data_weight_sorted) > 0.50) & (i_50 == 0):
            i_50 = 1
            per_50 = values[i]
            
        if (np.sum(data_weight_sorted) > 0.84) & (i_84 == 0):
            i_84 = 1
            per_84 = values[i]
        
        data_weight_sorted_sum.append(np.sum(data_weight_sorted))
    
    #self.per_16, self.per_50, self.mean_value, self.per_84, self.data_value_sorted, self.data_weight_sorted_sum = per_16, per_50, mean_value, per_84, data_value_sorted, data_weight_sorted_sum
    return per_16, per_50, mean_value, per_84, data_value_sorted, data_weight_sorted_sum
