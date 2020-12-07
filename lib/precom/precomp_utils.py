# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
import matplotlib.pyplot as plt


def find_bounding_box(faults):
    maxmin_pt_lon, maxmin_pt_lat = [], []
    for fi in range(nb_faults):
        maxmin_pt_lon.append([np.min([i[0] for i in faults[fi]['geometry']["coordinates"]]),
                             np.max([i[0] for i in faults[fi]['geometry']["coordinates"]])])
        maxmin_pt_lat.append([np.min([i[1] for i in faults[fi]['geometry']["coordinates"]]),
                             np.max([i[1] for i in faults[fi]['geometry']["coordinates"]])])
    return maxmin_pt_lon, maxmin_pt_lat
    
def find_possible_asso(maxmin_pt_lon,maxmin_pt_lat,plot_fig=False):
    d = 0.5
    assso_fault = []
    for lon_i,lat_i in zip(maxmin_pt_lon,maxmin_pt_lat):
        assso_fault_i = []
        j_fault = 0
        for lon_j,lat_j in zip(maxmin_pt_lon,maxmin_pt_lat):
            if lon_j[0] > lon_i[0]-d and lon_j[0] < lon_i[1]+d :
                if lat_j[0] > lat_i[0]-d and lat_j[0] < lat_i[1]+d :
                    assso_fault_i.append(j_fault)
                if lat_j[1] > lat_i[0]-d and lat_j[1] < lat_i[1]+d :
                    assso_fault_i.append(j_fault)
            if lon_j[1] > lon_i[0]-d and lon_j[1] < lon_i[1]+d :
                if lat_j[0] > lat_i[0]-d and lat_j[0] < lat_i[1]+d :
                    assso_fault_i.append(j_fault)
                if lat_j[1] > lat_i[0]-d and lat_j[1] < lat_i[1]+d :
                    assso_fault_i.append(j_fault)
            j_fault+=1
        assso_fault_i = list(set(assso_fault_i))
        assso_fault.append(assso_fault_i)
        
    if plot_fig == True :
        x =[]
        for i in assso_fault:
            x.append(len(i)-1-0.5)
        plt.hist(x)
        plt.xlabel("number of close faults to be considered for rupture jump")
        plt.ylabel("number of faults in this situation")
        plt.show()
    
    return assso_fault
