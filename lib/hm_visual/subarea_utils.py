#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

utils functions for dealing with subareas

@author: thomas chartier
"""


import numpy as np
import os
import geojson
import matplotlib.path as mplPath


def get_geom(sub_area_file,model):
    '''
    Extracts the names and geometries of the available sub-areas
    
    param :
    sub_area_file : str, path to the subarea file
    
    model : str, model for with the sub-area applies
    
    returns :
    
    sub_area_names : list, list of the name of the differente subareas
    
    sub_area_lon : list, list of lists of the longitudes of the subareas
    
    sub_area_lat : list, list of lists of the latitudes of the subareas
    
    '''
    sub_area_names = []
    sub_area_coord = []
    sub_area_lon = []
    sub_area_lat = []
    if not '.geojson' in sub_area_file:
        read_sub_area_file = open(sub_area_file,'rU')
        lines_sub_area = read_sub_area_file.readlines()
        for line in lines_sub_area:
            model_sub_area = line.split('\t')[0]
            if model == model_sub_area:
                if print_detail == True :
                    print("rates in sub area :",line.split('\t')[1])
                sub_area_names.append(line.split('\t')[1])
                sub_area_coord.append(line.split('\t')[2:])
                sub_area_lon_i = []
                sub_area_lat_i = []
                for sub_area_coord_i in line.split('\t')[2:]:
                    if not '\n' in sub_area_coord_i.split(','):
                        if not '' in sub_area_coord_i.split(','):
                            sub_area_lon_i.append(float(sub_area_coord_i.split(',')[0]))
                            sub_area_lat_i.append(float(sub_area_coord_i.split(',')[1]))
                sub_area_lon.append(sub_area_lon_i)
                sub_area_lat.append(sub_area_lat_i)
    else :
        with open(sub_area_file) as f:
            gj = geojson.load(f)
        areas = gj["features"]
        for area_i in areas:
            sub_area_names.append('sub_area_' + str(area_i["properties"]["Id"]))
            sub_area_lon_i = []
            sub_area_lat_i = []
            for pt in area_i["geometry"]["coordinates"][0][0]:
                sub_area_lon_i.append(pt[0])
                sub_area_lat_i.append(pt[1])
            sub_area_lon.append(sub_area_lon_i)
            sub_area_lat.append(sub_area_lat_i)
    return sub_area_names,sub_area_lon,sub_area_lat

def find_faults_in_sub(fault_names,Lon,Lat,sub_area_lon,sub_area_lat):
    '''
    Finds the faults that are inside the subare

    param :
    fault_names : list, list of all the fault names in the model

    Lon : list, list of lists of the longitudes of the faults

    Lat : list, list of lists of the latitudes of the faults
    
    sub_area_lon : list, list of longitudes of the subarea
    
    sub_area_lat : list, list of latitudes of the subarea

    returns :

    faults_in_sub_area : list, list of the name of the faults in the subarea

    '''
    Poly_sub = []
    for x1,y1 in zip(sub_area_lon,sub_area_lat): # creation du polygon de la zone
        Poly_sub.append((x1,y1))
    bbPath_sub_area_i = mplPath.Path(Poly_sub)

    faults_in_sub_area = []
    index_fault = 0
    for fault_name in fault_names:
        nb_point_in_sub_area = 0
        for lon_i,lat_i in zip(Lon[index_fault],Lat[index_fault]):
            if bbPath_sub_area_i.contains_point((lon_i,lat_i)) == 1: #test to know if the point is in the sub_area
                nb_point_in_sub_area += 1
        if nb_point_in_sub_area > len(Lon[index_fault])/2.: #if more than half the points of the trace are in the sub area
            faults_in_sub_area.append(fault_name)  #the fault is in the sub area
        index_fault +=1
    return faults_in_sub_area
