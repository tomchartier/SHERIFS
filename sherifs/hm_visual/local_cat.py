# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
import geojson

def read_geojson(file):
    '''
    file : str, path to file with the geojson with the polygon and the cat rate
    loc_dat : dict, contains for each zone the polygon and the rates of the local cat
    '''
    with open(file) as f:
        gj = geojson.load(f)
    zones = gj['features']
    loc_cat = {}
    for zone in zones :
        if "rate_cat" in zone["properties"].keys():
            cat_rates = zone["properties"]["rate_cat"]
        else :
            cat_rates = False

        loc_cat.update({zone["properties"]["id"]:{
        "cat_rates" : cat_rates,
        "poly":zone["geometry"]["coordinates"]}})
    return loc_cat

def get_model_rate(poly,OQ_entry_faults,OQ_entry_scenarios,pts_list,bin_mag,
param,faults_data,faults_names,index_faults_in_scenario,local_zone_mfd):
    '''
    poly : polypath , polygon of the zone
    OQ_entry_faults : list, rates and info on each single fault rupture
    OQ_entry_scenarios : list, rates and info on each complex fault rupture
    pts_list : dict, info and rates of the bg points
    param : dict, parameters used in this run
    bin_mag : list, bining in magnitude
    faults_data : dict, data of the faults
    faults_names : list , list of the names of the faults
    index_faults_in_scenario : list, list of the indexes of faults for each scenario
    local_zone_mfd : list, mfd of the bg seismicity if bg was (or is) a uiform zone

    '''
    rate_bg = np.zeros_like(bin_mag)
    smooth = np.zeros_like(bin_mag)
    if param["main"]["background"]["option_bg"]=="smooth" :
        txt_no_bg = ""
        for lon_lat in pts_list.keys() :
            lon = float(lon_lat.split("_")[0])
            lat = float(lon_lat.split("_")[1])
            if poly.contains_point((lon,lat)) == 1:
                scaled_mfd = pts_list[lon_lat]["scaled_mfd"]
                mfd_smooth = pts_list[lon_lat]["mfd_smooth"]
                i_mag = 0
                for r in scaled_mfd :
                    rate_bg[i_mag]+=r
                    i_mag += 1
                i_mag = 0
                for r in mfd_smooth :
                    smooth[i_mag]+=r
                    i_mag += 1

    elif param["main"]["background"]["option_bg"]=="zone" :
        txt_no_bg = ""
        rate_bg = local_zone_mfd
    else :
        txt_no_bg = " BG NOT INCLUDED"

    rate_faults = np.zeros_like(bin_mag)
    #loop on simple faults
    for index_fault,fault_name in zip(range(len(faults_names)),faults_names):
        f_in = False
        stop_loop = False
        while stop_loop == False :
            lons = faults_data[index_fault]['lon']
            lats = faults_data[index_fault]['lat']
            for lon,lat in zip(lons,lats):
                if poly.contains_point((lon,lat)) == 1:
                    f_in = True
                    stop_loop = True
            stop_loop = True
        if f_in == True :
            i_MFD = np.where(np.array(faults_names) == fault_name)[0][0]
            MFD = OQ_entry_faults[i_MFD]
            i_mag = 0
            for r in MFD :
                rate_faults[i_mag]+=MFD[i_mag]
                i_mag += 1


    if len(index_faults_in_scenario) != 0 :
        index_scenario = 0
        for index_faults_in_sc in index_faults_in_scenario:
            faults_in_poly = 0
            for index_fault in index_faults_in_sc[0] :
                f_in = False
                stop_loop = False
                while stop_loop == False :
                    lons = faults_data[index_fault]['lon']
                    lats = faults_data[index_fault]['lat']
                    for lon,lat in zip(lons,lats):
                        if poly.contains_point((lon,lat)) == 1:
                            f_in = True
                            stop_loop = True
                    stop_loop = True
                if f_in == True :
                    faults_in_poly += 1

            if faults_in_poly != 0 :
                MFD = OQ_entry_scenarios[index_scenario]
                i_mag = 0
                for r in MFD :
                    rate_faults[i_mag]+=MFD[i_mag] * float(faults_in_poly)/float(len(index_faults_in_sc[0]))
                    i_mag += 1

            index_scenario += 1


    return txt_no_bg,rate_faults,rate_bg,smooth
