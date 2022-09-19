# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
import matplotlib.pyplot as plt

from math import pi, cos, radians , sin, asin, sqrt, atan2, degrees, acos, radians

from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.scalerel.wc1994 import WC1994

from openquake.hazardlib.source.kite_fault import get_discrete_dimensions

from kite_for_sherifs import *

def distance(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = list(map(radians, [lon1, lat1, lon2, lat2]))
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    km = 6367 * c
    return km

def fault_length(lons,lats):
    length = 0.
    for i in range(len(lons)-1):
        length += distance(lons[i], lats[i], lons[i+1], lats[i+1])
    return length

def find_bounding_box(faults):
    nb_faults = len(faults)
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

def get_fault_dimentions(lons,lats,faults,fi,param):

    if "lsd" in faults[fi]['properties'].keys() :
        lsd = faults[fi]['properties']["lsd"]
    elif "lo_s_d" in faults[fi]['properties'].keys():
        lsd = faults[fi]['properties']["lo_s_d"]

    if "usd" in faults[fi]['properties'].keys() :
        usd = faults[fi]['properties']["usd"]
    elif "up_s_d" in faults[fi]['properties'].keys():
        usd = faults[fi]['properties']["up_s_d"]

    dip = faults[fi]['properties']["dip"]

    if "dip_dir" in faults[fi]['properties'].keys() :
        oriented = faults[fi]['properties']["dip_dir"]
    elif "oriented" in faults[fi]['properties'].keys():
        oriented = faults[fi]['properties']["oriented"]

    profiles = get_profiles_from_trace(lons,lats,dip,usd,lsd,oriented)


    if "msr" in param["pre"].keys():
        msr = param["pre"]["msr"]
    else :
        msr =  "WC1994"

    if msr == "WC1994" :
        msr = WC1994()

    if "rupture_mesh_spacing" in param["pre"].keys():
        rupture_mesh_spacing = param["pre"]["rupture_mesh_spacing"]
    else :
        rupture_mesh_spacing =  3.

    rupture_aspect_ratio = param["pre"]["rupture_aspect_ratio"]

    if "rupture_aspect_ratio" in param["pre"].keys():
        rupture_aspect_ratio = param["pre"]["rupture_aspect_ratio"]
    else :
        rupture_aspect_ratio =  1.

    if "v_rl" in faults[fi]['properties'].keys():
        #using the format of the GEM GFDB
        v_rl = abs(faults[fi]['properties']["v_rl"])
        v_ex = faults[fi]['properties']["v_ex"]
        slip_rate_moy = (v_rl**2 + v_ex**2)**0.5
        err = (faults[fi]['properties']["e_rl"]**2 + faults[fi]['properties']["e_ex"]**2)**0.5
        slip_rate_min = slip_rate_moy - err
        slip_rate_max = slip_rate_moy + err
        if slip_rate_min < 0. :
            print("ERROR - slip rate min is negative !")

        if v_ex != 0. :
            if abs(v_ex) < v_rl :
                rake = degrees(acos(abs(v_ex) / v_rl))
            else :
                rake = degrees(acos(v_rl / abs(v_ex)))
            if v_ex < 0 :
                rake = -rake
        else :
            rake = 0.
    elif "sr" in faults[fi]['properties'].keys():
        #using slip-rate direclty
        slip_rate_moy = abs(faults[fi]['properties']["sr"])
        err = abs(faults[fi]['properties']["e_sr"])
        slip_rate_min = slip_rate_moy - err
        slip_rate_max = slip_rate_moy + err
        rake = abs(faults[fi]['properties']["rake"])
    elif "sr_mean" in faults[fi]['properties'].keys():
        #using slip-rate direclty
        slip_rate_moy = abs(faults[fi]['properties']["sr_mean"])
        slip_rate_min = abs(faults[fi]['properties']["sr_min"])
        slip_rate_max = abs(faults[fi]['properties']["sr_max"])
        rake = abs(faults[fi]['properties']["rake"])

    mfd = TruncatedGRMFD(a_val=0.5, b_val=1.0, min_mag=2., max_mag=3.,bin_width=0.1)

    source = create_kite_source(profiles,msr,rupture_mesh_spacing,rake,rupture_aspect_ratio)

    area = source.surface.get_area()

    length, width = get_discrete_dimensions(area, rupture_mesh_spacing, rupture_aspect_ratio,rupture_mesh_spacing)
    '''
    length = 0.
    for i in range(len(lons)-1):
        length += distance(lons[i],lats[i],lons[i+1],lats[i+1])

    try :
        width =((faults[fi]['properties']['lsd']-
                  faults[fi]['properties']['usd'])/
                  sin(radians(faults[fi]['properties']['dip'])))
    except :
        width =((faults[fi]['properties']['lo_s_d']-
                  faults[fi]['properties']['up_s_d'])/
                  sin(radians(faults[fi]['properties']['dip'])))
    '''

    return area, length, width

def calc_f_dims(faults,param,plt_fig=False):


    nb_faults = len(faults)
    f_lengths = []
    f_areas = []
    for fi in range(nb_faults):
        # create kite source
        lons = [i[0] for i in faults[fi]['geometry']["coordinates"]]
        lats = [i[1] for i in faults[fi]['geometry']["coordinates"]]


        area, length, width = get_fault_dimentions(lons,lats,faults,fi,param)
        
        # print(area,length*width)

        #length_i = fault_length(lons_i,lats_i)
        #f_lengths.append(length_i)
        #try :
        #    width_i =((faults[fi]['properties']['lsd']-
        #              faults[fi]['properties']['usd'])/
        #              sin(radians(faults[fi]['properties']['dip'])))
        #except :
        #    width_i =((faults[fi]['properties']['lo_s_d']-
        #              faults[fi]['properties']['up_s_d'])/
        #              sin(radians(faults[fi]['properties']['dip'])))


        #f_areas.append(length_i*width_i)

        f_areas.append(area)
        f_lengths.append(length)

    if plt_fig ==True:
        plt.hist(f_lengths)
        plt.xlabel("Lengths (km)")
        plt.ylabel("Nb faults")
        plt.show()

    print("In total, there are ",round(sum(f_lengths))," km of faults in the model.")
    return f_lengths, f_areas



def calculate_initial_compass_bearing(pointA, pointB):
    """
    Calculates the bearing between two points.

    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))

    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees

    :Returns:
      The bearing in degrees

    :Returns Type:
      float
    """
#     if (type(pointA) != tuple) or (type(pointB) != tuple):
#         raise TypeError("Only tuples are supported as arguments")

    lat1 = radians(pointA[0])
    lat2 = radians(pointB[0])

    diffLong = radians(pointB[1] - pointA[1])

    x = sin(diffLong) * cos(lat2)
    y = cos(lat1) * sin(lat2) - (sin(lat1)
            * cos(lat2) * cos(diffLong))

    initial_bearing = atan2(x, y)

    # Now we have the initial bearing but atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def wc1994_median_mag( area, rake):
    """
    Return magnitude (Mw) given the area and rake.

    Setting the rake to ``None`` causes their "All" rupture-types
    to be applied.

    :param area:
        Area in square km.
    :param rake:
        Rake angle (the rupture propagation direction) in degrees,
        from -180 to 180.
    """
    assert rake is None or -180 <= rake <= 180
    if rake is None:
        # their "All" case
        return 4.07 + 0.98 * np.log10(area)
    elif (-45 <= rake <= 45) or (rake > 135) or (rake < -135):
        # strike slip
        return 3.98 + 1.02 * np.log10(area)
    elif rake > 0:
        # thrust/reverse
        return 4.33 + 0.90 * np.log10(area)
    else:
        # normal
        return 3.93 + 1.02 * np.log10(area)


def mag_to_M0(mag):
    #returns Mo for a given mag
    M0 = 10. ** (1.5 * mag + 9.1)
    return M0
