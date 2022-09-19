# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

@author: Thomas Chartier
"""
import math
from openquake.hazardlib.geo import Line

from openquake.hazardlib.scalerel.wc1994 import WC1994
from openquake.hazardlib.source.kite_fault import KiteFaultSource
from openquake.hazardlib.mfd import TruncatedGRMFD

from geometry_tools import *



def get_profiles_from_trace(lons,lats,dip,usd,lsd,oriented):
    ColLon = lons
    ColLat = lats

    dip_tan = math.tan(math.radians(dip))
    hdist_top = usd / dip_tan
    hdist_bottom = lsd / dip_tan
    """
    # Similar to :meth:`from_fault_data`, we just don't resample edges
    if vertical_faults == False :
        dip_tan = math.tan(math.radians(dip))
        hdist_top = usd / dip_tan
        hdist_bottom = lsd / dip_tan
    else :
        hdist_top = 0.
        hdist_bottom = 0."""

    # orienting the arrays in order to respect OQ right hand rule
    compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))

    if str('N') in str(oriented):
        if compass_bearing < 180. :
            ColLon = reversed(ColLon)
            ColLat = reversed(ColLat)
    elif str('S') in str(oriented):
        if compass_bearing > 180. :
            ColLon = reversed(ColLon)
            ColLat = reversed(ColLat)
    elif str('E') in str(oriented):
        if compass_bearing > 90. and compass_bearing < 270. :
            ColLon = reversed(ColLon)
            ColLat = reversed(ColLat)
    elif str('W') in str(oriented):
        if compass_bearing < 90. or compass_bearing > 270. :
            ColLon = reversed(ColLon)
            ColLat = reversed(ColLat)
    ColLon = list(ColLon)
    ColLat = list(ColLat)


    strike = compass_bearing
    mean_azimuth = (strike + 90.0) % 360

    prof_chek = []

    profiles = []
    # loop on points
    i_pt = 0
    for x,y in zip(ColLon,ColLat):
        xtop, ytop = point_at(x, y, mean_azimuth, hdist_top)
        xbot, ybot = point_at(x, y, mean_azimuth, hdist_bottom)
        profile_i = Line.from_vectors([xtop,xbot], [ytop,ybot], deps=[usd,lsd])
        #profile_i = Line([Point(xtop, ytop, usd), Point(xbot, ybot, lsd)])
        if not str(profile_i.points) in prof_chek:
            profiles.append(profile_i)
            prof_chek.append(str(profile_i.points))

    return profiles


def create_kite_source(profiles,msr,rupture_mesh_spacing,rake,rupture_aspect_ratio):
    WC1994()

    #profiles = [Line([Point(0.0, 0.0, 0.0), Point(0.0, 0.01, 15.0)]),
    #                Line([Point(0.3, 0.0, 0.0), Point(0.3, 0.01, 15.0)])]
    mfd = TruncatedGRMFD(a_val=0.5, b_val=1.0, min_mag=6.2, max_mag=6.4,
                                 bin_width=0.1)

    # Create the source instance
    source = KiteFaultSource(None, None, None, mfd, rupture_mesh_spacing,
                          msr,
                          rupture_aspect_ratio, None, profiles, rake,
                          floating_x_step=rupture_mesh_spacing,
                          floating_y_step=rupture_mesh_spacing)

    return source
