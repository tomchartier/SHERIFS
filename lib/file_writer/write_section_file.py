# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

@author: Thomas Chartier
"""
import numpy as np

from geometry_tools import *
import math
from decimal import Decimal, getcontext
getcontext().prec = 10


def build(f,txt):
    '''
    f : path to file name
    txt : str containing the file info
    '''
    f = open(f,'w')
    f.write(txt)
    f.close()

def start():
    '''
    txt : str containing the file info
    '''
    txt = ''
    # Initiate the xml file
    txt += '<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
    txt += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
    txt += '\txmlns="http://openquake.org/xmlns/nrml/0.5">\n'
    txt += '\t<faultSectionCollection name="'+model_name+'" id="fs1">"'

    return txt

def wrt_kite_geo(txt,fault_name,faults_names,faults_data,do_resample=False):
    '''
    txt : str containing the file info
    faults_names : list of the fault names in order (serves for indexing)
    faults_data : dict containing all faults info
    do_resample : when resample=ing the coordinates is necessary to lower the
    number of points.
    '''

#   for Fault_name in faults_in_scenario :
    index_fault = faults_names.index(fault_name)
    ColLon = faults_data[index_fault]['lon']
    ColLat = faults_data[index_fault]['lat']
    Depth = faults_data[index_fault]['depth']


    scenario_mechanism.append(faults_data[index_fault]['rake'])

    if Depth and all(elem == 'sf' for elem in Depth):
        type_of_fault = 'sf'
    else :
        type_of_fault = 'cf'
        Depth = [float(i) for i in Depth]
    if type_of_fault == 'sf':
        usd = faults_data[index_fault]['upper_sismo_depth']
        lsd = faults_data[index_fault]['lower_sismo_depth']
        # Similar to :meth:`from_fault_data`, we just don't resample edges
        dip_tan = math.tan(math.radians(faults_data[index_fault]['dip']))
        hdist_top = usd / dip_tan
        hdist_bottom = lsd / dip_tan

        # orienting the arrays in order to respect OQ right hand rule
        compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))

        if str('N') in str(faults_data[index_fault]['oriented']):
            if compass_bearing < 180. :
                ColLon = reversed(ColLon)
                ColLat = reversed(ColLat)
        elif str('S') in str(faults_data[index_fault]['oriented']):
            if compass_bearing > 180. :
                ColLon = reversed(ColLon)
                ColLat = reversed(ColLat)
        elif str('E') in str(faults_data[index_fault]['oriented']):
            if compass_bearing > 90. and compass_bearing < 270. :
                ColLon = reversed(ColLon)
                ColLat = reversed(ColLat)
        elif str('W') in str(faults_data[index_fault]['oriented']):
            if compass_bearing < 90. or compass_bearing > 270. :
                ColLon = reversed(ColLon)
                ColLat = reversed(ColLat)
        ColLon = list(ColLon)
        ColLat = list(ColLat)

        # does a resampling to reduce the number of points
        do_resample = True
        if do_resample == True :
            # parameters
            # min distance between two points
            min_d = 5.
            # max distance between two points
            max_d = 15.
            # change in azimuth to consider
            az_d = 15.

            resampled_ColLon, resampled_ColLat = [ColLon[0]], [ColLat[0]]
            for i_pt in range(len(ColLon)-2):
                add_point = False
                i_pt += 1
                dist_before = distance(resampled_ColLon[-1], resampled_ColLat[-1],
                 ColLon[i_pt], ColLat[i_pt])
                dist_after = distance(resampled_ColLon[-1], resampled_ColLat[-1],
                 ColLon[i_pt+1], ColLat[i_pt+1])

                az_before = calculate_initial_compass_bearing((resampled_ColLon[-1], resampled_ColLat[-1]),
                (ColLon[i_pt], ColLat[i_pt]))
                az_after = calculate_initial_compass_bearing((ColLon[i_pt], ColLat[i_pt]),
                (ColLon[i_pt+1], ColLat[i_pt+1]))

                if dist_before > min_d and (dist_after > max_d or abs(az_before - az_after)>az_d):
                    add_point = True
                    #special rule for the antepenultimate point
                    if i_pt == len(ColLon)-2:
                        last_dist_after = distance(ColLon[i_pt], ColLat[i_pt],
                                            ColLon[-1], ColLat[-1])
                        if last_dist_after < min_d :
                            add_point = False
                if add_point == True:
                    #add the point
                    resampled_ColLon.append(ColLon[i_pt])
                    resampled_ColLat.append(ColLat[i_pt])

            resampled_ColLon.append(ColLon[-1])
            resampled_ColLat.append(ColLat[-1])
            ColLon, ColLat = resampled_ColLon, resampled_ColLat

            plot_stuf_detail = True
            if plot_stuf_detail == True:
                min_dist_tmp = min_d
                for i_pt in range(len(ColLon)-2):
                    dist_tmp = distance(resampled_ColLon[i_pt+1], resampled_ColLat[i_pt+1],
                                            resampled_ColLon[i_pt], resampled_ColLat[i_pt])
                    if dist_tmp < min_dist_tmp :
                        min_dist_tmp = dist_tmp
                if min_dist_tmp < min_d/2.:
                    print("min dist :",round(min_dist_tmp),"id: ",faults_data[index_fault]['name'])


        str_geom+='\t\t\t<kiteSurface>\n'

        #mean azimuth of the section
        compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
        strike = compass_bearing
        mean_azimuth = (strike + 90.0) % 360

        # loop on profiles
        i_pt = 0
        for x,y in zip(ColLon,ColLat):
            # local azimuth
            if [x,y] == [ColLon[0],ColLat[0]]:
                compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[1],ColLon[1]))
                strike = compass_bearing
                azimuth = (strike + 90.0) % 360
            elif [x,y] == [ColLon[-1],ColLat[-1]]:
                compass_bearing = calculate_initial_compass_bearing((ColLat[-2],ColLon[-2]),(ColLat[-1],ColLon[-1]))
                strike = compass_bearing
                azimuth = (strike + 90.0) % 360
            else :
                compass_bearing = calculate_initial_compass_bearing((ColLat[i_pt-1],ColLon[i_pt-1]),
                (ColLat[i_pt+1],ColLon[i_pt+1]))
                strike = compass_bearing
                azimuth = (strike + 90.0) % 360

            azimuth = ((mean_azimuth+azimuth)/2.) % 360

            str_geom+='\t\t\t<profile>\n'
            str_geom+='\t\t\t\t<gml:LineString>\n'
            str_geom+='\t\t\t\t\t<gml:posList>\n'
            xt, yt = point_at(x, y, azimuth, hdist_top)
            str_geom+='\t\t\t\t\t\t'+str(xt)+' '+str(yt)+' '+str(usd)+' '
            xb, yb = point_at(x, y, azimuth, hdist_bottom)
            str_geom+=str(xb)+' '+str(yb)+' '+str(lsd)+'\n'
            str_geom+='\t\t\t\t\t</gml:posList>\n'
            str_geom+='\t\t\t\t</gml:LineString>\n'
            str_geom+='\t\t\t</profile>\n'
            i_pt +=1

        str_geom+='\t\t\t</kiteSurface>\n'

    if type_of_fault == 'cf':
        str_geom+='\t\t\t<kiteSurface>\n'

        for depth_i in sorted(set(Depth)):
            # write the upper profile
            str_geom+='\t\t\t<profile>\n'
            str_geom+='\t\t\t\t<gml:LineString>\n'
            str_geom+='\t\t\t\t\t<gml:posList>\n'
            for x,y in zip(ColLon,ColLat):
                str_geom+='\t\t\t\t\t\t'+str(x)+' '+str(y)+' '+str(hdist_d)+'\n'
            str_geom+='\t\t\t\t\t</gml:posList>\n'
            str_geom+='\t\t\t\t</gml:LineString>\n'
            str_geom+='\t\t\t<\profile>\n'

        str_geom+='\t\t\t</kiteSurface>\n'


    return txt

def wrt_section(txt,geotype):
    '''
    txt : str containing the file info
    geotype : type of geometry of the section (kite, complex, or  simple)
    '''

    txt += '    <section name="'+section_name+'" id="'+section_id+'">\n'
    if geotype == "kite":
        out = wrt_kite_geo()
    if geotype == "complex":
        out = wrt_complex_geo()
    if geotype == "simple":
        out = wrt_simple_geo()
    txt += out
    txt += '        </section>\n'
    return txt

def end(txt):
    '''
    txt : str containing the file info
    '''
    txt += '    </faultSectionCollection>\n'
    txt += '</nrml>\n'
    return txt


    return txt
