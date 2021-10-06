# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

@author: Thomas Chartier
"""
import math
from sherifs.utils.geometry_tools import *


def build(f, txt):
    '''
    f : path to file name
    txt : str containing the file info
    '''
    f = open(f, 'w')
    f.write(txt)
    f.close()


def start(model_name):
    '''
    txt : str containing the file info
    '''
    txt = ''
    # Initiate the xml file
    txt += '<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
    txt += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
    txt += '\txmlns="http://openquake.org/xmlns/nrml/0.5">\n'
    txt += '\t<geometryModel name="fault_sections">\n'

    return txt


def wrt_kite_geo(fault_name,faults_names,faults_data,resample,vertical_faults):
    '''
    txt : str containing the file info
    faults_names : list of the fault names in order (serves for indexing)
    faults_data : dict containing all faults info
    resample : list with a boolean and the resample parameters if needed
                when resampling the coordinates is necessary to lower the
                number of points.
    vertical_faults : bool , option to write all faults vertical
    '''

#   for Fault_name in faults_in_scenario :
    index_fault = faults_names.index(fault_name)
    ColLon = faults_data[index_fault]['lon']
    ColLat = faults_data[index_fault]['lat']
    Depth = faults_data[index_fault]['depth']


    if Depth and all(elem == 'sf' for elem in Depth):
        type_of_fault = 'sf'
    else :
        type_of_fault = 'cf'
        Depth = [float(i) for i in Depth]
    if type_of_fault == 'sf':
        usd = faults_data[index_fault]['upper_sismo_depth']
        lsd = faults_data[index_fault]['lower_sismo_depth']
        # Similar to :meth:`from_fault_data`, we just don't resample edges
        if vertical_faults == False :
            dip_tan = math.tan(math.radians(faults_data[index_fault]['dip']))
            hdist_top = usd / dip_tan
            hdist_bottom = lsd / dip_tan
        else :
            hdist_top = 0.
            hdist_bottom = 0.

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

        do_resample = resample[0]
        # does a resampling to reduce the number of points
        if do_resample == True :
            # parameters
            # min distance between two points
            min_d = resample[1]
            # max distance between two points
            max_d = resample[2]
            # change in azimuth to consider
            az_d = resample[3]

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

            plot_stuf_detail = False
            if plot_stuf_detail == True:
                min_dist_tmp = min_d
                for i_pt in range(len(ColLon)-2):
                    dist_tmp = distance(resampled_ColLon[i_pt+1], resampled_ColLat[i_pt+1],
                                            resampled_ColLon[i_pt], resampled_ColLat[i_pt])
                    if dist_tmp < min_dist_tmp :
                        min_dist_tmp = dist_tmp
                if min_dist_tmp < min_d/2.:
                    print("min dist :",round(min_dist_tmp),"id: ",faults_data[index_fault]['name'])


        txt='\t\t\t<kiteSurface>\n'

        #mean azimuth of the section
        compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
        strike = compass_bearing
        mean_azimuth = (strike + 90.0) % 360
        use_local_azimuth = False

        # loop on profiles
        i_pt = 0
        for x,y in zip(ColLon,ColLat):
            if use_local_azimuth == True :
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
            else :
                azimuth = mean_azimuth

            txt+='\t\t\t<profile>\n'
            txt+='\t\t\t\t<gml:LineString>\n'
            txt+='\t\t\t\t\t<gml:posList>\n'
            xt, yt = point_at(x, y, azimuth, hdist_top)
            txt+='\t\t\t\t\t\t'+"%.5f" %xt+' '+"%.5f" %yt+' '+"%.2f" %usd+' '
            xb, yb = point_at(x, y, azimuth, hdist_bottom)
            txt+="%.5f" % xb+' '+"%.5f" %yb+' '+"%.2f" %lsd+'\n'
            txt+='\t\t\t\t\t</gml:posList>\n'
            txt+='\t\t\t\t</gml:LineString>\n'
            txt+='\t\t\t</profile>\n'
            i_pt +=1

        txt+='\t\t\t</kiteSurface>\n'

    return txt

def wrt_section(txt,section_id,faults_names,faults_data,geotype,resample,vertical_faults):
    '''
    txt : str containing the file info
    geotype : type of geometry of the section (kite, complex, or  simple)
    resample : list with a boolean and the resample parameters if needed
    vertical_faults : bool, option to force all faults vertical
    '''

    section_name = faults_names[section_id]
    txt += '    <section name="'+section_name+'" id="'+str(section_id)+'">\n'

    if geotype == "kite":
        out = wrt_kite_geo(section_name,faults_names,faults_data,
        resample,vertical_faults)
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
    txt += '    </geometryModel>\n'
    txt += '</nrml>\n'
    return txt


    return txt
