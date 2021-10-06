# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems
@author: Thomas Chartier
"""
import sys
import math
import numpy as np
from sherifs.utils.geometry_tools import *
from decimal import Decimal, getcontext
getcontext().prec = 10


def write_simple_fault(index_fault, fault_name, OQ_entry_faults, faults_names,
                       faults_data, Model_name, Domain_in_the_model, ScL_oq,
                       log_mdf_file, M_min, ID_number):

    if fault_name in faults_names:
        i_MFD = np.where(np.array(faults_names) == fault_name)[0][0]
        MFD = OQ_entry_faults[i_MFD]
        ID_number = ID_number + 1
        if not faults_data[index_fault]['domain'] in str(Domain_in_the_model):
            Domain_in_the_model.append(faults_data[index_fault]['domain'])

        ColLon = faults_data[index_fault]['lon']
        ColLat = faults_data[index_fault]['lat']
        Depth = faults_data[index_fault]['depth']

        test_ok = 0
        if Depth and all(elem == 'sf' for elem in Depth):
            type_of_fault = 'sf'
        else :
            type_of_fault = 'cf'
            Depth = [float(i) for i in Depth]

        if type_of_fault == 'sf':
            fault_name = Model_name + '_' + str(fault_name)
            line='\t\t<simpleFaultSource id="'+ str(ID_number
            ) +'" name="'+ str(fault_name
            ) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']
            ) + '">\n'
            test_ok += 1


            # orienting the arrays in order to respect OQ right hand rule
            compass_bearing = calculate_initial_compass_bearing((ColLat[0]
            ,ColLon[0])
            ,(ColLat[-1]
            ,ColLon[-1]))

            if str('N') in str(faults_data[index_fault]['oriented']):
                if compass_bearing < 180. :
                    ColLon = reversed(ColLon)
                    ColLat = reversed(ColLat)
                    #reveresed = 'yes'
            elif str('S') in str(faults_data[index_fault]['oriented']):
                if compass_bearing > 180. :
                    ColLon = reversed(ColLon)
                    ColLat = reversed(ColLat)
                    #reveresed = 'yes'
            elif str('E') in str(faults_data[index_fault]['oriented']):
                if compass_bearing > 90. and compass_bearing < 270. :
                    ColLon = reversed(ColLon)
                    ColLat = reversed(ColLat)
                    #reveresed = 'yes'
            elif str('W') in str(faults_data[index_fault]['oriented']):
                if compass_bearing < 90. or compass_bearing > 270. :
                    ColLon = reversed(ColLon)
                    ColLat = reversed(ColLat)
                    #reveresed = 'yes'

            use_kite = False
            if use_kite == True:
                usd = faults_data[index_fault]['upper_sismo_depth']
                lsd = faults_data[index_fault]['lower_sismo_depth']
                # Similar to :meth:`from_fault_data`, we just don't resample edges
                dip_tan = math.tan(math.radians(faults_data[index_fault]['dip']))
                hdist_top = usd / dip_tan
                hdist_bottom = lsd / dip_tan

                ColLon = list(ColLon)
                ColLat = list(ColLat)

                compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
                strike = compass_bearing
                azimuth = (strike + 90.0) % 360

                line+='\t\t\t<kiteSurface>\n'
                # loop on profiles
                for x,y in zip(ColLon,ColLat):
                    line+='\t\t\t<profile>\n'
                    line+='\t\t\t\t<gml:LineString>\n'
                    line+='\t\t\t\t\t<gml:posList>\n'
                    xt, yt = point_at(x, y, azimuth, hdist_top)
                    line+='\t\t\t\t\t\t'+str(xt)+' '+str(yt)+' '+str(usd)+' '
                    xb, yb = point_at(x, y, azimuth, hdist_bottom)
                    line+=str(xb)+' '+str(yb)+' '+str(lsd)+'\n'
                    line+='\t\t\t\t\t</gml:posList>\n'
                    line+='\t\t\t\t</gml:LineString>\n'
                    line+='\t\t\t</profile>\n'
                line+='\t\t\t</kiteSurface>\n'

            else :
                line+='\t\t\t<simpleFaultGeometry>\n'
                line+='\t\t\t\t<gml:LineString>\n'
                line+='\t\t\t\t\t<gml:posList>\n'
                for x,y in zip(ColLon,ColLat):
                    #polygon.append((x,y)) #ecriture du polygone de la zone
                    line+='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                line+='\t\t\t\t\t</gml:posList>\n'
                line+='\t\t\t\t</gml:LineString>\n'

                line+='\t\t\t\t<dip>'+ str(faults_data[index_fault]['dip']) +'</dip>\n'
                line+='\t\t\t\t<upperSeismoDepth>'+ str(faults_data[index_fault]['upper_sismo_depth']) +'</upperSeismoDepth>\n'
                line+='\t\t\t\t<lowerSeismoDepth>'+ str(faults_data[index_fault]['lower_sismo_depth']) +'</lowerSeismoDepth>\n'
                line+='\t\t\t</simpleFaultGeometry>\n'

        if type_of_fault == 'cf':
            fault_name = Model_name + '_' + str(fault_name)
            line='\t\t<complexFaultSource id="'+ str(ID_number) +'" name="'+ str(fault_name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'
            test_ok += 1
            line+='\t\t\t<complexFaultGeometry>\n'

            index_edge = 0
            for depth_i in sorted(set(Depth)):
                indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                if index_edge == 0:
                    line+='\t\t\t\t<faultTopEdge>\n'
                elif index_edge == len(set(Depth))-1:
                    line+='\t\t\t\t<faultBottomEdge>\n'
                else :
                    line+='\t\t\t\t<intermediateEdge>\n'

                line+='\t\t\t\t<gml:LineString>\n'
                line+='\t\t\t\t\t<gml:posList>\n'
                for index in indexes_for_edge_i:
                    line+='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'

                line+='\t\t\t\t\t</gml:posList>\n'
                line+='\t\t\t\t</gml:LineString>\n'
                if index_edge == 0:
                    line+='\t\t\t\t</faultTopEdge>\n'
                elif index_edge == len(set(Depth))-1:
                    line+='\t\t\t\t</faultBottomEdge>\n'
                else :
                    line+='\t\t\t\t</intermediateEdge>\n'
                index_edge+=1

            line+='\t\t\t</complexFaultGeometry>\n'


        if test_ok == 0:
            print('!!!!!!!!!! Problem with the fault Geometry, please check input file''')
            sys.exit()

        ########################################################
        #scaling law and aspect ratio
        ########################################################

        line+='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
        line+='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'

        ########################################################
        # seismicity of the fault            ♀
        ########################################################
        line += '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'
        if sum(MFD)!=0 :
            log_mdf_file.write(str(fault_name)+'\t'+str(M_min)+'\t'+' '.join(list(map(str, MFD)))+'\n')
            line += '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
        else:
            MFD[0] += 0.00000000001 #so it's not zero
            log_mdf_file.write(str(fault_name)+'\t'+str(M_min)+'\t'+' '.join(list(map(str, MFD)))+'\n')
            line += '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'

        line += '\t\t\t</incrementalMFD>\n'

        line +='\t\t\t<rake>'+str(faults_data[index_fault]['rake'])+'</rake>\n'
        if type_of_fault == 'sf':
            line+='\t\t</simpleFaultSource>\n'
        if type_of_fault == 'cf':
            line+='\t\t</complexFaultSource>\n'

    return line,Domain_in_the_model,ID_number

def write_characteristic_scenario(scenarios_names,OQ_entry_scenarios,index_faults_in_scenario,scenario,faults_names,Model_name,faults_data,log_mdf_file,M_min,ID_number):

    index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
    MFD = OQ_entry_scenarios[index_scenario]
    if sum(MFD)!=0:
        ID_number = ID_number + 1
        index_fault = faults_names.index(scenario[1]['f_1'][0])
#                    self.FaultProperties(scenario[1]['f_1'])
        scenar_name = '_'.join("{!s}={!r}".format(key,val) for (key,val) in scenario[1].items())
        Fault_Name = Model_name + '_scenario_' + str(scenar_name)
        line='\t\t<characteristicFaultSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'


        ########################################################
        #Part concerning the geometry
        ########################################################
        index_faults_in_scenario =  index_faults_in_scenario[index_scenario][0]
        faults_in_scenario = np.take(faults_names,index_faults_in_scenario)

        line=+'\t\t\t<surface>\n'

        scenario_mechanism = []
        for Fault_name in faults_in_scenario :
            index_fault = faults_names.index(Fault_name)
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
                line+='\t\t\t<simpleFaultGeometry>\n'
                line+='\t\t\t\t<gml:LineString>\n'
                line+='\t\t\t\t\t<gml:posList>\n'
                #polygon = []

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

                for x,y in zip(ColLon,ColLat):
                    #polygon.append((x,y)) #ecriture du polygone de la zone
                    line+='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                line+='\t\t\t\t\t</gml:posList>\n'
                line+='\t\t\t\t</gml:LineString>\n'
                line+='\t\t\t\t<dip>'+ str(faults_data[index_fault]['dip']) +'</dip>\n'
                line+='\t\t\t\t<upperSeismoDepth>'+ str(faults_data[index_fault]['upper_sismo_depth']) +'</upperSeismoDepth>\n'
                line+='\t\t\t\t<lowerSeismoDepth>'+ str(faults_data[index_fault]['lower_sismo_depth']) +'</lowerSeismoDepth>\n'
                line+='\t\t\t</simpleFaultGeometry>\n'

            if type_of_fault == 'cf':
                line+='\t\t\t<complexFaultGeometry>\n'

                index_edge = 0
                for depth_i in sorted(set(Depth)):
                    indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                    if index_edge == 0:
                        line+='\t\t\t\t<faultTopEdge>\n'
                    elif index_edge == len(set(Depth))-1:
                        line+='\t\t\t\t<faultBottomEdge>\n'
                    else :
                        line+='\t\t\t\t<intermediateEdge>\n'

                    line+='\t\t\t\t<gml:LineString>\n'
                    line+='\t\t\t\t\t<gml:posList>\n'
                    for index in indexes_for_edge_i:
                        line+='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'

                    line+='\t\t\t\t\t</gml:posList>\n'
                    line+='\t\t\t\t</gml:LineString>\n'
                    if index_edge == 0:
                        line+='\t\t\t\t</faultTopEdge>\n'
                    elif index_edge == len(set(Depth))-1:
                        line+='\t\t\t\t</faultBottomEdge>\n'
                    else :
                        line+='\t\t\t\t</intermediateEdge>\n'
                    index_edge+=1

                line+='\t\t\t</complexFaultGeometry>\n'



        line+='\t\t\t</surface>\n'

        ########################################################
        #scaling law and aspect ratio
        ########################################################

        line+='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
        line+='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'

        ########################################################
        # seismicity of the scenario
        ########################################################
#
        log_mdf_file.write(str(Fault_Name) + '\t' + str(M_min) + '\t' + ' '.join(list(map(str, MFD)))+'\n')
        line += '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'
        line += '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
        line += '\t\t\t</incrementalMFD>\n'

        #find the dominant kinematic of the scenario
        rake= np.mean(scenario_mechanism)
        line+='\t\t\t<rake>'+str(rake)+'</rake>\n'
        line='\t\t</characteristicFaultSource>\n'
    return line,ID_number

def write_non_parametric_one_fault(index_fault,fault_name,OQ_entry_faults,
faults_names,faults_data,Model_name,Domain_in_the_model,
ScL_oq,log_mdf_file,M_min,ID_number,explo_time):
    '''
    Write a non parametric source with fault trace defined as set of kite fault surfaces
    for now the rupture occurs once of not. several occurences are not supported yet.
    For one fault only
    '''

    # index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
    # MFD = OQ_entry_scenarios[index_scenario]
    i_MFD = np.where(np.array(faults_names) == fault_name)[0][0]
    MFD = OQ_entry_faults[i_MFD]

    line = ''

    if sum(MFD)!=0:
        ID_number = ID_number + 1
        # index_fault = faults_names.index(scenario[1]['f_1'][0])
        #
        # scenar_name = '_'.join("{!s}={!r}".format(key,val) for (key,val) in scenario[1].items())
        # Fault_Name = Model_name + '_scenario_' + str(scenar_name)
        Fault_Name = Model_name + '_' + str(fault_name)
        if not faults_data[index_fault]['domain'] in str(Domain_in_the_model):
            Domain_in_the_model.append(faults_data[index_fault]['domain'])

        line='\t\t<nonParametricSeismicSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'


        ########################################################
        #Part concerning the geometry
        ########################################################
        # index_faults_in_scenario =  index_faults_in_scenario[index_scenario][0]
        # faults_in_scenario = np.take(faults_names,index_faults_in_scenario)

        str_geom = ''


        bin_mag = np.linspace(M_min,M_min+0.1*len(MFD)+0.1,num = 2+len(MFD))
        for mag,i_mag in zip(bin_mag,range(len(bin_mag))):
            if i_mag <= len(MFD)-1:
                # annual_rate
                l = MFD[i_mag]
                if l!=0. :
                    t = explo_time
                    # Poisson probability
                    p_occur_1 = np.float32((l*t)**1*np.exp(-l*t)/(np.math.factorial(1)))
#                    p_occur_2 = np.float32((l*t)**2*np.exp(-l*t)/(np.math.factorial(2)))
#                    p_occur_3 = np.float32((l*t)**3*np.exp(-l*t)/(np.math.factorial(3)))
#                    p_occur_4 = np.float32((l*t)**4*np.exp(-l*t)/(np.math.factorial(4)))
#                    p_occur_5 = np.float32((l*t)**5*np.exp(-l*t)/(np.math.factorial(5)))

                    xxx1 = Decimal('{:.8f}'.format(np.float32(p_occur_1)))
#                    xxx2 = Decimal('{:.8f}'.format(np.float32(p_occur_2)))
#                    xxx3 = Decimal('{:.8f}'.format(np.float32(p_occur_3)))
#                    xxx4 = Decimal('{:.8f}'.format(np.float32(p_occur_4)))
#                    xxx5 = Decimal('{:.8f}'.format(np.float32(p_occur_5)))
#                    p_not_occur = (Decimal('1')-(xxx1+xxx2+xxx3+xxx4+xxx5))
                    p_not_occur = (Decimal('1')-xxx1)
                    p_not_occur = '{:.8f}'.format(p_not_occur)


#                    line+='\t\t\t<multiPlanesRupture probs_occur="'+ str(p_not_occur) +' '+ str(xxx1) +' '+ str(xxx2) +' '+ str(xxx3) +' '+ str(xxx4) +' '+ str(xxx5) +'">\n'
                    line+='\t\t\t<multiPlanesRupture probs_occur="'+ str(p_not_occur) +' '+ str(xxx1) +'">\n'
                    line+='\t\t\t<magnitude>'+ str(mag) +'</magnitude>\n'
                    if str_geom == '':
                        scenario_mechanism = []

#                        for Fault_name in faults_in_scenario :
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



                    hypo_depth = (faults_data[index_fault]['upper_sismo_depth']+faults_data[index_fault]['lower_sismo_depth'])/2.

                    str_geom+='\t\t\t<hypocenter depth="'+str(hypo_depth)+'" lat="'+str(np.mean(faults_data[index_fault]['lat']))+'" lon="'+str(np.mean(faults_data[index_fault]['lon']))+'"/>\n'

#                        str_geom+='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
#                        str_geom+='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'

                    rake= np.mean(scenario_mechanism)
                    str_geom+='\t\t\t<rake>'+str(rake)+'</rake>\n'

                    line+=str_geom
                    line+='\t\t\t</multiPlanesRupture>\n'

        line+='\t\t</nonParametricSeismicSource>\n'


    return line,Domain_in_the_model,ID_number


def write_non_parametric_source(scenario,scenarios_names,OQ_entry_scenarios,index_faults_in_scenario,faults_names,faults_data,Model_name,Domain_in_the_model,ScL_oq,log_mdf_file,explo_time,M_min,ID_number):
    '''
    Write a non parametric source with fault trace defined as set of kite fault surfaces
    for now the rupture occurs once of not. several occurences are not supported yet.
    '''

    index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
    MFD = OQ_entry_scenarios[index_scenario]

    line = ''

    if sum(MFD)!=0:
        ID_number = ID_number + 1
        index_fault = faults_names.index(scenario[1]['f_1'][0])

        scenar_name = '_'.join("{!s}={!r}".format(key,val) for (key,val) in scenario[1].items())
        Fault_Name = Model_name + '_scenario_' + str(scenar_name)

        line='\t\t<nonParametricSeismicSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'


        ########################################################
        #Part concerning the geometry
        ########################################################
        index_faults_in_scenario =  index_faults_in_scenario[index_scenario][0]
        faults_in_scenario = np.take(faults_names,index_faults_in_scenario)

        str_geom = ''


        bin_mag = np.linspace(M_min,M_min+0.1*len(MFD)+0.1,num = 2+len(MFD))
        for mag,i_mag in zip(bin_mag,range(len(bin_mag))):
            if i_mag <= len(MFD)-1:
                # annual_rate
                l = MFD[i_mag]
                if l!=0. :
                    t = explo_time
                    # Poisson probability
                    p_occur_1 = np.float32((l*t)**1*np.exp(-l*t)/(np.math.factorial(1)))
#                    p_occur_2 = np.float32((l*t)**2*np.exp(-l*t)/(np.math.factorial(2)))
#                    p_occur_3 = np.float32((l*t)**3*np.exp(-l*t)/(np.math.factorial(3)))
#                    p_occur_4 = np.float32((l*t)**4*np.exp(-l*t)/(np.math.factorial(4)))
#                    p_occur_5 = np.float32((l*t)**5*np.exp(-l*t)/(np.math.factorial(5)))

                    xxx1 = Decimal('{:.8f}'.format(np.float32(p_occur_1)))
#                    xxx2 = Decimal('{:.8f}'.format(np.float32(p_occur_2)))
#                    xxx3 = Decimal('{:.8f}'.format(np.float32(p_occur_3)))
#                    xxx4 = Decimal('{:.8f}'.format(np.float32(p_occur_4)))
#                    xxx5 = Decimal('{:.8f}'.format(np.float32(p_occur_5)))
#                    p_not_occur = (Decimal('1')-(xxx1+xxx2+xxx3+xxx4+xxx5))
                    p_not_occur = (Decimal('1')-xxx1)
                    p_not_occur = '{:.8f}'.format(p_not_occur)


#                    line+='\t\t\t<multiPlanesRupture probs_occur="'+ str(p_not_occur) +' '+ str(xxx1) +' '+ str(xxx2) +' '+ str(xxx3) +' '+ str(xxx4) +' '+ str(xxx5) +'">\n'
                    line+='\t\t\t<multiPlanesRupture probs_occur="'+ str(p_not_occur) +' '+ str(xxx1) +'">\n'
                    line+='\t\t\t<magnitude>'+ str(mag) +'</magnitude>\n'
                    if str_geom == '':
                        scenario_mechanism = []

                        for Fault_name in faults_in_scenario :
                            index_fault = faults_names.index(Fault_name)
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



                        hypo_depth = (faults_data[index_fault]['upper_sismo_depth']+faults_data[index_fault]['lower_sismo_depth'])/2.

                        str_geom+='\t\t\t<hypocenter depth="'+str(hypo_depth)+'" lat="'+str(np.mean(faults_data[index_fault]['lat']))+'" lon="'+str(np.mean(faults_data[index_fault]['lon']))+'"/>\n'

#                        str_geom+='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
#                        str_geom+='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'

                        rake= np.mean(scenario_mechanism)
                        str_geom+='\t\t\t<rake>'+str(rake)+'</rake>\n'

                    line+=str_geom
                    line+='\t\t\t</multiPlanesRupture>\n'

        line+='\t\t</nonParametricSeismicSource>\n'


    return line,ID_number



def write_non_parametric_source_multi_types(scenario,scenarios_names,OQ_entry_scenarios,index_faults_in_scenario,faults_names,faults_data,Model_name,Domain_in_the_model,ScL_oq,log_mdf_file,explo_time,M_min,ID_number):

    index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
    MFD = OQ_entry_scenarios[index_scenario]

    line = ''

    if sum(MFD)!=0:
        ID_number = ID_number + 1
        index_fault = faults_names.index(scenario[1]['f_1'][0])
#                    self.FaultProperties(scenario[1]['f_1'])
        scenar_name = '_'.join("{!s}={!r}".format(key,val) for (key,val) in scenario[1].items())
        Fault_Name = Model_name + '_scenario_' + str(scenar_name)

        line='\t\t<nonParametricSeismicSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'


        ########################################################
        #Part concerning the geometry
        ########################################################
        index_faults_in_scenario =  index_faults_in_scenario[index_scenario][0]
        faults_in_scenario = np.take(faults_names,index_faults_in_scenario)

        str_geom = ''


        bin_mag = np.linspace(M_min,M_min+0.1*len(MFD)+0.1,num = 2+len(MFD))
        for mag,i_mag in zip(bin_mag,range(len(bin_mag))):
            if i_mag <= len(MFD)-1:
                # annual_rate
                l = MFD[i_mag]
                if l!=0. :
                    t = explo_time
                    # Poisson probability
                    p_occur_1 = np.float32((l*t)**1*np.exp(-l*t)/(np.math.factorial(1)))
#                    p_occur_2 = np.float32((l*t)**2*np.exp(-l*t)/(np.math.factorial(2)))
#                    p_occur_3 = np.float32((l*t)**3*np.exp(-l*t)/(np.math.factorial(3)))
#                    p_occur_4 = np.float32((l*t)**4*np.exp(-l*t)/(np.math.factorial(4)))
#                    p_occur_5 = np.float32((l*t)**5*np.exp(-l*t)/(np.math.factorial(5)))

                    xxx1 = Decimal('{:.8f}'.format(np.float32(p_occur_1)))
#                    xxx2 = Decimal('{:.8f}'.format(np.float32(p_occur_2)))
#                    xxx3 = Decimal('{:.8f}'.format(np.float32(p_occur_3)))
#                    xxx4 = Decimal('{:.8f}'.format(np.float32(p_occur_4)))
#                    xxx5 = Decimal('{:.8f}'.format(np.float32(p_occur_5)))
#                    p_not_occur = (Decimal('1')-(xxx1+xxx2+xxx3+xxx4+xxx5))
                    p_not_occur = (Decimal('1')-xxx1)
                    p_not_occur = '{:.8f}'.format(p_not_occur)


#                    line+='\t\t\t<multiPlanesRupture probs_occur="'+ str(p_not_occur) +' '+ str(xxx1) +' '+ str(xxx2) +' '+ str(xxx3) +' '+ str(xxx4) +' '+ str(xxx5) +'">\n'
                    line+='\t\t\t<multiPlanesRupture probs_occur="'+ str(p_not_occur) +' '+ str(xxx1) +'">\n'
                    line+='\t\t\t<magnitude>'+ str(mag) +'</magnitude>\n'
                    if str_geom == '':
                        scenario_mechanism = []

                        for Fault_name in faults_in_scenario :
                            index_fault = faults_names.index(Fault_name)
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
                                str_geom+='\t\t\t<simpleFaultGeometry>\n'
                                str_geom+='\t\t\t\t<gml:LineString>\n'
                                str_geom+='\t\t\t\t\t<gml:posList>\n'
                                #polygon = []

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

                                for x,y in zip(ColLon,ColLat):
                                    #polygon.append((x,y)) #ecriture du polygone de la zone
                                    str_geom+='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                                str_geom+='\t\t\t\t\t</gml:posList>\n'
                                str_geom+='\t\t\t\t</gml:LineString>\n'
                                str_geom+='\t\t\t\t<dip>'+ str(faults_data[index_fault]['dip']) +'</dip>\n'
                                str_geom+='\t\t\t\t<upperSeismoDepth>'+ str(faults_data[index_fault]['upper_sismo_depth']) +'</upperSeismoDepth>\n'
                                str_geom+='\t\t\t\t<lowerSeismoDepth>'+ str(faults_data[index_fault]['lower_sismo_depth']) +'</lowerSeismoDepth>\n'
                                str_geom+='\t\t\t</simpleFaultGeometry>\n'

                            if type_of_fault == 'cf':
                                str_geom+='\t\t\t<complexFaultGeometry>\n'

                                index_edge = 0
                                for depth_i in sorted(set(Depth)):
                                    indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                                    if index_edge == 0:
                                        str_geom+='\t\t\t\t<faultTopEdge>\n'
                                    elif index_edge == len(set(Depth))-1:                    str_geom+='\t\t\t\t<faultBottomEdge>\n'
                                    else :
                                        str_geom+='\t\t\t\t<intermediateEdge>\n'

                                    str_geom+='\t\t\t\t<gml:LineString>\n'
                                    str_geom+='\t\t\t\t\t<gml:posList>\n'
                                    for index in indexes_for_edge_i:
                                        str_geom+='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'

                                    str_geom+='\t\t\t\t\t</gml:posList>\n'
                                    str_geom+='\t\t\t\t</gml:LineString>\n'
                                    if index_edge == 0:
                                        str_geom+='\t\t\t\t</faultTopEdge>\n'
                                    elif index_edge == len(set(Depth))-1:
                                        str_geom+='\t\t\t\t</faultBottomEdge>\n'
                                    else :
                                        str_geom+='\t\t\t\t</intermediateEdge>\n'
                                    index_edge+=1

                                str_geom+='\t\t\t</complexFaultGeometry>\n'

                        hypo_depth = (faults_data[index_fault]['upper_sismo_depth']+faults_data[index_fault]['lower_sismo_depth'])/2.

                        str_geom+='\t\t\t<hypocenter depth="'+str(hypo_depth)+'" lat="'+str(np.mean(faults_data[index_fault]['lat']))+'" lon="'+str(np.mean(faults_data[index_fault]['lon']))+'"/>\n'
                        ########################################################
                        #scaling law and aspect ratio
                        ########################################################

                        str_geom+='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
                        str_geom+='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'

                        ########################################################
                        # seismicity of the scenario
                        ########################################################
                    #
                    #        log_mdf_file.write(str(Fault_Name) + '\t' + str(M_min) + '\t' + ' '.join(list(map(str, MFD)))+'\n')
                    #        line += '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'
                    #        line += '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
                    #        line += '\t\t\t</incrementalMFD>\n'

                        #find the dominant kinematic of the scenario
                        rake= np.mean(scenario_mechanism)
                        str_geom+='\t\t\t<rake>'+str(rake)+'</rake>\n'

                    line+=str_geom
                    line+='\t\t\t</multiPlanesRupture>\n'

        line+='\t\t</nonParametricSeismicSource>\n'


    return line,ID_number
