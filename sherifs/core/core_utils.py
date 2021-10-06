# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

@author: Thomas Chartier
"""
import geojson
import numpy as np
from scipy.interpolate import interp1d
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from sherifs.core import rates, mfd_shape


def seconds_to_str(seconds):
    time_str = " "
    if seconds >= 60:
        minutes = seconds / 60.
        if seconds > 60*60 :
            hours = seconds / (60.*60.)
            if seconds > 60*60*24 :
                days = seconds / (60.*60.*24.)
                time_str += str(int(days))+" d "
                time_str += str(int(hours-days*24.))+" h "
                time_str += str(int(minutes-hours*60.))+" m "
                time_str += str(int(seconds-minutes*60.))+" s"
            else :
                time_str += str(int(hours))+" h "
                time_str += str(int(minutes-hours*60.))+" m "
                time_str += str(int(seconds-minutes*60.))+" s"
        else :
            time_str += str(int(minutes))+" m "
            time_str += str(int(seconds-minutes*60.))+" s"
    else :
        time_str += str(int(seconds))+" s"

    return time_str

def progress(model_MFD,calculation_log,ratio_done,print_percent,rup_rates,fault_prop,bin_mag):
    if ratio_done > 0.01 and ratio_done <= 0.25 and print_percent == True :
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        print("1%")
        calculation_log.write("\n1%")
        model_MFD.append(rate_in_model)
        print_percent = False
    if ratio_done > 0.25 and ratio_done <= 0.5 and print_percent == False :
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        print( "25%")
        calculation_log.write("\n25%")
        model_MFD.append(rate_in_model)
        print_percent = True
    if ratio_done > 0.5 and ratio_done <= 0.75 and print_percent == True :
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        print( "50%")
        calculation_log.write("\n50%")
        model_MFD.append(rate_in_model)
        print_percent = False
    if ratio_done > 0.75 and ratio_done <= 0.9 and print_percent == False :
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        print( "75%")
        calculation_log.write("\n75%")
        model_MFD.append(rate_in_model)
        print_percent = True
    if ratio_done > 0.9 and ratio_done <= 0.9999 and print_percent == True :
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        print( "90%")
        calculation_log.write("\n90%")
        model_MFD.append(rate_in_model)
        print_percent = False
    return model_MFD,calculation_log,print_percent

def weight_fault_sampling(picked_bin,rup_in_bin,faults_names,faults_slip_rates,slip_rate_use_per_fault,faults_alone,scenarios_names,faults_isolated,index_faults_in_scenario,rup_rates,empty_rups):

    #the faults that haven't been picked often are more likely to be picked
    weight_fault = []
    for i_rup in rup_in_bin[picked_bin]:
        if str(i_rup) in empty_rups : #the rup is empty
            weight_fault.append(0.)
        else :
            involved_faults = rup_rates.get(str(i_rup)).get('involved_faults')
            #index_fault = np.where(np.array(faults_names) == fault)[0]
            if len(involved_faults) == 0: #it's a fault
                sr0 = faults_slip_rates[involved_faults[0]]
                sr_used = slip_rate_use_per_fault[involved_faults[0]]

                # calculate the sr factor to help the faster moving faults more
                # the srfactor goes from 1 to 6 with the ratio of the slip rate to the max of the slip rates
                if float(sr0)/float(max(faults_slip_rates)) <= 0.2:
                    srfactor = 1.
                else:
                    srfactor = 5.*float(sr0)/float(max(faults_slip_rates))

                if 1. - float(sr_used)/float(sr0) >= 0.:
                    if involved_faults[0] in faults_alone: #we give a boost for faults that are alone so they can  break more
                        weight_i = 4. * 5.
                        #weight_i = 2. * (sr0 - sr_used)**2.
                        weight_fault.append(weight_i)
                    elif involved_faults[0] in faults_isolated: #we give a boost for faults that are isolated so they can  break more
                        if (float(sr_used)/float(sr0)) < 0.2:
                            weight_i = 4. * srfactor
                        else :
                            weight_i = (4.-4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor
                        if weight_i < 1.:
                            weight_i = 1.
                        weight_fault.append(weight_i)
                    else :
                        if (float(sr_used)/float(sr0)) < 0.2:
                            weight_i = 0.5 * srfactor
                        else :
                            weight_i = (0.5-4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor

                        if weight_i < 1.:
                            weight_i = 1.
                        weight_fault.append(weight_i)
                else:
                    weight_i = 0.
                    weight_fault.append(weight_i)
            else : #it's a scenario, the weight is based on the largest weight of the infolveld faults
                #index_scenario = np.where(np.array(scenarios_names) == fault)[0]
                ratio_w = 0.
                for index in involved_faults :
                    sr0 = faults_slip_rates[index]
                    sr_used = slip_rate_use_per_fault[index]
                    fault_in_scenario = faults_names[index]

                    # calculate the sr factor to help the faster moving faults more
                    # the srfactor goes from 1 to 6 with the ratio of the slip rate to the max of the slip rates
                    if float(sr0)/float(max(faults_slip_rates)) <= 0.2:
                        srfactor = 1.
                    else:
                        srfactor = 5.*float(sr0)/float(max(faults_slip_rates))

                    if fault_in_scenario in faults_isolated: #we give a boost for faults that are isolated so they can  break more
                        if (float(sr_used)/float(sr0)) < 0.2:
                            ratio_w_i = 4. * srfactor
                        else :
                           ratio_w_i = (4.-4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor
                        if ratio_w_i < 1.:
                            ratio_w_i = 1.
                    else :
                        if (float(sr_used)/float(sr0)) < 0.2:
                            ratio_w_i = 0.5 * srfactor
                        else :
                            ratio_w_i = (0.5 -4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor

                        if ratio_w_i < 1.:
                            ratio_w_i = 1.
                    if ratio_w_i > ratio_w:
                        ratio_w = ratio_w_i

                if ratio_w >= 0.:
                    weight_i = ratio_w
                    weight_fault.append(weight_i)
                else :
                    weight_i = 0.
                    weight_fault.append(weight_i)


    if sum(weight_fault) == 0.:
        bin_is_empty = True
        weight_fault = [0]
    else :
        bin_is_empty = False
        weight_fault = [i**2 for i in weight_fault]
        weight_fault = np.array(weight_fault)
        weight_fault /= weight_fault.sum()

    return weight_fault


def variable_spending(index_fault,M_slip_repartition,faults_budget,slip_rate_use_per_fault,size_of_increment,faults_slip_rates,picked_rup,faults_names,sum_fault_budget):
    # spends the slip rate correlated with the slip-rate of the faults involeved in the rupture
    nb_sdr_used = 0
    sr_involved = []
    for index in index_fault:
        sr_involved.append(faults_slip_rates[index])
    norm_involved = [int(round(i/min(sr_involved))) for i in sr_involved]
    for index, factor in zip(index_fault,norm_involved) :
        #M_slip_repartition[index].append(picked_rup)
        M_slip_repartition[str(faults_names[index])][str(picked_rup)] += 1
        faults_budget[index]+=-(1*factor)
        sum_fault_budget+=-(1*factor)
        slip_rate_use_per_fault[index] += size_of_increment
        nb_sdr_used+=(1*factor)
    return M_slip_repartition,faults_budget,slip_rate_use_per_fault,nb_sdr_used,sum_fault_budget


def link_rup_mfd_area(rup_rates,f_mfd_area,faults_lon,faults_lat,bin_mag,bg_ratio):
    '''
    Link each rupture with the area that describe the
    local mfd to respect (on top of the global mfd).

    param :
    rup_rates : dict, contains the ruptures and the associated info and rates.

    f_mfd_area : str, path to the geojson file that contains the areas where a local mfd should be respected.

    faults_lon, faults_lat : list, list of the coordinates of the fault points.

    bg_ratio : list, ratio of the seismicity occuring on the faults and not in the bachground for the whole model.

    returns :
    local_mfds : list, list of mfd shapes tobe followed locally. (right now only GR is possible, but it can be modified.

    associated_rup : list, ruptures associated with each area where a local mfd must be respected.

    associated_weight : list, ratio between 0 and one representing how much of a given rupture is in the area.

    '''

    with open(f_mfd_area) as f:
        gj = geojson.load(f)
    areas = gj["features"]
    local_mfds, associated_rup, associated_weight = [], [], []
    for area_i in areas:
        if 'b_value' in area_i["properties"].keys():
            if not area_i["properties"]["b_value"] == None :
                b_value = area_i["properties"]["b_value"]
                poly = []
                for pt in area_i["geometry"]["coordinates"][0][0]:
                    poly.append((pt[0],pt[1]))
                polygon = Polygon(poly)

                mfd_param = {"b_value" : b_value}

                p_MFD = mfd_shape.GR(mfd_param,bin_mag)
                p_MFD /= sum(p_MFD)

                # Find if a local background has to be applied for the shape
                apply_global_bg = True
                if 'bg' in area_i["properties"].keys():
                    if not area_i["properties"]["bg"] in [None,"global"] :
                        apply_global_bg = False
                if apply_global_bg == True :
                    bg_ratio_loc = bg_ratio
                else : # apply a local background
                    if type(area_i["properties"]["bg"]) == type([0,0]):
                        bg_ratio_loc = area_i["properties"]["bg"]
                    else :
                        print("ERROR !!! please verify the local background proportions")

                # Apply the background ratio to the local target shape
                bin_mag_fault_prop = [ 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.]
                fault_prop_inc = bg_ratio_loc
                bin_mag_fault_prop.append(10.)
                fault_prop_inc = np.append(np.array(fault_prop_inc),1.)
                fault_prop = interp1d(bin_mag_fault_prop,fault_prop_inc)
                p_MFD_faults = []
                index_mag = 0
                for mag in bin_mag :
                    p_MFD_faults.append(fault_prop(mag) * p_MFD[index_mag])
                    index_mag+=1

                # add to the list of local MFDs
                local_mfds.append(p_MFD_faults)

                # Find the ruptures in the area
                # The sum of their eq rates will need to follow the local shape
                associated_rup_i = []
                associated_weight_i = []
                for rup_i in rup_rates:

                    id_sections = rup_rates.get(rup_i).get("involved_faults")
                    nb_sections = len(id_sections)
                    nb_in = 0
                    for id_s in id_sections:
                        is_in = False
                        for lon_i,lat_i in zip(faults_lon[id_s],faults_lat[id_s]):
                            if is_in == False :
                                if polygon.contains(Point(lon_i, lat_i)):
                                    is_in = True
                        if is_in == True :
                           nb_in += 1
                    if nb_in != 0 :
                        associated_rup_i.append(rup_rates.get(rup_i).get("rup_id"))
                        associated_weight_i.append(float(nb_in)/float(nb_sections))
                associated_rup.append(associated_rup_i)
                associated_weight.append(associated_weight_i)

    return local_mfds, associated_rup, associated_weight

def check_local_mfd(rup_rates, rup_in_bin, picked_bin, bin_mag, local_mfds, associated_rup, associated_weight):
    '''
    Calculate the factor to apply to the weight of the rupture in order to help to respect the local MFD.

    param :
    rup_rates : dict, contains the ruptures and the associated info and rates.

    rup_in_bin : list, ruptures that can produce the magnitude picked.

    picked_bin : int, index of the picked magnitude bin.

    bin_mag : list, bining in magnitude.

    local_mfds : list, list of mfd shapes tobe followed locally.

    associated_rup : list, ruptures associated with each area where a local mfd must be respected.

    associated_weight : list, ratio between 0 and one representing how much of a given rupture is in the area.

    returns :
    factor_on_weight : list, factor to be appleid to the weight of picking the rupture.
    '''
    factor_on_weight = [1. for _ in rup_in_bin]

    for local_mfd,associated_rup_i,associated_weight_i in zip(local_mfds,associated_rup,associated_weight):
        rates = np.zeros(len(bin_mag))
        for id_rup, w_i in zip(associated_rup_i,associated_weight_i) :
            rates  += rup_rates.get(str(id_rup)).get('rates') * w_i
        if sum(rates) != 0. :
            p_rates = (rates) / sum(rates) # transform to a probability distribution.

#            wiggle_room = 0.01
#            # 1% of the first rate of bin as wiggle room
#            # this allow the range of accepttalbe rates to be largers as the magnitude increase
#            min_acceptable = local_mfd[picked_bin] - local_mfd[0] * wiggle_room
#            max_acceptable = local_mfd[picked_bin] + local_mfd[0] * wiggle_room

            wiggle_room = 0.2
            # Uniform +- 20% for all bins
            min_acceptable = local_mfd[picked_bin] - local_mfd[picked_bin] * wiggle_room
            max_acceptable = local_mfd[picked_bin] + local_mfd[picked_bin] * wiggle_room

            # the factor is to be applied on the weight in order to help fit in the acceptable range.
            factor = 100.
            if p_rates[picked_bin] < min_acceptable :
                for id_w in range(len(rup_in_bin)):
                    if rup_in_bin[id_w] in associated_rup_i:
                        factor_on_weight[id_w] =  factor_on_weight[id_w] * factor

            if p_rates[picked_bin] > max_acceptable :
                for id_w in range(len(rup_in_bin)):
                    if rup_in_bin[id_w] in associated_rup_i:
                        factor_on_weight[id_w] =  factor_on_weight[id_w] / factor


    return factor_on_weight
