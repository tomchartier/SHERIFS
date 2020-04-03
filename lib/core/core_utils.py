# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
 
def progress(rate_in_model,model_MFD,calculation_log,ratio_done,print_percent):
    if ratio_done > 0.01 and ratio_done <= 0.25 and print_percent == True :
        print("1%")
        calculation_log.write("\n1%")
        model_MFD.append(rate_in_model)
        print_percent = False
    if ratio_done > 0.25 and ratio_done <= 0.5 and print_percent == False :
        print( "25%")
        calculation_log.write("\n25%")
        model_MFD.append(rate_in_model)
        print_percent = True
    if ratio_done > 0.5 and ratio_done <= 0.75 and print_percent == True :
        print( "50%")
        calculation_log.write("\n50%")
        model_MFD.append(rate_in_model)
        print_percent = False
    if ratio_done > 0.75 and ratio_done <= 0.9 and print_percent == False :
        print( "75%")
        calculation_log.write("\n75%")
        model_MFD.append(rate_in_model)
        print_percent = True
    if ratio_done > 0.9 and ratio_done <= 0.9999 and print_percent == True :
        print( "90%")
        calculation_log.write("\n90%")
        model_MFD.append(rate_in_model)
        print_percent = False
    return model_MFD,calculation_log,print_percent

def weight_fault_sampling(picked_bin,rup_in_bin,faults_names,faults_slip_rates,slip_rate_use_per_fault,faults_alone,scenarios_names,faults_isolated,index_faults_in_scenario,rup_rates):

    #the faults that haven't been picked often are more likely to be picked
    weight_fault = []
    for i_rup in rup_in_bin[picked_bin]:
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
            
    weight_fault = [i**2 for i in weight_fault]
    weight_fault = np.array(weight_fault)
    weight_fault /= weight_fault.sum()
    
    return weight_fault


def variable_spending(index_fault,M_slip_repartition,faults_budget,slip_rate_use_per_fault,size_of_increment,faults_slip_rates,picked_rup):
    # spends the slip rate correlated with the slip-rate of the faults involeved in the rupture
    nb_sdr_used = 0
    sr_involved = []
    for index in index_fault:
        sr_involved.append(faults_slip_rates[index])
    norm_involved = [int(round(i/min(sr_involved))) for i in sr_involved]
    for index, factor in zip(index_fault,norm_involved) :
        M_slip_repartition[index].append(picked_rup)
        faults_budget[index]+=-(1*factor)
        slip_rate_use_per_fault[index] += size_of_increment
        nb_sdr_used+=(1*factor)
    return M_slip_repartition,faults_budget,slip_rate_use_per_fault,nb_sdr_used



#
#def weight_fault_sampling_old(picked_bin,rup_in_bin,faults_names,faults_slip_rates,slip_rate_use_per_fault,faults_alone,scenarios_names,faults_isolated,index_faults_in_scenario):
#
#    #the faults that haven't been picked often are more likely to be picked
#    weight_fault = []
#    for fault in rup_in_bin[picked_bin]:
#        index_fault = np.where(np.array(faults_names) == fault)[0]
#        if len(index_fault) != 0: #it's a fault
#            sr0 = faults_slip_rates[index_fault[0]]
#            sr_used = slip_rate_use_per_fault[index_fault[0]]
#
#            # calculate the sr factor to help the faster moving faults more
#            # the srfactor goes from 1 to 6 with the ratio of the slip rate to the max of the slip rates
#            if float(sr0)/float(max(faults_slip_rates)) <= 0.2:
#                srfactor = 1.
#            else:
#                srfactor = 5.*float(sr0)/float(max(faults_slip_rates))
#
#            if 1. - float(sr_used)/float(sr0) >= 0.:
#                if fault in faults_alone: #we give a boost for faults that are alone so they can  break more
#                    weight_i = 4. * 5.
#                    #weight_i = 2. * (sr0 - sr_used)**2.
#                    weight_fault.append(weight_i)
#                elif fault in faults_isolated: #we give a boost for faults that are isolated so they can  break more
#                    if (float(sr_used)/float(sr0)) < 0.2:
#                        weight_i = 4. * srfactor
#                    else :
#                        weight_i = (4.-4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor
#                    if weight_i < 1.:
#                        weight_i = 1.
#                    weight_fault.append(weight_i)
#                else :
#                    if (float(sr_used)/float(sr0)) < 0.2:
#                        weight_i = 0.5 * srfactor
#                    else :
#                        weight_i = (0.5-4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor
#
#                    if weight_i < 1.:
#                        weight_i = 1.
#                    weight_fault.append(weight_i)
#            else:
#                weight_i = 0.
#                weight_fault.append(weight_i)
#        else : #it's a scenario, the weight is based on the largest weight of the infolveld faults
#            index_scenario = np.where(np.array(scenarios_names) == fault)[0]
#            ratio_w = 0.
#            for index in index_faults_in_scenario[index_scenario[0]][0] :
#                sr0 = faults_slip_rates[index]
#                sr_used = slip_rate_use_per_fault[index]
#                fault_in_scenario = faults_names[index]
#
#                # calculate the sr factor to help the faster moving faults more
#                # the srfactor goes from 1 to 6 with the ratio of the slip rate to the max of the slip rates
#                if float(sr0)/float(max(faults_slip_rates)) <= 0.2:
#                    srfactor = 1.
#                else:
#                    srfactor = 5.*float(sr0)/float(max(faults_slip_rates))
#
#                if fault_in_scenario in faults_isolated: #we give a boost for faults that are isolated so they can  break more
#                    if (float(sr_used)/float(sr0)) < 0.2:
#                        ratio_w_i = 4. * srfactor
#                    else :
#                       ratio_w_i = (4.-4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor
#                    if ratio_w_i < 1.:
#                        ratio_w_i = 1.
#                else :
#                    if (float(sr_used)/float(sr0)) < 0.2:
#                        ratio_w_i = 0.5 * srfactor
#                    else :
#                        ratio_w_i = (0.5 -4.*(float(sr_used)/float(sr0)+0.3)**6.)*srfactor
#
#                    if ratio_w_i < 1.:
#                        ratio_w_i = 1.
#                if ratio_w_i > ratio_w:
#                    ratio_w = ratio_w_i
#            if ratio_w >= 0.:
#                weight_i = ratio_w
#                weight_fault.append(weight_i)
#            else :
#                weight_i = 0.
#                weight_fault.append(weight_i)
#
#    weight_fault = [i**2 for i in weight_fault]
#    weight_fault = np.array(weight_fault)
#    weight_fault /= weight_fault.sum()
#
#    return weight_fault
