# -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

Populates the magnitude bins with the faults and scenarios that can generate these magnitudes.

@author: Thomas Chartier
"""
import numpy as np
import pickle

def pop(bin_mag,index_rup,rup_rates,M_min,re_use,f_bin_pop):
    if re_use == False :
        print("\t - Populating magnitude bins with ruptures")
        
        # Find ruptures that are smaller and in included in rupt i
        rupt_included = []
        for i_rup in index_rup:
            rupt_included_i = []
            involved_faults = rup_rates.get(str(i_rup)).get('involved_faults')
            for j_rup in index_rup: #check if a smaller rupture included can generate the mag
                if set(rup_rates.get(str(j_rup)).get('involved_faults'))<=set(involved_faults):
                    mmax_j = rup_rates.get(str(j_rup)).get('Mmax')
                    rupt_included_i.append([j_rup,mmax_j])
            rupt_included.append(rupt_included_i)
        
        rup_in_bin = [] #for each bin, list the fault and scenarios in it.
        for mag in bin_mag:
            rup_in_bin_i = []
            for i_rup in index_rup:
                involved_faults = rup_rates.get(str(i_rup)).get('involved_faults')
                Mmax = rup_rates.get(str(i_rup)).get('Mmax')
                if len(involved_faults) == 1 : #its a fault
                    if Mmax >= M_min and mag <= Mmax :
                        rup_in_bin_i.append(i_rup)
                else : #its a scenario
                    if mag <= Mmax :
                        add_scenario_to_bin =True
                        for j_rup in rupt_included[i_rup]:
                            if mag < j_rup[1] and  (j_rup[0] != i_rup):
                                add_scenario_to_bin =False
                        if add_scenario_to_bin == True:
                            rup_in_bin_i.append(i_rup)
            rup_in_bin.append(rup_in_bin_i)
        
        print("\t\t -> Bins populated.")
        with open(f_bin_pop, 'wb') as f:
            pickle.dump(rup_in_bin, f)
        
    else :
        print('Reloading bin pop from data file')
        with open(f_bin_pop, 'rb') as f:
            rup_in_bin = pickle.load(f)
                
    return rup_in_bin

#def pop_old(bin_mag,faults_names,Mmax_faults,M_min,scenarios_names,index_faults_in_scenario,Mmax_scenario):
#    rup_in_bin = [] #for each bin, list the fault and scenarios in it.
#    for mag in bin_mag:
#        rup_in_bin_i = []
#        index_fault = 0
#        for fault in faults_names :
#            if Mmax_faults[index_fault] >= M_min and mag <= Mmax_faults[index_fault] :
#                rup_in_bin_i.append(fault)
#            index_fault += 1
#        rup_in_bin_ii = []
#        index_scenario = 0
#        for scenario in scenarios_names :
#            nb_of_faults_too_small = 0
#            for index_fault in index_faults_in_scenario[index_scenario][0] :
#                if mag > Mmax_faults[index_fault] :
#                    nb_of_faults_too_small += 1
#            # if one of the fault is too small for the bin, the scenario can be used
#            if nb_of_faults_too_small >= (len(index_faults_in_scenario[index_scenario][0]) - 1) and mag <= Mmax_scenario[index_scenario] :
#                rup_in_bin_ii.append(scenario)
#            index_scenario += 1
#        #cleaning so large scenario don't participate on small magnitude earthquake if smaller scenarios are available with the same faults
#        index_scenario =0
#        for scenario_i in rup_in_bin_ii:
#            if len(np.where(np.array(faults_names) == scenario_i)[0]) == 0:
#                index_scenario_i = np.where(np.array(scenarios_names) == scenario_i)[0]
#                index_scenario_i = index_faults_in_scenario[index_scenario_i[0]][0]
#                tagged_faults= []
#                for scenario_j in rup_in_bin_ii:
#                    if len(np.where(np.array(faults_names) == scenario_j)[0]) == 0:
#                        index_scenario_j = np.where(np.array(scenarios_names) == scenario_j)[0]
#                        index_scenario_j = index_faults_in_scenario[index_scenario_j[0]][0]
#                        fault_in = 0
#                        if len(index_scenario_i)>len(index_scenario_j):
#                            for j in index_scenario_j :
#                                if j in index_scenario_i:
#                                    fault_in +=1
#                                    if not j in tagged_faults:
#                                        tagged_faults.append(j)
#                if len(tagged_faults) == len(index_scenario_i):
#                    if mag == Mmax_scenario[index_scenario] :
#                        rup_in_bin_i.append(scenario_i)
#                else :
#                    rup_in_bin_i.append(scenario_i)
#            index_scenario+=1
#
#        rup_in_bin.append(rup_in_bin_i)
#    return rup_in_bin
