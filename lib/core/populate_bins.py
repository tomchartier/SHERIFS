# -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

Populates the magnitude bins with the faults and scenarios that can generate these magnitudes.

@author: Thomas Chartier
"""
import numpy as np

def pop(bin_mag,index_rup,rup_rates,M_min):
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
                add_scenario_to_bin =True
                for j_rup in index_rup: #check if a smaller rupture included can generate the mag
                    
                    if mag < rup_rates.get(str(j_rup)).get('Mmax') and  set(rup_rates.get(str(j_rup)).get('involved_faults'))<=set(involved_faults) and  (j_rup != i_rup):
                        add_scenario_to_bin =False
                if add_scenario_to_bin == True:
                    rup_in_bin_i.append(i_rup)
        rup_in_bin.append(rup_in_bin_i)
    return rup_in_bin

def pop_old(bin_mag,faults_names,Mmax_faults,M_min,scenarios_names,index_faults_in_scenario,Mmax_scenario):
    rup_in_bin = [] #for each bin, list the fault and scenarios in it.
    for mag in bin_mag:
        rup_in_bin_i = []
        index_fault = 0
        for fault in faults_names :
            if Mmax_faults[index_fault] >= M_min and mag <= Mmax_faults[index_fault] :
                rup_in_bin_i.append(fault)
            index_fault += 1
        rup_in_bin_ii = []
        index_scenario = 0
        for scenario in scenarios_names :
            nb_of_faults_too_small = 0
            for index_fault in index_faults_in_scenario[index_scenario][0] :
                if mag > Mmax_faults[index_fault] :
                    nb_of_faults_too_small += 1
            # if one of the fault is too small for the bin, the scenario can be used
            if nb_of_faults_too_small >= (len(index_faults_in_scenario[index_scenario][0]) - 1) and mag <= Mmax_scenario[index_scenario] :
                rup_in_bin_ii.append(scenario)
            index_scenario += 1
        #cleaning so large scenario don't participate on small magnitude earthquake if smaller scenarios are available with the same faults
        index_scenario =0
        for scenario_i in rup_in_bin_ii:
            if len(np.where(np.array(faults_names) == scenario_i)[0]) == 0:
                index_scenario_i = np.where(np.array(scenarios_names) == scenario_i)[0]
                index_scenario_i = index_faults_in_scenario[index_scenario_i[0]][0]
                tagged_faults= []
                for scenario_j in rup_in_bin_ii:
                    if len(np.where(np.array(faults_names) == scenario_j)[0]) == 0:
                        index_scenario_j = np.where(np.array(scenarios_names) == scenario_j)[0]
                        index_scenario_j = index_faults_in_scenario[index_scenario_j[0]][0]
                        fault_in = 0
                        if len(index_scenario_i)>len(index_scenario_j):
                            for j in index_scenario_j :
                                if j in index_scenario_i:
                                    fault_in +=1
                                    if not j in tagged_faults:
                                        tagged_faults.append(j)
                if len(tagged_faults) == len(index_scenario_i):
                    if mag == Mmax_scenario[index_scenario] :
                        rup_in_bin_i.append(scenario_i)
                else :
                    rup_in_bin_i.append(scenario_i)
            index_scenario+=1
        
        rup_in_bin.append(rup_in_bin_i)
    return rup_in_bin
