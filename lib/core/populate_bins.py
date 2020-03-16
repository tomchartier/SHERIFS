# -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

Populates the magnitude bins with the faults and scenarios that can generate these magnitudes.

@author: Thomas Chartier
"""
import numpy as np

def pop(bin_mag,faults_names,Mmax_faults,M_min,scenarios_names,index_faults_in_scenario,Mmax_scenario):
    fault_n_scenario_in_bin = [] #for each bin, list the fault and scenarios in it.
    for mag in bin_mag:
        fault_n_scenario_in_bin_i = []
        index_fault = 0
        for fault in faults_names :
            if Mmax_faults[index_fault] >= M_min and mag <= Mmax_faults[index_fault] :
                fault_n_scenario_in_bin_i.append(fault)
            index_fault += 1
        fault_n_scenario_in_bin_ii = []
        index_scenario = 0
        for scenario in scenarios_names :
            nb_of_faults_too_small = 0
            for index_fault in index_faults_in_scenario[index_scenario][0] :
                if mag > Mmax_faults[index_fault] :
                    nb_of_faults_too_small += 1
            # if one of the fault is too small for the bin, the scenario can be used
            if nb_of_faults_too_small >= (len(index_faults_in_scenario[index_scenario][0]) - 1) and mag <= Mmax_scenario[index_scenario] :
                fault_n_scenario_in_bin_ii.append(scenario)
            index_scenario += 1
        #cleaning so large scenario don't participate on small magnitude earthquake if smaller scenarios are available with the same faults
        index_scenario =0
        for scenario_i in fault_n_scenario_in_bin_ii:
            if len(np.where(np.array(faults_names) == scenario_i)[0]) == 0:
                index_scenario_i = np.where(np.array(scenarios_names) == scenario_i)[0]
                index_scenario_i = index_faults_in_scenario[index_scenario_i[0]][0]
                tagged_faults= []
                for scenario_j in fault_n_scenario_in_bin_ii:
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
                        fault_n_scenario_in_bin_i.append(scenario_i)
                else :
                    fault_n_scenario_in_bin_i.append(scenario_i)
            index_scenario+=1
        
        fault_n_scenario_in_bin.append(fault_n_scenario_in_bin_i)
    return fault_n_scenario_in_bin
