# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np


def select(sr_values,sample,index_fault,linked_level,list_quater_picked):
    sr_min = sr_values[0]
    sr_mean = sr_values[1]
    sr_max = sr_values[2]
    #linked level (:is the fault very connected with others?)
    
    ########################################################
    #ramdom sampling the slip rate  . uniform sampling
    ########################################################
    if sample == 1 :
        slip_rate = sr_mean
    else :
        if sum(linked_level)==0 : #the fault is alone
            slip_rate_inf = np.random.uniform(sr_min,sr_mean)
            slip_rate_sup = np.random.uniform(sr_mean,sr_max)
            slip_rate = np.random.choice([slip_rate_inf,slip_rate_sup])
        else: #the fault is connected
            value_lvl = 10
            quarters_picked = []
            for index_c in range(len(linked_level)):
                if linked_level[index_c] != 0:
                    if list_quater_picked[index_c]!=0:
                        if linked_level[index_c] < value_lvl :
                            quarters_picked = []
                            value_lvl = linked_level[index_c]
                            #print faults_names[index_fault],faults_names[index_c]
                        if linked_level[index_c] == value_lvl : quarters_picked.append(list_quater_picked[index_c])
            if quarters_picked == []: #none of the faults have bin picked yet
                slip_rate_inf = np.random.uniform(sr_min,sr_mean)
                slip_rate_sup = np.random.uniform(sr_mean,sr_max)
                slip_rate = np.random.choice([slip_rate_inf,slip_rate_sup])
                if slip_rate < (sr_min + 1./2. * (sr_mean-sr_min)):
                    quarter_to_pick = 1
                elif slip_rate < (sr_mean):
                    quarter_to_pick = 2
                elif slip_rate < (sr_mean + 1./2. * (sr_max-sr_mean)):
                    quarter_to_pick = 3
                else:
                    quarter_to_pick = 4
                    
            else :
                quarter_to_pick = max(set(quarters_picked),key=quarters_picked.count)
                if quarter_to_pick == 1 :
                    slip_rate = np.random.uniform(sr_min,sr_min + 1./2. * (sr_mean-sr_min) +0.000001)
                elif quarter_to_pick == 2 :
                    slip_rate = np.random.uniform(sr_min + 1./2. * (sr_mean-sr_min),sr_mean+0.000001)
                elif quarter_to_pick == 3 :
                    slip_rate = np.random.uniform(sr_mean,sr_mean + 1./2. * (sr_max-sr_mean)+0.000001)
                elif quarter_to_pick == 4 :
                    slip_rate = np.random.uniform(sr_mean + 1./2. * (sr_max-sr_mean),sr_max+0.000001)
                
                    
            list_quater_picked[index_fault] = quarter_to_pick
                    
    return slip_rate
