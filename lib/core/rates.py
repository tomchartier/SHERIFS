# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np

def get_rate_model(rup_rates,fault_prop,bin_mag):
    # calculates the earthquake rate in the model (faults+background)
    rate_f_in_model = []
    for dict in rup_rates:
        rate_f_in_model.append(rup_rates.get(dict).get('rates'))
    rate_f_in_model = np.sum(rate_f_in_model,axis = 0)
    rate_bg_in_model=np.zeros(len(bin_mag)) #rate of EQ in the BG for each magnitude
    for index_mag in range(len(bin_mag)) :
        rate_bg_in_model[index_mag] += ((1-fault_prop(bin_mag[index_mag])) * rate_f_in_model[index_mag]) / fault_prop(bin_mag[index_mag]) #earthquake incremental rate of the BG in the model at the time of the calculation
    rate_in_model = rate_f_in_model + rate_bg_in_model
    return rate_in_model

def get_rate_faults_n_bg(rup_rates,fault_prop,bin_mag):
    # calculates the earthquake rate in the model (faults+background)
    rate_f_in_model = []
    for dict in rup_rates:
        rate_f_in_model.append(rup_rates.get(dict).get('rates'))
    rate_f_in_model = np.sum(rate_f_in_model,axis = 0)
    rate_bg_in_model=np.zeros(len(bin_mag)) #rate of EQ in the BG for each magnitude
    for index_mag in range(len(bin_mag)) :
        rate_bg_in_model[index_mag] += ((1-fault_prop(bin_mag[index_mag])) * rate_f_in_model[index_mag]) / fault_prop(bin_mag[index_mag]) #earthquake incremental rate of the BG in the model at the time of the calculation
    return rate_f_in_model, rate_bg_in_model
