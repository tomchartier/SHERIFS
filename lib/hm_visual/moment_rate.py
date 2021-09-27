#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import os
import numpy as np
import matplotlib.pyplot as plt

def box_plot(data_for_boxplot,label_for_boxplot,title_for_boxplot,path_for_boxplot):
           
    data = data_for_boxplot
    label = label_for_boxplot
    
    fs = 9  # fontsize
   
    plt.boxplot(data, labels=label, showmeans=True)
    plt.title(title_for_boxplot, fontsize=fs)
    plt.grid()
    plt.savefig(path_for_boxplot,dpi = 180, transparent=True)
    #plt.show()
    plt.close()
        
        
        
def moment_rate(Run_name,plot_moment_rate,geologic_moment_rate_no_as,geologic_moment_rate,
                seismological_moment_rate,scenarios_names_list,total_list_scenario_name,
                MFD_type_list,total_list_MFD_type):

    
    '''#############################
    ###############################
    #     comparison of the moment rate
    #   geologic, geodetic, sismologic
    ###############################
    ##############################'''
    if plot_moment_rate == True:
        name_list = ['Potential moment rate\nof the fault system',
        'Modelled moment rate\n(NMS) taken into account',
        'Moment rate calculated\nfrom the earthquake catalog']
        data_for_boxplot = []
        label_for_boxplot = name_list
        data_for_boxplot.append(geologic_moment_rate_no_as)
        data_for_boxplot.append(geologic_moment_rate)
        data_for_boxplot.append(seismological_moment_rate)
#            data_for_boxplot.append(geodetic_moment_rate)
#            data_for_boxplot.append(geodetic_moment_rate_as)
    

        if not os.path.exists(str(Run_name) + '/analysis/figures/compare_moment_rate'):
            os.makedirs(str(Run_name) + '/analysis/figures/compare_moment_rate')

        path_for_boxplot = str(Run_name) + '/analysis/figures/compare_moment_rate/M0.png'     
        title_for_boxplot = 'Compare M0'
        box_plot(data_for_boxplot,label_for_boxplot,title_for_boxplot,path_for_boxplot)    
        
        
    '''####################################
    #plot moment rate for each set of rupture scenario used
    #######################################'''   
    if plot_moment_rate == True:
        data_for_boxplot = []
        label_for_boxplot = []
        for scenario_name in scenarios_names_list :
            label_for_boxplot.append(scenario_name)
            label_for_boxplot.append(scenario_name + '_no_NMS')
            
            M0_scenario_name = []
            rows = np.where(np.array(total_list_scenario_name) == scenario_name)[0]
            M0_no_as_scenario_name = []
            for index in rows :
                M0_scenario_name.append(geologic_moment_rate[index])
                M0_no_as_scenario_name.append(geologic_moment_rate_no_as[index])
            data_for_boxplot.append(M0_scenario_name)
            data_for_boxplot.append(M0_no_as_scenario_name)
            
        path_for_boxplot = str(Run_name) + '/analysis/figures/compare_moment_rate/M0_each_model_sce_rup.png'     
        title_for_boxplot = 'Compare M0 - each model \nof rupture scenario'
        box_plot(data_for_boxplot,label_for_boxplot,title_for_boxplot,path_for_boxplot) 
     
        
        
    '''####################################
    #plot moment rate for each shape of MFD
    #######################################'''   
    if plot_moment_rate == True:
        data_for_boxplot = []
        label_for_boxplot = []
        for mfd in MFD_type_list :
            label_for_boxplot.append(mfd)
            label_for_boxplot.append(mfd + '_no_NMS')
            
            M0_mfd = []
            rows = np.where(np.array(total_list_MFD_type) == mfd)[0]
            M0_no_as_mfd = []
            for index in rows :
                M0_mfd.append(geologic_moment_rate[index])
                M0_no_as_mfd.append(geologic_moment_rate_no_as[index])
            data_for_boxplot.append(M0_mfd)
            data_for_boxplot.append(M0_no_as_mfd)
            
        path_for_boxplot = str(Run_name) + '/analysis/figures/compare_moment_rate/M0_each_mfd.png'     
        title_for_boxplot = 'Compare M0 - each model \nof MFD'
        box_plot(data_for_boxplot,label_for_boxplot,title_for_boxplot,path_for_boxplot) 
     
