#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import xml.etree.ElementTree as ET
import pandas as pd

import matplotlib.path as mplPath
from math import pi, cos, radians

#''' display parameters for pandas '''
#np.set_printoptions(threshold=np.nan)
##pd.set_option('display.multi_sparse', False)
#pd.set_option('display.height', 1000000)
#pd.set_option('display.max_rows', 500000)
#pd.set_option('display.max_columns', 500000)
#pd.set_option('display.width', 1000000)


def Geom_bg(Model_name,File_bg):
    Lon_bg = []
    Lat_bg = []

    # manually defined  in the file Background geometry
    geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
    
    column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
    index_model = np.where(np.array(column_model) == Model_name)[0]
    Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
    Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
    
    return Lon_bg, Lat_bg
    
def FaultGeometry(File_geom,model):       
    #CritereDistance = 3. 
    NomFichier_InfosZonage = File_geom
    InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8')],skip_header = 1)
    Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
    index_model = np.where(np.array(Column_model_name) == model)
    Column_Fault_name = list(map(lambda i : InfosZonage[i][1],index_model[0]))
    Longitudes = list(map(lambda i : InfosZonage[i][2],index_model[0]))
    Latitudes = list(map(lambda i : InfosZonage[i][3],index_model[0]))
    return Column_Fault_name,Longitudes,Latitudes
        
def reproject(latitude, longitude):
    """Returns the x & y coordinates in meters using a sinusoidal projection"""
    earth_radius = 6371009 # in meters
    lat_dist = pi * earth_radius / 180.0
    y = [lat * lat_dist for lat in latitude]
    x = [long * lat_dist * cos(radians(lat)) 
                for lat, long in zip(latitude, longitude)]
    return x, y
    
def area_of_polygon(x, y):
    """Calculates the area of an arbitrary polygon given its verticies"""
    area = 0.0
    for i in range(-1, len(x)-1):
        area += x[i] * (y[i+1] - y[i-1])
    return abs(area) / 2.0
    


   
def do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                     path_for_boxplot,title_for_boxplot,data_on_fault_available,
                                     data_M,data_sig_M,data_rate,data_sig_rate,data_type,sub_area_flag): 
    
    
    for i in range(len(mfd_X)):
        plt.scatter(mega_bining_in_mag,mfd_X[i], c='darkcyan', s=50, edgecolor='',marker = '_',alpha = 0.5)
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    for index_mag in range(len(mega_bining_in_mag)): 
        rate_plus = np.percentile(mfd_X,84,axis=0)[index_mag]
        rate_minus = np.percentile(mfd_X,16,axis=0)[index_mag]
        mag = mega_bining_in_mag[index_mag]
        mag_plus = mag+0.05
        mag_minus = mag-0.05
        verts = [(mag_minus, rate_minus ),
                 (mag_minus, rate_plus),
                 (mag_plus, rate_plus),
                 (mag_plus, rate_minus),
                 (mag_minus, rate_minus)]
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY]
                 
        path_poly = Path(verts, codes)
        
        patch = patches.PathPatch(path_poly,facecolor = '#598556', lw = 0., alpha = 0.15)
        axes.add_patch(patch)
                    
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,50,axis=0),
                c='darkgreen', s=25, edgecolor='',marker = 'o',alpha = 0.8)
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,16,axis=0),
            c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,84,axis=0),
        c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
    plt.plot(mega_bining_in_mag,np.array(mfd_X).mean(axis=0),
                color='darkgreen', linewidth = 2)
    plt.grid()
    
    
    #draw the data on the plot if they are availabla and the plot is cumulative
    if data_on_fault_available == True and sub_area_flag == False:
        for i in range(len(data_M)):
            verts = [(data_M[i] - data_sig_M[i]-0.02, data_rate[i] - data_sig_rate[i]-0.01 * data_rate[i]),
                     (data_M[i] - data_sig_M[i]-0.02, data_rate[i] + data_sig_rate[i]+0.01 * data_rate[i]),
                     (data_M[i] + data_sig_M[i]+0.02, data_rate[i] + data_sig_rate[i]+0.01 * data_rate[i]),
                     (data_M[i] + data_sig_M[i]+0.02, data_rate[i] - data_sig_rate[i]-0.01 * data_rate[i]),
                     (data_M[i] - data_sig_M[i]-0.02, data_rate[i] - data_sig_rate[i]-0.01 * data_rate[i])]
            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY]
                     
            path = Path(verts, codes)
            
            if data_type[i] == 'cat' :
                patch = patches.PathPatch(path,facecolor = 'red', lw = 0.3, alpha = 0.2)
                plt.scatter(data_M[i],data_rate[i],color='red',s=4,marker = 'o',alpha=0.6)
            if data_type[i] == 'pal' :
                patch = patches.PathPatch(path,facecolor = 'indigo', lw = 0.3, alpha = 0.2)
                plt.scatter(data_M[i],data_rate[i],color='indigo',s=4,marker = 'o',alpha=0.6)
            axes.add_patch(patch)
    
    plt.yscale('log')
    plt.title(title_for_boxplot)
    plt.savefig(path_for_boxplot ,dpi = 100, transparent=True)
    plt.close()
        
    file_eq_rate = open(path_for_boxplot[:-3]+'txt','w')
    index_mag=0
    for mag in mega_bining_in_mag:
        file_eq_rate.write(str(mag)+'\t'+str(np.percentile(mfd_X,16,axis=0)[index_mag])+'\t'
                           +str(np.percentile(mfd_X,50,axis=0)[index_mag])+'\t'
                            +str(np.array(mfd_X).mean(axis=0)[index_mag])+'\t'
                            +str(np.percentile(mfd_X,84,axis=0)[index_mag])+'\n')
        index_mag+=1
    file_eq_rate.close()      
    
    
    
def plt_EQ_rates(Run_name,mega_MFD,df_mega_MFD, scenarios_names_list, ScL_complet_list, ScL_list, Model_list,BG_hyp_list,
             dimension_used_list,faults_name_list,sample_list,b_value_list,MFD_type_list,m_Mmax,
             mega_bining_in_mag,a_s_model,b_sample,sm_sample,Mt_sample,plot_mfd,plot_as_rep,plot_Mmax,xmin,xmax,ymin,ymax,
             file_faults_data,File_bg,File_geom,sub_area_file):
    #extract the faults data
    faults_data = np.genfromtxt(file_faults_data,dtype=[('U100000'),('U100000'),('U100000'),('f8'),('f8'),('f8'),('f8')]
                                                             ,delimiter = '\t',skip_header = 1)
    #print faults_data
    data_model = list(map(lambda i : faults_data[i][0], range(len(faults_data))))
    data_fault_name =list( map(lambda i : faults_data[i][1], range(len(faults_data))))
    data_type =list( map(lambda i : faults_data[i][2], range(len(faults_data))))
    data_M =list( map(lambda i : float(faults_data[i][3]), range(len(faults_data))))
    data_sig_M =list( map(lambda i : float(faults_data[i][4]), range(len(faults_data))))
    data_rate = list(map(lambda i : float(faults_data[i][5]), range(len(faults_data))))
    data_sig_rate =list( map(lambda i : float(faults_data[i][6]), range(len(faults_data))))
    
    magnitude_groups = np.linspace(xmin,xmax,int((xmax-xmin)*10)+1) 
    #print df_mega_MFD
    '''############################################
    ###############################################
    #   for each model extract the data
    ###############################################
    ############################################'''
    
    for model in Model_list:
        i_mfd = 0
        while mega_MFD[i_mfd][3] != model :
            i_mfd += 1
            
        
        input_faults_file = (str(Run_name) + '/' + str(mega_MFD[i_mfd][3]) + '/' + 'bg_' + str(mega_MFD[i_mfd][4]) 
        + '/' + str(mega_MFD[i_mfd][0]) + '_' + str(mega_MFD[i_mfd][1]) + '_' + str(mega_MFD[i_mfd][2]) 
        + '/sc_' + str(mega_MFD[i_mfd][8]) + '/faults_n_scenarios.txt')
        
        
        fault_names = np.genfromtxt(input_faults_file,dtype=[('S100000')],delimiter = '\n') #extract from the text file
        if np.size(fault_names) == 1 :
            list_fault_names = str(fault_names)[2::]
            list_fault_names = list_fault_names[:-3]
            list_fault_names = list_fault_names.split(' ')
        else :
            list_fault_names = str(np.array(fault_names[0]))[2::]
            list_fault_names = list_fault_names[:-3]
            list_fault_names = list_fault_names.split(' ') #adapt format to be usable (there is probably a better way to do that)
        
            
        for fault_name in list_fault_names:
            label_for_boxplot = []
            #data_for_boxplot = []
            data_for_boxplot_cum = []

            #find if there is data conserning that fault
            self_data_on_fault_available = False
            self_data_type = []
            self_data_M = []
            self_data_sig_M = []
            self_data_rate = []
            self_data_sig_rate = []
            if fault_name in data_fault_name:

                index_fault_in_data = np.where(np.array(data_fault_name)==fault_name)[0]
                for index_i in index_fault_in_data:
                    if data_model[index_i] == model :
                        self_data_on_fault_available = True
                        self_data_type.append(data_type[index_i])
                        self_data_M.append(data_M[index_i]) 
                        self_data_sig_M.append(data_sig_M[index_i])
                        self_data_rate.append(data_rate[index_i])
                        self_data_sig_rate.append(data_sig_rate[index_i]) 
            
            df_fault_mfd = df_mega_MFD[(df_mega_MFD.Model == model) & (df_mega_MFD.source.str.contains(fault_name))]
            if df_fault_mfd.empty == False:
                df_fault_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                            'b_min','b_max','MFD_type','scenario_set','sample','source',
                                            '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                            '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                            '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                            '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                            '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
                grouped_df_mfd = df_fault_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
                'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
                
                #print 'grouped_df_mfd',grouped_df_mfd
                index_mag = 0    
                for group in magnitude_groups: #loop on the magnitudes
                     
                    reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
                    rec_cum_fault_in_model = grouped_df_mfd[str(round(group,1))].tolist()
                    reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                        
                    
                    if str(group)[-1] == '0' or str(group)[-1] == '5' :
                        label_for_boxplot.append(str(group))
                    else :
                        label_for_boxplot.append(' ')
                    
                    data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                    
                    
                index_mag += 1    
                if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum'):
                    os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum')
                if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name):
                    os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name)
      
                #data_for_boxplot = data_for_boxplot_cum
                path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/'  + 'all_' + fault_name +'.png'     
                title_for_boxplot = 'Frequency of rupture '+ model +' ' +'all_ ' +fault_name+' cumulative rate'
#                box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
#                         self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot)
                grouped_df_mfd=grouped_df_mfd.drop('source',1)
                grouped_df_mfd=grouped_df_mfd.drop('Model',1)
                mfd_X =   grouped_df_mfd.values
                mfd_X = np.array(mfd_X)
                do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                     path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,False)
                
        #for the background    
        label_for_boxplot = []
        #data_for_boxplot = []
        data_for_boxplot_cum = []

        #find if there is data conserning that fault
        self_data_on_fault_available = False
        self_data_type = []
        self_data_M = []
        self_data_sig_M = []
        self_data_rate = []
        self_data_sig_rate = []    
        df_fault_mfd = df_mega_MFD[(df_mega_MFD.Model == model) & (df_mega_MFD.source.str.contains('Background') )]
        if df_fault_mfd.empty == False:
            df_fault_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                        'b_min','b_max','MFD_type','scenario_set','sample','source',
                                        '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                        '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                        '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                        '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                        '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
            grouped_df_mfd = df_fault_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
            'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
            
            #print 'grouped_df_mfd',grouped_df_mfd
            index_mag = 0    
            for group in magnitude_groups: #loop on the magnitudes
                 
                reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
                rec_cum_fault_in_model = grouped_df_mfd[str(round(group,1))].tolist()
                
                
                reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                    
                
                if str(group)[-1] == '0' or str(group)[-1] == '5' :
                    label_for_boxplot.append(str(group))
                else :
                    label_for_boxplot.append(' ')
                
                data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                
                
            index_mag += 1    
            if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum'):
                os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum')
            if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/Background'):
                os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/Background')
  
            #data_for_boxplot = data_for_boxplot_cum
            path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/Background/'  + 'all_Background.png'     
            title_for_boxplot = 'Frequency of rupture '+ model +' ' +'all_ Background cumulative rate'
#            box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
#                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot)
            grouped_df_mfd=grouped_df_mfd.drop('source',1)
            grouped_df_mfd=grouped_df_mfd.drop('Model',1)
            mfd_X =   grouped_df_mfd.values
            mfd_X = np.array(mfd_X)
            do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                     path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,False)

    '''###############################
    ##################################
#    #   for each scenario
#    ##################################
#    ###############################'''


    for model in Model_list:
        for MFD_type in MFD_type_list :  
            for scenario in scenarios_names_list : 
                for fault_name in list_fault_names:
                    if fault_name in data_fault_name:
                        label_for_boxplot = []
                        #data_for_boxplot = []
                        data_for_boxplot_cum = []
            
                        #find if there is data conserning that fault
                        self_data_on_fault_available = False
                        self_data_type = []
                        self_data_M = []
                        self_data_sig_M = []
                        self_data_rate = []
                        self_data_sig_rate = []
        
                        index_fault_in_data = np.where(np.array(data_fault_name)==fault_name)[0]
                        for index_i in index_fault_in_data:
                            if data_model[index_i] == model :
                                self_data_on_fault_available = True
                                self_data_type.append(data_type[index_i])
                                self_data_M.append(data_M[index_i]) 
                                self_data_sig_M.append(data_sig_M[index_i])
                                self_data_rate.append(data_rate[index_i])
                                self_data_sig_rate.append(data_sig_rate[index_i]) 
                    
                        df_fault_mfd = df_mega_MFD[(df_mega_MFD.Model == model)
                        & (df_mega_MFD.source.str.contains(fault_name))
                        & (df_mega_MFD.scenario_set.str.contains(scenario))
                        & (df_mega_MFD.MFD_type.str.contains(MFD_type))]
                        df_fault_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                                    'b_min','b_max','MFD_type','scenario_set','sample','source',
                                                    '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                                    '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                                    '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                                    '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                                    '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
                        grouped_df_mfd = df_fault_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
                        'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
                                    
                        index_mag = 0    
                        for group in magnitude_groups: #loop on the magnitudes
                             
                            reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
                            rec_cum_fault_in_model = grouped_df_mfd[str(round(group,1))].tolist()
                            
                            
                            reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                                
                            
                            if str(group)[-1] == '0' or str(group)[-1] == '5' :
                                label_for_boxplot.append(str(group))
                            else :
                                label_for_boxplot.append(' ')
                            
                            data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                            
                            
                        index_mag += 1    
                        if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/scenario_set'):
                            os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/scenario_set')
              
                        #data_for_boxplot = data_for_boxplot_cum
                        path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/scenario_set/'+'MFD_'+MFD_type +'_'  + scenario + '_' + fault_name +'.png'     
                        title_for_boxplot = 'Frequency of rupture '+ model +' ' +MFD_type+' '+ scenario+' ' +fault_name+' cumulative rate'
    #                    box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
    #                                 self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot)
                        grouped_df_mfd=grouped_df_mfd.drop('source',1)
                        grouped_df_mfd=grouped_df_mfd.drop('Model',1)
                        mfd_X =   grouped_df_mfd.values
                        mfd_X = np.array(mfd_X)
                        do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                         path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                         self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,False)

    '''###############################
    ##################################
    #   for each BG
    ##################################
    ###############################'''

    for model in Model_list:
        if len(BG_hyp_list) > 1:
            for BG_hyp in BG_hyp_list : 
                for fault_name in list_fault_names:
                    if fault_name in data_fault_name:
                        label_for_boxplot = []
                        #data_for_boxplot = []
                        data_for_boxplot_cum = []
            
                        #find if there is data conserning that fault
                        self_data_on_fault_available = False
                        self_data_type = []
                        self_data_M = []
                        self_data_sig_M = []
                        self_data_rate = []
                        self_data_sig_rate = []
            
                        index_fault_in_data = np.where(np.array(data_fault_name)==fault_name)[0]
                        for index_i in index_fault_in_data:
                            if data_model[index_i] == model :
                                self_data_on_fault_available = True
                                self_data_type.append(data_type[index_i])
                                self_data_M.append(data_M[index_i]) 
                                self_data_sig_M.append(data_sig_M[index_i])
                                self_data_rate.append(data_rate[index_i])
                                self_data_sig_rate.append(data_sig_rate[index_i]) 
                    
                        df_fault_mfd = df_mega_MFD[(df_mega_MFD.Model == model)
                        & (df_mega_MFD.source.str.contains(fault_name))
                        & (df_mega_MFD.BG_hyp.str.contains(BG_hyp))]
        
                        df_fault_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                                    'b_min','b_max','MFD_type','scenario_set','sample','source',
                                                    '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                                    '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                                    '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                                    '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                                    '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
                        grouped_df_mfd = df_fault_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
                        'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
                                    
                        index_mag = 0    
                        for group in magnitude_groups: #loop on the magnitudes
                             
                            reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
                            rec_cum_fault_in_model = grouped_df_mfd[str(round(group,1))].tolist()
                            
                            
                            reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                                
                            
                            if str(group)[-1] == '0' or str(group)[-1] == '5' :
                                label_for_boxplot.append(str(group))
                            else :
                                label_for_boxplot.append(' ')
                            
                            data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                            
                            
                        index_mag += 1    
                        if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/BG'):
                            os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/BG')
              
                        #data_for_boxplot = data_for_boxplot_cum
                        path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/' + fault_name + '/BG/'  + BG_hyp + '_' + fault_name +'.png'     
                        title_for_boxplot = 'Frequency of rupture '+ model +' ' + BG_hyp+' ' +fault_name+' cumulative rate'
#                        box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
#                                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot)  
                        grouped_df_mfd=grouped_df_mfd.drop('source',1)
                        grouped_df_mfd=grouped_df_mfd.drop('Model',1)
                        mfd_X =   grouped_df_mfd.values
                        mfd_X = np.array(mfd_X)
                        do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                     path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,False)   
                           
                df_fault_mfd = df_mega_MFD[(df_mega_MFD.Model == model)
                & (df_mega_MFD.source.str.contains('Background'))
                & (df_mega_MFD.BG_hyp.str.contains(BG_hyp))]
                label_for_boxplot = []
                #data_for_boxplot = []
                data_for_boxplot_cum = []
                
                #find if there is data conserning that fault
                self_data_on_fault_available = False
                self_data_type = []
                self_data_M = []
                self_data_sig_M = []
                self_data_rate = []
                self_data_sig_rate = []    
                if df_fault_mfd.empty == False:
                    df_fault_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                                'b_min','b_max','MFD_type','scenario_set','sample','source',
                                                '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                                '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                                '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                                '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                                '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
                    grouped_df_mfd = df_fault_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
                    'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
                    
                    #print 'grouped_df_mfd',grouped_df_mfd
                    index_mag = 0    
                    for group in magnitude_groups: #loop on the magnitudes
                         
                        reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
                        rec_cum_fault_in_model = grouped_df_mfd[str(round(group,1))].tolist()
                        
                        
                        reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                            
                        
                        if str(group)[-1] == '0' or str(group)[-1] == '5' :
                            label_for_boxplot.append(str(group))
                        else :
                            label_for_boxplot.append(' ')
                        
                        data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                        
                        
                    index_mag += 1    
                    if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum'):
                        os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum')
                    if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/Background'+ '/BG'):
                        os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/Background'+ '/BG')
              
                    #data_for_boxplot = data_for_boxplot_cum
                    path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/Background'+ '/BG/'  + BG_hyp + '_' +'_Background.png'     
                    title_for_boxplot = 'Frequency of rupture '+ model+' ' + BG_hyp+' ' +'_ Background cumulative rate'
#                    box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
#                             self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot)
                    grouped_df_mfd=grouped_df_mfd.drop('source',1)
                    grouped_df_mfd=grouped_df_mfd.drop('Model',1)
                    mfd_X =   grouped_df_mfd.values
                    mfd_X = np.array(mfd_X)
                    do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                     path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,False)
                

    '''###############################
    ##################################
    #   plot for the sub areas
    ##################################
    ###############################'''
    for model in Model_list:
        #extract the name and geometry of the faults
        Column_Fault_name,Longitudes,Latitudes = FaultGeometry(File_geom,model)
        fault_names = []
        Lon = []
        Lat = []
        for fault_name in list_fault_names:
            fault_names.append(fault_name)
            index_fault = np.where(np.array(Column_Fault_name)== fault_name)[0]
            Lon.append(np.take(Longitudes,index_fault))
            Lat.append(np.take(Latitudes,index_fault))
        
        # extract the geometry of the zone ( geometry of the background)
        Lon_bg, Lat_bg = Geom_bg(model,File_bg)
        
        #calculate the area of the background
        x,y = reproject(Lat_bg,Lon_bg) 
        area_of_the_bg = area_of_polygon( x,y)  #to be verified!!
        #extract the name of the sources
        df_sources_names =  df_mega_MFD[(df_mega_MFD.Model == model)]
        df_sources_names.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                    'b_min','b_max','MFD_type','scenario_set','sample','source',
                                    '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                    '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                    '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                    '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                    '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
        source_names = np.array(df_sources_names.source.unique())
        
        
        bbPath_sub_areas = []
        if os.path.exists(sub_area_file):
            read_sub_area_file = open(sub_area_file,'rU')
            lines_sub_area = read_sub_area_file.readlines()
            sub_area_names = []
            sub_area_coord = []
            sub_area_lon = []
            sub_area_lat = []
            for line in lines_sub_area:
                model_sub_area = line.split('\t')[0]
                if model == model_sub_area:
                    sub_area_names.append(line.split('\t')[1])
                    sub_area_coord.append(line.split('\t')[2:])
                    sub_area_lon_i = []
                    sub_area_lat_i = []
                    for sub_area_coord_i in line.split('\t')[2:]:
                        if not '\n' in sub_area_coord_i.split(','):
                            if not '' in sub_area_coord_i.split(','):
                                sub_area_lon_i.append(float(sub_area_coord_i.split(',')[1]))
                                sub_area_lat_i.append(float(sub_area_coord_i.split(',')[0]))
                    sub_area_lon.append(sub_area_lon_i)
                    sub_area_lat.append(sub_area_lat_i)
                    if not os.path.exists(str(Run_name) + '/analysis/figures/catalogue/sub_area'):
                        os.makedirs(str(Run_name) + '/analysis/figures/catalogue/sub_area')  
                                 
                    Poly_sub = []   
                    for x1,y1 in zip(sub_area_lon_i,sub_area_lat_i): # creation du polygon de la zone
                        Poly_sub.append((x1,y1))    
                    bbPath_sub_area_i = mplPath.Path(Poly_sub)
                    bbPath_sub_areas.append(bbPath_sub_area_i)
                    # calculate the area of the sub_area
                    x,y = reproject(sub_area_lat_i,sub_area_lon_i) 
                    area_of_the_sub_area = area_of_polygon( x,y)  #to be verified!!
                    
                    faults_in_sub_area = []
                    index_fault = 0
                    
                    for fault_name in fault_names:
                        nb_point_in_sub_area = 0
                        for lon_i,lat_i in zip(Lon[index_fault],Lat[index_fault]):
                            if bbPath_sub_area_i.contains_point((lon_i,lat_i)) == 1: #test to know if the point is in the sub_area
                                nb_point_in_sub_area += 1
                        if nb_point_in_sub_area > len(Lon[index_fault])/2.: #if more than half the points of the trace are in the sub area
                            faults_in_sub_area.append(fault_name)  #the fault is in the sub area
                             
                        index_fault +=1
                    sources_in_sub_area = []
                    ratio_in_subarea = []
                    for source_name_i in source_names:
                        #print 'source_name_i',source_name_i
                        nb_faults_in_source_n_area = 0.
                        for fault_name in faults_in_sub_area:
                            if fault_name in source_name_i:
                                nb_faults_in_source_n_area += 1.
                        if not ']_f_' in source_name_i : #it's a single fault
                            if nb_faults_in_source_n_area >= 1.:
                                sources_in_sub_area.append(source_name_i)
                                ratio_in_subarea.append(1.)
                                #print source_name_i
                        else : 
                            nb_faults_in_source = len(source_name_i.split(']_f_'))
                            if nb_faults_in_source_n_area >= 1.:
                                sources_in_sub_area.append(source_name_i)
                                ratio_in_subarea.append(nb_faults_in_source_n_area/nb_faults_in_source)
                                #print source_name_i,nb_faults_in_source_n_area/nb_faults_in_source
                    df_subarea_mfd = pd.DataFrame(columns=['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                        'b_min','b_max','MFD_type','scenario_set','sample','source',
                                        '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                        '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                        '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                        '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                        '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0'], index = range(len(sources_in_sub_area)*10000))
                    index_source = 0
                    #print sources_in_sub_area
                    for source in sources_in_sub_area:
                        df_source_i_mfd = df_mega_MFD[(df_mega_MFD.Model == model) & (df_mega_MFD.source == source)]
                        if df_source_i_mfd.empty == False:
                            df_source_i_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                                        'b_min','b_max','MFD_type','scenario_set','sample','source',
                                                        '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                                        '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                                        '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                                        '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                                        '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']
#                            print
#                            print
#                            print source,df_source_i_mfd['5.0']
#                            print 'ratio',ratio_in_subarea[index_source]
                            for magnitude in ['4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                            '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                            '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                            '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                            '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']:
                                df_source_i_mfd[magnitude] = df_source_i_mfd[magnitude].astype(float)*ratio_in_subarea[index_source]
                            #print source,df_source_i_mfd['5.0']
                            df_subarea_mfd = pd.concat([df_subarea_mfd,df_source_i_mfd])  
                    
                        index_source+=1
                    #extract the background    
                    df_source_i_mfd = df_mega_MFD[(df_mega_MFD.Model == model) & (df_mega_MFD.source == 'Background')]
                    if df_source_i_mfd.empty == False:
                        df_source_i_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                                    'b_min','b_max','MFD_type','scenario_set','sample','source',
                                                    '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                                    '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                                    '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                                    '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                                    '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0'] 
                        for magnitude in ['4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                        '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                        '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                        '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                        '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0']:
                            df_source_i_mfd[magnitude] = df_source_i_mfd[magnitude].astype(float)*float(area_of_the_sub_area)/float(area_of_the_bg)
                        #print source,df_source_i_mfd['5.0']
                        df_subarea_mfd = pd.concat([df_subarea_mfd,df_source_i_mfd])         
                    grouped_df_mfd = df_subarea_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
                    'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
                    
#                    label_for_boxplot = []
#                    data_for_boxplot = []
#                    data_for_boxplot_cum = []
#                    for group in magnitude_groups: #loop on the magnitudes
#                         
#                        reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
#                        rec_cum_fault_in_model = grouped_df_mfd[str(round(group,1))].tolist()
#                        
#                        
#                        reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                            
                        
#                        if str(group)[-1] == '0' or str(group)[-1] == '5' :
#                            label_for_boxplot.append(str(group))
#                        else :
#                            label_for_boxplot.append(' ')
#                        
#                        data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                        
                    #find if there is data conserning that fault
#                    self_data_on_fault_available = False
#                    self_data_type = []
#                    self_data_M = []
#                    self_data_sig_M = []
#                    self_data_rate = []
#                    self_data_sig_rate = []    
                    if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum'):
                        os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum')
                    if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]):
                        os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1])
          
#                    data_for_boxplot = data_for_boxplot_cum
                    path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]+'/'  + sub_area_names[-1]+'.png'     
                    title_for_boxplot = 'Frequency of rupture '+ model +' ' +'all_ '  + sub_area_names[-1]+' cumulative rate'
#                    box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
#                             self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot,xmin,xmax,ymin,ymax)
                    mfd_X =   grouped_df_mfd.values
                    do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                     path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                     self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,True)
                    
                    
                    #for a bit more detail
                    for MFD_type in MFD_type_list :  
                        for scenario in scenarios_names_list :
                            if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]+ '/' + MFD_type):
                                os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]+ '/' + MFD_type)
                            if not os.path.exists(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]+ '/' + MFD_type+ '/' +scenario):
                                os.makedirs(str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]+ '/' + MFD_type+ '/' +scenario)
                                
                            
                            df_source_i_mfd = df_subarea_mfd[(df_subarea_mfd.scenario_set == scenario) & (df_subarea_mfd.MFD_type == MFD_type)]
                            if df_source_i_mfd.empty == False:
                                df_source_i_mfd.columns = ['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                                            'b_min','b_max','MFD_type','scenario_set','sample','source',
                                                            '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                                            '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                                            '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                                            '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                                            '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0'] 
       
                            grouped_df_mfd = df_source_i_mfd.groupby(['selected_ScL','dim_used','str_all_data','BG_hyp',
                            'b_min','b_max','MFD_type','scenario_set','sample']).sum() 
                        
                        
        #                    data_for_boxplot = data_for_boxplot_cum
                            path_for_boxplot = str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/'+sub_area_names[-1]+ '/' + MFD_type+ '/' +scenario+'/'  + sub_area_names[-1]+'.png'     
                            title_for_boxplot = 'Frequency of rupture '+ model +' ' +scenario+' ' +MFD_type+' '  + sub_area_names[-1]+' cumulative rate'
        #                    box_plot_log(data_for_boxplot,label_for_boxplot,title_for_boxplot,self_data_on_fault_available,
        #                             self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,path_for_boxplot,xmin,xmax,ymin,ymax)
                            mfd_X =   grouped_df_mfd.values
                            do_the_plots(mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                                             path_for_boxplot,title_for_boxplot,self_data_on_fault_available,
                                             self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate,self_data_type,True)
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        