#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import numpy as np
import os

def rup_freq(Run_Name,scenarios_names_list,mega_MFD):    
    magnitude_groups = [[5.,5.4],[5.5,5.9],[6.,6.4],[6.5,6.9],[7.,7.4]]   
    
    for group in magnitude_groups:
        for scenario_name in scenarios_names_list :
            reccurence_source = []
            list_sources = []
            for mega_mfd_i in mega_MFD :
                if mega_mfd_i[9] == scenario_name :
                    index_init = int(group[0]*10.) - 40 + 12
                    index_final = int(group[1]*10.) - 40 + 12
                    if len(list_sources) != 0:  
                        if mega_mfd_i[11] == list_sources[0]:
                            index = 0
                    if not mega_mfd_i[11] in list_sources :
                        list_sources.append(mega_mfd_i[11])
                        reccurence_source_i = [np.sum(mega_mfd_i[index_init::index_final].astype(np.float))]  #indice 12 = M4 ;  indice 33 = M6
                        reccurence_source.append(reccurence_source_i)  #does the sum of the incremental rates for this source                          
                    else :
                        reccurence_source[index].append(np.sum(mega_mfd_i[index_init::index_final].astype(np.float)))                        
                    index +=1
                    
            if not os.path.exists(str(self.Run_name) + '/analysis/figures/rupture_rate'):
                os.makedirs(str(self.Run_name) + '/analysis/figures/rupture_rate')
                        
            file_legend = str(self.Run_name) + '/analysis/figures/rupture_rate/'+ scenario_name+'_mag_'+str(group[0])+'-'+str(group[1])+'.txt'        
            File_legend = open(file_legend,'w')
        
            label_for_boxplot = []
            data_for_boxplot = []
            index = 0
            for source in list_sources :
                File_legend.write(str(index) + '\t' + str(source) + '\t' +
                str(len(reccurence_source[index])) + '\t' +
                str(np.percentile(reccurence_source[index],16)) + '\t' +
                str(np.percentile(reccurence_source[index],50)) + '\t' +
                str(np.mean(reccurence_source[index])) + '\t' +
                str(np.percentile(reccurence_source[index],84)) + '\n')
                
                self.label_for_boxplot.append(str(index))
                self.data_for_boxplot.append(reccurence_source[index])
                index += 1
                
    
            self.path_for_boxplot = str(self.Run_name) + '/analysis/figures/rupture_rate/'+ scenario_name+'_mag_'+str(group[0])+'-'+str(group[1])+'.png'     
            self.title_for_boxplot = 'Frequency of rupture'+ scenario_name+'_mag_'+str(group[0])+'-'+str(group[1])
            box_plot_log()  
            File_legend.close()

'''####################################
# rupture rate fault by fault
#######################################'''
    
if plot_rup_freq == True :
    magnitude_groups = [[5.,5.],[5.1,5.1],[5.2,5.2],[5.3,5.3],[5.4,5.4],
                        [5.5,5.5],[5.6,5.6],[5.7,5.7],[5.8,5.8],[5.9,5.9],
                        [6.,6.],[6.1,6.1],[6.2,6.2],[6.3,6.3],[6.4,6.4],
                        [6.5,6.5],[6.6,6.6],[6.7,6.7],[6.8,6.8],[6.9,6.9],
                        [7.,7.],[7.1,7.1],[7.2,7.2],[7.3,7.3],[7.4,7.4]]
    
    for scenario_name in scenarios_names_list :
        
        #extract the list of faults in this model
        i_mfd = 0
        while mega_MFD[i_mfd][9] != scenario_name :
            i_mfd += 1
            
        input_faults_file = (str(self.Run_name) + '/' + str(mega_MFD[i_mfd][0]) + '_' + str(mega_MFD[i_mfd][1]) + '_' + str(mega_MFD[i_mfd][2]) 
        + '/' + str(mega_MFD[i_mfd][3]) + '/' + str(scenario_name) + '/faults_n_scenarios.txt')
        
        data_fault_names = np.genfromtxt(input_faults_file,dtype=[('S10000')],delimiter = '\n') #extract from the text file
        if np.size(data_fault_names) == 1 :
            list_fault_names = str(data_fault_names)[2::]
            list_fault_names = list_fault_names[:-3]
            list_fault_names = list_fault_names.split(' ')
        else :
            list_fault_names = str(np.array(data_fault_names[0]))[2::]
            list_fault_names = list_fault_names[:-3]
            list_fault_names = list_fault_names.split(' ') #adapt format to be usable (there is probably a better way to do that)
        
        for fault_name in list_fault_names:
            self.label_for_boxplot = []
            self.data_for_boxplot = []
            data_for_boxplot_cum = []
            for group in magnitude_groups: #loop on the magnitudes
                reccurence_fault_mag = []  #frequency this fault produce this magnitude
                reccurence_cum_fault_mag = []  #frequency this fault produce this magnitude cumulative
                list_sources = [] # list of the sources the fault is involved in (to avoid double counting)
                for mega_mfd_i in mega_MFD : # loop on the MFD of all the sources of all the model
                    if len(list_sources) != 0:  
                        if mega_mfd_i[11] == list_sources[0]:
                            if rec_fault_in_model != 0. :
                                reccurence_fault_mag.append(rec_fault_in_model)
                            if rec_cum_fault_in_model != 0. :
                                reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                            rec_fault_in_model = 0.  #rec incremental
                            rec_cum_fault_in_model = 0. #rec_cumulative
                            
                    else :
                        rec_fault_in_model = 0.
                        rec_cum_fault_in_model = 0.
                        
                    if mega_mfd_i[9] == scenario_name :
                        index_init = int(group[0]*10.) - 40 + 12
                        index_final = int(group[1]*10.) - 40 + 12
                        if fault_name in str(mega_mfd_i[11]): 
                            if len(list_sources) != 0:  
                                if mega_mfd_i[11] == list_sources[0]:
                                    index = 0
                            if not mega_mfd_i[11] in list_sources :
                                list_sources.append(mega_mfd_i[11])
                                reccurence_source_i = [np.sum(mega_mfd_i[index_init::index_final].astype(np.float))]  #indice 12 = M4 ;  indice 33 = M6

                                rec_fault_in_model += reccurence_source_i[0]

                                reccurence_cum_source_i = [np.sum(mega_mfd_i[index_init:].astype(np.float))]  #indice 12 = M4 ;  indice 33 = M6
                                rec_cum_fault_in_model += reccurence_cum_source_i[0] 
                            else :
                                rec_fault_in_model += np.sum(mega_mfd_i[index_init::index_final].astype(np.float)) 
                                rec_cum_fault_in_model += np.sum(mega_mfd_i[index_init:].astype(np.float)) 
                            
                            index +=1
                            
                reccurence_fault_mag.append(rec_fault_in_model)
                
                reccurence_cum_fault_mag.append(rec_cum_fault_in_model)
                    
                
                if str(group[-1])[-1] == '0' or str(group[-1])[-1] == '5' :
                    self.label_for_boxplot.append(str(group[0]) + '\n' + str(len(reccurence_fault_mag)))
                else :
                    self.label_for_boxplot.append(' ')
                
                #feeding the matrice for the boxplot
                self.data_for_boxplot.append(reccurence_fault_mag)
                data_for_boxplot_cum.append(reccurence_cum_fault_mag)
                
            if not os.path.exists(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_inc'):
                os.makedirs(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_inc')
            if not os.path.exists(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_inc/' + fault_name):
                os.makedirs(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_inc/' + fault_name)
                
            self.path_for_boxplot = str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_inc/' + fault_name + '/' + scenario_name + '_' + fault_name +'.png'     
            self.title_for_boxplot = 'Frequency of rupture'+ scenario_name+' ' +fault_name+' incremental rate'
            self.box_plot_log()
                
            if not os.path.exists(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum'):
                os.makedirs(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum')
            if not os.path.exists(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + fault_name):
                os.makedirs(str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + fault_name)
  
            self.data_for_boxplot = data_for_boxplot_cum
            self.path_for_boxplot = str(self.Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + fault_name + '/' + scenario_name + '_' + fault_name +'.png'     
            self.title_for_boxplot = 'Frequency of rupture'+ scenario_name+' ' +fault_name+' cumulative rate'
            self.box_plot_log()
        
    '''####################################
    # use of the slip rate per fault
    #######################################'''  
    time_i = time.time()
    input_file_name = 'init'
    slip_rep_faults_all_data = open(self.Run_name + '/analysis/txt_files/slip_rep_on_faults_all_data.txt','w') #file with all the information about slip repartition
    if plot_rup_freq == True :
        for scenario_name in scenarios_names_list :
            for mega_mfd_i in mega_MFD :
                if mega_mfd_i[9] == scenario_name:
                    input_file_name_i = (str(self.Run_name) + '/' + str(mega_mfd_i[0]) + '_' + str(mega_mfd_i[1]) + '_' + str(mega_mfd_i[2]) 
                    + '/' + str(mega_mfd_i[3]) + '/' + str(scenario_name)  + '/Mt_' + str(mega_mfd_i[4]) + '/bmin_' + str(mega_mfd_i[5]) + '_bmax_' + str(mega_mfd_i[6]) 
                    + '/sm_' + str(mega_mfd_i[7]) + '/as_' +str(mega_mfd_i[8]) + '/Log/sliprep_sample_' + str(mega_mfd_i[10]) + '.txt')
                    if input_file_name_i != input_file_name :
                        input_file_name = input_file_name_i
                        input_file = np.genfromtxt(input_file_name,dtype = 'S1000000', delimiter = '], [')
                        for index_fault in range(len(input_file)):
                            if index_fault == 0 :
                                fault_name = input_file[index_fault][3:(input_file[index_fault].find(',')-1)]
                            else :
                                fault_name = input_file[index_fault][1:(input_file[index_fault].find(',')-1)]

                            total_number = input_file[index_fault].count('aseismic_slip') + input_file[index_fault].count(fault_name)
                            n_FtF_6 = input_file[index_fault].count('f_6')
                            n_FtF_5 = input_file[index_fault].count('f_5') - input_file[index_fault].count('f_6')
                            n_FtF_4 = input_file[index_fault].count('f_4') - input_file[index_fault].count('f_5')
                            n_FtF_3 = input_file[index_fault].count('f_3') - input_file[index_fault].count('f_4')
                            n_FtF_2 = input_file[index_fault].count('f_2') - input_file[index_fault].count('f_3')
                            n_fault_alone = input_file[index_fault].count(fault_name) - input_file[index_fault].count('f_2')
                            n_NMS = input_file[index_fault].count('aseismic_slip')
                            
                            p_FtF_6 = round(float(n_FtF_6) / float(total_number) * 100., 1)
                            p_FtF_5 = round(float(n_FtF_5) / float(total_number) * 100., 1)
                            p_FtF_4 = round(float(n_FtF_4) / float(total_number) * 100., 1)
                            p_FtF_3 = round(float(n_FtF_3) / float(total_number) * 100., 1)
                            p_FtF_2 = round(float(n_FtF_2) / float(total_number) * 100., 1)
                            p_fault_alone = round(float(n_fault_alone) / float(total_number) * 100., 1)
                            p_NMS = round(float(n_NMS) / float(total_number) * 100., 1)
                            
                            line = (str(mega_mfd_i[0]) + '_' + str(mega_mfd_i[1]) + '_' + str(mega_mfd_i[2]) 
                            + '\t' + str(mega_mfd_i[3]) + '\tMt_' + str(mega_mfd_i[4]) + '\tbmin_' + str(mega_mfd_i[5]) + '_bmax_' + str(mega_mfd_i[6]) 
                            + '\tsm_' + str(mega_mfd_i[7]) + '\tas_' +str(mega_mfd_i[8]) + '\t' + str(scenario_name)  + '\tsample_' + str(mega_mfd_i[10])
                            + '\t' + fault_name + '\t' + str(p_fault_alone) + '\t' + str(p_FtF_2) + '\t' + str(p_FtF_3) + '\t' + str(p_FtF_4) + '\t' 
                            + str(p_FtF_5) + '\t' + str(p_FtF_6) + '\t' + str(p_NMS) )
                            slip_rep_faults_all_data.write(line+'\n')
    slip_rep_faults_all_data.close()
                        
    if plot_rup_freq == True :
        slip_rep_data = np.genfromtxt(self.Run_name + '/analysis/txt_files/slip_rep_on_faults_all_data.txt',
                                      dtype = [('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),
                                               ('S100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t') 
        scenario_names_sli_rep = map(lambda i : slip_rep_data[i][6], range(len(slip_rep_data)))
        p_fault_alone = map(lambda i : slip_rep_data[i][9], range(len(slip_rep_data)))
        p_FtF_2 = map(lambda i : slip_rep_data[i][10], range(len(slip_rep_data)))
        p_FtF_3 = map(lambda i : slip_rep_data[i][11], range(len(slip_rep_data)))
        p_FtF_4 = map(lambda i : slip_rep_data[i][12], range(len(slip_rep_data)))
        p_FtF_5 = map(lambda i : slip_rep_data[i][13], range(len(slip_rep_data)))
        p_FtF_6 = map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data)))
        p_NMS = map(lambda i : slip_rep_data[i][15], range(len(slip_rep_data)))
        list_faults_slip_rep = map(lambda i : slip_rep_data[i][8], range(len(slip_rep_data)))
                    
        for scenario_name in scenarios_names_list :
            index = np.where(np.array(scenario_names_sli_rep) == scenario_name)[0]
            p_fault_alone_i = np.take(p_fault_alone,index)
            p_FtF_2_i = np.take(p_FtF_2,index)
            p_FtF_3_i = np.take(p_FtF_3,index)
            p_FtF_4_i = np.take(p_FtF_4,index)
            p_FtF_5_i = np.take(p_FtF_5,index)
            p_FtF_6_i = np.take(p_FtF_6,index)
            p_NMS_i = np.take(p_NMS,index)
            list_faults_slip_rep_i = np.take(list_faults_slip_rep,index)
                
            
            list_fault_i = []
            for fault in list_faults_slip_rep_i :
                if fault not in list_fault_i :
                    list_fault_i.append(fault)
                    
            slip_rep_faults_mean = open(self.Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_' + scenario_name +'.txt','w') #file with the mean slip repartition for each fault in this scenario

            for fault in list_fault_i :
                index = np.where(list_faults_slip_rep_i == fault)[0]
                p_fault_alone_j = np.take(p_fault_alone_i,index)
                p_FtF_2_j = np.take(p_FtF_2_i,index)
                p_FtF_3_j = np.take(p_FtF_3_i,index)
                p_FtF_4_j = np.take(p_FtF_4_i,index)
                p_FtF_5_j = np.take(p_FtF_5_i,index)
                p_FtF_6_j = np.take(p_FtF_6_i,index)
                p_NMS_j = np.take(p_NMS_i,index)
                slip_rep_faults_mean.write(fault + '\t' + str(np.mean(p_fault_alone_j)) + '\t' + str(np.mean(p_FtF_2_j)) 
                + '\t' + str(np.mean(p_FtF_3_j)) + '\t' + str(np.mean(p_FtF_4_j)) + '\t' + str(np.mean(p_FtF_5_j)) 
                + '\t' + str(np.mean(p_FtF_6_j)) + '\t' + str(np.mean(p_NMS_j)) +'\n')
            
            
            slip_rep_faults_mean.close()

    print '\nTime to see how the slip rate in distributed : ' + str(time.time() - time_i) +' s.\n'