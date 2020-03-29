#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import numpy as np
        
        
def sr_rate(Run_name,scenarios_names_list,mega_MFD,Model_list,MFD_type_list):

    '''####################################
    # use of the slip rate per fault
    #######################################'''  
    input_file_name = 'init'
    slip_rep_faults_all_data = open(Run_name + '/analysis/txt_files/slip_rep_on_faults_all_data.txt','w') #file with all the information about slip repartition

    for model in Model_list:
        for scenario_name in scenarios_names_list :
            for mega_mfd_i in mega_MFD :
                if mega_mfd_i[8] == scenario_name and mega_mfd_i[3] == model :
                    i_mfd =  0
                    input_file_name_i = (str(Run_name) + '/' + str(mega_mfd_i[3]) + '/' + 'bg_' + str(mega_mfd_i[4]) 
                                        + '/' + str(mega_mfd_i[0]) + '_' + str(mega_mfd_i[1]) + '_' + str(mega_mfd_i[2]) 
                                        + '/sc_' + str(mega_mfd_i[8])  + '/bmin_' + str(mega_mfd_i[5]) + '_bmax_' + str(mega_mfd_i[6]) 
                                        + '/MFD_' + str(mega_mfd_i[7]) + '/Log/sliprep_sample_' + str(mega_mfd_i[9]) + '.txt')
                                        
                    rup_file = (str(Run_name) + '/' + str(mega_mfd_i[3]) + '/' + 'bg_' + str(mega_mfd_i[4])
                    + '/' + str(mega_mfd_i[0]) + '_' + str(mega_mfd_i[1]) + '_' + str(mega_mfd_i[2])
                    + '/sc_' + str(mega_mfd_i[8])  + '/bmin_' + str(mega_mfd_i[5]) + '_bmax_' + str(mega_mfd_i[6])
                    + '/MFD_' + str(mega_mfd_i[7]) + '/Log/ruptures.txt')
                    rup_id = []
                    rup_length = []
                    for line in open(rup_file):
                        if not 'rup_id' in line :
                            rup_i = line.split('\t')[0]
                            rup_id.append(rup_i)
                            rup_length_i = line.split('\t')[1]
                            rup_length.append(len(rup_length_i.split(' '))-1)
                        
                    if input_file_name_i != input_file_name :
                        input_file_name = input_file_name_i
                        with open(input_file_name,'r') as f:
                            for line in f:
                                content = line.split(' ')
                                fault_name = content[0]
                                total_number = len(content)-1
                                n_fault_alone= 0.
                                n_FtF_2= 0.
                                n_FtF_3= 0.
                                n_FtF_4= 0.
                                n_FtF_5= 0.
                                n_FtF_6= 0.
                                n_FtF_7= 0.
                                n_FtF_8= 0.
                                n_FtF_9= 0.
                                n_FtF_10= 0.
                                n_FtF_11= 0.
                                n_FtF_12= 0.
                                n_FtF_13 = 0.
                                for rup_i,rup_length_i in zip(rup_id,rup_length):
                                    c = content[1:].count(str(rup_i))
                                    if rup_length_i >= 13 :
                                        n_FtF_13 += c
                                    if rup_length_i == 12 :
                                        n_FtF_12 += c
                                    if rup_length_i == 11 :
                                        n_FtF_11 += c
                                    if rup_length_i == 10 :
                                        n_FtF_10 += c
                                    if rup_length_i == 9 :
                                        n_FtF_9 += c
                                    if rup_length_i == 8 :
                                        n_FtF_8 += c
                                    if rup_length_i == 7 :
                                        n_FtF_7 += c
                                    if rup_length_i == 6 :
                                        n_FtF_6 += c
                                    if rup_length_i == 5 :
                                        n_FtF_5 += c
                                    if rup_length_i == 4 :
                                        n_FtF_4 += c
                                    if rup_length_i == 3 :
                                        n_FtF_3 += c
                                    if rup_length_i == 2 :
                                        n_FtF_2 += c
                                    if rup_length_i == 1 :
                                        n_fault_alone += c
                                n_NMS = content[1:].count('NMS')
                                
                                p_FtF_13 = round(float(n_FtF_13) / float(total_number) * 100., 1)
                                p_FtF_12 = round(float(n_FtF_12) / float(total_number) * 100., 1)
                                p_FtF_11 = round(float(n_FtF_11) / float(total_number) * 100., 1)
                                p_FtF_10 = round(float(n_FtF_10) / float(total_number) * 100., 1)
                                p_FtF_9 = round(float(n_FtF_9) / float(total_number) * 100., 1)
                                p_FtF_8 = round(float(n_FtF_8) / float(total_number) * 100., 1)
                                p_FtF_7 = round(float(n_FtF_7) / float(total_number) * 100., 1)
                                p_FtF_6 = round(float(n_FtF_6) / float(total_number) * 100., 1)
                                p_FtF_5 = round(float(n_FtF_5) / float(total_number) * 100., 1)
                                p_FtF_4 = round(float(n_FtF_4) / float(total_number) * 100., 1)
                                p_FtF_3 = round(float(n_FtF_3) / float(total_number) * 100., 1)
                                p_FtF_2 = round(float(n_FtF_2) / float(total_number) * 100., 1)
                                p_fault_alone = round(float(n_fault_alone) / float(total_number) * 100., 1)
                                p_NMS = round(float(n_NMS) / float(total_number) * 100., 1)
                                
                                line = (str(mega_mfd_i[0]) + '_' + str(mega_mfd_i[1]) + '_' + str(mega_mfd_i[2])
                                + '\t' + str(model) + '\tbg_' + str(mega_mfd_i[4]) + '\tbmin_' + str(mega_mfd_i[5]) + '_bmax_' + str(mega_mfd_i[6])
                                + '\tMFD_' + str(mega_mfd_i[7]) + '\t' + str(scenario_name)  + '\tsample_' + str(mega_mfd_i[9])
                                + '\t' + fault_name + '\t' + str(p_fault_alone) + '\t' + str(p_FtF_2) + '\t' + str(p_FtF_3) + '\t' + str(p_FtF_4) + '\t'
                                + str(p_FtF_5) + '\t' + str(p_FtF_6) + '\t' + str(p_FtF_7) + '\t'+ str(p_FtF_8) + '\t'+ str(p_FtF_9) + '\t'+ str(p_FtF_10) + '\t'
                                + str(p_FtF_11) + '\t'+ str(p_FtF_12) + '\t'+ str(p_FtF_13) + '\t'+ str(p_NMS) )
                                slip_rep_faults_all_data.write(line+'\n')
    slip_rep_faults_all_data.close()
                        
    slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_all_data.txt',
                                  dtype = [('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),
                                           ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t') 
    model_slip_rep = list(map(lambda i : slip_rep_data[i][1].decode("utf-8"), range(len(slip_rep_data))))
    MFD_type_sli_rep = list(map(lambda i : slip_rep_data[i][4].decode("utf-8"), range(len(slip_rep_data))))
    scenario_names_sli_rep = list(map(lambda i : slip_rep_data[i][5].decode("utf-8"), range(len(slip_rep_data))))
    p_fault_alone = list(map(lambda i : slip_rep_data[i][8], range(len(slip_rep_data))))
    p_FtF_2 = list(map(lambda i : slip_rep_data[i][9], range(len(slip_rep_data))))
    p_FtF_3 = list(map(lambda i : slip_rep_data[i][10], range(len(slip_rep_data))))
    p_FtF_4 = list(map(lambda i : slip_rep_data[i][11], range(len(slip_rep_data))))
    p_FtF_5 = list(map(lambda i : slip_rep_data[i][12], range(len(slip_rep_data))))
    p_FtF_6 = list(map(lambda i : slip_rep_data[i][13], range(len(slip_rep_data))))
    p_FtF_7 = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))
    p_FtF_8 = list(map(lambda i : slip_rep_data[i][15], range(len(slip_rep_data))))
    p_FtF_9 = list(map(lambda i : slip_rep_data[i][16], range(len(slip_rep_data))))
    p_FtF_10 = list(map(lambda i : slip_rep_data[i][17], range(len(slip_rep_data))))
    p_FtF_11 = list(map(lambda i : slip_rep_data[i][18], range(len(slip_rep_data))))
    p_FtF_12 = list(map(lambda i : slip_rep_data[i][19], range(len(slip_rep_data))))
    p_FtF_13 = list(map(lambda i : slip_rep_data[i][20], range(len(slip_rep_data))))
    p_NMS = list(map(lambda i : slip_rep_data[i][21], range(len(slip_rep_data))))
    list_faults_slip_rep = list(map(lambda i : slip_rep_data[i][7].decode("utf-8"), range(len(slip_rep_data))))
    
    for model in Model_list:            
        for scenario_name in scenarios_names_list :
            for MFD_type in MFD_type_list:
                index_sc = np.where(np.array(scenario_names_sli_rep) == scenario_name)[0]
                index_mfd_type = np.where(np.array(MFD_type_sli_rep) == ('MFD_'+MFD_type))[0]
                index_m = np.where(np.array(model_slip_rep) == model)[0]
                index = list(set(index_sc).intersection(index_m))
                index = list(set(index).intersection(index_mfd_type))
                p_fault_alone_i = np.take(p_fault_alone,index)
                p_FtF_2_i = np.take(p_FtF_2,index)
                p_FtF_3_i = np.take(p_FtF_3,index)
                p_FtF_4_i = np.take(p_FtF_4,index)
                p_FtF_5_i = np.take(p_FtF_5,index)
                p_FtF_6_i = np.take(p_FtF_6,index)
                p_FtF_7_i = np.take(p_FtF_7,index)
                p_FtF_8_i = np.take(p_FtF_8,index)
                p_FtF_9_i = np.take(p_FtF_9,index)
                p_FtF_10_i = np.take(p_FtF_10,index)
                p_FtF_11_i = np.take(p_FtF_11,index)
                p_FtF_12_i = np.take(p_FtF_12,index)
                p_FtF_13_i = np.take(p_FtF_13,index)
                p_NMS_i = np.take(p_NMS,index)
                list_faults_slip_rep_i = np.take(list_faults_slip_rep,index)
                    
                
                list_fault_i = []
                for fault in list_faults_slip_rep_i :
                    if fault not in list_fault_i :
                        list_fault_i.append(fault)
                        
                slip_rep_faults_mean = open(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+model+'_' + MFD_type +'_' + scenario_name +'.txt','w') #file with the mean slip repartition for each fault in this scenario
    
                for fault in list_fault_i :
                    index = np.where(list_faults_slip_rep_i == fault)[0]
                    p_fault_alone_j = np.take(p_fault_alone_i,index)
                    p_FtF_2_j = np.take(p_FtF_2_i,index)
                    p_FtF_3_j = np.take(p_FtF_3_i,index)
                    p_FtF_4_j = np.take(p_FtF_4_i,index)
                    p_FtF_5_j = np.take(p_FtF_5_i,index)
                    p_FtF_6_j = np.take(p_FtF_6_i,index)
                    p_FtF_7_j = np.take(p_FtF_7_i,index)
                    p_FtF_8_j = np.take(p_FtF_8_i,index)
                    p_FtF_9_j = np.take(p_FtF_9_i,index)
                    p_FtF_10_j = np.take(p_FtF_10_i,index)
                    p_FtF_11_j = np.take(p_FtF_11_i,index)
                    p_FtF_12_j = np.take(p_FtF_12_i,index)
                    p_FtF_13_j = np.take(p_FtF_13_i,index)
                    p_NMS_j = np.take(p_NMS_i,index)
                    slip_rep_faults_mean.write(fault + '\t' + str(np.mean(p_fault_alone_j)) + '\t' + str(np.mean(p_FtF_2_j)) 
                    + '\t' + str(np.mean(p_FtF_3_j)) + '\t' + str(np.mean(p_FtF_4_j)) + '\t' + str(np.mean(p_FtF_5_j)) 
                    + '\t' + str(np.mean(p_FtF_6_j)) + '\t' + str(np.mean(p_FtF_7_j)) + '\t' + str(np.mean(p_FtF_8_j)) 
                    + '\t' + str(np.mean(p_FtF_9_j)) + '\t' + str(np.mean(p_FtF_10_j)) + '\t' + str(np.mean(p_FtF_11_j)) 
                    + '\t' + str(np.mean(p_FtF_12_j)) + '\t' + str(np.mean(p_FtF_13_j)) + '\t' + str(np.mean(p_NMS_j)) +'\n')
                
                
                slip_rep_faults_mean.close()
