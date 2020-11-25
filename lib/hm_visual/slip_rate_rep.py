#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import numpy as np
import pickle
import os
import geojson
from subarea_utils import find_faults_in_sub,get_geom
from Participation_rates import FaultGeometry
from Read_file import read_file_Mmax_log
import time


def get_shearmod(File_prop,Name_of_fault,Model):
    # TODO Clean up this and put it in the read file stuff
    if not ".geojson" in File_prop:
        FileName_Prop = File_prop
        Prop = np.genfromtxt(FileName_Prop,
                                   dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
                                          ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
        Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
        Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
        index_model = np.where(np.array(Column_model_name) == Model)[0]
        
        Prop = np.take(Prop,index_model)
        index_fault = np.where(np.array(Column_fault_name[index_model[0]:index_model[-1]+1]) == Name_of_fault)
        Indexfault_final = index_fault[0]
        shear_mod = Prop[Indexfault_final][0][11]
            
    else : #it's a geojson file
        with open(File_prop) as f:
            gj = geojson.load(f)
        faults = gj['features']
        
        for fi in range(len(faults)):
            if str(faults[fi]['properties']['si']) == Name_of_fault :
                if faults[fi]['properties']['model'] == Model :
                    shear_mod = faults[fi]['properties']['shear_modulus']
    return shear_mod
        
def sr_rate(Run_name,scenarios_names_list,mega_MFD,
Model_list,MFD_type_list,sub_area_file,File_geom,File_prop):

    '''####################################
    # use of the slip rate per fault
    #######################################'''  
    input_file_name = 'init'
    slip_rep_faults_all_data = open(Run_name + '/analysis/txt_files/slip_rep_on_faults_all_data.txt','w') #file with all the information about slip repartition
    
    time_tmp = time.time()
    
    for model in Model_list:
        for scenario_name in scenarios_names_list :
        
            print("t",round(time.time()-time_tmp))
            time_tmp = time.time()
            
            rup_file = (str(Run_name) + '/' + model + '/' + 'bg_' + str(mega_MFD[0][4])
            + '/' + str(mega_MFD[0][0]) + '_' + str(mega_MFD[0][1]) + '_' + str(mega_MFD[0][2])
            + '/sc_' + scenario_name  + '/bmin_' + str(mega_MFD[0][5]) + '_bmax_' + str(mega_MFD[0][6])
            + '/MFD_' + str(mega_MFD[0][7]) + '/Log/ruptures.txt')
            rup_id = []
            rup_length = []
            for line in open(rup_file):
                if not 'rup_id' in line :
                    rup_i = line.split('\t')[0]
                    rup_id.append(rup_i)
                    rup_length_i = line.split('\t')[1]
                    rup_length.append(len(rup_length_i.split(' '))-1)
                    
            for mega_mfd_i in mega_MFD :
                if mega_mfd_i[8] == scenario_name and mega_mfd_i[3] == model :
                    i_mfd =  0
                    input_file_name_i = (str(Run_name) + '/' + str(mega_mfd_i[3]) + '/' + 'bg_' + str(mega_mfd_i[4])
                                        + '/' + str(mega_mfd_i[0]) + '_' + str(mega_mfd_i[1]) + '_' + str(mega_mfd_i[2])
                                        + '/sc_' + str(mega_mfd_i[8])  + '/bmin_' + str(mega_mfd_i[5]) + '_bmax_' + str(mega_mfd_i[6])
                                        + '/MFD_' + str(mega_mfd_i[7]) + '/Log/sliprep_sample_' + str(mega_mfd_i[9]) + '.pkl')
                    if input_file_name_i != input_file_name :
                        print('reading pkl')
                        input_file_name = input_file_name_i
                        with open(input_file_name, 'rb') as f:
                            data_sr_rep = pickle.load(f)
                                            
                        print('reading pkl-done')
                        print()
                        for f_i in data_sr_rep :
                            fault_name = f_i
                            total_number = 0
                            n_fault_alone = 0.
                            n_FtF_2 = 0.
                            n_FtF_3 = 0.
                            n_FtF_4 = 0.
                            n_FtF_5 = 0.
                            n_FtF_6 = 0.
                            n_FtF_7 = 0.
                            n_FtF_8 = 0.
                            n_FtF_9 = 0.
                            n_FtF_10 = 0.
                            n_FtF_11 = 0.
                            n_FtF_12 = 0.
                            n_FtF_13 = 0.
                            for rup_i,rup_length_i in zip(rup_id,rup_length):
                                if str(rup_i) in data_sr_rep[f_i]:
                                    c = data_sr_rep[f_i][str(rup_i)]
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
                                    total_number += c
                            n_NMS = data_sr_rep[f_i]['NMS']
                            total_number += n_NMS
                            
                                
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
                            + str(p_FtF_11) + '\t'+ str(p_FtF_12) + '\t'+ str(p_FtF_13) + '\t'+ str(p_NMS)+'\n' )
                            slip_rep_faults_all_data.write(line)
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
                index_tmp = list(set(index_sc).intersection(index_m))
                index = list(set(index_tmp).intersection(index_mfd_type))
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
                
                #################
                # group up things for subareas
                # Calculate the NMS (in terms of moment rate)
                #################
                if os.path.exists(sub_area_file):
                    if not os.path.exists(str(Run_name) + '/analysis/txt_files/sub_area'):
                        os.makedirs(str(Run_name) + '/analysis/txt_files/sub_area')
                        
                    f_subarea_nms = open(Run_name + '/analysis/txt_files/sub_area/NMS_M0_'+model+'_' + MFD_type +'_' + scenario_name +'.txt','w')
                    f_subarea_nms.write('subarea_name\tmoment_rate_potential\tactual_moment_rate_faults\tNMS\n')
                    
                    # Get geometries of the faults and subarea and associate the two
                    sub_area_names,sub_area_lon,sub_area_lat = get_geom(sub_area_file,model)
                    Column_Fault_name,Longitudes,Latitudes = FaultGeometry(File_geom,model)
                    Lon, Lat = [] , []
#                    i_f = 0
#                    for fault_name in Column_Fault_name:
#                        if fault_name in list_fault_i:
#                            Lon.append(Longitudes[i_f])
#                            Lat.append(Latitudes[i_f])
#                        i_f +=1
#                    print(len(list_fault_i),len(Lon))
#                    print(Lon[0])
                    for fault_name in list_fault_i:
                        index_fault = np.where(np.array(Column_Fault_name)== fault_name)[0]
                        Lon.append(np.take(Longitudes,index_fault))
                        Lat.append(np.take(Latitudes,index_fault))

                    # Fault slip-rate
                    mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
                                  dtype = [('U100'),('U100'),('U1000'),('f8'),('f8')], delimiter = '\t')
                    model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
#                    set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
                    fault_name_mean_param =list( map(lambda i : mean_param[i][2], range(len(mean_param))))
                    sr_mean =list( map(lambda i : mean_param[i][3], range(len(mean_param))))
#                    Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
                    index_model = np.where(np.array(model_mean_param)==model)[0]
#                    set_mean_param_model = np.take(set_mean_param,index_model)
                    fault_name_mean_param_model = np.take(fault_name_mean_param,index_model)
                    sr_mean_model = np.take(sr_mean,index_model)
#                    Mmax_mean_model = np.take(Mmax_mean,index_model)
                    
                    # Faults areas
#                    log_Mmax_file = (str(Run_name) + '/' + str(mega_MFD[0][3]) + '/' + 'bg_' + str(mega_MFD[0][4])
#                    + '/' + str(mega_MFD[0][0]) + '_' + str(mega_MFD[0][1]) + '_' + str(mega_MFD[0][2])
#                    + '/sc_' + str(scenario_name)  + '/bmin_' + str(mega_MFD[0][5]) + '_bmax_' + str(mega_MFD[0][6])
#                    + '/MFD_' + str(MFD_type) + '/Log/Mmax_sample_1.txt')
                    log_Mmax_file = (str(Run_name)  + '/' + str(mega_MFD[0][3]) + '/Log/' + 'Mmax_sample_'+str(mega_MFD[0][0]) + '_' + str(mega_MFD[0][1]) + '_' + str(mega_MFD[0][2])+'_sc_' +  str(scenario_name)+ '_1.txt')
                    sources_names,_,_,sources_Areas = read_file_Mmax_log(log_Mmax_file)
                    
                    for i_sub_area in range(len(sub_area_names)):
                        sub_area_names_i = sub_area_names[i_sub_area]
                        sub_area_lon_i = sub_area_lon[i_sub_area]
                        sub_area_lat_i = sub_area_lat[i_sub_area]
                        faults_in_sub_area = find_faults_in_sub(list_fault_i,Lon,Lat,sub_area_lon_i,sub_area_lat_i)
                        
                        
                        sub_area_potential_M0 = 0.
                        sub_area_actual_M0 = 0.
                        for fault in list_fault_i :
                            if fault in faults_in_sub_area:
                                index = np.where(list_faults_slip_rep_i == fault)[0]
                                p_NMS_k = float(np.mean(np.take(p_NMS_i,index)))/100.
                                
                                index_fault = np.where(np.array(fault_name_mean_param_model)==fault)[0][0]
                                sr = sr_mean_model[index_fault]
                        
                                index_fault = np.where(np.array(sources_names)==fault)[0][0]
                                area = sources_Areas[index_fault]
                                
                                shear_mod = get_shearmod(File_prop,fault_name,model)
                                shear_mod = float(shear_mod)*10**9
                                
                                # Moment rate in m.yr
                                # converting the area in m2 and the sliprate in m/yr
                                M0 = float(area * 1000000. * shear_mod * sr / 1000.)
                                
                                # adding the sum for the subarea
                                sub_area_potential_M0 += M0
                                sub_area_actual_M0 += M0 * (1.-p_NMS_k)
                                
                        # NMS for the subarea
                        NMS = 1. - sub_area_actual_M0/sub_area_potential_M0
                        
                        print(sub_area_names_i,round(sub_area_potential_M0),round(sub_area_actual_M0),round(NMS,3))
                        f_subarea_nms.write(sub_area_names_i+'\t'+
                        str(round(sub_area_potential_M0))+'\t'+
                        str(round(sub_area_actual_M0))+'\t'+str(round(NMS,3)))
                    f_subarea_nms.close()
