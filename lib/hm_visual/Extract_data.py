#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import os, sys
import numpy as np
import xml.etree.ElementTree as ET
import Read_file
import pandas as pd


def extract(Run_name):
    if not os.path.exists(str(Run_name) + '/analysis'):
        os.makedirs(str(Run_name) + '/analysis')
    if not os.path.exists(str(Run_name) + '/analysis/txt_files'):
        os.makedirs(str(Run_name) + '/analysis/txt_files')
    
    ## defining the binning in magnitude
    mega_bining_in_mag = np.linspace(4,9.9,num =60)
    mega_MFD = []
    
    f_sr = []
    m_Mmax = [] #model Mmax
    a_s_model = []
    faults_names = []
    scenarios_names = []
    
    #list of the general parameters of each model
    b_sample = []
    Mt_sample = []
    sm_sample = []
    
    faults_name_list = []
    faults_id_list = []
    fault_id = 0
    scenarios_names_list = []
    sample_list = []
    BG_hyp_list = []
    ScL_complet_list = [] #ScL + dimension used + all data
    ScL_list = []
    dimension_used_list = []
    Model_list = []
    b_value_list = []
    MFD_type_list = []
    
    boolean_mfd = True
    
    LT_file = str(Run_name)+'/Sources_Logic_tree.xml'

    tree = ET.parse(LT_file)        
    nrml = tree.getroot()            
    Branch_names = []
    branch_path = []
    general_weight = []                   
    for logicTree in nrml:
        for logicTreeBranchLevel in logicTree:
            for logicTreeBranchSet in logicTreeBranchLevel:
                for logicTreeBranch in logicTreeBranchSet:
                    Branch_names.append(logicTreeBranch.attrib['branchID'])
                    path_i = logicTreeBranch[0].text[:-4]
                    branch_path.append(path_i.split('/'))
                    general_weight.append(logicTreeBranch[1].text)
    
    df_mega_MFD = pd.DataFrame(columns=['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                        'b_min','b_max','MFD_type','scenario_set','sample','source',
                                        '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                        '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                        '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                        '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                        '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9',
                                        '9.0','9.1','9.2','9.3','9.4','9.5','9.6','9.7','9.8','9.9'], index = range(len(Branch_names)*10000)) 

    index_df = 0
               
    slip_rate_sampling = open(Run_name + '/analysis/txt_files/slip_rate_sampling.txt','w') 
    mean_parameters_faults = open(Run_name + '/analysis/txt_files/mean_parameters_faults.txt','w') 
    
    for branch in Branch_names:
        branch = branch.split('-')


        Model = branch[0]
        BG_hyp = branch[1][3::]
        selected_ScL = branch[2]
        dim_used = branch[3][0]
        str_all_data = branch[4]
        scenario_set= branch[5][3::]
        b_value = branch[6]
        b_min = float(b_value.split('_')[1])
        b_max = float(b_value.split('_')[3])
        MFD_type =branch[7][4:]
        sample = (branch[8].split('_')[1])  
        
              
        
        # extract the slip-rates and the faults names
        log_sr_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/' 
                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(MFD_type)
                            + '/Log/slip_rate_sample_' + str(sample) + '.txt') 
        
           
        faults_names_i, faults_slip_rates__i = Read_file.read_file_sr_log(log_sr_file) #read the log of slip rateassociated with the model
                    
        f_sr.append(faults_slip_rates__i)
        faults_names.append(faults_names_i)
        for fault_name,sr in zip(faults_names_i,faults_slip_rates__i):
            slip_rate_sampling.write(str(Run_name)  +'\t'+ str(Model) +'\t'+ str(BG_hyp) +'\t'+ str(selected_ScL) + '_' 
                           + str(dim_used) + '_' + str_all_data +'\t'+ str(scenario_set) +'\t'
                        + 'bmin_' + str(b_min) + '_bmax_' + str(b_max)+'\t'+str(MFD_type)+'\t'+str(sample)+'\t'+str(fault_name)+'\t'+str(sr)+'\n')
        
        
        #extract the final aseismic slip in the model 
        log_as_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/' 
                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(MFD_type)
                            + '/Log/calculation_sample_' + str(sample) + '.txt')    
        a_s_i = Read_file.read_file_as_log(log_as_file) #read the log of aseismic slip associated with the model
                    
        a_s_model.append(a_s_i)  
        
        #extract the general parameters of the model
        log_general_param_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/' 
                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(MFD_type)
                            + '/Log/general_parameters_sample_' + str(sample) + '.txt')                
        M_trunc_model,b_value_model = Read_file.read_file_general_param_log(log_general_param_file) #read the log of the general parameters associated with the model
        
        b_sample.append(b_value_model)
        Mt_sample.append(M_trunc_model)
        
        
        # extract the Mmax of the faults and the scenarios
        log_Mmax_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/' 
                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(MFD_type)
                            + '/Log/Mmax_sample_' + str(sample) + '.txt')                
        sources_names,sources_Mmax,sources_Lengths,sources_Areas = Read_file.read_file_Mmax_log(log_Mmax_file) #read the log of Mmax associated with the model
        
        m_Mmax.append(max(sources_Mmax))
        
                
       
        #build the array containing the name of each fault
        
        for name in faults_names :
            if not name in faults_name_list :
                faults_name_list.append(name)
                faults_id_list.append(fault_id)
                fault_id += 1
                
        #build the array containing the list of the samples
        if not sample in sample_list :
                sample_list.append(sample)
                
        #build the array containing the list of the b value
        if not b_value in b_value_list :
                b_value_list.append(b_value)
                
        #build the array containing the list of the b value
        if not MFD_type in MFD_type_list :
                MFD_type_list.append(MFD_type)
                
        #build the array containing the list of the dimention used
        if not dim_used in dimension_used_list :
                dimension_used_list.append(dim_used)
                
        scenarios_names.append(scenario_set)
        
        if not scenario_set in scenarios_names_list :
            scenarios_names_list.append(scenario_set)
            
                
        if not BG_hyp in BG_hyp_list :
            BG_hyp_list.append(BG_hyp)
            
            
        if not (str(selected_ScL) + '_' + str(dim_used) + '_' + str(str_all_data)) in ScL_complet_list :
            ScL_complet_list.append((str(selected_ScL) + '_' + str(dim_used) + '_' + str(str_all_data)))
            
        if not selected_ScL in ScL_list :
            ScL_list.append(selected_ScL)
            
        if not Model in Model_list :
            Model_list.append(Model)
            
            
        #gettin gthe mean value of slip_rate and Mmax(sample 1)
        for fault_name,sr in zip(faults_names_i,faults_slip_rates__i):
            if sample == '1' and BG_hyp==BG_hyp_list[0] and (str(selected_ScL) + '_' + str(dim_used) + '_' + str(str_all_data)) == ScL_complet_list[0]:
                if b_value == b_value_list[0] and MFD_type==MFD_type_list[0]:
                    Mmax_fault =0.
                    for source,Mmax_i in zip(sources_names,sources_Mmax):
                        if fault_name in source:
                            if float(Mmax_i) > Mmax_fault :
                                Mmax_fault = Mmax_i
                    mean_parameters_faults.write(str(Model) +'\t'+  str(scenario_set) +'\t'+str(fault_name)+'\t'+str(sr)+'\t'+str(Mmax_fault)+'\n')
         
                
                
        if boolean_mfd == True : #there is probably a smateer way to buit the data frame so it's lighter...
            log_mfd_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/' 
                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(MFD_type)
                            + '/Log/mdf_sample_' + str(sample) + '.txt')
            sources_names,Mmin,MFD = Read_file.read_file_mfd_log(log_mfd_file) #read the log of mfd associated with the model
                
            ### storing in the MFD super matrix
            index_Mmin = np.where(np.array(mega_bining_in_mag) == min(Mmin))[0]
            
            index_source = 0
            for source in sources_names :
                '''
                #the source is replaced by a list of id of fault as a tring to gain memory place
                list_id  =[]
                for i in source.split("=['"): 
                    try :
                        iii = i.index("'")
                        fault_id = np.where(np.array(faults_name_list[0])==i[:iii])[0][0]
                        list_id.append(fault_id)
                    except:
                        i=i
                if len(source.split("=['")) == 1:
                    if source == 'Background':
                        list_id = '[Background]'
                    else:
                        #print(source[len(Model)+1:])
                        fault_id = np.where(np.array(faults_name_list[0])==source[len(Model)+1:])[0][0]
                        list_id.append(fault_id)'''
                
                mega_mfd_i = []
#                print(source,sys.getsizeof(source),str(list_id)[1:-1],sys.getsizeof(str(list_id)[1:-1]))
#                print()
                mega_mfd_i.append(selected_ScL)
                mega_mfd_i.append(dim_used)
                mega_mfd_i.append(str_all_data)
                mega_mfd_i.append(Model)
                mega_mfd_i.append(BG_hyp)
                mega_mfd_i.append(str(b_min))
                mega_mfd_i.append(str(b_max))
                mega_mfd_i.append(MFD_type)
                mega_mfd_i.append(scenario_set)
                mega_mfd_i.append(sample)
#                try: #shorten the string from useless stuff
#                    for i_str in range(10):
#                        source= source.replace("_f_"+str(i_str)+"=['","")
#                        source= source.replace("scenario","")
#                        source= source[len(Model)+1:]
#                except:
#                    source = source
                mega_mfd_i.append(source)
                #mega_mfd_i.append(str(list_id)[1:-1])
                
                for index_mag in range(len(mega_bining_in_mag)):
                    if (index_mag < (len(MFD[index_source]) + index_Mmin)) and (index_mag >= index_Mmin):
                        try :
                            mega_mfd_i.append(float(MFD[index_source][int(index_mag-index_Mmin)]))
                        except TypeError:
                            print('!!!!!!!!!!!!!!!!\n\n\n'+
                            'There is a problem with a log file'+
                            'Delete the .xml corresponding to this file'+
                            'Then re-run SHERIFS with overwrite_files = False'+
                            'file with a problem : ')
                            print(log_mfd_file)
                    else :
                        mega_mfd_i.append(0.)
                mega_MFD.append(np.array(mega_mfd_i))
                dict_df_mfd_i = {'selected_ScL':selected_ScL,
                                 'dim_used':dim_used,
                                 'str_all_data':str_all_data,
                                 'Model':Model,
                                 'BG_hyp':BG_hyp,
                                 'b_min':str(b_min),
                                 'b_max':str(b_max),
                                 'MFD_type':MFD_type,
                                 'scenario_set':scenario_set,
                                 'sample':sample,
                                 'source':source,
                                 '4.0':np.sum(mega_mfd_i[11:]),
                                 '4.1':np.sum(mega_mfd_i[12:]),
                                 '4.2':np.sum(mega_mfd_i[13:]),
                                 '4.3':np.sum(mega_mfd_i[14:]),
                                 '4.4':np.sum(mega_mfd_i[15:]),
                                 '4.5':np.sum(mega_mfd_i[16:]),
                                 '4.6':np.sum(mega_mfd_i[17:]),
                                 '4.7':np.sum(mega_mfd_i[18:]),
                                 '4.8':np.sum(mega_mfd_i[19:]),
                                 '4.9':np.sum(mega_mfd_i[20:]),
                                 '5.0':np.sum(mega_mfd_i[21:]),
                                 '5.1':np.sum(mega_mfd_i[22:]),
                                 '5.2':np.sum(mega_mfd_i[23:]),
                                 '5.3':np.sum(mega_mfd_i[24:]),
                                 '5.4':np.sum(mega_mfd_i[25:]),
                                 '5.5':np.sum(mega_mfd_i[26:]),
                                 '5.6':np.sum(mega_mfd_i[27:]),
                                 '5.7':np.sum(mega_mfd_i[28:]),
                                 '5.8':np.sum(mega_mfd_i[29:]),
                                 '5.9':np.sum(mega_mfd_i[30:]),
                                 '6.0':np.sum(mega_mfd_i[31:]),
                                 '6.1':np.sum(mega_mfd_i[32:]),
                                 '6.2':np.sum(mega_mfd_i[33:]),
                                 '6.3':np.sum(mega_mfd_i[34:]),
                                 '6.4':np.sum(mega_mfd_i[35:]),
                                 '6.5':np.sum(mega_mfd_i[36:]),
                                 '6.6':np.sum(mega_mfd_i[37:]),
                                 '6.7':np.sum(mega_mfd_i[38:]),
                                 '6.8':np.sum(mega_mfd_i[39:]),
                                 '6.9':np.sum(mega_mfd_i[40:]),
                                 '7.0':np.sum(mega_mfd_i[41:]),
                                 '7.1':np.sum(mega_mfd_i[42:]),
                                 '7.2':np.sum(mega_mfd_i[43:]),
                                 '7.3':np.sum(mega_mfd_i[44:]),
                                 '7.4':np.sum(mega_mfd_i[45:]),
                                 '7.5':np.sum(mega_mfd_i[46:]),
                                 '7.6':np.sum(mega_mfd_i[47:]),
                                 '7.7':np.sum(mega_mfd_i[48:]),
                                 '7.8':np.sum(mega_mfd_i[49:]),
                                 '7.9':np.sum(mega_mfd_i[50:]),
                                 '8.0':np.sum(mega_mfd_i[51:]),
                                 '8.1':np.sum(mega_mfd_i[52:]),
                                 '8.2':np.sum(mega_mfd_i[53:]),
                                 '8.3':np.sum(mega_mfd_i[54:]),
                                 '8.4':np.sum(mega_mfd_i[55:]),
                                 '8.5':np.sum(mega_mfd_i[56:]),
                                 '8.6':np.sum(mega_mfd_i[57:]),
                                 '8.7':np.sum(mega_mfd_i[58:]),
                                 '8.8':np.sum(mega_mfd_i[59:]),
                                 '8.9':np.sum(mega_mfd_i[60:]),
                                 '9.0':np.sum(mega_mfd_i[61:]),
                                 '9.1':np.sum(mega_mfd_i[62:]),
                                 '9.2':np.sum(mega_mfd_i[63:]),
                                 '9.3':np.sum(mega_mfd_i[64:]),
                                 '9.4':np.sum(mega_mfd_i[65:]),
                                 '9.5':np.sum(mega_mfd_i[66:]),
                                 '9.6':np.sum(mega_mfd_i[67:]),
                                 '9.7':np.sum(mega_mfd_i[68:]),
                                 '9.8':np.sum(mega_mfd_i[69:]),
                                 '9.9':np.sum(mega_mfd_i[70:])}
                df_mega_MFD.loc[index_df] = pd.Series(dict_df_mfd_i)
                
                '''df_mega_MFD.loc[index_df] = mega_mfd_i
                ,columns=['selected_ScL','dim_used','str_all_data','Model','BG_hyp',
                                        'b_min','b_max','MFD_type','scenario_set','sample','source',
                                        '4.0','4.1','4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9',
                                        '5.0','5.1','5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9',
                                        '6.0','6.1','6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9',
                                        '7.0','7.1','7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9',
                                        '8.0','8.1','8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9',
                                        '9.0','9.1','9.2','9.3','9.4','9.5','9.6','9.7','9.8','9.9'])'''
                index_df +=1
                index_source += 1
        
            
            
    boolean_mfd = True        
    
    file_MFD_name = str(Run_name) + '/analysis/txt_files/faults_MFD.txt'
    file_MFD = open(file_MFD_name,'w')
    file_MFD.write(str(mega_MFD))
    file_MFD.close()
    slip_rate_sampling.close()
    mean_parameters_faults.close()
    
    
    df_mega_MFD = df_mega_MFD.dropna(how = 'all')
#    df_mega_MFD.info(memory_usage = 'deep')
    
    return (mega_MFD,df_mega_MFD, scenarios_names_list, ScL_complet_list, ScL_list, Model_list,BG_hyp_list,
            dimension_used_list,faults_name_list,sample_list,b_value_list,MFD_type_list,m_Mmax,
            mega_bining_in_mag,a_s_model,b_sample,sm_sample,Mt_sample,sources_Lengths,sources_Areas)
