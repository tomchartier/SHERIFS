# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0

This code creates a simple earthquake source model and
exports it as  OQ file.

@author: Thomas Chartier
"""
import os
import numpy as np
import sys
path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/..'
sys.path.append(path_lib)
path_f = path_lib + '/core'
sys.path.append(path_f)
path_f = path_lib + '/utils'
sys.path.append(path_f)
path_f = path_lib + '/sources'
sys.path.append(path_f)
import EQ_on_faults
import Geometry_scenario
import select_sr
from background import bg
from geometry_tools import *
import host_model
import math
import fault_source
import pickle
import write_section_file as wsf
import write_multifault_source_file as wmfs
import xml.etree.ElementTree as ET
import matplotlib.path as mplPath
from seismic_moment import mag_to_M0

import matplotlib.pyplot as plt

class Source_Model_Creator:
    def __init__(self,path,pathlog,param,Model_name,rupture_set,sample,
                 Domain_in_model,selected_ScL,dimention_used,use_all_ScL_data,
                 mfd_param,mfd_hyp,bg_ratio,calculation_log_file,faults_names,
                 scenarios_names,faults_data,faults_lon,faults_lat,
                 list_fbg,fbgpath,branch):

        self.Run_Name = param["Run_Name"]
        self.param = param
        self.Domain_in_the_model = Domain_in_model

        self.rupture_set = rupture_set

        #a envoyer dans logic tree
        self.Model_name = Model_name
        # self.b_value = float(b_value)
        self.mfd_param = mfd_param

        self.mfd_hyp = mfd_hyp

        self.sample = sample

        self.path = path
        self.pathlog = pathlog

        np.random.seed = param["main"]["parameters"]["Random_seed"]
        self.Mmin = param["main"]["parameters"]["Mmin"]

        self.selected_ScL = selected_ScL
        self.dimention_used = dimention_used
        self.use_all_ScL_data = use_all_ScL_data
        self.bg_ratio = bg_ratio

        self.sr_correl = param["main"]["parameters"]["SR_correl"]
        self.size_of_increment = param["main"]["parameters"]["dsr"]
        self.fit_quality = param["main"]["parameters"]["fit_quality"] / 100.
        self.Mmax_range = param["main"]["parameters"]["Mmax_range"]

        self.calculation_log_file = calculation_log_file

        # background option and parameters
        if not param["main"]["background"]["option_bg"] in ["None","none"]:
            self.File_bg = param["main"]["background"]["File_bg"]

        if param["main"]["background"]["option_bg"] == "zone":
            if param["main"]["background"]["use_host_model"] in ["true","True"]:
                self.use_host_model = True
                self.host_model_file = param["main"]["background"]["host_model_file"]
                self.file_prop_bg =  param["main"]["background"]["file_prop_bg"]

        self.faults_names = faults_names
        self.scenarios_names = scenarios_names
        self.faults_data=faults_data
        self.faults_lon,self.faults_lat= faults_lon,faults_lat

        self.list_fbg = list_fbg
        self.fbgpath = fbgpath

        self.branch = branch

        self.initialize()

    def initialize(self):

        explo_time = self.param["main"]["parameters"]["explo_time"]

        # for the use of multifault source
        if self.param["main"]["parameters"]["use_multiF"] in ["True",'true']:
            use_multiF = True
        else :
            use_multiF == False

        #maxi number of sources per input files
        # cutting the files in smaller chunks allow OQ to read the sources in
        # parallel
        max_rup_per_file = 2000

        # init
        faults_names = self.faults_names
        scenarios_names = self.scenarios_names
        faults_data = self.faults_data
        faults_lon,faults_lat= self.faults_lon,self.faults_lat

        # list files for branch
        list_src_files = []

        # create log files and repo
        log_sr_file=open(self.pathlog+'/slip_rate_sample_' + str(self.sample) + '.txt','w')
        log_mdf_file=open(self.pathlog+'/mdf_sample_' + str(self.sample) + '.txt','w')

        # test if the model is large enough to deserve to be cut

        cut_sm_file = False
        if len(faults_names) + len(scenarios_names) > max_rup_per_file :
            cut_sm_file = True

        if cut_sm_file == False :
            XMLfile=open(self.path +'/Source_model_' + str(self.sample) + '.xml','w')
        else :
            # find number of subparts for simple faults sources
            n_cut_sf = 1
            while len(faults_names)/n_cut_sf > max_rup_per_file:
                n_cut_sf +=1
            # find number of subparts for multi faults sources
            n_cut_mf = 1
            while len(scenarios_names)/n_cut_mf > max_rup_per_file:
                n_cut_mf +=1

            if use_multiF == False  :
                " create the file list"
                sf_files, sf_counter = [], []

                for i in range(n_cut_sf):
                    f = self.path +'/sm_' + str(self.sample) + '_sf_'+str(i+1)+'.xml'
                    sf_files.append(open(f,'w'))
                    sf_counter.append(0)
                    list_src_files.append(f)
                mf_files, mf_counter = [], []
                for i in range(n_cut_mf):
                    f = self.path +'/sm_' + str(self.sample) + '_mf_'+str(i+1)+'.xml'
                    mf_files.append(open(f,'w'))
                    mf_counter.append(0)
                    list_src_files.append(f)

            if use_multiF == True  :
                n_cut_f = n_cut_sf + n_cut_mf
                s_files = []
                for i in range(n_cut_mf):
                    f = self.path +'/sm_'+str(i)+'.xml'
                    s_files.append(f)
                    list_src_files.append(f)

                dict_txt, dict_n_rup = {}, {}
                for f in s_files  :
                    dict_txt.update({f:''})
                    dict_n_rup.update({f:0})




        if use_multiF == True  :
            '''
            #######################################
            #######################################
            USING MULTIFAULT TYPOLOGY IN OQ
            #######################################
            #######################################
            '''
            # write the section files
            sections_xml = self.Run_Name+'/ssm/'+self.Model_name+'_sections.xml'
            if not os.path.exists(sections_xml):
                txt = wsf.start(self.Model_name)
                for section_id in range(len(faults_names)):
                    txt = wsf.wrt_section(txt,section_id,faults_names,
                    faults_data,"kite")
                txt = wsf.end(txt)
                wsf.build(sections_xml,txt)


            #initialisation of the general parameters (M_min, shear modulus and b value)
            log_general_parameters_file = open(self.pathlog +'/general_parameters_sample_' + str(self.sample) + '.txt','w')
            M_min = float(self.Mmin)
            log_general_parameters_file.write('M_tronc\t'+str(M_min)+'\n')
            # if self.sample == 1 :
            #     b_value = (self.b_min + self.b_max)/2.
            # else :
            #     b_value = np.random.triangular(self.b_min,(self.b_min + self.b_max)/2.,self.b_max)
            # b_value = self.b_value
            # mfd_param = {}
            # mfd_param.update({'b_value' : b_value})
            mfd_param = self.mfd_param
            log_general_parameters_file.write('b_value\t'+str(mfd_param["b_value"])+'\n')
            log_general_parameters_file.close()

            #initialization of array for faults information
            faults_area = []
            faults_length = []
            faults_width = []
            faults_slip_rates = []
            faults_mecanism = []
            faults_shear_mod = []

            '''########################################################
            # Random sampling of the fault slip-rate
            ########################################################'''
            # find which faults interact a lot with each other
            print("Picking slip-rates...")
            if self.sample == 1 :
                self.sr_correl = False

            if self.sr_correl ==True :
                M_faults_correl = []
                for Fault_name in faults_names:
                    M_faults_correl_i =  []
                    sc_with_fault = []
                    for index_scenario in range(len(scenarios_names)) :
                        scenario_i = self.rupture_set[index_scenario]
                        if Fault_name in scenario_i :
                            sc_with_fault += scenario_i
                    for Fault_name_i in faults_names:
                        nb_join_ruptures = 0
                        if Fault_name == Fault_name_i :
                            nb_join_ruptures = 0
                        elif Fault_name_i in sc_with_fault:
                            nb_join_ruptures = sc_with_fault.count(Fault_name_i)
                        else :
                            nb_join_ruptures = 0
                        M_faults_correl_i.append(nb_join_ruptures)
                    M_faults_correl.append(M_faults_correl_i)

                #calculating the linked level (up to five)
                M_linked_lvl = np.zeros_like(M_faults_correl)
                for lvl in [1,2,3,4]:
                    for i in range(len(faults_names)):
                        for j in range(len(faults_names)):
                            if M_linked_lvl[i][j] == 0:
                                if M_faults_correl[i][j] != 0:
                                    M_linked_lvl[i][j]=lvl
                                else :
                                    for jj in  range(len(M_faults_correl[j])): #for all faults connecting to a fault connected to i
                                        if M_faults_correl[j][jj]!=0:
                                            if M_linked_lvl[i][j] == 0:
                                                M_linked_lvl[i][j]=lvl+1
                #print(M_linked_lvl)

                max_correl = np.max(M_faults_correl)

                list_quater_picked = list(np.zeros(len(faults_names)))
                # pick the faults slip rates
                index_fault = 0
                for Fault_name in faults_names:
                    faults_length.append(faults_data[index_fault]['length'])
                    faults_area.append(faults_data[index_fault]['area'])
                    faults_width.append(faults_data[index_fault]['width'])
                    faults_mecanism.append(faults_data[index_fault]['mecanism'])
                    sr_values = [faults_data[index_fault]['slip_rate_min'],
                    faults_data[index_fault]['slip_rate_moy'],
                    faults_data[index_fault]['slip_rate_max']]

                    #selects randomly the slip rate
                    slip_rate = select_sr.select(sr_values,
                                                self.sample,
                                                index_fault,
                                                M_linked_lvl[index_fault],
                                                list_quater_picked)

                    log_line = str(Fault_name) + '\t' + str(slip_rate) + '\n' #writting in the log file
                    log_sr_file.write(log_line)

                    faults_slip_rates.append(slip_rate * 0.001) # transform from mm to m

    #                faults_shear_mod.append(float(self.shear_mod)*10**9 )
                    faults_shear_mod.append(faults_data[index_fault]['shear_mod'])
                    index_fault += 1

            else: # the option of correlation bewteen corroption fault not used
                index_fault = 0
                for Fault_name in faults_names:
                    faults_length.append(faults_data[index_fault]['length'])
                    faults_area.append(faults_data[index_fault]['area'])
                    faults_width.append(faults_data[index_fault]['width'])
                    faults_mecanism.append(faults_data[index_fault]['mecanism'])
                    #ramdom sampling the slip rate  . uniform sampling
                    slip_rate_min = faults_data[index_fault]['slip_rate_min']
                    slip_rate_moy = faults_data[index_fault]['slip_rate_moy']
                    slip_rate_max = faults_data[index_fault]['slip_rate_max']

                    if self.sample == 1 :
                        slip_rate = slip_rate_moy
                    else :
                        slip_rate_inf = np.random.uniform(slip_rate_min,slip_rate_moy)
                        slip_rate_sup = np.random.uniform(slip_rate_moy,slip_rate_max)
                        slip_rate = np.random.choice([slip_rate_inf,slip_rate_sup])

                    log_line = str(Fault_name) + '\t' + str(slip_rate) + '\n' #writting in the log file
                    log_sr_file.write(log_line)
                    faults_slip_rates.append(slip_rate * 0.001) # transform from mm to m
                    faults_shear_mod.append(faults_data[index_fault]['shear_mod'])

                    index_fault += 1
            print("\t\tslip-rates picked.")

            ratio_test = 0.5
            count_reruns = 1 #used to divide the sr increment if the fit is not good
            count_mfd90 = 1

            f_pkl_mdf = self.pathlog +'/mfd_' + str(self.sample) + '.pkl'
            re_use_mfd_pkl= True
            if not os.path.isfile(f_pkl_mdf):
                re_use_mfd_pkl = False
            if re_use_mfd_pkl == False:
                while abs(ratio_test-1) >self.fit_quality or math.isnan(ratio_test) == True:
                    MFDs = EQ_on_faults.EQ_on_faults_from_sr(self.Run_Name,
                                                            M_min,
                                                            mfd_param,
                                                            faults_names,
                                                            faults_area,
                                                            faults_length,
                                                            faults_width,
                                                            faults_slip_rates,
                                                            scenarios_names,
                                                            faults_shear_mod,
                                                            self.path,
                                                            self.pathlog,
                                                            self.sample,
                                                            self.selected_ScL,
                                                            self.dimention_used,
                                                            self.use_all_ScL_data,
                                                            faults_mecanism,
                                                            self.bg_ratio,
                                                            self.size_of_increment,
                                                            self.mfd_hyp,
                                                            count_reruns,
                                                            faults_lon,
                                                            faults_lat,
                                                            self.Mmax_range,
                                                            self.calculation_log_file,
                                                            self.branch,
                                                            self.param)

                    ratio_test = MFDs.ratio_test
                    if abs(ratio_test-1) > self.fit_quality:
                        print('bad sampling => re-run')
                        count_reruns += 1
                        if MFDs.ratio_NMS> 90.:
                            count_reruns -= 1
                            count_mfd90 +=1
                        else :
                            count_mfd90 =1
                    else:
                        if MFDs.ratio_NMS> 90.:
                            ratio_test = 0.5
                            print('bad sampling => re-run')
                            count_mfd90 +=1
                        else :
                            count_mfd90 =1

                    if math.isnan(ratio_test) == True:
                        print('bad sampling => re-run')
                        count_reruns = 1
                    if count_reruns > 3 or count_mfd90 > 3:
                        print('\n\n\n!!!!!! maybe there is a problem!!!')
                        ratio_test = 1.

                with open(f_pkl_mdf, 'wb') as f:
                    MFDs_to_pkl = [
                                    MFDs.EQ_rate_BG,
                                    MFDs.faults_names,
                                    MFDs.OQ_entry_faults,
                                    MFDs.scenarios_names,
                                    MFDs.OQ_entry_scenarios,
                                    MFDs.index_faults_in_scenario
                                    ]
                    pickle.dump(MFDs_to_pkl, f)

            else :
                print('Reloading MFDs from previous run')
                with open(f_pkl_mdf, 'rb') as f:
                    MFDs_to_pkl = pickle.load(f)

            EQ_rate_BG = MFDs_to_pkl[0]
            faults_names = MFDs_to_pkl[1]
            OQ_entry_faults = MFDs_to_pkl[2]
            scenarios_names = MFDs_to_pkl[3]
            OQ_entry_scenarios = MFDs_to_pkl[4]
            index_faults_in_scenario = MFDs_to_pkl[5]

            # sclaling law as called by openquake
            if self.selected_ScL == 'Le2010' :
                ScL_oq = 'Leonard2014_SCR'
            if self.selected_ScL == 'WC1994' :
                ScL_oq = 'WC1994'
            else : #default Scaling relationship for opanquake
                ScL_oq = 'WC1994'

            # Extract the background geometry
            Lon_bg, Lat_bg  = bg.geom(self.Model_name,self.File_bg )

            #Create the BG polygon
            Poly_bg = []
            for x1,y1 in zip(Lon_bg,Lat_bg):
                Poly_bg.append((x1,y1))
            bbPath_BG = mplPath.Path(Poly_bg)

            # find the faults in the bg
            f_in_bg = []
            for index_fault,fault_name in zip(range(len(faults_names)),faults_names):
                for lon,lat in zip(faults_data[index_fault]['lon'],faults_data[index_fault]['lat']):
                    if bbPath_BG.contains_point((lon,lat)):
                        f_in_bg.append(index_fault)


            source_id = 0

            # write source file
            if cut_sm_file == False :
                source_xml = self.path +'/Source_model_' + str(self.sample) + '.xml'
                trt = faults_data[0]['domain']
                txt = wmfs.start(explo_time,trt)
                name = "multifaultsource"
                txt = wmfs.start_multifault_source(txt,name,trt,
                sections_xml,source_id)
            else :
                for f in s_files  :
                    trt = faults_data[0]['domain']
                    dict_txt[f] = wmfs.start(explo_time,trt)
                    name = "multifaultsource_"+str(source_id)
                    dict_txt[f] = wmfs.start_multifault_source(dict_txt[f],
                    name,trt,sections_xml,source_id)
                    source_id+=1
                # working file
                f = list(dict_txt.keys())[0]

            #loop on simple faults
            for index_fault,fault_name in zip(range(len(faults_names)),faults_names):
                if index_fault in f_in_bg:
                    # write the single fault ruptures in the file
                    i_MFD = np.where(np.array(faults_names) == fault_name)[0][0]
                    MFD = OQ_entry_faults[i_MFD]

                    if sum(MFD)!=0:
                        source_id = source_id + 1
                        Fault_Name = self.Model_name + '_' + str(fault_name)
                        fault_trt = faults_data[index_fault]['domain']
                        if not fault_trt in str(self.Domain_in_the_model):
                            self.Domain_in_the_model.append(fault_trt)

                        trt = str(fault_trt)
                        rake = faults_data[index_fault]['rake']

                        name = "single_fault_"+fault_name
                        if cut_sm_file == False :
                            txt = wmfs.wrt_multifault_source(txt,MFD,
                                self.Mmin,explo_time,rake,[index_fault])
                        else :
                            if not dict_n_rup[f] < max_rup_per_file :
                                i = 0
                                while dict_n_rup[f] > max_rup_per_file:
                                    i+=1
                                    f = list(dict_txt.keys())[i]
                            txt = dict_txt[f]
                            dict_txt[f] = wmfs.wrt_multifault_source(txt,MFD,
                                self.Mmin,explo_time,rake,[index_fault])
                            dict_n_rup[f] += 1



            # loop on the multifault_rup
            if len(self.rupture_set) != 0 :
                for scenario in enumerate(scenarios_names):
                    # check if the scenario has at least one fault in the bg
                    sc_in = False
                    index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
                    index_faults_in_sc =  index_faults_in_scenario[index_scenario][0]

                    for i in index_faults_in_sc:
                        if i in f_in_bg:
                            sc_in = True

                    if sc_in == True:

                        index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
                        MFD = OQ_entry_scenarios[index_scenario]

                        if sum(MFD)!=0:
                            faults_in_scenario = np.take(faults_names,index_faults_in_sc)

                            source_id = source_id + 1

                            scenar_name = '_'.join("{!s}={!r}".format(key,val)
                             for (key,val) in scenario[1].items())
                            source_name = self.Model_name + '_scenario_'
                            source_name += str(scenar_name)

                            list_trt = []
                            scenario_mechanism =[]
                            for fname in faults_in_scenario :
                                i = faults_names.index(fname)
                                scenario_mechanism.append(faults_data[i]['rake'])
                                fault_trt = faults_data[i]['domain']
                                list_trt.append(fault_trt)
                                if not fault_trt in str(self.Domain_in_the_model):
                                    self.Domain_in_the_model.append(fault_trt)
                            fault_trt = max(list_trt, key = list_trt.count)
                            trt = str(fault_trt)


                            rake= np.mean(scenario_mechanism)

                            if cut_sm_file == False :
                                txt = wmfs.wrt_multifault_source(txt,MFD,
                                    self.Mmin,explo_time,rake,index_faults_in_sc)
                            else :
                                if not dict_n_rup[f] < max_rup_per_file :
                                    i = 0
                                    while dict_n_rup[f] > max_rup_per_file:
                                        i+=1
                                        f = list(dict_txt.keys())[i]
                                txt = dict_txt[f]
                                dict_txt[f] = wmfs.wrt_multifault_source(txt,MFD,
                                    self.Mmin,explo_time,rake,index_faults_in_sc)
                                dict_n_rup[f] += 1

            # completing the file
            if cut_sm_file == False :
                txt = wmfs.end_multifault_source(txt)
                txt = wmfs.end(txt)
                wmfs.build(source_xml,txt)
            else :
                for f in s_files  :
                    if dict_n_rup[f] > 0 :
                        dict_txt[f] = wmfs.end_multifault_source(dict_txt[f])
                        txt = wmfs.end(dict_txt[f])
                        wmfs.build(f,txt)


        else :
            '''
            ###########
            ############
            NOT USING MULTIFAULTSOURCE TYPOLOGY in OQ

            '''
            # Initiate the xml file
            line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
            line+='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
            line+='\txmlns="http://openquake.org/xmlns/nrml/0.5">\n'
            line+='\t<sourceModel name="Hazard Model"'

            line+=' investigation_time="'+str(round(explo_time,1))+'">\n'
            if cut_sm_file == False :
                XMLfile.write(line)
            else :
                for f in sf_files + mf_files  :
                    f.write(line)

            line='\t\t<sourceGroup\nname="group 1"\nrup_interdep="indep"\nsrc_interdep="indep"\ntectonicRegion="Active Shallow Crust"\n>\n'
            if cut_sm_file == False :
                XMLfile.write(line) # TODO change the tectonic Region
            else :
                for f in sf_files + mf_files  :
                    f.write(line)

            #initialisation of the general parameters (M_min, shear modulus and b value)
            log_general_parameters_file = open(self.pathlog +'/Log/general_parameters_sample_' + str(self.sample) + '.txt','w')

            M_min = float(self.Mmin)
            log_general_parameters_file.write('M_tronc\t'+str(M_min)+'\n')
            # if self.sample == 1 :
            #     b_value = (self.b_min + self.b_max)/2.
            # else :
            #     b_value = np.random.triangular(self.b_min,(self.b_min + self.b_max)/2.,self.b_max)
            mfd_param = {}
            mfd_param.update({'b_value' : self.b_value})

            log_general_parameters_file.write('b_value\t'+str(self.b_value)+'\n')
            if self.mfd_hyp == "double_GR" :
                mfd_param.update({'Mrupt' : self.})
                log_general_parameters_file.write('Mrupt\t'+str(Mrupt)+'\n')

            if self.mfd_hyp == "YC_modified" :
                Mf = mfd_param['Mf'] #Mmax with the best a priori fit, value for the Marmara sea
                size_of_bump = mfd_param['size_of_bump'] #in order to modify the respective size of the two parts of the
                log_general_parameters_file.write('Mf\t'+str(Mf)+'\n')
                log_general_parameters_file.write('size_of_bump\t'+str(size_of_bump)+'\n')

            log_general_parameters_file.close()

            #initialization of array for faults information
            faults_area = []
            faults_length = []
            faults_width = []
            faults_slip_rates = []
            faults_mecanism = []
            faults_shear_mod = []

            '''########################################################
            # Random sampling of the fault slip-rate
            ########################################################'''

            ########################################################
            # find which faults interact a lot with each other
            ########################################################
            print("Picking slip-rates...")
            if self.sample == 1 :
                self.sr_correl = False

            if self.sr_correl ==True :
                M_faults_correl = []
                for Fault_name in faults_names:
                    M_faults_correl_i =  []
                    sc_with_fault = []
                    for index_scenario in range(len(scenarios_names)) :
                        scenario_i = self.rupture_set[index_scenario]
                        if Fault_name in scenario_i :
                            sc_with_fault += scenario_i
                    for Fault_name_i in faults_names:
                        nb_join_ruptures = 0
                        if Fault_name == Fault_name_i :
                            nb_join_ruptures = 0
                        elif Fault_name_i in sc_with_fault:
                            nb_join_ruptures = sc_with_fault.count(Fault_name_i)
                        else :
                            nb_join_ruptures = 0
                        M_faults_correl_i.append(nb_join_ruptures)
                    M_faults_correl.append(M_faults_correl_i)

                #calculating the linked level (up to five)
                M_linked_lvl = np.zeros_like(M_faults_correl)
                for lvl in [1,2,3,4]:
                    for i in range(len(faults_names)):
                        for j in range(len(faults_names)):
                            if M_linked_lvl[i][j] == 0:
                                if M_faults_correl[i][j] != 0:
                                    M_linked_lvl[i][j]=lvl
                                else :
                                    for jj in  range(len(M_faults_correl[j])): #for all faults connecting to a fault connected to i
                                        if M_faults_correl[j][jj]!=0:
                                            if M_linked_lvl[i][j] == 0:
                                                M_linked_lvl[i][j]=lvl+1
                #print(M_linked_lvl)

                max_correl = np.max(M_faults_correl)

                list_quater_picked = list(np.zeros(len(faults_names)))
                # pick the faults slip rates
                index_fault = 0
                for Fault_name in faults_names:
                    faults_length.append(faults_data[index_fault]['length'])
                    faults_area.append(faults_data[index_fault]['area'])
                    faults_width.append(faults_data[index_fault]['width'])
                    faults_mecanism.append(faults_data[index_fault]['mecanism'])
                    sr_values = [faults_data[index_fault]['slip_rate_min'],
                    faults_data[index_fault]['slip_rate_moy'],
                    faults_data[index_fault]['slip_rate_max']]

                    #selects randomly the slip rate
                    slip_rate = select_sr.select(sr_values,self.sample,index_fault,M_linked_lvl[index_fault],list_quater_picked)

                    log_line = str(Fault_name) + '\t' + str(slip_rate) + '\n' #writting in the log file
                    log_sr_file.write(log_line)

                    faults_slip_rates.append(slip_rate * 0.001) # transform from mm to m

    #                faults_shear_mod.append(float(self.shear_mod)*10**9 )
                    faults_shear_mod.append(faults_data[index_fault]['shear_mod'])
                    index_fault += 1

            else: # the option of correlation bewteen corroption fault not used
                index_fault = 0
                for Fault_name in faults_names:
                    faults_length.append(faults_data[index_fault]['length'])
                    faults_area.append(faults_data[index_fault]['area'])
                    faults_width.append(faults_data[index_fault]['width'])
                    faults_mecanism.append(faults_data[index_fault]['mecanism'])
                    ########################################################
                    #ramdom sampling the slip rate  . uniform sampling
                    ########################################################
                    slip_rate_min = faults_data[index_fault]['slip_rate_min']
                    slip_rate_moy = faults_data[index_fault]['slip_rate_moy']
                    slip_rate_max = faults_data[index_fault]['slip_rate_max']

                    if self.sample == 1 :
                        slip_rate = slip_rate_moy
                    else :
                        slip_rate_inf = np.random.uniform(slip_rate_min,slip_rate_moy)
                        slip_rate_sup = np.random.uniform(slip_rate_moy,slip_rate_max)
                        slip_rate = np.random.choice([slip_rate_inf,slip_rate_sup])

                    log_line = str(Fault_name) + '\t' + str(slip_rate) + '\n' #writting in the log file
                    log_sr_file.write(log_line)

                    faults_slip_rates.append(slip_rate * 0.001) # transform from mm to m

    #                faults_shear_mod.append(float(self.shear_mod)*10**9 )
                    faults_shear_mod.append(faults_data[index_fault]['shear_mod'])

                    index_fault += 1
            print("\t\tslip-rates picked.")

            ratio_test = 0.5
            count_reruns = 1 #used to divide the sr increment if the fit is not good
            count_mfd90 = 1

            f_pkl_mdf = self.pathlog +'/mfd_' + str(self.sample) + '.pkl'
            re_use_mfd_pkl= True
            if not os.path.isfile(f_pkl_mdf):
                re_use_mfd_pkl = False
            if re_use_mfd_pkl == False:
                while abs(ratio_test-1) >self.fit_quality or math.isnan(ratio_test) == True:
                    MFDs = EQ_on_faults.EQ_on_faults_from_sr(M_min,
                                                            mfd_param,
                                                            faults_names,
                                                            faults_area,
                                                            faults_length,
                                                            faults_width,
                                                            faults_slip_rates,
                                                            scenarios_names,
                                                            faults_shear_mod,
                                                            self.path,
                                                            self.sample,
                                                            self.selected_ScL,
                                                            self.dimention_used,
                                                            self.use_all_ScL_data,
                                                            faults_mecanism,
                                                            self.bg_ratio,
                                                            self.size_of_increment,
                                                            self.mfd_hyp,
                                                            count_reruns,
                                                            faults_lon,
                                                            faults_lat,
                                                            self.Mmax_range,
                                                            self.calculation_log_file,
                                                            self.branch,
                                                            self.param)
                    ratio_test = MFDs.ratio_test
                    if abs(ratio_test-1) > self.fit_quality:
                        print('bad sampling => re-run')
                        count_reruns += 1
                        if MFDs.ratio_NMS> 90.:
                            count_reruns -= 1
                            count_mfd90 +=1
                        else :
                            count_mfd90 =1
                    else:
                        if MFDs.ratio_NMS> 90.:
                            ratio_test = 0.5
                            print('bad sampling => re-run')
                            count_mfd90 +=1
                        else :
                            count_mfd90 =1

                    if math.isnan(ratio_test) == True:
                        print('bad sampling => re-run')
                        count_reruns = 1
                    if count_reruns > 3 or count_mfd90 > 3:
                        print('\n\n\n!!!!!! maybe there is a problem!!!')
                        ratio_test = 1.

                with open(f_pkl_mdf, 'wb') as f:
                    MFDs_to_pkl = [
                                    MFDs.EQ_rate_BG,
                                    MFDs.faults_names,
                                    MFDs.OQ_entry_faults,
                                    MFDs.scenarios_names,
                                    MFDs.OQ_entry_scenarios,
                                    MFDs.index_faults_in_scenario
                                    ]
                    pickle.dump(MFDs_to_pkl, f)

            else :
                print('Reloading MFDs from previous run')
                with open(f_pkl_mdf, 'rb') as f:
                    MFDs_to_pkl = pickle.load(f)

            EQ_rate_BG = MFDs_to_pkl[0]
            faults_names = MFDs_to_pkl[1]
            OQ_entry_faults = MFDs_to_pkl[2]
            scenarios_names = MFDs_to_pkl[3]
            OQ_entry_scenarios = MFDs_to_pkl[4]
            index_faults_in_scenario = MFDs_to_pkl[5]

            # sclaling law as called by openquake
            if self.selected_ScL == 'Le2010' :
                ScL_oq = 'Leonard2014_SCR'
            if self.selected_ScL == 'WC1994' :
                ScL_oq = 'WC1994'
            else : #default Scaling relationship for opanquake
                ScL_oq = 'WC1994'

            # Extract the background geometry
            Lon_bg, Lat_bg  = bg.geom(self.Model_name,self.File_bg )

            #Create the BG polygon
            Poly_bg = []
            for x1,y1 in zip(Lon_bg,Lat_bg):
                Poly_bg.append((x1,y1))
            bbPath_BG = mplPath.Path(Poly_bg)

            # find the faults in the bg
            f_in_bg = []
            for index_fault,fault_name in zip(range(len(faults_names)),faults_names):
                for lon,lat in zip(faults_data[index_fault]['lon'],faults_data[index_fault]['lat']):
                    if bbPath_BG.contains_point((lon,lat)):
                        f_in_bg.append(index_fault)

            ########################################################
            #loop on simple faults
            ########################################################
            ID_number = 0
            for index_fault,fault_name in zip(range(len(faults_names)),faults_names):
                if index_fault in f_in_bg:
                    use_simple_faults = False
                    if use_simple_faults == True :
                        line,self.Domain_in_the_model,ID_number =fault_source.write_simple_fault(index_fault,fault_name,
                        OQ_entry_faults,faults_names,faults_data,self.Model_name,
                        self.Domain_in_the_model,ScL_oq,log_mdf_file,M_min,ID_number)
                    else :
                        line,self.Domain_in_the_model,ID_number =fault_source.write_non_parametric_one_fault(index_fault,fault_name,
                        OQ_entry_faults,faults_names,faults_data,self.Model_name,
                        self.Domain_in_the_model,ScL_oq,log_mdf_file,M_min,ID_number,explo_time)


                    if cut_sm_file == False :
                        XMLfile.write(line)
                    else :
                        #find index to write in
                        i_w = 0
                        while sf_counter[i_w] > 500 :
                            i_w += 1
                        sf_files[i_w].write(line)
                        sf_counter[i_w]+=1

            if len(self.rupture_set) != 0 :
                for scenario in enumerate(scenarios_names):

                    # check if the scenario has at least one fault in the bg
                    sc_in = False
                    index_scenario = np.where(np.array(scenarios_names) == scenario[1])[0][0]
                    for i in index_faults_in_scenario[index_scenario][0]:
                        if i in f_in_bg:
                            sc_in = True

                    if sc_in == True:
                        use_non_param = True
                        if use_non_param == False:
                            line,ID_number = fault_source. write_characteristic_scenario(scenarios_names,
                                                                                        OQ_entry_scenarios,
                                                                                        index_faults_in_scenario,
                                                                                        scenario,
                                                                                        faults_names,
                                                                                        self.Model_name,
                                                                                        faults_data,
                                                                                        log_mdf_file,
                                                                                        M_min,
                                                                                        ID_number)
                        else :

                            line,ID_number = fault_source. write_non_parametric_source(scenario,
                                                                                        scenarios_names,
                                                                                        OQ_entry_scenarios,
                                                                                        index_faults_in_scenario,
                                                                                        faults_names,
                                                                                        faults_data,
                                                                                        self.Model_name,
                                                                                        self.Domain_in_the_model,
                                                                                        ScL_oq,
                                                                                        log_mdf_file,
                                                                                        explo_time,
                                                                                        M_min,
                                                                                        ID_number)
                        if cut_sm_file == False :
                            XMLfile.write(line)
                        else :
                            #find index to write in
                            i_w = 0
                            while mf_counter[i_w] > 500 :
                                i_w += 1
                            mf_files[i_w].write(line)
                            mf_counter[i_w]+=1

        '''#########################
        # Defining the background seismicity
        #########################'''
        MFD = EQ_rate_BG
        # if self.param["main"]["background"]["option_bg"] in ["None","none"]:
        #     do_bg_in_SHERIFS = False
        # else :
        #     do_bg_in_SHERIFS = True

        if sum(MFD) != 0. and self.param["main"]["background"]["option_bg"] == "zone":
            bg_file = open(self.path +'/Source_model_bg_' + str(self.sample) + '.xml','w')
            list_src_files.append(bg_file)
            # Initiate the xml file
            line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
            line+='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
            line+='\txmlns="http://openquake.org/xmlns/nrml/0.5">\n'
            line+='\t<sourceModel name="Hazard Model">\n'

            upperSeismoDepth, lowerSeismoDepth, ruptAspectRatio, nodalPlanes, hypoDepths = bg.prop(self.Model_name,self.file_prop_bg)
            line+='\t\t<sourceGroup\nname="group 2"\nrup_interdep="indep"\nsrc_interdep="indep"\ntectonicRegion="' + str(self.Domain_in_the_model[0]) + '"\n>\n'
            line+='\t\t<areaSource id="'+ str(source_id + 1 ) +'" name="Background" tectonicRegion="' + str(self.Domain_in_the_model[0]) + '">\n'
            line+='\t\t\t<areaGeometry>\n'
            line+='\t\t\t\t<gml:Polygon>\n'
            line+='\t\t\t\t\t<gml:exterior>\n'
            line+='\t\t\t\t\t\t<gml:LinearRing>\n'
            line+='\t\t\t\t\t\t\t<gml:posList>\n'
            for x,y in zip(Lon_bg,Lat_bg):
                #polygon.append((x,y)) #ecriture du polygone de la zone
                line+='\t\t\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
            line+='\t\t\t\t\t\t\t</gml:posList>\n'
            line+='\t\t\t\t\t\t</gml:LinearRing>\n'
            line+='\t\t\t\t\t</gml:exterior>\n'
            line+='\t\t\t\t</gml:Polygon>\n'
            line+='\t\t\t\t<upperSeismoDepth>' + str(upperSeismoDepth) + '</upperSeismoDepth>\n'
            line+='\t\t\t\t<lowerSeismoDepth>' + str(lowerSeismoDepth) + '</lowerSeismoDepth>\n'
            line+='\t\t\t</areaGeometry>\n'
            line+='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
            line+='\t\t\t<ruptAspectRatio>' + str(ruptAspectRatio) + '</ruptAspectRatio>\n'
            log_mdf_file.write('Background' + '\t' + str(M_min) + '\t' + ' '.join(list(map(str, MFD)))+'\n')
            line+='\t\t\t<incrementalMFD binWidth=\"0.10\" minMag="'+ str(M_min)+'">\n'
            line+='\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
            line+='\t\t\t</incrementalMFD>\n'
            line+='\t\t\t<nodalPlaneDist>\n'
            for i in range(len(nodalPlanes)) :
                line+=('\t\t\t\t<nodalPlane probability="' + str(nodalPlanes[i][0]) + '" strike="' + str(nodalPlanes[i][1]) + '" dip="'
                      + str(nodalPlanes[i][2]) + '" rake="' + str(nodalPlanes[i][3]) + '" />\n')
            line+='\t\t\t</nodalPlaneDist>\n'
            line+='\t\t\t<hypoDepthDist>\n'
            for i in range(len(hypoDepths)) :
                line+='\t\t\t\t<hypoDepth probability="' + str(hypoDepths[i][0]) + '" depth="' + str(hypoDepths[i][1]) + '" />\n'
            line+='\t\t\t</hypoDepthDist>\n'
            line+='\t\t</areaSource>\n'
            line+='\t</sourceGroup>\n'
            line+='\t</sourceModel>\n'
            line+='</nrml>\n'

            bg_file.write(line)
            bg_file.close()

        elif sum(MFD) != 0. and self.param["main"]["background"]["option_bg"]=="smooth":
            Mmin_checked = False
            # read the xml and stores the list of aValues

            list_bg_xml = self.list_fbg
            if os.path.isdir(self.fbgpath):
                list_bg_xml = [self.fbgpath+"/"+i for i in list_bg_xml]
            pts_list = {}
            sum_rates = 0.

            for fbg in list_bg_xml:
                tree = ET.parse(fbg)
                ET.register_namespace('', "http://openquake.org/xmlns/nrml/0.5")
                nrml = tree.getroot()
                i_point =0
                for pointSource in nrml[0][0]:
                    pt_in_BG = False
                    i_child = 0
                    for child in pointSource.getchildren():
                        if "pointGeometry" in str(child) :
                            s_tmp = pointSource[i_child][0][0].text
                            s_tmp=s_tmp.replace('\n','')
                            s_tmp=[float(i) for i in s_tmp.split(' ') if i != '']
                            if bbPath_BG.contains_point((s_tmp[0],s_tmp[1])) == 1:
                                pt_in_BG = True
                            str_loc = str(s_tmp[0])+'_'+str(s_tmp[1])
                        if "truncGutenbergRichterMFD" in str(child) :
                            aValue = nrml[0][0][i_point][i_child].get('aValue')
                            minMag = nrml[0][0][i_point][i_child].get("minMag")
                            i_trGR = i_child
                        i_child+=1

                    if pt_in_BG == True :
                        aValue = nrml[0][0][i_point][i_trGR].get('aValue')
                        minMag = nrml[0][0][i_point][i_trGR].get('minMag')
                        bValue = nrml[0][0][i_point][i_trGR].get('bValue')
                        maxMag = nrml[0][0][i_point][i_trGR].get('maxMag')
                        pts_list.update({str_loc:{"aValue":aValue,
                        "bValue":bValue,
                        "maxMag":maxMag,
                        "minMag":minMag}})
                        sum_rates += float(10.**float(aValue))

                        if Mmin_checked == False :
                            if float(minMag) < M_min :
                                print("!!!!")
                                print("WARNING : BG has a small Mmin than the SHERIFS input")
                                print("!!!!")
                                Mmin_checked = True
                    i_point +=1

            # Normalize to fit the BG MFD
            Mmax = M_min+len(EQ_rate_BG)*0.1-1
            mags = np.linspace(M_min,Mmax,len(EQ_rate_BG))

            i_bg = 0
            for fbg in list_bg_xml:
                tree = ET.parse(fbg)
                ET.register_namespace('', "http://openquake.org/xmlns/nrml/0.5")
                nrml = tree.getroot()
                i_point =0
                for pointSource in nrml[0][0]:
                    pt_in_BG = False
                    i_child = 0
                    for child in pointSource.getchildren():
                        if "magScaleRel" in str(child) :
                            nrml[0][0][i_point][i_child].text = ScL_oq
                        if "pointGeometry" in str(child) :
                            s_tmp = pointSource[i_child][0][0].text
                            s_tmp=s_tmp.replace('\n','')
                            s_tmp=[float(i) for i in s_tmp.split(' ') if i != '']
                            if bbPath_BG.contains_point((s_tmp[0],s_tmp[1])) == 1:
                                pt_in_BG = True
                            str_loc = str(s_tmp[0])+'_'+str(s_tmp[1])
                        if "truncGutenbergRichterMFD" in str(child) :
                            i_trGR = i_child
                        i_child +=1
                    if pt_in_BG == True :
                        b_value = float(pts_list[str_loc]["bValue"])
                        a_value = float(pts_list[str_loc]["aValue"])

                        attrib = {"minMag":str(M_min), "binWidth":"0.10"}
                        element = nrml[0][0][i_point].makeelement('incrementalMFD', attrib)
                        nrml[0][0][i_point].append(element)
                        element = nrml[0][0][i_point][-1].makeelement('occurRates', {})
                        nrml[0][0][i_point][-1].append(element)
                        str_tmp = " "
                        i_mag = 0
                        for mag in mags:
                            mag_lo = mag - 0.05
                            mag_hi = mag + 0.05
                            r = (10 ** (a_value - b_value * mag_lo)
                                - 10 ** (a_value - b_value * mag_hi))
                            norm_r = r * ((10.**a_value) / sum_rates)
                            str_tmp += str(EQ_rate_BG[i_mag]*norm_r)
                            str_tmp += " "

                        nrml[0][0][i_point][-1][0].text =str_tmp
                        nrml[0][0][i_point].remove(nrml[0][0][i_point][i_trGR])
                    i_point+=1
                i_bg += 1
                fbg_out = self.path + '/bg_'+str(self.sample)+'_'+str(i_bg)+'.xml'
                tree.write(fbg_out)
                list_src_files.append(fbg_out)



        '''#############################
        ### defining the other sources based on the host model
        ##############################'''
        if self.param["main"]["background"]["use_host_model"] == True and cut_sm_file == False:
            host_model.build(XMLfile,self.host_model_file,Lon_bg,Lat_bg)
        if self.param["main"]["background"]["use_host_model"] == True and cut_sm_file == True:
            print("WARNING : can't use host model and cut files yet !")

        if use_multiF == False :
            #end of the file
            line='\t\t\t</sourceGroup>\n'
            line+='\t</sourceModel>\n'
            line+='</nrml>\n'

            if cut_sm_file == False :
                XMLfile.write(line)
                XMLfile.close()
            else :
                for f in sf_files + mf_files  :
                    f.write(line)
                    f.close()

        log_sr_file.close()

        log_mdf_file.close()

        self.list_src_files = list_src_files

        '''#######################
        ### some figures
        ######################'''

    if self.param["figures"]["print"] in ["true","True"]:
        make_figures = True

    if make_figure == True :
        if self.param["figures"]["model_mfd"] in ["true","True"]:
            plt_model_mfd = True
