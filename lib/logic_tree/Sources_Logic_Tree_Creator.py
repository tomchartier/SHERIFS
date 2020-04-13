# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

This code create a source logic tree 

@author: Thomas Chartier
"""
import xml.etree.ElementTree as ET
import numpy as np
import os   
import sys
from shutil import copyfile

path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/lib'
sys.path.append(path_lib)
path_f = path_lib + '/file_writer'
sys.path.append(path_f)
path_f = path_lib + '/utils'
sys.path.append(path_f)

from Source_Model import *
#import Source_LT
import BG_ratio
import faults_n_scenarios
from OQ_job_Creator import OQ_job_Creator
import read_input


class Sources_Logic_Tree_Creator:
    def __init__(self,Run_Name,File_geom,File_prop,File_bg,file_prop_bg,Domain_in_model,
                 nb_random_sampling,seed,Mmin,sr_correl,size_of_increment,Mmax_range,
                 overwrite_files,fit_quality,calculation_log_file,use_host_model,host_model_file):
        self.Run_Name = Run_Name
        self.File_geom = File_geom
        self.File_prop = File_prop
        self.File_bg = File_bg
        self.file_prop_bg = file_prop_bg
        
        
        self.Domain_in_model = Domain_in_model
        
        self.nb_random_sampling = nb_random_sampling
        
        
        self.seed = seed
        np.random.seed = seed
        
        self.Mmin = Mmin
        self.sr_correl = sr_correl
        self.size_of_increment = size_of_increment
        self.Mmax_range = Mmax_range
        
        self.overwrite = overwrite_files
        self.fit_quality = fit_quality
        
        self.calculation_log_file = calculation_log_file
        self.use_host_model = use_host_model
        self.host_model_file = host_model_file
        
        self.initialize()

    def initialize(self):
        LT_file = str(self.Run_Name)+'/Sources_Logic_tree.xml'
        LT_log_name  =  str(self.Run_Name)+'/LT_log.txt' 
#        reading_file = False
        if not os.path.exists(LT_file) :
            print('ERROR : Please provide a LT_log.txt file \n See the user amnual for guidelines and the example for file setup example.')
            exit()
        
        else : #get from the xml file
            
#            reading_file = True
    
            file_log_LT = open(LT_log_name,'r')
            self.log_LT = file_log_LT.readlines()
            
            OQ_job = OQ_job_Creator(self.Run_Name) # ask the info about the run and create the job.ini file
            
            nb_random_sampling = OQ_job.nb_sample               
            
            selected_Model = self.log_LT[1].split('\t')
            if '\r\n' in selected_Model:
                selected_Model.remove('\r\n')
            if '\n' in selected_Model:
                selected_Model.remove('\n')
                
            selected_ScL = self.log_LT[3].split('\t')
            if '\r\n' in selected_ScL:
                selected_ScL.remove('\r\n')
            if '\n' in selected_ScL:
                selected_ScL.remove('\n')
            
                
            index_advance = 0
            mfd_hyps = []
            b_values_hyps = []
            while self.log_LT[5+index_advance][0:3] == 'MFD':
                mfd_hyps.append(self.log_LT[5+index_advance].split('\t')[0])
                b_values_hyps_i = []
                for b_hyp in self.log_LT[5+index_advance].split('\t')[1:]:
                    b_values_hyps_i.append(b_hyp)
                if '\r\n' in b_values_hyps_i:
                    b_values_hyps_i.remove('\r\n')
                if '\n' in b_values_hyps_i:
                    b_values_hyps_i.remove('\n')
                b_values_hyps.append(b_values_hyps_i)
                index_advance += 1
                
            # Background
            bg_names = self.log_LT[6+index_advance].split('\t')
            if '\r\n' in bg_names:
                bg_names.remove('\r\n')
            if '\n' in bg_names:
                bg_names.remove('\n')
            if '' in bg_names:
                bg_names.remove('')
                
            available_bg = read_input.extract_bg_input('input/'+self.Run_Name+'/bg_seismicity.txt')
            
            # Scenario Set
            sc_names = self.log_LT[8+index_advance].split('\t')
            if '\r\n' in sc_names:
                sc_names.remove('\r\n')
            if '\n' in sc_names:
                sc_names.remove('\n')
            if '' in sc_names:
                sc_names.remove('')
                
            available_sets = read_input.extract_sc_input('input/'+self.Run_Name+'/ruptures.txt')
            
            # Build branches
            branches = []
            for model_i in selected_Model: 
                index_mfd = 0
                for mfd_i in mfd_hyps:               
                    for bvalue in b_values_hyps[index_mfd]:            
                        for bg_hyp_i in bg_names:             
                            for sc_name in sc_names:  
                                for ScL_i in selected_ScL :  
                                    ScL_i = ScL_i.split(' ')
                                    ScL_name_i= ScL_i[0]
                                    use_all_i = ScL_i[2][0]
                                    dim_i = ScL_i[1]
                    
                                    branch_i = [model_i,ScL_name_i,use_all_i,dim_i,mfd_i,bvalue,bg_hyp_i,sc_name]
                                    branches.append(branch_i)
                    index_mfd += 1
                            
#
#
#        if reading_file == False :     # create a new XML file
#            line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
#            line+='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
#            line+='\txmlns="http://openquake.org/xmlns/nrml/0.4">\n'
#            line+='\t<logicTree logicTreeID="lt1">\n'
#
#            str_all_data = []
#            id_number = 1
#            scenario_done = []
#            scenario_path = []
#
#            line+='\t\t<logicTreeBranchingLevel branchingLevelID="bl' + str(id_number) + '">\n'
#            line+='\t\t\t<logicTreeBranchSet uncertaintyType="sourceModel"\n'
#            line+='\t\t\t\t\t\t\tbranchSetID="bs' + str(id_number) + '">\n'
#            for branch in branches :
#                Model = branch[0]
#                selected_ScL = branch[1]
#                dim_used = branch[3][0]
#
#                if branch[2] == True :
#                    str_all_data = 'a' # ' a ' is for 'all data is used'
#                else :
#                    str_all_data = 'm' # ' m ' is for 'mechanic specific data only'
#
#                b_min = float(branch[5].split('_')[0])
#                b_max = float(branch[5].split('_')[1])
#                mfd_hyp = str(branch[4])
#                BG_hyp = branch[6]
#                scenario_set = branch[7]
#
#                for sample in range(1,self.nb_random_sampling+1):
#                    print
#                    print (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
#                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/'
#                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(mfd_hyp) + ' sample : ' + str(sample)) # name of the branch
#
#
#                    path = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
#                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/'
#                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(mfd_hyp)) # path to the source file
#
#                    if not os.path.exists(path):
#                        os.makedirs(path)
#
#
#                    # setting the ration of seismicity that is in the background
#                    bg_ratio_file = self.Run_Name + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/bg_ratio.txt'
#                    if not os.path.exists(bg_ratio_file):
#                        BG_r = BG_ratio.BG_ratio(self.Run_Name,Model,BG_hyp,bg_ratio_file)
#                        bg_ratio = BG_r.bg_ratio
#                    else :
#                        bg_ratio = np.genfromtxt(bg_ratio_file)
#
#
#
#                    File_faults_n_scenarios = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
#                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/faults_n_scenarios.txt')
#                    # file containing the faults included in the model
#
#                    ''' If needed, launch the grafical interface to select the faults and scenarios in the model'''
#                    if not str('model_' + str(Model) +'_sc_' +  str(scenario_set)) in scenario_done :
#                        f_n_s = faults_n_scenarios.selecFAULT_tk(self.Run_Name,Model,self.File_geom,File_faults_n_scenarios,scenario_set)
#                        path_scenario_set = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
#                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set))
#                        scenario_done.append(str('model_' + str(Model) +'_sc_' +  str(scenario_set)))
#                        scenario_path.append(path_scenario_set)
#                    else :
#                        if not os.path.exists(File_faults_n_scenarios) :
#                            index = np.where(np.array(scenario_done) == str('model_' + str(Model) +'_sc_' +  str(scenario_set)))[0]
#                            copyfile(scenario_path[index]+ '/faults_n_scenarios.txt',File_faults_n_scenarios)
#
#                    line+=('\t\t\t\t<logicTreeBranch branchID="' + str(Model) + '-' + 'bg_' + str(BG_hyp) + '-' + str(selected_ScL) + '-'
#                           + str(dim_used) + '-' + str_all_data + '-sc_' +  str(scenario_set) + '-'
#                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '-' + 'MFD_'+ str(mfd_hyp)  + '-s_' + str(sample) + '">\n')
#
#                    line+=('\t\t\t\t\t<uncertaintyModel>' + (str(Model) + '/' + 'bg_' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
#                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set) + '/'
#                            + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'sMFD_'+ str(mfd_hyp)) + '/Source_model_'
#                    + str(sample) + '.xml</uncertaintyModel>\n')
#
#                    line+='\t\t\t\t\t<uncertaintyWeight>' + str(1./ float(nb_model_in_LT)) + '</uncertaintyWeight>\n'
#
#
#                    Source_model = Source_Model_Creator(path,self.Run_Name,
#                                                        Model,self.File_geom,self.File_prop,self.File_bg,self.file_prop_bg,File_faults_n_scenarios,
#                                                        self.Domain_in_model,sample,self.seed,self.Mmin,selected_ScL,
#                                                        dim_used,branch[2],b_min,b_max,mfd_hyp,bg_ratio,self.sr_correl,self.size_of_increment,self.fit_quality,
#                                                        self.Mmax_range,self.calculation_log_file,self.use_host_model,self.host_model_file)  #create the source modele
#
#                    self.Domain_in_model = Source_model.Domain_in_the_model
#
#                    line+='\t\t\t\t</logicTreeBranch>\n'
#
#                    id_number += 1
#
#            line+='\t\t\t</logicTreeBranchSet>\n'
#            line+='\t\t</logicTreeBranchingLevel>\n'
#            line+='\t</logicTree>\n'
#            line+='</nrml>\n'
#
#            XMLfile=open(LT_file,'w')
#            XMLfile.write(line)
#            XMLfile.close()
#
#
#
#        else : #read the logic tree log file
#
        str_all_data = []
        id_number = 1
        scenario_done = []
        scenario_path = []
        # writting the xml file for the logic tree
        line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
        line+='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
        line+='\txmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        line+='\t<logicTree logicTreeID="lt1">\n'
        line+='\t\t<logicTreeBranchingLevel branchingLevelID="bl' + str(id_number) + '">\n'
        line+='\t\t\t<logicTreeBranchSet uncertaintyType="sourceModel"\n'
        line+='\t\t\t\t\t\t\tbranchSetID="bs' + str(id_number) + '">\n'
        
        for branch in branches :
            # Branch info
            Model = branch[0]
            selected_ScL = branch[1]
            dim_used = branch[3][0]
            str_all_data = branch[2]
            bvalue = branch[5]
            mfd_hyp = str(branch[4])
            BG_hyp = branch[6]
            scenario_set = branch[7]
            b_min = bvalue.split('_')[1]
            b_max = bvalue.split('_')[3]

                
            if not len(Model)==0 or len(BG_hyp[3::])==0 or len(scenario_set[3::])==0 or len(mfd_hyp[4::])==0 or len(selected_ScL)==0  or len(dim_used)==0  or len(str_all_data)==0 :
                path = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp[3::]) + '/' + str(selected_ScL) + '_'
                       + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set[3::]) + '/'
                        + 'bmin_' + str(b_min) + '_bmax_' + str(b_max) + '/' + 'MFD_'+ str(mfd_hyp[4::])) # path to the source file
                
                for sample in range(1,nb_random_sampling+1):
                    rerun_the_files = False
                    if not os.path.exists(path):
                        os.makedirs(path)
                    path_file = (str(self.Run_Name) + '/' +(str(Model) + '/' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
                           + str(dim_used) + '_' + str_all_data + '/' +  str(scenario_set) + '/'
                            + str(bvalue)+ '/' + str(mfd_hyp)) + '/Source_model_'+ str(sample) + '.xml')
                    if not os.path.isfile(path_file):
                        rerun_the_files = True
                    if self.overwrite == True:
                        rerun_the_files = True
    
                    line+=('\t\t\t\t<logicTreeBranch branchID="' + str(Model) + '-' + str(BG_hyp) + '-' + str(selected_ScL) + '-'
                           + str(dim_used) + '-' + str_all_data + '-' +  str(scenario_set) + '-'
                            + str(bvalue)+ '-' + str(mfd_hyp)  + '-s_' + str(sample) + '">\n')
                    print('\nRunning logic tree branch:')
                    
                    self.calculation_log_file.write('\n\nRunning logic tree branch:')
                    print(str(Model) + '-' + str(BG_hyp) + '-' + str(selected_ScL) + '-'
                           + str(dim_used) + '-' + str_all_data + '-' +  str(scenario_set) + '-'
                            + str(bvalue)+ '-' + str(mfd_hyp)  + '-s_' + str(sample))
                    self.calculation_log_file.write('\n'+str(Model) + '-' + str(BG_hyp) + '-' + str(selected_ScL) + '-'
                           + str(dim_used) + '-' + str_all_data + '-' +  str(scenario_set) + '-'
                            + str(bvalue)+ '-' + str(mfd_hyp)  + '-s_' + str(sample))
                    line+=('\t\t\t\t\t<uncertaintyModel>' + (str(Model) + '/' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
                           + str(dim_used) + '_' + str_all_data + '/' +  str(scenario_set) + '/'
                            + str(bvalue)+ '/' + str(mfd_hyp)) + '/Source_model_'
                    + str(sample) + '.xml</uncertaintyModel>\n')
                    
                    line+='\t\t\t\t\t<uncertaintyWeight>' + str(round(float(len(branches)*nb_random_sampling),5)) + '</uncertaintyWeight>\n'
                    line+='\t\t\t\t</logicTreeBranch>\n'
                    
                    path = (str(self.Run_Name) + '/' + str(Model) + '/' + str(BG_hyp) + '/' + str(selected_ScL) + '_'
                           + str(dim_used) + '_' + str_all_data + '/' +  str(scenario_set) + '/'
                            + str(bvalue)+ '/' + str(mfd_hyp)) # path to the source file
                    #print path
#                    bg_ratio_file = self.Run_Name + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp[3::]) + '/bg_ratio.txt'
#                    if not os.path.exists(bg_ratio_file):
#                        BG_r = BG_ratio.BG_ratio(self.Run_Name,Model,BG_hyp[3::],bg_ratio_file)
#                        bg_ratio = BG_r.bg_ratio
#                    else :
#                        bg_ratio = np.genfromtxt(bg_ratio_file)


                    # setting the ration of seismicity that is in the background
                    bg_ratio = available_bg[BG_hyp]
                    
                    # extracting the complexe multi fault ruptures
                    rupture_set = available_sets[scenario_set]
                    
#                    File_faults_n_scenarios = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp[3::]) + '/' + str(selected_ScL) + '_'
#                           + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set[3::]) + '/faults_n_scenarios.txt')
#                    # file containing the faults included in the model
#
#                    ''' If needed, launch the grafical interface to select the faults and scenarios in the model'''
#                    if not os.path.exists(File_faults_n_scenarios) :
#                        if not str('model_' + str(Model) +'_sc_' +  str(scenario_set[3::])) in scenario_done :
#                            f_n_s = faults_n_scenarios.selecFAULT_tk(self.Run_Name,Model,self.File_geom,File_faults_n_scenarios,scenario_set[3::])
#                            path_scenario_set = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp[3::]) + '/' + str(selected_ScL) + '_'
#                               + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set[3::]))
#                            scenario_done.append(str('model_' + str(Model) +'_sc_' +  str(scenario_set[3::])))
#                            scenario_path.append(path_scenario_set)
#                        else :
#                            if not os.path.exists(File_faults_n_scenarios) :
#                                index = np.where(np.array(scenario_done) == str('model_' + str(Model) +'_sc_' +  str(scenario_set[3::])))[0]
#                                copyfile(scenario_path[index]+ '/faults_n_scenarios.txt',File_faults_n_scenarios)
#                    else :
#                        if not str('model_' + str(Model) +'_sc_' +  str(scenario_set[3::])) in scenario_done :
#                            path_scenario_set = (str(self.Run_Name) + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp[3::]) + '/' + str(selected_ScL) + '_'
#                               + str(dim_used) + '_' + str_all_data + '/sc_' +  str(scenario_set[3::]))
#                            scenario_done.append(str('model_' + str(Model) +'_sc_' +  str(scenario_set[3::])))
#                            scenario_path.append(path_scenario_set)
    
                    if str_all_data == 'a' :
                        use_all_ScL_data = True
                    elif str_all_data == 'm' :
                        use_all_ScL_data = False
                        
                    if rerun_the_files == True :
#                        Source_model = Source_Model_Creator(path,self.Run_Name,Model,
#                                                            self.File_geom,self.File_prop,self.File_bg,self.file_prop_bg,File_faults_n_scenarios,self.Domain_in_model,
#                                                            sample,self.seed,self.Mmin,selected_ScL,dim_used,
#                                                            use_all_ScL_data,b_min,b_max,mfd_hyp[4::],bg_ratio,self.sr_correl,self.size_of_increment,self.fit_quality,
#                                                            self.Mmax_range,self.calculation_log_file,self.use_host_model,self.host_model_file)  #create the source model
#
                        Source_model = Source_Model_Creator(path,self.Run_Name,Model,
                                                            self.File_geom,self.File_prop,self.File_bg,self.file_prop_bg,rupture_set,self.Domain_in_model,
                                                            sample,self.seed,self.Mmin,selected_ScL,dim_used,
                                                            use_all_ScL_data,b_min,b_max,mfd_hyp[4::],bg_ratio,self.sr_correl,self.size_of_increment,self.fit_quality,
                                                            self.Mmax_range,self.calculation_log_file,self.use_host_model,self.host_model_file)  #create the source model
                        
                        
                        self.Domain_in_model = Source_model.Domain_in_the_model
                    
                    id_number += 1

        line+='\t\t\t</logicTreeBranchSet>\n'
        line+='\t\t</logicTreeBranchingLevel>\n'
        line+='\t</logicTree>\n'
        line+='</nrml>\n'
        
        LT_file = str(self.Run_Name)+'/Sources_Logic_tree.xml'
        XMLfile=open(LT_file,'w')
        XMLfile.write(line)
        XMLfile.close()
     
