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
        LT_log_name  =  'input/'+str(self.Run_Name)+'/LT_log.txt'
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

            # setting the ration of seismicity that is in the background
            try:
                available_bg = read_input.extract_bg_input('input/'+self.Run_Name+'/bg_seismicity.txt')
            except:
                print('Error related to the background file \n'+
                'Please make sure input/run_name/bg_seismicity.txt is correctly set up')
                    
            
            # Scenario Set
            sc_names = self.log_LT[8+index_advance].split('\t')
            if '\r\n' in sc_names:
                sc_names.remove('\r\n')
            if '\n' in sc_names:
                sc_names.remove('\n')
            if '' in sc_names:
                sc_names.remove('')

            # extracting the complexe multi fault ruptures
            try:
                available_sets = read_input.extract_sc_input('input/'+self.Run_Name+'/ruptures.txt')
            except:
                print('Error related to the rupture scenario set file \n'+
                'Please make sure input/run_name/ruptures.txt is correctly set up')
                
            
            
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

                    # setting the ration of seismicity that is in the background
                    bg_ratio = available_bg[BG_hyp]
                    
                    # extracting the complexe multi fault ruptures
                    rupture_set = available_sets[scenario_set]
                    
    
                    if str_all_data == 'a' :
                        use_all_ScL_data = True
                    elif str_all_data == 'm' :
                        use_all_ScL_data = False
                        
                    if rerun_the_files == True :
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
     
