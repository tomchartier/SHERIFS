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
        
        last_bg = 'impossible_name'
        last_set = 'impossible_name'
        last_model = 'impossible_name'
        
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

                    if last_bg != BG_hyp :
                        # setting the ration of seismicity that is in the background
                        bg_ratio = available_bg[BG_hyp]
                        last_bg = BG_hyp
                    if last_set != scenario_set :
                        # extracting the complexe multi fault ruptures
                        rupture_set = available_sets[scenario_set]
                        
                        index_scenario = 0
                        scenarios_names = []
                        if np.size(rupture_set) == 0 :
                            scenarios_names = []
                        else :
                            for index_scenario in range(len(rupture_set)):
                                faults_in_scenario = rupture_set[index_scenario]
                                if len(faults_in_scenario) > 1:
                                    scenario = {}
                                    faults_done = []
                                    for i in range(len(faults_in_scenario)):
                                        if not str(faults_in_scenario[i]).replace('\r','') in faults_done:
                                            scenario["f_%s" % str(i+1)] = [str(faults_in_scenario[i]).replace('\r','').replace('\t','').replace('\n','')]
                                            faults_done.append(str(faults_in_scenario[i]).replace('\r','').replace('\t','').replace('\n',''))
                                    if len(scenario)!=0:
                                        scenarios_names.append(scenario)
                                index_scenario += 1
                                
                        last_set = scenario_set
                        
                    if last_model != Model :
                        # Extraction of the fault geometry and properties
                        last_model = Model
                        print("Importing fault data")

                        #Extraction of the faults and scenarios present in the model from the text file
                        Prop = np.genfromtxt(self.File_prop,
                                                   dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
                                                          ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
                        Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
                        Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
                        index_model = np.where(np.array(Column_model_name) == Model)[0]
                        Prop = np.take(Prop,index_model)
                        faults_names = np.array(Column_fault_name[index_model[0]:index_model[-1]+1])
                        faults_names = list(faults_names)
                        


                    ########################################################
                    #Extraction of the properties and geometries of faults
                    ########################################################
                    faults_data = {}
                    index_fault = 0
                    #extractions of the geometries of the faults
                    geom_scenar = Geometry_scenario.Geom_scenar(faults_names,scenarios_names,self.File_geom,Model)
                    faults_lon = geom_scenar.faults_lon
                    faults_lat = geom_scenar.faults_lat
                    
                    self.FaultGeometry(Model)  #extract the geometries from the geometry file
                    
                    for Fault_name in faults_names:
                        # extract depth
                        i_d = np.where(np.array(self.Column_Fault_name) == Fault_name)
                        depth = list(map(lambda i : self.Depths[i],i_d[0]))
                        #extractions of the properties of the fault
                        self.FaultProperties(Fault_name,Model)
                        dip = self.dip
                        upper_sismo_depth = self.upper_sismo_depth
                        lower_sismo_depth = self.lower_sismo_depth
                        width = (lower_sismo_depth - upper_sismo_depth) / math.sin(math.radians(dip))
                        length = geom_scenar.length[index_fault] * 1000.
                        area = length * width * 1000.
                        
                        if self.rake> -135. and self.rake< -45:
                            mecanism = 'N'
                        elif self.rake< 135. and self.rake> 45:
                            mecanism = 'R'
                        else :
                            mecanism = 'S'
                            
                        slip_rate_min = self.slip_rate_min
                        slip_rate_moy = self.slip_rate_moy
                        slip_rate_max = self.slip_rate_max

                        faults_data.update({index_fault:{'name':Fault_name,
                        'dip':dip,
                        'oriented':self.oriented,
                        'upper_sismo_depth':upper_sismo_depth,
                        'lower_sismo_depth':lower_sismo_depth,
                        'width':width,'length':length,'area':area,
                        'mecanism':mecanism,'rake':self.rake,
                        'slip_rate_min':slip_rate_min,
                        'slip_rate_moy':slip_rate_moy,
                        'slip_rate_max':slip_rate_max,
                        'shear_mod':float(self.shear_mod)*10**9,
                        'domain':self.Domain,
                        'lon':faults_lon[index_fault],
                        'lat':faults_lat[index_fault],
                        'depth':depth}})
                        index_fault += 1
    
    
    
                    if str_all_data == 'a' :
                        use_all_ScL_data = True
                    elif str_all_data == 'm' :
                        use_all_ScL_data = False
                        
                    if rerun_the_files == True :
                        # Create the source model
                        Source_model = Source_Model_Creator(path,self.Run_Name,
                                                            Model,
                                                            self.File_geom,
                                                            self.File_prop,
                                                            self.File_bg,
                                                            self.file_prop_bg,
                                                            rupture_set,
                                                            self.Domain_in_model,
                                                            sample,self.seed,
                                                            self.Mmin,
                                                            selected_ScL,
                                                            dim_used,
                                                            use_all_ScL_data,
                                                            b_min,
                                                            b_max,
                                                            mfd_hyp[4::],
                                                            bg_ratio,
                                                            self.sr_correl,
                                                            self.size_of_increment,
                                                            self.fit_quality,
                                                            self.Mmax_range,
                                                            self.calculation_log_file,
                                                            self.use_host_model,
                                                            self.host_model_file,
                                                            faults_names,
                                                            scenarios_names,
                                                            faults_data,
                                                            faults_lon,faults_lat)
                        
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
     

    def FaultProperties(self,Name_of_fault,Model):
        FileName_Prop = self.File_prop
        Prop = np.genfromtxt(FileName_Prop,
                                   dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
                                          ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
        Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
        Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
        index_model = np.where(np.array(Column_model_name) == Model)[0]
        
        Prop = np.take(Prop,index_model)
        index_fault = np.where(np.array(Column_fault_name[index_model[0]:index_model[-1]+1]) == Name_of_fault)
        Indexfault_final = index_fault[0]

        self.dip = Prop[Indexfault_final][0][2]
        self.oriented = Prop[Indexfault_final][0][3]
        self.rake = Prop[Indexfault_final][0][4]
        self.upper_sismo_depth = Prop[Indexfault_final][0][5]
        self.lower_sismo_depth = Prop[Indexfault_final][0][6]
        
        self.slip_rate_min = Prop[Indexfault_final][0][7]
        self.slip_rate_moy = Prop[Indexfault_final][0][8]
        self.slip_rate_max = Prop[Indexfault_final][0][9]
        self.Domain = Prop[Indexfault_final][0][10]
        self.shear_mod = Prop[Indexfault_final][0][11]


        if self.rake == 'N' :
            self.rake = -90.00
        if self.rake == 'S' :
            self.rake = 00.00
        if self.rake == 'SS' :
            self.rake = 00.00
        if self.rake == 'R' :
            self.rake = 90.00
        self.rake = float(self.rake)

        if len(str(self.dip)) == 0:
            print('\nError!!! please verify your input file for fault parameters\n')


    def FaultGeometry(self,Model):
        NomFichier_InfosZonage = self.File_geom
        InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8'),('U100')],skip_header = 1)
        Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
        index_model = np.where(np.array(Column_model_name) == Model)
        self.Column_Fault_name = list(map(lambda i : InfosZonage[i][1],index_model[0]))
        self.Longitudes = list(map(lambda i : InfosZonage[i][2],index_model[0]))
        self.Latitudes = list(map(lambda i : InfosZonage[i][3],index_model[0]))
        self.Depths = list(map(lambda i : InfosZonage[i][4],index_model[0]))

        ZoneSelec = self.Column_Fault_name
        DicoZone = dict([(k,ZoneSelec.count(k)) for k in set(ZoneSelec)])
        Longitudes = []
        Latitudes = []
        Depths = []
        Column_Fault_name = []
        for cle in DicoZone.keys():
            indices_ZonesSelec = np.where(np.array(self.Column_Fault_name) == cle)
            ColonneNomZone_inter = np.take(self.Column_Fault_name,indices_ZonesSelec)
            Longitudes_inter = np.take(self.Longitudes,indices_ZonesSelec)
            Latitudes_inter = np.take(self.Latitudes,indices_ZonesSelec)
            depth_inter = np.take(self.Depths,indices_ZonesSelec)

            Longitudes_inter = Longitudes_inter[0].tolist()
            Latitudes_inter = Latitudes_inter[0].tolist()
            depth_inter = depth_inter[0].tolist()
            ColonneNomZone_inter = ColonneNomZone_inter[0].tolist()
            compt = 0
            for xx,yy,nn,dd in zip(Longitudes_inter,Latitudes_inter,ColonneNomZone_inter,depth_inter):
                compt+=1
                Longitudes.append(xx)
                Latitudes.append(yy)
                Depths.append(dd)
                Column_Fault_name.append(nn)

        self.Longitudes =Longitudes
        self.Latitudes =Latitudes
        self.Depths =Depths
        self.Column_Fault_name = Column_Fault_name
        self.Nb_data_per_zone = dict([(k,self.Column_Fault_name.count(k)) for k in set(self.Column_Fault_name)])
        self.Fault_Names = sorted(self.Nb_data_per_zone.keys())
