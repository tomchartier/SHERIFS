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

import matplotlib.pyplot as plt

class Source_Model_Creator:
    def __init__(self,path,Run_Name,Model_name,File_geom,File_prop,File_bg,file_prop_bg,rupture_set,
                 Domain_in_model,sample,seed,Mmin,selected_ScL,dimention_used,use_all_ScL_data,
                 b_min,b_max,mfd_hyp,bg_ratio,sr_correl,size_of_increment,fit_quality,Mmax_range,calculation_log_file,use_host_model,host_model_file,faults_names,scenarios_names,faults_data,faults_lon,faults_lat):
        self.Run_Name = Run_Name
        self.Domain_in_the_model = Domain_in_model #list of the different domain included in the model
        #a envoyer dans job
        self.File_geom = File_geom
        self.File_prop = File_prop
        self.File_bg = File_bg
        self.file_prop_bg = file_prop_bg
        self.rupture_set = rupture_set
        
        #a envoyer dans logic tree
        self.Model_name = Model_name
        self.b_min = float(b_min)
        self.b_max = float(b_max)
        self.mfd_hyp = mfd_hyp
        
        self.sample = sample
        
        self.path = path
        
        np.random.seed = seed
        self.Mmin = Mmin
        self.selected_ScL = selected_ScL        
        self.dimention_used = dimention_used
        self.use_all_ScL_data = use_all_ScL_data
        self.bg_ratio = bg_ratio
        self.sr_correl = sr_correl
        self.size_of_increment = size_of_increment
        self.fit_quality = fit_quality
        self.Mmax_range = Mmax_range
        
        self.calculation_log_file = calculation_log_file
        self.use_host_model = use_host_model
        self.host_model_file = host_model_file
        
        self.faults_names = faults_names
        self.scenarios_names = scenarios_names
        self.faults_data=faults_data
        self.faults_lon,self.faults_lat= faults_lon,faults_lat

        self.initialize()

    def initialize(self):
        
        faults_names = self.faults_names
        scenarios_names = self.scenarios_names
        faults_data = self.faults_data
        faults_lon,faults_lat= self.faults_lon,self.faults_lat
        #self.FaultGeometry()  #extract the geometries from the geometry file
        
        XMLfile=open(self.path +'/Source_model_' + str(self.sample) + '.xml','w')
        
        if not os.path.exists(self.path +'/Log'):
            os.makedirs(self.path +'/Log')
        log_sr_file=open(self.path +'/Log/slip_rate_sample_' + str(self.sample) + '.txt','w')
        log_mdf_file=open(self.path +'/Log/mdf_sample_' + str(self.sample) + '.txt','w')
        
        
        
        line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
        line+='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
        line+='\txmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        line+='\t<sourceModel name="Hazard Model">\n'
        line+='\n'
        XMLfile.write(line)

        #initialisation of the general parameters (M_min, shear modulus and b value)
        log_general_parameters_file = open(self.path +'/Log/general_parameters_sample_' + str(self.sample) + '.txt','w')
        
        M_min = float(self.Mmin)
        log_general_parameters_file.write('M_tronc\t'+str(M_min)+'\n')
        if self.sample == 1 :
            b_value = (self.b_min + self.b_max)/2.
        else :
            b_value = np.random.triangular(self.b_min,(self.b_min + self.b_max)/2.,self.b_max)
        mfd_param = {}
        mfd_param.update({'b_value' : b_value})
        log_general_parameters_file.write('b_value\t'+str(b_value)+'\n')
        
        log_general_parameters_file.close()
        
        #initialization of array for faults information
        faults_area = []
        faults_length = []
        faults_width = []
        faults_slip_rates = []
        faults_mecanism = []
        faults_shear_mod = []
                
#        ########################################################
#        #Extraction of the faults and scenarios present in the model from the text file
#        ########################################################
#
#        Prop = np.genfromtxt(self.File_prop,
#                                   dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
#                                          ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
#        Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
#        Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
#        index_model = np.where(np.array(Column_model_name) == self.Model_name)[0]
#        Prop = np.take(Prop,index_model)
#        faults_names = np.array(Column_fault_name[index_model[0]:index_model[-1]+1])
#        faults_names = list(faults_names)
        
        
#        lines_FtF = [line.rstrip('\n') for line in open(self.File_faults_n_scenarios)]
#        line_fauts = lines_FtF[0]
#        faults_names = line_fauts.split(' ')

#        faults_names[-1] = "\n".join(faults_names[-1].splitlines())

#        index_scenario = 0
#        scenarios_names = []
#        if np.size(self.rupture_set) == 0 :
#            scenarios_names = []
#        else :
#            for index_scenario in range(len(self.rupture_set)):
#                faults_in_scenario = self.rupture_set[index_scenario]
#                if len(faults_in_scenario) == 1.:
#                    print(faults_in_scenario)
#                    print('!!!!! ERROR !!!!!')
#                    print('A scenario have to be build with at least two faults')
#                    print('!!!!!!!!!!!!!!!!!')
#                else :
#                    scenario = {}
#                    faults_done = []
#                    for i in range(len(faults_in_scenario)):
#                        if not str(faults_in_scenario[i]).replace('\r','') in faults_done:
#                            scenario["f_%s" % str(i+1)] = [str(faults_in_scenario[i]).replace('\r','').replace('\t','').replace('\n','')]
#                            faults_done.append(str(faults_in_scenario[i]).replace('\r','').replace('\t','').replace('\n',''))
#                    if len(scenario)!=0:
#                        scenarios_names.append(scenario)
#                index_scenario += 1
#
#        ########################################################
#        #extractions of the geometries of the scenarios
#        ########################################################
#        geom_scenar = Geometry_scenario.Geom_scenar(faults_names,scenarios_names,self.File_geom,self.Model_name)
#        faults_lon = geom_scenar.faults_lon
#        faults_lat = geom_scenar.faults_lat
#
#
                
        '''########################################################
        # Random sampling of the fault slip-rate
        ########################################################'''
        
        ########################################################
        # find which faults interact a lot with each other
        ########################################################
        if self.sr_correl ==True :
            M_faults_correl = []
            for Fault_name in faults_names:
                M_faults_correl_i =  []
                for Fault_name_i in faults_names:
                    nb_join_ruptures = 0
                    if Fault_name == Fault_name_i :
                        nb_join_ruptures = 0
                    else:
                        for scenario_name in scenarios_names :
                            if Fault_name in str(scenario_name) :
                                if Fault_name_i in str(scenario_name) :
                                    nb_join_ruptures += 1
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
#                self.FaultProperties(Fault_name)
#                dip = self.dip
#                upper_sismo_depth = self.upper_sismo_depth
#                lower_sismo_depth = self.lower_sismo_depth
#                width = (lower_sismo_depth - upper_sismo_depth) / math.sin(math.radians(dip))
#                area = geom_scenar.length[index_fault] * 1000. * width * 1000.
#                faults_length.append(geom_scenar.length[index_fault] * 1000.)
#                faults_area.append(area)
#                faults_width.append(width)
                faults_length.append(faults_data[index_fault]['length'])
                faults_area.append(faults_data[index_fault]['area'])
                faults_width.append(faults_data[index_fault]['width'])
                
#                if self.rake> -135. and self.rake< -45:
#                    mecanism = 'N'
#                elif self.rake< 135. and self.rake> 45:
#                    mecanism = 'R'
#                else :
#                    mecanism = 'S'
#
#                faults_mecanism.append(mecanism)
                faults_mecanism.append(faults_data[index_fault]['mecanism'])

#                slip_rate_min = self.slip_rate_min
#                slip_rate_moy = self.slip_rate_moy
#                slip_rate_max = self.slip_rate_max
#                sr_values = [slip_rate_min,slip_rate_moy,slip_rate_max]
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
#                self.FaultProperties(Fault_name)
#                dip = self.dip
#                upper_sismo_depth = self.upper_sismo_depth
#                lower_sismo_depth = self.lower_sismo_depth
#                width = (lower_sismo_depth - upper_sismo_depth) / math.sin(math.radians(dip))
#                area = geom_scenar.length[index_fault] * 1000. * width * 1000.
#                faults_length.append(geom_scenar.length[index_fault] * 1000.)
#                faults_area.append(area)
#                faults_width.append(width)
#
#                if self.rake> -135. and self.rake< -45:
#                    mecanism = 'N'
#                elif self.rake< 135. and self.rake> 45:
#                    mecanism = 'R'
#                else :
#                    mecanism = 'S'
#
#
#                faults_mecanism.append(mecanism)
                    

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
            
        ratio_test = 0.5
        count_reruns = 1 #used to divide the sr increment if the fit is not good
        count_mfd90 = 1
        while abs(ratio_test-1) >self.fit_quality or math.isnan(ratio_test) == True:
            MFDs = EQ_on_faults.EQ_on_faults_from_sr(M_min,mfd_param,faults_names,faults_area,faults_length,faults_width,faults_slip_rates,scenarios_names,
                                                     faults_shear_mod,self.path,self.sample,self.selected_ScL,
                                                     self.dimention_used,self.use_all_ScL_data,faults_mecanism,self.bg_ratio,self.size_of_increment,self.mfd_hyp,count_reruns,
                                                     faults_lon,faults_lat,self.Mmax_range,self.calculation_log_file)
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
                
        # sclaling law as called by openquake
        if self.selected_ScL == 'Le2010' :
            ScL_oq = 'Leonard2014_SCR'
        if self.selected_ScL == 'WC1994' :
            ScL_oq = 'WC1994'
        else : #default Scaling relationship for opanquake
            ScL_oq = 'WC1994'
        
        ########################################################
        ########################################################
        #loop on simple faults  
        ########################################################
        ########################################################
        ID_number = 0
#        for Fault_name in enumerate(self.Fault_Names):
        for index_fault,fault_name in zip(range(len(faults_names)),faults_names):
            if fault_name in MFDs.faults_names :
                i_MFD = np.where(np.array(MFDs.faults_names) == fault_name)[0][0]
                MFD = MFDs.OQ_entry_faults[i_MFD]
                ID_number = ID_number + 1
                
#                self.FaultProperties(Fault_name[1])
                
#                if not self.Domain in str(self.Domain_in_the_model):
#                    self.Domain_in_the_model.append(self.Domain)
                if not faults_data[index_fault]['domain'] in str(self.Domain_in_the_model):
                    self.Domain_in_the_model.append(faults_data[index_fault]['domain'])
                    
                    
                ########################################################
                #Part concerning the geometry            
                ########################################################
                
#                index_zones = np.where(np.array(self.Column_Fault_name) == Fault_name[1])
#
#                ColLon = list(map(lambda i : self.Longitudes[i],index_zones[0]))
#                ColLat = list(map(lambda i : self.Latitudes[i],index_zones[0]))
#
#                ColLon = list(map(lambda i : self.Longitudes[i],index_zones[0]))
#                ColLat = list(map(lambda i : self.Latitudes[i],index_zones[0]))
#                Depth = list(map(lambda i : self.Depths[i],index_zones[0]))
                ColLon = faults_data[index_fault]['lon']
                ColLat = faults_data[index_fault]['lat']
                Depth = faults_data[index_fault]['depth']
                
                test_ok = 0
                if Depth and all(elem == 'sf' for elem in Depth):
                    type_of_fault = 'sf'
                else :
                    type_of_fault = 'cf'
                    Depth = [float(i) for i in Depth]
                    
                if type_of_fault == 'sf':
                    fault_name = self.Model_name + '_' + str(fault_name)
                    line='\t\t<simpleFaultSource id="'+ str(ID_number) +'" name="'+ str(fault_name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'
                    test_ok += 1
                    line+='\t\t\t<simpleFaultGeometry>\n'
                    line+='\t\t\t\t<gml:LineString>\n'
                    line+='\t\t\t\t\t<gml:posList>\n'
                    
                    #polygon = []
                    
                    # orienting the arrays in order to respect OQ right hand rule
                    compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
                    
                    if str('N') in str(faults_data[index_fault]['oriented']):
                        if compass_bearing < 180. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat) 
                            #reveresed = 'yes'
                    if str('S') in str(faults_data[index_fault]['oriented']):
                        if compass_bearing > 180. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat) 
                            #reveresed = 'yes'  
                    if str('E') in str(faults_data[index_fault]['oriented']):
                        if compass_bearing > 90. and compass_bearing < 270. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat)
                            #reveresed = 'yes'
                    if str('W') in str(faults_data[index_fault]['oriented']):
                        if compass_bearing < 90. or compass_bearing > 270. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat)
                            #reveresed = 'yes'
                        
                    for x,y in zip(ColLon,ColLat):
                        #polygon.append((x,y)) #ecriture du polygone de la zone
                        line+='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                    line+='\t\t\t\t\t</gml:posList>\n'
                    line+='\t\t\t\t</gml:LineString>\n'
                    line+='\t\t\t\t<dip>'+ str(faults_data[index_fault]['dip']) +'</dip>\n'
                    line+='\t\t\t\t<upperSeismoDepth>'+ str(faults_data[index_fault]['upper_sismo_depth']) +'</upperSeismoDepth>\n'
                    line+='\t\t\t\t<lowerSeismoDepth>'+ str(faults_data[index_fault]['lower_sismo_depth']) +'</lowerSeismoDepth>\n'
                    line+='\t\t\t</simpleFaultGeometry>\n'
                    XMLfile.write(line)
                    
                if type_of_fault == 'cf':
                    fault_name = self.Model_name + '_' + str(fault_name)
                    #line='\t\t<complexFaultSource id="'+ str(self.Model_name) + '_' + str(Fault_Name)+ '_' + str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    line='\t\t<complexFaultSource id="'+ str(ID_number) +'" name="'+ str(fault_name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'
                    test_ok += 1
                    line+='\t\t\t<complexFaultGeometry>\n'
                    
                    index_edge = 0
                    for depth_i in sorted(set(Depth)):
                        indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                        if index_edge == 0:
                            line+='\t\t\t\t<faultTopEdge>\n'
                        elif index_edge == len(set(Depth))-1:
                            line+='\t\t\t\t<faultBottomEdge>\n'
                        else :
                            line+='\t\t\t\t<intermediateEdge>\n'
                            
                        line+='\t\t\t\t<gml:LineString>\n'
                        line+='\t\t\t\t\t<gml:posList>\n'
                        for index in indexes_for_edge_i:
                            line+='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'
                            
                        line+='\t\t\t\t\t</gml:posList>\n'
                        line+='\t\t\t\t</gml:LineString>\n'
                        if index_edge == 0:
                            line+='\t\t\t\t</faultTopEdge>\n'
                        elif index_edge == len(set(Depth))-1:
                            line+='\t\t\t\t</faultBottomEdge>\n'
                        else :
                            line+='\t\t\t\t</intermediateEdge>\n'
                        index_edge+=1   
                        
                    line+='\t\t\t</complexFaultGeometry>\n'
                    XMLfile.write(line)
                    
                    
                if test_ok == 0:
                    print('!!!!!!!!!! Problem with the fault Geometry, please check input file''')
                    sys.exit()
                    
                    
                ########################################################
                #scaling law and aspect ratio
                ########################################################
                
                line='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
                XMLfile.write(line)
                line='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'
                XMLfile.write(line)
                
                ########################################################
                # seismicity of the fault            ♀
                ########################################################
                line = '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'
                XMLfile.write(line)
                #MFD_model += MFD
                if sum(MFD)!=0 :
                    log_mdf_file.write(str(fault_name)+'\t'+str(M_min)+'\t'+' '.join(list(map(str, MFD)))+'\n')
                    line = '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
                    XMLfile.write(line)
                else:
                    MFD[0] += 0.00000000001 #so it's not zero
                    log_mdf_file.write(str(fault_name)+'\t'+str(M_min)+'\t'+' '.join(list(map(str, MFD)))+'\n')
                    line = '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
                    XMLfile.write(line)
                    
                line = '\t\t\t</incrementalMFD>\n'
                XMLfile.write(line)
                
                line='\t\t\t<rake>'+str(faults_data[index_fault]['rake'])+'</rake>\n'
                XMLfile.write(line)
                if type_of_fault == 'sf':
                    line='\t\t</simpleFaultSource>\n'
                    XMLfile.write(line)
                if type_of_fault == 'cf':
                    line='\t\t</complexFaultSource>\n'
                    XMLfile.write(line)
        
        if len(self.rupture_set) != 0 :
            index_scenario = 0    
            
            for scenario in enumerate(MFDs.scenarios_names):
                index_scenario = np.where(np.array(MFDs.scenarios_names) == scenario[1])[0][0]
                MFD = MFDs.OQ_entry_scenarios[index_scenario]
                if sum(MFD)!=0:         
                    ID_number = ID_number + 1
                    index_fault = faults_names.index(scenario[1]['f_1'][0])
#                    self.FaultProperties(scenario[1]['f_1'])
                    scenar_name = '_'.join("{!s}={!r}".format(key,val) for (key,val) in scenario[1].items())
                    Fault_Name = self.Model_name + '_scenario_' + str(scenar_name)
                    line='\t\t<characteristicFaultSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(faults_data[index_fault]['domain']) + '">\n'
                    XMLfile.write(line)

                    
                    ########################################################
                    #Part concerning the geometry            
                    ########################################################
                    index_faults_in_scenario =  MFDs.index_faults_in_scenario[index_scenario][0]
                    faults_in_scenario = np.take(faults_names,index_faults_in_scenario)
                    
                    line='\t\t\t<surface>\n'
                    XMLfile.write(line)
                    
                    scenario_mechanism = []
                    for Fault_name in faults_in_scenario :
                        index_fault = faults_names.index(Fault_name)
                        
#                        self.FaultProperties(Fault_name)
                        
#                        index_zones = np.where(np.array(self.Column_Fault_name) == Fault_name)
                                       
#                        ColLon = list(map(lambda i : self.Longitudes[i],index_zones[0]))
#                        ColLat = list(map(lambda i : self.Latitudes[i],index_zones[0]))
#                        Depth = list(map(lambda i : self.Depths[i],index_zones[0]))
                        ColLon = faults_data[index_fault]['lon']
                        ColLat = faults_data[index_fault]['lat']
                        Depth = faults_data[index_fault]['depth']

                        scenario_mechanism.append(faults_data[index_fault]['rake'])
                        
                        if Depth and all(elem == 'sf' for elem in Depth):
                            type_of_fault = 'sf'
                        else :
                            type_of_fault = 'cf'
                            Depth = [float(i) for i in Depth]
                        if type_of_fault == 'sf':
                            line='\t\t\t<simpleFaultGeometry>\n'
                            XMLfile.write(line)
                            line='\t\t\t\t<gml:LineString>\n'
                            XMLfile.write(line)
                            line='\t\t\t\t\t<gml:posList>\n'
                            XMLfile.write(line)
                            
                            #polygon = []
                            
                            # orienting the arrays in order to respect OQ right hand rule
                            compass_bearing = calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
                            
                            if str('N') in str(faults_data[index_fault]['oriented']):
                                if compass_bearing < 180. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat)
                            if str('S') in str(faults_data[index_fault]['oriented']):
                                if compass_bearing > 180. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat)
                            if str('E') in str(faults_data[index_fault]['oriented']):
                                if compass_bearing > 90. and compass_bearing < 270. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat)
                            if str('W') in str(faults_data[index_fault]['oriented']):
                                if compass_bearing < 90. or compass_bearing > 270. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat)
                                
                            for x,y in zip(ColLon,ColLat):
                                #polygon.append((x,y)) #ecriture du polygone de la zone
                                line='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                                XMLfile.write(line)
                            line='\t\t\t\t\t</gml:posList>\n'
                            XMLfile.write(line)
                            line='\t\t\t\t</gml:LineString>\n'
                            XMLfile.write(line)
                            line='\t\t\t\t<dip>'+ str(faults_data[index_fault]['dip']) +'</dip>\n'
                            XMLfile.write(line)
                            line='\t\t\t\t<upperSeismoDepth>'+ str(faults_data[index_fault]['upper_sismo_depth']) +'</upperSeismoDepth>\n'
                            XMLfile.write(line)
                            line='\t\t\t\t<lowerSeismoDepth>'+ str(faults_data[index_fault]['lower_sismo_depth']) +'</lowerSeismoDepth>\n'
                            XMLfile.write(line)
                            line='\t\t\t</simpleFaultGeometry>\n'
                            XMLfile.write(line)
                            
                        if type_of_fault == 'cf':
                            line='\t\t\t<complexFaultGeometry>\n'
                            XMLfile.write(line)
                                    
                            index_edge = 0
                            for depth_i in sorted(set(Depth)):
                                indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                                if index_edge == 0:
                                    line='\t\t\t\t<faultTopEdge>\n'
                                    XMLfile.write(line)
                                elif index_edge == len(set(Depth))-1:
                                    line='\t\t\t\t<faultBottomEdge>\n'
                                    XMLfile.write(line)
                                else :
                                    line='\t\t\t\t<intermediateEdge>\n'
                                    XMLfile.write(line)
                                    
                                line='\t\t\t\t<gml:LineString>\n'
                                XMLfile.write(line)
                                line='\t\t\t\t\t<gml:posList>\n'
                                XMLfile.write(line)
                                for index in indexes_for_edge_i:
                                    line='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'
                                    XMLfile.write(line)
                                    
                                line='\t\t\t\t\t</gml:posList>\n'
                                XMLfile.write(line)
                                line='\t\t\t\t</gml:LineString>\n'
                                XMLfile.write(line)
                                if index_edge == 0:
                                    line='\t\t\t\t</faultTopEdge>\n'
                                    XMLfile.write(line)
                                elif index_edge == len(set(Depth))-1:
                                    line='\t\t\t\t</faultBottomEdge>\n'
                                    XMLfile.write(line)
                                else :
                                    line='\t\t\t\t</intermediateEdge>\n'
                                    XMLfile.write(line)
                                index_edge+=1   
                                
                            line='\t\t\t</complexFaultGeometry>\n'
                            XMLfile.write(line)
                            
                    
                    
                    line='\t\t\t</surface>\n'
                    XMLfile.write(line)
                    
                    ########################################################
                    #scaling law and aspect ratio
                    ########################################################
                    
                    line='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
                    XMLfile.write(line)
                    line='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'
                    XMLfile.write(line)
                    
                    ########################################################
                    # seismicity of the scenario           
                    ########################################################
        #            
                    log_mdf_file.write(str(Fault_Name) + '\t' + str(M_min) + '\t' + ' '.join(list(map(str, MFD)))+'\n')
                    #while MFD[i_Mmin] == 0:
                    #    i_Mmin += 1
                    line = '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'# + 0.1 * i_Mmin
                    XMLfile.write(line)
                    #MFD_model += MFD
                    line = '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'#[-(len(MFD) - i_Mmin):]
                    XMLfile.write(line)
                    line = '\t\t\t</incrementalMFD>\n'
                    XMLfile.write(line)
                    
                    #find the dominant kinematic of the scenario
                    rake= np.mean(scenario_mechanism)
                    line='\t\t\t<rake>'+str(rake)+'</rake>\n'
                    XMLfile.write(line)
                    line='\t\t</characteristicFaultSource>\n'
                    XMLfile.write(line)
                index_scenario += 1
                
                
        '''#########################
        # Defining the background
        #########################'''
        Lon_bg, Lat_bg  = bg.geom(self.Model_name,self.File_bg )
        upperSeismoDepth, lowerSeismoDepth, ruptAspectRatio, nodalPlanes, hypoDepths = bg.prop(self.Model_name,self.file_prop_bg)
        
        MFD = MFDs.EQ_rate_BG   
        if sum(MFD) != 0. :
            line='\t\t<areaSource id="'+ str(ID_number + 1 ) +'" name="Background" tectonicRegion="' + str(self.Domain) + '">\n'
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
            XMLfile.write(line)
        
        '''#############################
        ### defining the other sources based on the host model
        ##############################'''
        if self.use_host_model == True :
            host_model.build(XMLfile,self.host_model_file,Lon_bg,Lat_bg)

        #end of the file
        line='\t</sourceModel>\n'
        XMLfile.write(line)
        line='</nrml>\n'
        XMLfile.write(line)
        XMLfile.close()
        log_sr_file.close()
        log_mdf_file.close()
                
#    def FaultGeometry(self):
#        NomFichier_InfosZonage = self.File_geom
#        InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8'),('U100')],skip_header = 1)
#        Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
#        index_model = np.where(np.array(Column_model_name) == self.Model_name)
#        self.Column_Fault_name = list(map(lambda i : InfosZonage[i][1],index_model[0]))
#        self.Longitudes = list(map(lambda i : InfosZonage[i][2],index_model[0]))
#        self.Latitudes = list(map(lambda i : InfosZonage[i][3],index_model[0]))
#        self.Depths = list(map(lambda i : InfosZonage[i][4],index_model[0]))
#
#        ZoneSelec = self.Column_Fault_name
#        DicoZone = dict([(k,ZoneSelec.count(k)) for k in set(ZoneSelec)])
#        Longitudes = []
#        Latitudes = []
#        Depths = []
#        Column_Fault_name = []
#        for cle in DicoZone.keys():
#            indices_ZonesSelec = np.where(np.array(self.Column_Fault_name) == cle)
#            ColonneNomZone_inter = np.take(self.Column_Fault_name,indices_ZonesSelec)
#            Longitudes_inter = np.take(self.Longitudes,indices_ZonesSelec)
#            Latitudes_inter = np.take(self.Latitudes,indices_ZonesSelec)
#            depth_inter = np.take(self.Depths,indices_ZonesSelec)
#
#            Longitudes_inter = Longitudes_inter[0].tolist()
#            Latitudes_inter = Latitudes_inter[0].tolist()
#            depth_inter = depth_inter[0].tolist()
#            ColonneNomZone_inter = ColonneNomZone_inter[0].tolist()
#            compt = 0
#            for xx,yy,nn,dd in zip(Longitudes_inter,Latitudes_inter,ColonneNomZone_inter,depth_inter):
#                compt+=1
#                Longitudes.append(xx)
#                Latitudes.append(yy)
#                Depths.append(dd)
#                Column_Fault_name.append(nn)
#
#        self.Longitudes =Longitudes
#        self.Latitudes =Latitudes
#        self.Depths =Depths
#        self.Column_Fault_name = Column_Fault_name
#        self.Nb_data_per_zone = dict([(k,self.Column_Fault_name.count(k)) for k in set(self.Column_Fault_name)])
#        self.Fault_Names = sorted(self.Nb_data_per_zone.keys())
#
#    def FaultProperties(self,Name_of_fault):
#        FileName_Prop = self.File_prop
#        Prop = np.genfromtxt(FileName_Prop,
#                                   dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
#                                          ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
#        Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
#        Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
#        index_model = np.where(np.array(Column_model_name) == self.Model_name)[0]
#
#        Prop = np.take(Prop,index_model)
#        index_fault = np.where(np.array(Column_fault_name[index_model[0]:index_model[-1]+1]) == Name_of_fault)
#        Indexfault_final = index_fault[0]
#
#        self.dip = Prop[Indexfault_final][0][2]
#        self.oriented = Prop[Indexfault_final][0][3]
#        self.rake = Prop[Indexfault_final][0][4]
#        self.upper_sismo_depth = Prop[Indexfault_final][0][5]
#        self.lower_sismo_depth = Prop[Indexfault_final][0][6]
#
#        self.slip_rate_min = Prop[Indexfault_final][0][7]
#        self.slip_rate_moy = Prop[Indexfault_final][0][8]
#        self.slip_rate_max = Prop[Indexfault_final][0][9]
#        self.Domain = Prop[Indexfault_final][0][10]
#        self.shear_mod = Prop[Indexfault_final][0][11]
#
#
#        if self.rake == 'N' :
#            self.rake = -90.00
#        if self.rake == 'S' :
#            self.rake = 00.00
#        if self.rake == 'SS' :
#            self.rake = 00.00
#        if self.rake == 'R' :
#            self.rake = 90.00
#        self.rake = float(self.rake)
#
#        if len(str(self.dip)) == 0:
#            print('\nError!!! please verify your input file for fault parameters\n')


