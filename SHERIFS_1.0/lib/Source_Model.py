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
path_dossier = path_actuel + '/lib'
sys.path.append(path_dossier)
#import EQ_on_faults
import EQ_on_faults
import Geometry_scenario
import math
from math import pi, cos, radians , sin, asin, sqrt
import matplotlib.path as mplPath

import matplotlib.pyplot as plt

def reproject(latitude, longitude):
    """Returns the x & y coordinates in meters using a sinusoidal projection"""
    earth_radius = 6371009 # in meters
    lat_dist = pi * earth_radius / 180.0
    y = [lat * lat_dist for lat in latitude]
    x = [long * lat_dist * cos(radians(lat)) 
                for lat, long in zip(latitude, longitude)]
    return x, y

def points_aligned(a, b, c):
    crossproduct = (c[1] - a[1]) * (b[0] -a[0]) - (c[0] -a[0]) * (b[1] - a[1])
    epsilon = 10000000.
    if abs(crossproduct) > epsilon:
        return False
    dotproduct = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1])*(b[1] - a[1])
    if dotproduct < 0:
        return False
    
    squaredlengthba = (b[0]-a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1]-a[1])
    if dotproduct > squaredlengthba:
        return False
    
    return True
    
    
def area_of_polygon(x, y):
    """Calculates the area of an arbitrary polygon given its verticies"""
    # first, does the list of vertices
    x_vertices = []
    y_vertices = []
    inn = []
    for i in range(len(x)):
        if i == 0:
            if points_aligned([x[-1],y[-1]],[x[0],y[0]], [x[1],y[1]]) == False:
                 x_vertices.append(x[i])
                 y_vertices.append(y[i])
                 inn.append(1)
            else:
                inn.append(0)
        elif i == len(x)-1:
            if points_aligned([x[-2],y[-2]], [x[-1],y[-1]],[x[0],y[0]]) == False:
                 x_vertices.append(x[i])
                 y_vertices.append(y[i])
                 inn.append(1)
            else:
                 inn.append(0)
        else:
             if points_aligned([x[i-1],y[i-1]], [x[i],y[i]],[x[i+1],y[i+1]]) == False:
                 x_vertices.append(x[i])
                 y_vertices.append(y[i])
                 inn.append(1)
             else:
                 inn.append(0)
    #print len(x),len(x_vertices),inn
    
    area = 0.0
    for i in range(-1, len(x_vertices)-1):
        area += x_vertices[i] * (y_vertices[i+1] - y_vertices[i-1])
    return abs(area) / 2.0
#
#def PolyArea(x, y):
#    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1]) 

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       #raise Exception('lines do not intersect')
       x = 'no_intesection'
       y = 'no_intesection'
    else: 
        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
    return x, y

class Source_Model_Creator:
    def __init__(self,path,Run_Name,Model_name,File_geom,File_prop,File_bg,file_prop_bg,File_faults_n_scenarios,
                 Domain_in_model,sample,seed,Mmin,selected_ScL,dimention_used,use_all_ScL_data,
                 b_min,b_max,mfd_hyp,bg_ratio,sr_correl,size_of_increment,fit_quality,Mmax_range,calculation_log_file,use_host_model,host_model_file):
        self.Run_Name = Run_Name
#        self.Latitude = Latitude
#        self.Longitude = Longitude
        self.Domain_in_the_model = Domain_in_model #list of the different domain included in the model
        #a envoyer dans job
        self.File_geom = File_geom
        self.File_prop = File_prop
        self.File_bg = File_bg
        self.file_prop_bg = file_prop_bg
        self.File_faults_n_scenarios = File_faults_n_scenarios
        
        #a envoyer dans logic tree
        self.Model_name = Model_name 
        #self.M_trunc = M_trunc
        self.b_min = float(b_min)
        self.b_max = float(b_max)
        self.mfd_hyp = mfd_hyp
        #self.a_s = a_s
        
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

        self.initialize()

    def initialize(self):

        self.FaultGeometry()  #extract the geometries from the geometry file
        
        XMLfile=open(self.path +'/Source_model_' + str(self.sample) + '.xml','w')
        
        if not os.path.exists(self.path +'/Log'):
            os.makedirs(self.path +'/Log')
        log_sr_file=open(self.path +'/Log/slip_rate_sample_' + str(self.sample) + '.txt','w')
        log_mdf_file=open(self.path +'/Log/mdf_sample_' + str(self.sample) + '.txt','w')
        
        
        
        Line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
        XMLfile.write(Line)
        Line='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
        XMLfile.write(Line)
        Line='\txmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        XMLfile.write(Line)
        Line='\t<sourceModel name="Hazard Model">\n'
        XMLfile.write(Line)
        Line='\n'
        XMLfile.write(Line)

        #initialisation of the general parameters (M_min, shear modulus and b value)
        log_general_parameters_file = open(self.path +'/Log/general_parameters_sample_' + str(self.sample) + '.txt','w')
        
        M_min = float(self.Mmin)
        log_general_parameters_file.write('M_tronc\t'+str(M_min)+'\n')

#        shear_modulus = float(self.shear_mod)*10**9 #convertion from GPa in Newton.m²
#        log_general_parameters_file.write('shear_modulus\t'+str(self.shear_mod)+'\n')
        
        if self.sample == 1 :
            b_value = (self.b_min + self.b_max)/2.
        else :
            b_value = np.random.triangular(self.b_min,(self.b_min + self.b_max)/2.,self.b_max)
        log_general_parameters_file.write('b_value\t'+str(b_value)+'\n')
        
        log_general_parameters_file.close()
        
        #initialization of array for faults information
        faults_area = []
        faults_length = []
        faults_width = []
        faults_slip_rates = []
        faults_mecanism = []
        faults_shear_mod = []
                
        ########################################################
        #Extraction of the faults and scenarios present in the model from the text file
        ########################################################
        
        lines_FtF = [line.rstrip('\n') for line in open(self.File_faults_n_scenarios)]
        line_fauts = lines_FtF[0] 
        faults_names = line_fauts.split(' ')
            
        faults_names[-1] = "\n".join(faults_names[-1].splitlines())

        index_scenario = 0
        scenarios_names = []  
        if np.size(lines_FtF) == 1 :
            scenarios_names = [] 
        else :
            for index_scenario in range(len(lines_FtF)-1):
                line_scenario = str(lines_FtF[index_scenario+1])
                faults_in_scenario = list(line_scenario.split(' '))
                if len(faults_in_scenario) == 1.:
                    print(faults_in_scenario)
                    print('!!!!! ERROR !!!!!')
                    print('One scenario have to be build with at least two faults')
                    print('!!!!!!!!!!!!!!!!!')
                else :
                    scenario = {}
                    faults_done = []
                    for i in range(len(faults_in_scenario)):
                        if not str(faults_in_scenario[i]).replace('\r','') in faults_done:
                            scenario["f_%s" % str(i+1)] = [str(faults_in_scenario[i]).replace('\r','').replace('\n','').replace('\t','')]
                            faults_done.append(str(faults_in_scenario[i]).replace('\r','').replace('\n','').replace('\t',''))
                    if len(scenario)!=0:
                        scenarios_names.append(scenario)
                index_scenario += 1
                
        ########################################################
        #extractions of the geometries of the scenarios and of the MFDs (Magnitude Frequency Distribution)
        ########################################################
        geom_scenar = Geometry_scenario.Geom_scenar(faults_names,scenarios_names,self.File_geom,self.Model_name)
        faults_lon = geom_scenar.faults_lon
        faults_lat = geom_scenar.faults_lat
        
        
                
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
            
            slip_rates_min = []
            slip_rates_moy = []
            slip_rates_max = []
            list_quater_picked = list(np.zeros(len(faults_names)))
            # pick the faults slip rates
            index_fault = 0
            for Fault_name in faults_names:
                self.FaultProperties(Fault_name)
                dip = self.dip
                upper_sismo_depth = self.upper_sismo_depth
                lower_sismo_depth = self.lower_sismo_depth
                width = (lower_sismo_depth - upper_sismo_depth) / math.sin(math.radians(dip))
                area = geom_scenar.length[index_fault] * 1000. * width * 1000.
                faults_length.append(geom_scenar.length[index_fault] * 1000.)
                faults_area.append(area)
                faults_width.append(width)
                
                if self.rake> -135. and self.rake< -45:
                    mecanism = 'N'
                elif self.rake< 135. and self.rake> 45: 
                    mecanism = 'R'
                else :
                    mecanism = 'S'
                    
                    
                faults_mecanism.append(mecanism)
                    
                ########################################################
                #ramdom sampling the slip rate  . uniform sampling
                ########################################################            
                slip_rate_min = self.slip_rate_min
                slip_rate_moy = self.slip_rate_moy
                slip_rate_max = self.slip_rate_max
                
                slip_rates_min.append(slip_rate_min)
                slip_rates_moy.append(slip_rates_moy)
                slip_rates_max.append(slip_rates_max)
                
                
                if self.sample == 1 :
                    slip_rate = self.slip_rate_moy
                    #print Fault_name,self.slip_rate_moy
#                elif self.sample == 2 :
#                    slip_rate = self.slip_rate_min
#                    #print Fault_name,self.slip_rate_moy
#                elif self.sample == 3 :
#                    slip_rate = self.slip_rate_max
#                    #print Fault_name,self.slip_rate_moy
                else :
#                    #find the fault most correleted with
#                    M_faults_correl_i = M_faults_correl[index_fault]
#                    index_max_corr = np.argmax(M_faults_correl_i)
#                    while index_max_corr > len(faults_slip_rates) and np.sum(M_faults_correl_i) != 0 : #the carrelate fault has been done already
#                        M_faults_correl_i[index_max_corr] = 0
#                        index_max_corr = np.argmax(M_faults_correl_i)
#                    if np.sum(M_faults_correl_i) != 0 :
#                        if M_faults_correl_i[index_max_corr] >= 2 and M_faults_correl_i[index_max_corr] > 2./3. * max_correl :
#                            # in a first order correlation bewteen the faults
#                            if faults_slip_rates[index_max_corr] >= slip_rates_moy[index_max_corr] :
#                                slip_rate = np.random.uniform(slip_rate_moy - 0.00001,slip_rate_max)
#                            else :
#                                slip_rate = np.random.uniform(slip_rate_min,slip_rate_moy+ 0.00001)
#                        else : # it's a second order correlation between the faults
#                            if faults_slip_rates[index_max_corr] >= slip_rates_moy[index_max_corr] :
#                                slip_rate = np.random.uniform(slip_rate_min + 2./3. * (slip_rate_moy - slip_rate_min)
#                                ,slip_rate_max)
#                            else :
#                                slip_rate = np.random.uniform(slip_rate_min,slip_rate_max - 2./3. * (slip_rate_max - slip_rate_moy))
#                    else :  # no interaction yet       
#                         slip_rate = np.random.uniform(slip_rate_min,slip_rate_max)

                    if sum(M_linked_lvl[index_fault])==0 : #the fault is alone
                        slip_rate = np.random.uniform(slip_rate_min,slip_rate_max)
                    else: #the fault is connected
                        value_lvl = 10
                        quarters_picked = []
                        for index_c in range(len(M_linked_lvl[index_fault])):
                            if M_linked_lvl[index_fault][index_c] != 0:
                                if list_quater_picked[index_c]!=0:
                                    if M_linked_lvl[index_fault][index_c] < value_lvl :
                                        quarters_picked = []
                                        value_lvl = M_linked_lvl[index_fault][index_c]
                                        #print faults_names[index_fault],faults_names[index_c]
                                    if M_linked_lvl[index_fault][index_c] == value_lvl :
                                            quarters_picked.append(list_quater_picked[index_c])
                        if quarters_picked == []: #none of the faults have bin picked yet
                            slip_rate = np.random.uniform(slip_rate_min,slip_rate_max)
                            if slip_rate < (slip_rate_min + 1./2. * (slip_rate_moy-slip_rate_min)):
                                quarter_to_pick = 1
                            elif slip_rate < (slip_rate_moy):
                                quarter_to_pick = 2
                            elif slip_rate < (slip_rate_moy + 1./2. * (slip_rate_max-slip_rate_moy)):
                                quarter_to_pick = 3
                            else:
                                quarter_to_pick = 4
                                
                        else :
                            quarter_to_pick = max(set(quarters_picked),key=quarters_picked.count)
                            if quarter_to_pick == 1 :
                                slip_rate = np.random.uniform(slip_rate_min,slip_rate_min + 1./2. * (slip_rate_moy-slip_rate_min) +0.000001)
                            elif quarter_to_pick == 2 :
                                slip_rate = np.random.uniform(slip_rate_min + 1./2. * (slip_rate_moy-slip_rate_min),slip_rate_moy+0.000001)
                            elif quarter_to_pick == 3 :
                                slip_rate = np.random.uniform(slip_rate_moy,slip_rate_moy + 1./2. * (slip_rate_max-slip_rate_moy)+0.000001)
                            elif quarter_to_pick == 4 :
                                slip_rate = np.random.uniform(slip_rate_moy + 1./2. * (slip_rate_max-slip_rate_moy),slip_rate_max+0.000001)
                            
                                
                        list_quater_picked[index_fault] = quarter_to_pick
                                
                        #print(faults_names[index_fault],quarter_to_pick)
                
            
                log_line = str(Fault_name) + '\t' + str(slip_rate) + '\n' #writting in the log file
                log_sr_file.write(log_line)                        
                
                faults_slip_rates.append(slip_rate * 0.001) # transform from mm to m
                
                faults_shear_mod.append(float(self.shear_mod)*10**9 )
                index_fault += 1  
        
        else: # the option of correlation bewteen corroption fault not used
            index_fault = 0
            for Fault_name in faults_names:
                self.FaultProperties(Fault_name)
                dip = self.dip
                upper_sismo_depth = self.upper_sismo_depth
                lower_sismo_depth = self.lower_sismo_depth
                width = (lower_sismo_depth - upper_sismo_depth) / math.sin(math.radians(dip))
                area = geom_scenar.length[index_fault] * 1000. * width * 1000.
                faults_length.append(geom_scenar.length[index_fault] * 1000.)
                faults_area.append(area)
                faults_width.append(width)
                
                if self.rake> -135. and self.rake< -45:
                    mecanism = 'N'
                elif self.rake< 135. and self.rake> 45: 
                    mecanism = 'R'
                else :
                    mecanism = 'S'
                    
                    
                faults_mecanism.append(mecanism)
                    
                ########################################################
                #ramdom sampling the slip rate  . uniform sampling
                ########################################################            
                slip_rate_min = self.slip_rate_min
                slip_rate_moy = self.slip_rate_moy
                slip_rate_max = self.slip_rate_max
                
                if self.sample == 1 :
                    slip_rate = self.slip_rate_moy
                else :
                    slip_rate = np.random.uniform(slip_rate_min,slip_rate_max)
                
                log_line = str(Fault_name) + '\t' + str(slip_rate) + '\n' #writting in the log file
                log_sr_file.write(log_line)                        
                
                faults_slip_rates.append(slip_rate * 0.001) # transform from mm to m
                
                faults_shear_mod.append(float(self.shear_mod)*10**9 )
                
                index_fault += 1            
            
        ratio_test = 0.5
        count_reruns = 1 #used to divide the sr increment if the fit is not good
        count_mfd90 = 1
        while abs(ratio_test-1) >self.fit_quality or math.isnan(ratio_test) == True:
            MFDs = EQ_on_faults.EQ_on_faults_from_sr(M_min,b_value,faults_names,faults_area,faults_length,faults_width,faults_slip_rates,scenarios_names,
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
        for Fault_name in enumerate(self.Fault_Names):
            if Fault_name[1] in MFDs.faults_names :
                MFD = MFDs.OQ_entry_faults[np.where(np.array(MFDs.faults_names) == Fault_name[1])][0]
                ID_number = ID_number + 1
                
                self.FaultProperties(Fault_name[1])
                
                if not self.Domain in str(self.Domain_in_the_model):
                    self.Domain_in_the_model.append(self.Domain)
                ########################################################
                #Part concerning the geometry            
                ########################################################
                
                index_zones = np.where(np.array(self.Column_Fault_name) == Fault_name[1])
                               
                ColLon = list(map(lambda i : self.Longitudes[i],index_zones[0]))
                ColLat = list(map(lambda i : self.Latitudes[i],index_zones[0]))
                          
                ColLon = list(map(lambda i : self.Longitudes[i],index_zones[0]))
                ColLat = list(map(lambda i : self.Latitudes[i],index_zones[0]))
                Depth = list(map(lambda i : self.Depths[i],index_zones[0]))
                test_ok = 0
                if Depth and all(elem == 'sf' for elem in Depth):
                    type_of_fault = 'sf'
                else :
                    type_of_fault = 'cf'
                    Depth = [float(i) for i in Depth]
                    
                if type_of_fault == 'sf':
                    Fault_Name = self.Model_name + '_' + str(Fault_name[1])
                    #Line='\t\t<simpleFaultSource id="'+ str(self.Model_name) + '_' + str(Fault_Name)+ '_' + str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    Line='\t\t<simpleFaultSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    XMLfile.write(Line)
                    test_ok += 1
                    Line='\t\t\t<simpleFaultGeometry>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t\t<gml:LineString>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t\t\t<gml:posList>\n'
                    XMLfile.write(Line)
                    
                    #polygon = []
                    
                    # orienting the arrays in order to respect OQ right hand rule
                    compass_bearing = self.calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
                    
                    if str('N') in str(self.oriented):
                        if compass_bearing < 180. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat) 
                            #reveresed = 'yes'
                    if str('S') in str(self.oriented):
                        if compass_bearing > 180. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat) 
                            #reveresed = 'yes'  
                    if str('E') in str(self.oriented):
                        if compass_bearing > 90. and compass_bearing < 270. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat)
                            #reveresed = 'yes'
                    if str('W') in str(self.oriented):
                        if compass_bearing < 90. or compass_bearing > 270. :
                            ColLon = reversed(ColLon)
                            ColLat = reversed(ColLat)
                            #reveresed = 'yes'
                        
                    for x,y in zip(ColLon,ColLat):
                        #polygon.append((x,y)) #ecriture du polygone de la zone
                        Line='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                        XMLfile.write(Line)
                    Line='\t\t\t\t\t</gml:posList>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t\t</gml:LineString>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t\t<dip>'+ str(self.dip) +'</dip>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t\t<upperSeismoDepth>'+ str(self.upper_sismo_depth) +'</upperSeismoDepth>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t\t<lowerSeismoDepth>'+ str(self.lower_sismo_depth) +'</lowerSeismoDepth>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t</simpleFaultGeometry>\n'
                    XMLfile.write(Line)
                    
                if type_of_fault == 'cf':
                    Fault_Name = self.Model_name + '_' + str(Fault_name[1])
                    #Line='\t\t<complexFaultSource id="'+ str(self.Model_name) + '_' + str(Fault_Name)+ '_' + str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    Line='\t\t<complexFaultSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    XMLfile.write(Line)
                    test_ok += 1
                    Line='\t\t\t<complexFaultGeometry>\n'
                    XMLfile.write(Line)
#                    
#                    depth_i = Depth[0]
#                    indexes_for_edges= []
#                    indexes_for_edge_i = []
#                    for index in range(len(Depth)):
#                        if Depth[index] == depth_i:
#                            print depth_i
#                            indexes_for_edge_i.append(index)
#                            depth_i = Depth[index]
#                        else : 
#                            print depth_i
#                            indexes_for_edges.append(indexes_for_edge_i)
#                            depth_i = Depth[index]
#                    indexes_for_edges.append(indexes_for_edge_i)
#                    print indexes_for_edges
                    
                    index_edge = 0
                    for depth_i in sorted(set(Depth)):
                        indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                        if index_edge == 0:
                            Line='\t\t\t\t<faultTopEdge>\n'
                            XMLfile.write(Line)
                        elif index_edge == len(set(Depth))-1:
                            Line='\t\t\t\t<faultBottomEdge>\n'
                            XMLfile.write(Line)
                        else :
                            Line='\t\t\t\t<intermediateEdge>\n'
                            XMLfile.write(Line)
                            
                        Line='\t\t\t\t<gml:LineString>\n'
                        XMLfile.write(Line)
                        Line='\t\t\t\t\t<gml:posList>\n'
                        XMLfile.write(Line)
                        for index in indexes_for_edge_i:
                            Line='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'
                            XMLfile.write(Line)
                            
                        Line='\t\t\t\t\t</gml:posList>\n'
                        XMLfile.write(Line)
                        Line='\t\t\t\t</gml:LineString>\n'
                        XMLfile.write(Line)
                        if index_edge == 0:
                            Line='\t\t\t\t</faultTopEdge>\n'
                            XMLfile.write(Line)
                        elif index_edge == len(set(Depth))-1:
                            Line='\t\t\t\t</faultBottomEdge>\n'
                            XMLfile.write(Line)
                        else :
                            Line='\t\t\t\t</intermediateEdge>\n'
                            XMLfile.write(Line)
                        index_edge+=1   
                        
                    Line='\t\t\t</complexFaultGeometry>\n'
                    XMLfile.write(Line)
                    
                    
                if test_ok == 0:
                    print('!!!!!!!!!! Problem with the fault Geometry, please check input file''')
                    
                    
                ########################################################
                #scaling law and aspect ratio
                ########################################################
                
                Line='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
                XMLfile.write(Line)
                Line='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'
                XMLfile.write(Line)
                
                ########################################################
                # seismicity of the fault            ♀
                ########################################################
                #Line='\t\t\t<truncGutenbergRichterMFD aValue="'+ str(self.a_value) +'" bValue="'+ str(self.b_value) +'" minMag="6.0" maxMag="'+ str(self.Mmax_moy) +'" />\n'
                Line = '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'
                XMLfile.write(Line)
                #MFD_model += MFD
                if sum(MFD)!=0 :
                    log_mdf_file.write(str(Fault_Name)+'\t'+str(M_min)+'\t'+' '.join(list(map(str, MFD)))+'\n')
                    Line = '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
                    XMLfile.write(Line)
                else:
                    MFD[0] += 0.00000000001 #so it's not zero
                    log_mdf_file.write(str(Fault_Name)+'\t'+str(M_min)+'\t'+' '.join(list(map(str, MFD)))+'\n')
                    Line = '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
                    XMLfile.write(Line)
                    
                Line = '\t\t\t</incrementalMFD>\n'
                XMLfile.write(Line)
                
                Line='\t\t\t<rake>'+str(self.rake)+'</rake>\n'  
                XMLfile.write(Line)
                if type_of_fault == 'sf':
                    Line='\t\t</simpleFaultSource>\n'
                    XMLfile.write(Line)
                if type_of_fault == 'cf':
                    Line='\t\t</complexFaultSource>\n'
                    XMLfile.write(Line)
        
        if np.size(lines_FtF) != 1 :
                    
            index_scenario = 0    
            
            for scenario in enumerate(MFDs.scenarios_names):
                MFD = MFDs.OQ_entry_scenarios[np.where(np.array(MFDs.scenarios_names) == scenario[1])][0]
                if sum(MFD)!=0:         
                    ID_number = ID_number + 1
                    self.FaultProperties(scenario[1]['f_1'])
                    scenar_name = '_'.join("{!s}={!r}".format(key,val) for (key,val) in scenario[1].items())
                    Fault_Name = self.Model_name + '_scenario_' + str(scenar_name)
                    #Line='\t\t<characteristicFaultSource id="'+ str(self.Model_name) + '_' + str(Fault_Name)+ '_' + str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    Line='\t\t<characteristicFaultSource id="'+ str(ID_number) +'" name="'+ str(Fault_Name) +'" tectonicRegion="' + str(self.Domain) + '">\n'
                    XMLfile.write(Line)

                    
                    ########################################################
                    #Part concerning the geometry            
                    ########################################################
                    index_faults_in_scenario =  MFDs.index_faults_in_scenario[index_scenario][0]
                    faults_in_scenario = np.take(faults_names,index_faults_in_scenario)
                    
                    Line='\t\t\t<surface>\n'
                    XMLfile.write(Line)
                    
                    scenario_mechanism = []
                    for Fault_name in faults_in_scenario :
                        
                        self.FaultProperties(Fault_name)
                        
                        index_zones = np.where(np.array(self.Column_Fault_name) == Fault_name)
                                       
                        ColLon = list(map(lambda i : self.Longitudes[i],index_zones[0]))
                        ColLat = list(map(lambda i : self.Latitudes[i],index_zones[0]))
                        Depth = list(map(lambda i : self.Depths[i],index_zones[0]))
                        scenario_mechanism.append(self.rake)
                        
                        if Depth and all(elem == 'sf' for elem in Depth):
                            type_of_fault = 'sf'
                        else :
                            type_of_fault = 'cf'
                            Depth = [float(i) for i in Depth]
                        if type_of_fault == 'sf':
                            Line='\t\t\t<simpleFaultGeometry>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t\t<gml:LineString>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t\t\t<gml:posList>\n'
                            XMLfile.write(Line)
                            
                            #polygon = []
                            
                            # orienting the arrays in order to respect OQ right hand rule
                            compass_bearing = self.calculate_initial_compass_bearing((ColLat[0],ColLon[0]),(ColLat[-1],ColLon[-1]))
                            
                            if str('N') in str(self.oriented):
                                if compass_bearing < 180. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat) 
                                    #reveresed = 'yes'
                            if str('S') in str(self.oriented):
                                if compass_bearing > 180. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat) 
                                    #reveresed = 'yes'  
                            if str('E') in str(self.oriented):
                                if compass_bearing > 90. and compass_bearing < 270. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat)
                                    #reveresed = 'yes'
                            if str('W') in str(self.oriented):
                                if compass_bearing < 90. or compass_bearing > 270. :
                                    ColLon = reversed(ColLon)
                                    ColLat = reversed(ColLat)
                                    #reveresed = 'yes'
                                
                            for x,y in zip(ColLon,ColLat):
                                #polygon.append((x,y)) #ecriture du polygone de la zone
                                Line='\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                                XMLfile.write(Line)
                            Line='\t\t\t\t\t</gml:posList>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t\t</gml:LineString>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t\t<dip>'+ str(self.dip) +'</dip>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t\t<upperSeismoDepth>'+ str(self.upper_sismo_depth) +'</upperSeismoDepth>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t\t<lowerSeismoDepth>'+ str(self.lower_sismo_depth) +'</lowerSeismoDepth>\n'
                            XMLfile.write(Line)
                            Line='\t\t\t</simpleFaultGeometry>\n'
                            XMLfile.write(Line)
                            
                        if type_of_fault == 'cf':
                            Line='\t\t\t<complexFaultGeometry>\n'
                            XMLfile.write(Line)
                                    
        #                    depth_i = Depth[0]
        #                    indexes_for_edges= []
        #                    indexes_for_edge_i = []
        #                    for index in range(len(Depth)):
        #                        if Depth[index] == depth_i:
        #                            print depth_i
        #                            indexes_for_edge_i.append(index)
        #                            depth_i = Depth[index]
        #                        else : 
        #                            print depth_i
        #                            indexes_for_edges.append(indexes_for_edge_i)
        #                            depth_i = Depth[index]
        #                    indexes_for_edges.append(indexes_for_edge_i)
        #                    print indexes_for_edges
                            
                            index_edge = 0
                            for depth_i in sorted(set(Depth)):
                                indexes_for_edge_i = np.where(np.array(Depth)==depth_i)[0]
                                if index_edge == 0:
                                    Line='\t\t\t\t<faultTopEdge>\n'
                                    XMLfile.write(Line)
                                elif index_edge == len(set(Depth))-1:
                                    Line='\t\t\t\t<faultBottomEdge>\n'
                                    XMLfile.write(Line)
                                else :
                                    Line='\t\t\t\t<intermediateEdge>\n'
                                    XMLfile.write(Line)
                                    
                                Line='\t\t\t\t<gml:LineString>\n'
                                XMLfile.write(Line)
                                Line='\t\t\t\t\t<gml:posList>\n'
                                XMLfile.write(Line)
                                for index in indexes_for_edge_i:
                                    Line='\t\t\t\t\t\t' + str(ColLon[index]) + ' ' + str(ColLat[index])  + ' ' + str(Depth[index]) + '\n'
                                    XMLfile.write(Line)
                                    
                                Line='\t\t\t\t\t</gml:posList>\n'
                                XMLfile.write(Line)
                                Line='\t\t\t\t</gml:LineString>\n'
                                XMLfile.write(Line)
                                if index_edge == 0:
                                    Line='\t\t\t\t</faultTopEdge>\n'
                                    XMLfile.write(Line)
                                elif index_edge == len(set(Depth))-1:
                                    Line='\t\t\t\t</faultBottomEdge>\n'
                                    XMLfile.write(Line)
                                else :
                                    Line='\t\t\t\t</intermediateEdge>\n'
                                    XMLfile.write(Line)
                                index_edge+=1   
                                
                            Line='\t\t\t</complexFaultGeometry>\n'
                            XMLfile.write(Line)
                            
                    
                    
                    Line='\t\t\t</surface>\n'
                    XMLfile.write(Line)
                    
                    ########################################################
                    #scaling law and aspect ratio
                    ########################################################
                    
                    Line='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
                    XMLfile.write(Line)
                    Line='\t\t\t<ruptAspectRatio>1.0</ruptAspectRatio>\n'
                    XMLfile.write(Line)
                    
                    ########################################################
                    # seismicity of the scenario           
                    ########################################################
        #            
                    log_mdf_file.write(str(Fault_Name) + '\t' + str(M_min) + '\t' + ' '.join(list(map(str, MFD)))+'\n')
                    #while MFD[i_Mmin] == 0:
                    #    i_Mmin += 1
                    Line = '\t\t\t<incrementalMFD minMag=\"'+ str(M_min)+'\" binWidth=\"0.10\">\n'# + 0.1 * i_Mmin
                    XMLfile.write(Line)
                    #MFD_model += MFD
                    Line = '\t\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'#[-(len(MFD) - i_Mmin):]
                    XMLfile.write(Line)
                    Line = '\t\t\t</incrementalMFD>\n'
                    XMLfile.write(Line)
                    
                    #find the dominant kinematic of the scenario
#                    unique,pos = np.unique(scenario_mechanism,return_inverse=True) #Finds all unique elements and their positions
#                    counts = np.bincount(pos)                     #Count the number of each unique element
#                    maxpos = counts.argmax()                      #Finds the positions of the maximum count
#                    rake = unique[maxpos]
                    rake= np.mean(scenario_mechanism)
                    Line='\t\t\t<rake>'+str(rake)+'</rake>\n'
                    XMLfile.write(Line)
                    Line='\t\t</characteristicFaultSource>\n'
                    XMLfile.write(Line)
                index_scenario += 1
                
                
        '''#########################
        # Defining the background
        #########################'''
        self.Geom_bg()
        self.Prop_bg()
        
        MFD = MFDs.EQ_rate_BG   
        if sum(MFD) != 0. :
            #Line='\t\t<areaSource id="'+ str(self.Model_name) + '_bg_' + str(ID_number + 1 ) +'" name="Background" tectonicRegion="' + str(self.Domain) + '">\n'
            Line='\t\t<areaSource id="'+ str(ID_number + 1 ) +'" name="Background" tectonicRegion="' + str(self.Domain) + '">\n'
            XMLfile.write(Line)
            Line='\t\t\t<areaGeometry>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t<gml:Polygon>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t\t<gml:exterior>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t\t\t<gml:LinearRing>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t\t\t\t<gml:posList>\n'  
            XMLfile.write(Line)
            for x,y in zip(self.Lon_bg,self.Lat_bg):
                #polygon.append((x,y)) #ecriture du polygone de la zone
                Line='\t\t\t\t\t\t\t\t' + str(x) + ' ' + str(y) + '\n'
                XMLfile.write(Line)
            Line='\t\t\t\t\t\t\t</gml:posList>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t\t\t</gml:LinearRing>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t\t</gml:exterior>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t</gml:Polygon>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t<upperSeismoDepth>' + str(self.upperSeismoDepth) + '</upperSeismoDepth>\n'
            XMLfile.write(Line)
            Line='\t\t\t\t<lowerSeismoDepth>' + str(self.lowerSeismoDepth) + '</lowerSeismoDepth>\n'
            XMLfile.write(Line)
            Line='\t\t\t</areaGeometry>\n'
            XMLfile.write(Line)
            Line='\t\t\t<magScaleRel>'+ ScL_oq +'</magScaleRel>\n'
            XMLfile.write(Line)
            Line='\t\t\t<ruptAspectRatio>' + str(self.ruptAspectRatio) + '</ruptAspectRatio>\n'
            XMLfile.write(Line)         
            log_mdf_file.write('Background' + '\t' + str(M_min) + '\t' + ' '.join(list(map(str, MFD)))+'\n')
            Line='\t\t\t<incrementalMFD binWidth=\"0.10\" minMag="'+ str(M_min)+'">\n'
            XMLfile.write(Line)
            #MFD_model += MFD
            Line='\t\t\t<occurRates> ' + ' '.join(list(map(str, MFD))) + '</occurRates>\n'
            XMLfile.write(Line)
            Line = '\t\t\t</incrementalMFD>\n'
            XMLfile.write(Line)
            Line='\t\t\t<nodalPlaneDist>\n'
            XMLfile.write(Line)
            for i in range(len(self.nodalPlanes)) :
                Line=('\t\t\t\t<nodalPlane probability="' + str(self.nodalPlanes[i][0]) + '" strike="' + str(self.nodalPlanes[i][1]) + '" dip="' 
                      + str(self.nodalPlanes[i][2]) + '" rake="' + str(self.nodalPlanes[i][3]) + '" />\n')
                XMLfile.write(Line)
            Line='\t\t\t</nodalPlaneDist>\n'
            XMLfile.write(Line)
            Line='\t\t\t<hypoDepthDist>\n'
            XMLfile.write(Line)
            for i in range(len(self.hypoDepths)) :
                Line='\t\t\t\t<hypoDepth probability="' + str(self.hypoDepths[i][0]) + '" depth="' + str(self.hypoDepths[i][1]) + '" />\n'
                XMLfile.write(Line)
            Line='\t\t\t</hypoDepthDist>\n'
            XMLfile.write(Line)
            Line='\t\t</areaSource>\n'
            XMLfile.write(Line)
        
        '''#############################
        ### defining the other sources based on the host model
        ##############################'''
        if self.use_host_model == True :
            # extract the geometry of the zone ( geometry of the background)
            self.Geom_bg()                      
            Poly = []   
            for x1,y1 in zip(self.Lon_bg,self.Lat_bg ): # creation du polygon de la zone
                Poly.append((x1,y1))    
            bbPath = mplPath.Path(Poly)
        
            read_host_file = open(self.host_model_file,'r')
            lines_of_the_host_file = read_host_file.readlines()
            lines_of_the_host_file = [x.strip('L\n') for x in lines_of_the_host_file] 
            lines_of_the_host_file = [x.strip('\r\n') for x in lines_of_the_host_file] 
            lines_of_the_host_file = [x.strip('\n') for x in lines_of_the_host_file] 
            line_number = 0
            source_read = False
            simple_fault = False
            complex_fault = False
            area_source = False
            point_source = False
            subduction_source = False
            
            for line in lines_of_the_host_file:
                if '<simpleFaultSource' in line :
                    line_start = line_number
                    index_id = line.find('id="')+4
                    source_read = False
                    simple_fault = True
                    complex_fault = False
                    area_source = False
                    point_source = False
                    subduction_source = False
                    if 'Subduction' in line or 'subduction' in line:
                        subduction_source = True
                if '<areaSource' in line :
                    Xing_bg = False
                    type_increment = False
                    type_mfd = False
                    zone_defined = False
                    line_start = line_number
                    index_id = line.find('id="')+4
                    source_read = False
                    simple_fault = False
                    complex_fault = False
                    area_source = True
                    point_source = False
                    subduction_source = False
                    if 'Subduction' in line or 'subduction' in line:
                        subduction_source = True
                if '<complexFaultSource' in line :
                    line_start = line_number
                    index_id = line.find('id="')+4
                    source_read = False
                    simple_fault = False
                    complex_fault = True
                    area_source = False
                    point_source = False
                    subduction_source = False
                    if 'Subduction' in line or 'subduction' in line:
                        subduction_source = True
                if '<pointSource' in line :
                    line_start = line_number
                    index_id = line.find('id="')+4
                    source_read = False
                    simple_fault = False
                    complex_fault = False
                    area_source = False
                    point_source = True
                    subduction_source = False
                    if 'Subduction' in line or 'subduction' in line:
                        subduction_source = True
                if '</sourceModel' in line :
                    source_read = False
                    simple_fault = False
                    complex_fault = False
                    area_source = False
                    point_source = False
                    subduction_source = False
                    
                    
                if simple_fault == True or complex_fault == True:
                    print_source = True
                    if '<gml:posList>' in line :
                        line_start_lonlat = line_number 
                    if '</gml:posList>' in line :
                        line_stop_lonlat = line_number 
                        lon_lat = ''
                        for line_lon_lat in lines_of_the_host_file[line_start_lonlat:line_stop_lonlat+1]:
                            line_lon_lat = line_lon_lat.replace('<gml:posList>','')
                            line_lon_lat = line_lon_lat.replace('</gml:posList>','')
                            lon_lat += ' ' + line_lon_lat
                        lon_lat = lon_lat.replace('  ',' ')
                        lon_lat = lon_lat.replace('  ',' ')
                        lon_lat = lon_lat.split(' ')
                        points = []
                        for i in range(len(lon_lat)-1):
                            if lon_lat[i] !=  '':
                                if len(points)!= 0:
                                    if float(lon_lat[i]) != points[-1][1]:
                                        points.append([float(lon_lat[i]), float(lon_lat[i+1])])
                                else :
                                    points.append([float(lon_lat[i]), float(lon_lat[i+1])])
                        for point in points :
                            if bbPath.contains_point((point[0],point[1])) == True:
                                print_source = False
                    if '<\simpleFaultSource' in line or '</complexFaultSource' in line:
                        line_end = line_number
                        source_read = True
                    if print_source == True or subduction_source == True:
                        if source_read == True :
                            line_to_print = lines_of_the_host_file[line_start][:index_id]+'1111'+lines_of_the_host_file[line_start][index_id:]
                            XMLfile.write(line_to_print+'\n')
                            for line_to_print in lines_of_the_host_file[line_start+1:line_end+1] :
                                XMLfile.write(line_to_print+'\n')
                            
                if point_source == True:
                    print_source = True
                    if '<gml:posList>' in line :
                        line_start_lonlat = line_number 
                    if '</gml:posList>' in line :
                        line_stop_lonlat = line_number 
                        lon_lat = ''
                        for line_lon_lat in lines_of_the_host_file[line_start_lonlat:line_stop_lonlat+1]:
                            line_lon_lat = line_lon_lat.replace('<gml:posList>','')
                            line_lon_lat = line_lon_lat.replace('</gml:posList>','')
                            lon_lat += ' ' + line_lon_lat
                        lon_lat = lon_lat.replace('  ',' ')
                        lon_lat = lon_lat.replace('  ',' ')
                        lon_lat = lon_lat.split(' ')
                        points = []
                        for i in range(len(lon_lat)-1):
                            if lon_lat[i] !=  '':
                                if len(points)!= 0:
                                    if float(lon_lat[i]) != points[-1][1]:
                                        points.append([float(lon_lat[i]), float(lon_lat[i+1])])
                                else :
                                    points.append([float(lon_lat[i]), float(lon_lat[i+1])])
                        for point in points :
                            if bbPath.contains_point((point[0],point[1])) == True:
                                print_source = False
                    if '<\pointSource' in line:
                        line_end = line_number
                        source_read = True
                    if print_source == True or subduction_source == True:
                        if source_read == True :
                            line_to_print = lines_of_the_host_file[line_start][:index_id]+'1111'+lines_of_the_host_file[line_start][index_id:]
                            XMLfile.write(line_to_print+'\n')
                            for line_to_print in lines_of_the_host_file[line_start+1:line_end+1] :
                                XMLfile.write(line_to_print+'\n')
                        
                if area_source == True:
                    if '<gml:posList>' in line :
                        line_start_lonlat = line_number 
                    if '</gml:posList>' in line :
                        line_stop_lonlat = line_number 
                        lon_lat = ''
                        for line_lon_lat in lines_of_the_host_file[line_start_lonlat:line_stop_lonlat+1]:
                            line_lon_lat = line_lon_lat.replace('<gml:posList>','')
                            line_lon_lat = line_lon_lat.replace('</gml:posList>','')
                            lon_lat += ' ' + line_lon_lat
                        lon_lat = lon_lat.replace('  ',' ')
                        lon_lat = lon_lat.replace('  ',' ')
                        lon_lat = lon_lat.split(' ')
                        points_zone = []
                        for i in range(len(lon_lat)-1):
                            if lon_lat[i] !=  '':
                                if len(points_zone)!= 0:
                                    if float(lon_lat[i]) != points_zone[-1][1]:
                                        points_zone.append([float(lon_lat[i]), float(lon_lat[i+1])])
                                else :
                                    points_zone.append([float(lon_lat[i]), float(lon_lat[i+1])])
                        ColX = []
                        ColY = []
                        for point in points_zone :
                            ColX.append(point[0])
                            ColY.append(point[1]) 
                            if bbPath.contains_point((point[0],point[1])) == True:
                                Xing_bg = True
                        Poly = []   
                        for x1,y1 in zip(ColX,ColY): 
                            Poly.append((x1,y1))    
                        bbPath_zone = mplPath.Path(Poly)
                        for lon,lat in zip(self.Lon_bg,self.Lat_bg):
                            if bbPath_zone.contains_point((lon,lat)) == True:
                                Xing_bg = True
                    if '</areaSource>' in line:
                        #print lines_of_the_host_file[line_start][-9:],Xing_bg,subduction_source
                        line_end = line_number
                        source_read = True
                    if Xing_bg == False or subduction_source == True: #doesn't change anything, just write the source down
                        if source_read == True :
                            line_to_print = lines_of_the_host_file[line_start][:index_id]+'1111'+lines_of_the_host_file[line_start][index_id:]
                            XMLfile.write(line_to_print+'\n')
                            for line_to_print in lines_of_the_host_file[line_start+1:line_end+1] :
                                XMLfile.write(line_to_print+'\n')    
                    elif Xing_bg ==True and subduction_source == False:
                        if zone_defined == False :
                            zone_defined = True
                            #number_of_Xing=0
                            listpoint_in_bg = []
                            ColX = []
                            ColY = []
                            for point in points_zone :
                                # find the points that are in the zone and the couple of point that crosses the other zone limit
    #                            Xsing_couple_zone = []
                                ColX.append(point[0])
                                ColY.append(point[1])          
                                if bbPath.contains_point((point[0],point[1])) == True:
                                    listpoint_in_bg.append(1)
                                else:
                                    listpoint_in_bg.append(0)
                                
                            #does the same for the background points
    #                        Xsing_points_bg =[]
                            bg_point_inzone = []            
                            Poly = []   
                            for x1,y1 in zip(ColX,ColY): 
                                Poly.append((x1,y1))    
                            bbPath_zone = mplPath.Path(Poly)
                            for lon,lat in zip(self.Lon_bg,self.Lat_bg):
                                if bbPath_zone.contains_point((lon,lat)) == True:
                                    bg_point_inzone.append(1)
                                else :
                                    bg_point_inzone.append(0)
                            #find the number of crossing
                            nb_Xing_zone = 0
                            for index_pt_z in range(len(listpoint_in_bg)):
                                if index_pt_z != 0:
                                    if listpoint_in_bg[index_pt_z]!=listpoint_in_bg[index_pt_z-1]:
                                        nb_Xing_zone+=1
                                else:
                                    if listpoint_in_bg[-1]!=listpoint_in_bg[index_pt_z]:
                                        nb_Xing_zone+=1
                                    
                            nb_Xing_bg = 0
                            for index_pt_bg in range(len(bg_point_inzone)):
                                if index_pt_bg != 0:
                                    if bg_point_inzone[index_pt_bg]!=bg_point_inzone[index_pt_bg-1]:
                                        nb_Xing_bg+=1
                                else:
                                    if bg_point_inzone[-1]!=bg_point_inzone[index_pt_bg]:
                                        nb_Xing_bg+=1
                            number_of_Xing = max([nb_Xing_zone,nb_Xing_bg])
                                
                            
                            if sum(bg_point_inzone) == len(bg_point_inzone) and sum(listpoint_in_bg)==0: # if the background if completely included in the zone
                                lon_zone_modif = []
                                lat_zone_modif = []
                                for point in points_zone :
                                    lon_zone = point[0]
                                    lat_zone = point[1]
                                    lon_zone_modif.append(lon_zone)
                                    lat_zone_modif.append(lat_zone)
                                distances = []
                                for lon_bg,lat_bg in zip(self.Lon_bg,self.Lat_bg):
                                    distances.append(self.distance(lon_bg,lat_bg,points_zone[-1][0],points_zone[-1][1]))
                                index_dist_min = np.argmin(distances)
                                lon_bg_modif = self.Lon_bg[index_dist_min:]+self.Lon_bg[:index_dist_min]
                                lat_bg_modif = self.Lat_bg[index_dist_min:]+self.Lat_bg[:index_dist_min]
                                if (self.distance(lon_bg_modif[-1],lat_bg_modif[-1],points_zone[0][0],points_zone[0][1])
                                >self.distance(lon_bg_modif[0],lat_bg_modif[0],points_zone[0][0],points_zone[0][1])):
                                    lon_bg_modif = list(reversed(lon_bg_modif))
                                    lat_bg_modif = list(reversed(lat_bg_modif))
                                for lon_bg,lat_bg in zip(lon_bg_modif,lat_bg_modif):
                                    lon_zone_modif.append(lon_bg)
                                    lat_zone_modif.append(lat_bg)
                                    
                                    
#                                lon_zone_modif.append(lon_bg_modif[0])
#                                lat_zone_modif.append(lat_bg_modif[0])
#                                lon_zone_modif.append(points_zone[-1][0])
#                                lat_zone_modif.append(points_zone[-1][1])
                                
                                #avoid intersections of the bg and moves the point slightly so OQ doesn't freak out
                                line1 = [[points_zone[-1][0],points_zone[-1][1]],[lon_bg_modif[0],lat_bg_modif[0]]]
                                line2 = [[np.mean([lon_bg_modif[0],points_zone[-1][0]])+0.0001,np.mean([lat_bg_modif[0],points_zone[-1][1]])],
                                          [np.mean([lon_bg_modif[0],points_zone[-1][0]]),np.mean([lat_bg_modif[-1],points_zone[-1][1]])]]                                
                                x,y = line_intersection(line1, line2)
                                if x != 'no_intesection':
                                    if (points_aligned([np.mean([lon_bg_modif[0],points_zone[-1][0]])+0.0001,np.mean([lat_bg_modif[0],points_zone[-1][1]])],
                                                        [np.mean([lon_bg_modif[0],points_zone[-1][0]]),np.mean([lat_bg_modif[-1],points_zone[-1][1]])],
                                                         [x,y]) == False):
                                        #print 'to the right'
                                        lon_zone_modif.append(lon_bg_modif[0]+0.0001)
                                        lat_zone_modif.append(lat_bg_modif[0])
                                        lon_zone_modif.append(points_zone[-1][0]+0.0001)
                                        lat_zone_modif.append(points_zone[-1][1])
                                    else:
                                        #print 'to the left!'
                                        lon_zone_modif.append(lon_bg_modif[0]-0.0001)
                                        lat_zone_modif.append(lat_bg_modif[0])
                                        lon_zone_modif.append(points_zone[-1][0]-0.0001)
                                        lat_zone_modif.append(points_zone[-1][1])
#                                plt.plot(lon_zone_modif,lat_zone_modif,':k')
#                                plt.scatter(points_zone[-1][0],points_zone[-1][1],c='grey')
#                                plt.scatter(points_zone[-2][0],points_zone[-2][1],c='k')
#                                plt.scatter(lon_bg_modif[0],lat_bg_modif[0],c='cyan')
#                                plt.scatter(lon_bg_modif[1],lat_bg_modif[1],c='b',linewidth=0.01)
#                                plt.scatter(lon_bg_modif[-1],lat_bg_modif[-1],c='b')
#                                plt.scatter(x,y,c='purple')
#                                plt.scatter(np.mean([lon_bg_modif[0],points_zone[-1][0]])+0.0001,np.mean([lat_bg_modif[0],points_zone[-1][1]]),c='r')
#                                plt.scatter(np.mean([lon_bg_modif[-1],points_zone[-1][0]]),np.mean([lat_bg_modif[-1],points_zone[-1][1]]),c='pink')
#                                plt.show()
#                                print(lon_zone_modif[-5:])
  
                                
                            else:
                                lon_zone_modif = []
                                lat_zone_modif = []
                                index_point_z = 0
                                #loop on the points of the zone
                                for point in points_zone :
                                    lon_zone = point[0]
                                    lat_zone = point[1]
                                    if listpoint_in_bg[index_point_z] == 0 : #if the point is not in the bg, we add it
                                        lon_zone_modif.append(lon_zone)
                                        lat_zone_modif.append(lat_zone)
                                        
                                    index_bg_intercept = None    
                                    if index_point_z != len(points_zone)-1: #it's not the last point
                                        if (listpoint_in_bg[index_point_z] == 0 and listpoint_in_bg[index_point_z+1] == 1
                                            )or(listpoint_in_bg[index_point_z] == 1 and listpoint_in_bg[index_point_z+1] == 0): #it's a crossing point
                                            
                                            index_point_bg = 0
                                            for lon_bg,lat_bg in zip(self.Lon_bg,self.Lat_bg):
                                                if index_point_bg != len(bg_point_inzone)-1: #it's not the last point
                                                        #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]], [x,y])==True):
                                                            #print('Xing1')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[index_point_bg] == 1:
                                                                index_bg_intercept = index_point_bg
                                                            else:
                                                                index_bg_intercept = index_point_bg+1
                                                                
                                                else: #it's the last point
                                                    #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[0],self.Lat_bg[0]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[0],self.Lat_bg[0]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]], [x,y])==True):
                                                            #print('Xing2')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[-1] == 1:
                                                                index_bg_intercept = -1
                                                            else:
                                                                index_bg_intercept = 0
                                                            
                                                index_point_bg +=1    
                                        else : #check if a bg point intersects but no zone points are in the background
                                            index_point_bg = 0
                                            for lon_bg,lat_bg in zip(self.Lon_bg,self.Lat_bg):
                                                if index_point_bg != len(bg_point_inzone)-1: #it's not the last point
                                                        #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]], [x,y])==True):
                                                            #print('Xing1')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[index_point_bg] == 1:
                                                                index_bg_intercept = index_point_bg
                                                            else:
                                                                index_bg_intercept = index_point_bg+1
                                                                
                                                else: #it's the last point
                                                    #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[0],self.Lat_bg[0]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[0],self.Lat_bg[0]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[index_point_z+1][0],points_zone[index_point_z+1][1]], [x,y])==True):
                                                            #print('Xing2')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[-1] == 1:
                                                                index_bg_intercept = -1
                                                            else:
                                                                index_bg_intercept = 0
                                                            
                                                index_point_bg +=1    
                                                
                                    else : #it's the last point
                                        if (listpoint_in_bg[index_point_z] == 0 and listpoint_in_bg[0] == 1
                                            )or(listpoint_in_bg[index_point_z] == 1 and listpoint_in_bg[0] == 0): #it's a crossing point
                                            index_point_bg = 0
                                            for lon_bg,lat_bg in zip(self.Lon_bg,self.Lat_bg):
                                                if index_point_bg != len(bg_point_inzone)-1: #it's not the last point
                                                        #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[0][0],points_zone[0][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[0][0],points_zone[0][1]], [x,y])==True):
                                                            #print('Xing3')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[index_point_bg] == 1:
                                                                index_bg_intercept = index_point_bg
                                                            else:
                                                                index_bg_intercept = index_point_bg+1
                                                else: #it's the last point
                                                    #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[0][0],points_zone[0][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[0],self.Lat_bg[0]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[0],self.Lat_bg[0]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[0][0],points_zone[0][1]], [x,y])==True):
                                                            #print('Xing4')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[-1] == 1:
                                                                index_bg_intercept = -1
                                                            else:
                                                                index_bg_intercept = 0
                                                index_point_bg +=1
                                        else:
                                            index_point_bg = 0
                                            for lon_bg,lat_bg in zip(self.Lon_bg,self.Lat_bg):
                                                if index_point_bg != len(bg_point_inzone)-1: #it's not the last point
                                                        #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[0][0],points_zone[0][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[index_point_bg+1],self.Lat_bg[index_point_bg+1]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[0][0],points_zone[0][1]], [x,y])==True):
                                                            #print('Xing3')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[index_point_bg] == 1:
                                                                index_bg_intercept = index_point_bg
                                                            else:
                                                                index_bg_intercept = index_point_bg+1
                                                else: #it's the last point
                                                    #find the intersecting point
                                                    line1 = [[lon_zone,lat_zone],[points_zone[0][0],points_zone[0][1]]]
                                                    line2 = [[lon_bg,lat_bg],[self.Lon_bg[0],self.Lat_bg[0]]]
                                                    x,y = line_intersection(line1, line2)
                                                    if x != 'no_intesection':
                                                        if (points_aligned([lon_bg,lat_bg], [self.Lon_bg[0],self.Lat_bg[0]], [x,y]) == True
                                                            )and(points_aligned([lon_zone,lat_zone], [points_zone[0][0],points_zone[0][1]], [x,y])==True):
                                                            #print('Xing4')
                                                            lon_zone_modif.append(float(x))
                                                            lat_zone_modif.append(float(y))
                                                            if bg_point_inzone[-1] == 1:
                                                                index_bg_intercept = -1
                                                            else:
                                                                index_bg_intercept = 0
                                                index_point_bg +=1
                                                
                                        
                                    #integrating the bg points in the modified source
                                    if listpoint_in_bg[index_point_z] == 0 and index_bg_intercept != None: #only need to write bg points if the zone point is outside of the bg
                                        lon_bg_modif = self.Lon_bg[index_bg_intercept:]+self.Lon_bg[:index_bg_intercept]
                                        lat_bg_modif = self.Lat_bg[index_bg_intercept:]+self.Lat_bg[:index_bg_intercept]
                                        bg_point_inzone_modif = bg_point_inzone[index_bg_intercept:]+bg_point_inzone[:index_bg_intercept]
                                        if index_bg_intercept != 0 and bg_point_inzone_modif[-1]==1:
                                            lon_bg_modif = list(reversed(lon_bg_modif))
                                            lat_bg_modif = list(reversed(lat_bg_modif))
                                            bg_point_inzone_modif = list(reversed(bg_point_inzone_modif))
                                        #print(bg_point_inzone_modif)
                                        i = 0
                                        while bg_point_inzone_modif[i] == 1:
                                            lon_zone_modif.append(lon_bg_modif[i])
                                            lat_zone_modif.append(lat_bg_modif[i])
                                            i+=1
                                        
                                    index_point_z +=1 
                                
                            
                            
                            #caltulating the ratio of surfaces between the old zone and the new one
                            # calculate the area of the sub_area
                            x,y = reproject(ColY,ColX) 
                            area_of_the_zone = area_of_polygon( x,y)  #to be verified!!
                            # calculate the area of the sub_area
                            x,y = reproject(lat_zone_modif,lon_zone_modif) 
                            area_of_the_zone_modified = area_of_polygon( x,y)  #to be verified!!
                            ratio_areas = area_of_the_zone_modified/area_of_the_zone
                            
                            if number_of_Xing >=3 :
                                from itertools import groupby
                                #check if the zones needs to be cut in different small zone for OQ to handle it
                                #if one of the points from the subzone is included in the zone without these points 
                                #print bg_point_inzone
                                indexes_for_test_init = [list(j) for i,j in groupby(listpoint_in_bg)]
                                #print(indexes_for_test_init)
                                indexes_for_test = []
                                index = 0
                                for i in indexes_for_test_init:
                                    indexes_for_test_i=[]
                                    for ii in i:
                                        if ii==1:
                                            indexes_for_test_i.append(index)
                                        index+=1
                                    if ii==1:
                                        indexes_for_test.append(indexes_for_test_i)
                                #print(indexes_for_test)
                                indexes_for_test_modif = []
                                for indexes in list(indexes_for_test):
                                    indexes_modif = []
                                    if len(indexes)> 1:
                                        for index in indexes :
                                            i=0
                                            for lon,lat in zip(lon_zone_modif,lat_zone_modif):
#                                                print indexes
#                                                print lon,lat
#                                                print ColX[index],len(ColX),ColY[index],len(ColY)
                                                if lon == ColX[index] and lat == ColY[index]:
                                                    indexes_modif.append(i)
                                                i+=1
                                    else :
                                        i=0
                                        for lon,lat in zip(lon_zone_modif,lat_zone_modif):
#                                            print indexes
#                                            print lon,lat
#                                            print ColX[indexes[0]],ColY[indexes[0]]
                                            if lon == ColX[indexes[0]] and lat == ColY[indexes[0]]:
                                                indexes_modif.append(i)
                                            i+=1
                                    indexes_for_test_modif.append(indexes_modif)
                                
                                for indexes in list(indexes_for_test_modif):
                                    if len(indexes)> 1:
                                        indexes = sorted(indexes)
                                        lon_to_test=lon_zone_modif
                                        lat_to_test=lat_zone_modif
                                        i=0
                                        #print(indexes)
                                        for index in indexes :
                                            #print(len(lon_to_test))
                                            del lon_to_test[index-i]
                                            del lat_to_test[index-i]
                                            i+=1
                                        x,y = reproject(lat_to_test,lon_to_test) 
                                        area_of_the_zone_to_test = area_of_polygon( x,y)  #to be verified!!
                                        ratio_areas_test = area_of_the_zone_to_test/area_of_the_zone
                                    else :
                                        lon_to_test=lon_zone_modif
                                        lat_to_test=lat_zone_modif
                                        #print(indexes,len(lon_to_test))
                                        if len(indexes)!=0:
                                            del lon_to_test[indexes[0]]
                                            del lat_to_test[indexes[0]]
                                        x,y = reproject(lat_to_test,lon_to_test) 
                                        area_of_the_zone_to_test = area_of_polygon( x,y)  #to be verified!!
                                        ratio_areas_test = area_of_the_zone_to_test/area_of_the_zone
                                        
                                    if ratio_areas_test > 1. :
                                        print('included')
                                    else :
                                        print('excluded')
                        
                        if sum(bg_point_inzone) == len(bg_point_inzone) and sum(listpoint_in_bg)==0:
                            x,y = reproject(lat_bg_modif,lon_bg_modif)
                            area_of_the_bg = area_of_polygon( x,y)  #to be verified!!
                            ratio_areas = (area_of_the_zone-area_of_the_bg)/area_of_the_zone
                        
                        #reading the seismic actvity parameters        
                        if 'hterMFD' in line:
                            line_mfd_start_number = line_number
                            type_mfd = True
                        if type_mfd == True and source_read ==True:
                            index_a = line.find('aValue="')+8
                            if line.find('aValue="') == -1:
                                index_a = line.find('aValue= "')+9
                            if line.find('aValue="') == -1:
                                print('PROBLEM with reading the host file to find the a value')
                            a_str = line[index_a:]
                            i = a_str.find('"')
                            a = float(a_str[:i])
                            a_modif = a * ratio_areas
                            line_mfd_modif = line[:index_a]+str(a_modif)+line[line.find('"',index_a+1):]
                        if  '<occurRates'  in line:
                            line_mfd_start_number = line_number
                            type_increment = True
                        if  '/occurRates'  in line:
                            line_mfd_stop_number = line_number
                        
                        if type_increment == True and source_read ==True:
                            mfd_str = ''
                            for line_mfd_str in lines_of_the_host_file[line_mfd_start_number:line_mfd_stop_number+1]:
                                mfd_str += ' ' + line_mfd_str
                            mfd_str = mfd_str.replace('<occurRates>','')
                            mfd_str = mfd_str.replace('</occurRates>','')
                            mfd_str = mfd_str.split(' ')
                            mfd_modif = []
                            for value in mfd_str:
                                if value != '':
                                    mfd_modif.append(float(value)*ratio_areas)
                            line_mfd_modif = '<occurRates>'
                            for value in mfd_modif:
                                line_mfd_modif += str(value)+' '
                            line_mfd_modif+='</occurRates>'
                            
                        if source_read == True :
                            line_to_print = lines_of_the_host_file[line_start][:index_id]+'1111'+lines_of_the_host_file[line_start][index_id:]
                            XMLfile.write(line_to_print+'\n')
                            for line_to_print in lines_of_the_host_file[line_start+1:line_start_lonlat] :
                                XMLfile.write(line_to_print+'\n')
                                
                            line_geom = '<gml:posList> '
                            for lon,lat in zip(lon_zone_modif,lat_zone_modif):
                                line_geom += str(lon) + ' ' + str(lat) + ' '
                            line_geom += '</gml:posList> '
                            XMLfile.write(line_geom+'\n')
                            
                            if number_of_Xing >= 3:
                                #print('please check if the host model is incorporate correctly, problems might have occured!!!')
                                #print(lines_of_the_host_file[line_start][-9:-2],'number_of_Xing=',number_of_Xing)
                                #print('ratio_areas',ratio_areas)
                                plt.scatter(ColX,ColY,c='b',alpha=0.2)
                                plt.scatter(self.Lon_bg,self.Lat_bg,c='r',alpha=0.2)
                                plt.scatter(lon_zone_modif,lat_zone_modif,c='k',alpha=0.2,marker='s')
                                plt.plot(lon_zone_modif,lat_zone_modif,':k')
                                plt.xlim(min(self.Lon_bg)-0.5,max(self.Lon_bg)+0.5)
                                plt.ylim(min(self.Lat_bg)-0.5,max(self.Lat_bg)+0.5)
                                plt.show()
                            
                            for line_to_print in lines_of_the_host_file[line_stop_lonlat+1:line_mfd_start_number] :
                                XMLfile.write(line_to_print+'\n')
                            XMLfile.write(line_mfd_modif+'\n')
                            for line_to_print in lines_of_the_host_file[line_mfd_stop_number+1:line_end+1] :
                                XMLfile.write(line_to_print+'\n')  
                line_number +=1
            
        #end of the file
        Line='\t</sourceModel>\n'
        XMLfile.write(Line)
        Line='</nrml>\n'
        XMLfile.write(Line)
        XMLfile.close()
        log_sr_file.close()
        log_mdf_file.close()
        
        


                
    def FaultGeometry(self): 
        #CritereDistance = 3. 
        NomFichier_InfosZonage = self.File_geom
        InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8'),('U100')],skip_header = 1)
        Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
        index_model = np.where(np.array(Column_model_name) == self.Model_name)
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

    def FaultProperties(self,Name_of_fault):   
        FileName_Prop = self.File_prop
        Prop = np.genfromtxt(FileName_Prop,
                                   dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
                                          ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
        Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
        Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
        index_model = np.where(np.array(Column_model_name) == self.Model_name)[0]
        
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

    def Prop_bg(self):
        prop_bg = open(self.file_prop_bg,'r').readlines()
        # background general parameters read from the input file
        self.nodalPlanes = []
        self.hypoDepths = []
        for line in prop_bg :
            if line.split('\t')[0] == self.Model_name:
                if line.split('\t')[1] == 'upperSeismoDepth':
                    self.upperSeismoDepth = float(line.split('\t')[2])
                if line.split('\t')[1] == 'lowerSeismoDepth':
                    self.lowerSeismoDepth = float(line.split('\t')[2])
                if line.split('\t')[1] == 'ruptAspectRatio':
                    self.ruptAspectRatio = float(line.split('\t')[2])
                if line.split('\t')[1] == 'nodalPlane':
                    self.nodalPlanes.append([float(line.split('\t')[2]),
                                                   float(line.split('\t')[3]),
                                                         float(line.split('\t')[4]),
                                                               float(line.split('\t')[5])])
                if line.split('\t')[1] == 'hypoDepth':
                    self.hypoDepths.append([float(line.split('\t')[2]),
                                                   float(line.split('\t')[3])])
        if len(str(self.nodalPlanes))==0:
            print('Error!! Verify your Background parameters file')                               
    def calculate_initial_compass_bearing(self,pointA, pointB):
        """
        Calculates the bearing between two points.
    
        The formulae used is the following:
            θ = atan2(sin(Δlong).cos(lat2),
                      cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    
        :Parameters:
          - `pointA: The tuple representing the latitude/longitude for the
            first point. Latitude and longitude must be in decimal degrees
          - `pointB: The tuple representing the latitude/longitude for the
            second point. Latitude and longitude must be in decimal degrees
    
        :Returns:
          The bearing in degrees
    
        :Returns Type:
          float
        """
        if (type(pointA) != tuple) or (type(pointB) != tuple):
            raise TypeError("Only tuples are supported as arguments")
    
        lat1 = math.radians(pointA[0])
        lat2 = math.radians(pointB[0])
    
        diffLong = math.radians(pointB[1] - pointA[1])
    
        x = math.sin(diffLong) * math.cos(lat2)
        y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
                * math.cos(lat2) * math.cos(diffLong))
    
        initial_bearing = math.atan2(x, y)
    
        # Now we have the initial bearing but math.atan2 return values
        # from -180° to + 180° which is not what we want for a compass bearing
        # The solution is to normalize the initial bearing as shown below
        initial_bearing = math.degrees(initial_bearing)
        compass_bearing = (initial_bearing + 360) % 360
    
        return compass_bearing

    def Geom_bg(self):
        self.Lon_bg = []
        self.Lat_bg = []

        # manually defined  in the file Background geometry
        geom_bg = np.genfromtxt(self.File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
        
        column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
        index_model = np.where(np.array(column_model) == self.Model_name)[0]
        self.Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
        self.Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
        
        
        if self.Lon_bg == 0:
            print('Error!! Check your input background geometry')
            
    def distance(self,lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points 
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians 
        lon1, lat1, lon2, lat2 = list(map(radians, [lon1, lat1, lon2, lat2]))
        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 
        km = 6367 * c
        return km    
