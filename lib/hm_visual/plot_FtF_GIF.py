#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""


import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import xml.etree.ElementTree as ET
import numpy as np
import math
import pylab as pl
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import Read_file as Read_file
import maps.maps_utils as maps_utils
import maps.geom as geom


import utils.read_input as read_input


def map_faults(Run_name,Model_list,scenarios_names_list,
               ScL_complet_list, BG_hyp_list,
               sample_list,b_value_list,MFD_type_list,
               llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,File_bg,
               FileName_Prop,plot_sr_use,visual_FtF,sub_area_file):
    nb_on_maps = False
    
    #exctracting all the scenario set information
    available_sets = read_input.extract_sc_input('input/'+Run_name+'/ruptures.txt')
    
    for Model in Model_list :
        for scenario_set in scenarios_names_list:
            file_names = []
            if not os.path.exists(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set):
                os.makedirs(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set)
            file_source = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp_list[0]) + '/' 
                            + str(ScL_complet_list[0]) + '/sc_' +  str(scenario_set) + '/' 
                            + str(b_value_list[0]) + '/' + 'MFD_'+ str(MFD_type_list[0])
                            + '/Source_model_1.xml') 
            
            # extract the Mmax of the faults and the scenarios
            log_Mmax_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp_list[0]) + '/' 
                            + str(ScL_complet_list[0]) + '/sc_' +  str(scenario_set) + '/' 
                            + str(b_value_list[0]) + '/' + 'MFD_'+ str(MFD_type_list[0])
                            +  '/Log/Mmax_sample_1.txt')                
            sources_names,sources_Mmax,sources_Lengths,sources_Areas = Read_file.read_file_Mmax_log(log_Mmax_file) #read the log of Mmax associated with the model
            
            ########################################################
            #Extraction of the faults and scenarios present in the model from the text file
            ########################################################
#            File_faults_n_scenarios = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp_list[0]) + '/'
#                            + str(ScL_complet_list[0]) + '/sc_' +  str(scenario_set)
#                            + '/faults_n_scenarios.txt')
#            Infosfaults_n_scenarios = np.genfromtxt(File_faults_n_scenarios,dtype=[('U1000')],delimiter = '\t')
    
            # extracting the complexe multi fault ruptures
            rupture_set = available_sets['sc_'+scenario_set]
                                 
            ############
            ##    map the faults
            ############                     
             
            tree = ET.parse(file_source)
            nrml = tree.getroot()
            #print
            
            
            #getting the info for each source
            source_name = []
            Lon = []
            Lat = []
            Dip = []
            U_sism_depth = []
            L_sism_depth = []
            Lon_bg = []
            fault_type = []
            
            nb_sources = 0
            for sourceModel in nrml:
                for Source in sourceModel:
                    if 'simpleFaultSource' in str(Source):
                        source_name_i = Source.attrib['name']
                        u_sism_depth = float(Source[0][2].text)
                        l_sism_depth = float(Source[0][3].text)
                        dip = float(Source[0][1].text)
                        geom = str(Source[0][0][0].text)
                        geom = np.array(geom.split(),float)
                        lon = geom[::2]       
                        lat = geom[1::2]

                        if not '_f_1' in source_name_i :
                            source_name.append(source_name_i)
                            Dip.append(dip)
                            U_sism_depth.append(u_sism_depth)
                            L_sism_depth.append(l_sism_depth)
                            Lon.append(lon.ravel())
                            Lat.append(lat.ravel())
                            fault_type.append('sf')
                            
                            nb_sources += 1
                            
                    if 'complexFaultSource' in str(Source):
                        source_name_i = Source.attrib['name']
                        dip = 0.
                        u_sism_depth = 0.
                        l_sism_depth = 0.
                        geom = str(Source[0][0][0][0].text)
                        geom = np.array(geom.split(),float)
                        lon_top = geom[::3]       
                        lat_top = geom[1::3]  
                        geom = str(Source[0][-1][0][0].text)
                        geom = np.array(geom.split(),float) 
                        lon_bottom = geom[::3]       
                        lat_bottom = geom[1::3]    

                        if not '_f_1' in source_name_i :
                            source_name.append(source_name_i)
                            Dip.append(dip)
                            U_sism_depth.append(u_sism_depth)
                            L_sism_depth.append(l_sism_depth)
                            Lon.append(np.concatenate([lon_top.ravel(),np.flipud(lon_bottom.ravel())]))
                            Lat.append(np.concatenate([lat_top.ravel(),np.flipud(lat_bottom.ravel())]))
                            fault_type.append('cf')
                            
                            nb_sources += 1
                        
                    if 'areaSource' in str(Source): #fetching for the background
                        geom = str(Source[0][0][0][0][0].text)
                        geom = np.array(geom.split(),float)
                        lon = geom[::2]       
                        lat = geom[1::2]
                        Lon_bg = (lon.ravel())
                        Lat_bg = (lat.ravel())
 
                        
            #'''
        
            if visual_FtF == True and '{' in str(sources_names[-1]):
                index_scenario = 0
                index_Mmax_0 = 0
                i=0
                while index_Mmax_0 == 0:
                    if '{' in str(sources_names[i]):
                        index_Mmax_0=i
                    i+=1
#                if not np.size(Infosfaults_n_scenarios) == 1 :
#                    for index_scenario in range(len(Infosfaults_n_scenarios)-1):
#                        line_scenario = str(Infosfaults_n_scenarios[index_scenario+1][0])
                if not rupture_set == [] :
                    for index_scenario in range(len(rupture_set)):
                        #line_scenario = str(Infosfaults_n_scenarios[index_scenario+1][0])
                        #faults_in_scenario = list(line_scenario.split(' '))
                        faults_in_scenario = rupture_set[index_scenario]
                        m = Basemap(projection='mill',
                                      llcrnrlon=llcrnrlon, 
                                      llcrnrlat=llcrnrlat, 
                                      urcrnrlon=urcrnrlon, 
                                      urcrnrlat=urcrnrlat,resolution='l')
                        
                        if len(Lon_bg) != 0 : #draw the background zone
                            maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.2, 0.5, 'k')
                        Mmax = sources_Mmax[index_Mmax_0+index_scenario]
                        #for each fault  
                        
                        for index_source in range(nb_sources):
                            source_name_i = source_name[index_source].replace(Model+'_','')
                            if fault_type[index_source] == 'sf':
                                dip_tan = math.tan(math.radians(Dip[index_source]))
                                
                                hdist_u = U_sism_depth[index_source] / dip_tan
                                hdist_d = L_sism_depth[index_source] / dip_tan
                                
                                azimuth = [0.,0.,0.,0.]
                                
                                lon_u = []
                                lat_u = []
                                lon_d = []
                                lat_d = []
                                if len(Lon[index_source]) < 9 :
                                    azimuth = []
                                    indexx = range(len(Lon[index_source]))
                                    for i in indexx:
                                        strike = abs(math.degrees(math.atan((
                                        Lon[index_source][-1]-Lon[index_source][0])/(
                                        Lat[index_source][-1]-Lat[index_source][0]))))
                                        if strike > 180. :
                                            strike = strike - 180.
                                        azimuth.append((strike + 90.0) % 360)
                                        
                                    
                                else:
                                    indexx = range(len(Lon[index_source]))
                                    for i in indexx[4:-4]:
                                        strike = abs(math.degrees(math.atan((
                                        Lon[index_source][i+4]-Lon[index_source][i-4])/(
                                        Lat[index_source][i+4]-Lat[index_source][i-4]))))
                                        if strike > 180. :
                                            strike = strike - 180.
                                        azimuth.append((strike + 90.0) % 360)
                                        
                                    azimuth[0] = azimuth[4]
                                    azimuth[1] = azimuth[4]
                                    azimuth[2] = azimuth[4]
                                    azimuth[3] = azimuth[4]
                                    azimuth.append(azimuth[-4])
                                    azimuth.append(azimuth[-4])
                                    azimuth.append(azimuth[-4])
                                    azimuth.append(azimuth[-4])
                                
                                for i in indexx:
                                    if (Lon[index_source][0])>(Lon[index_source][-1]):
                                        x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                                        y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
                                        x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                                        y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
                                        lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                        lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
                                        lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                        lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
                                    if (Lon[index_source][0])<(Lon[index_source][-1]):
                                        x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                                        y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
                                        x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                                        y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
                                        lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                        lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
                                        lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                        lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
                                    
                            
                                
                                source_name_i = source_name[index_source].replace(Model+'_','')
                                
                                x, y = m(Lon[index_source], Lat[index_source]) 
                                if source_name_i in faults_in_scenario:
                                    m.plot(x, y, 'D-', markersize=2.5, linewidth=0.001, color='r', markerfacecolor='r', markeredgewidth = 0.001)
                                else :
                                    m.plot(x, y, 'D-', markersize=0.2, linewidth=0.001, color='k', markerfacecolor='k')
                                
                            
                                poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                                poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                                
                                if len(Lon[index_source]) < 3 :
                                    poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                                    poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                                
                                if source_name_i in faults_in_scenario:
                                    maps_utils.draw_screen_poly(poly_lons, poly_lats,  m ,'r' , 0.5, 0.05, 'r')
                                else :
                                    maps_utils.draw_screen_poly(poly_lons, poly_lats,  m ,'k' , 0.2, 0.05, 'r')
                                    
                            if fault_type[index_source] == 'cf':
                                if source_name_i in faults_in_scenario:
                                    maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m ,'r' , 0.5, 0.05, 'r')
                                else :
                                    maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m ,'k' , 0.2, 0.05, 'k')
                                                                                            
                            
                        Lon_bg = []
                        Lat_bg = []
                    
                        # manually defined  in the file Background geometry
                        geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
                        
                        column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
                        index_model = np.where(np.array(column_model) == Model)[0]
                        Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
                        Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
                        if len(Lon_bg) != 0 : #draw the background zone
                            maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.05, 0.05, 'k')
                            
                        m.drawcoastlines(linewidth=0.2)
                        m.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
                        title = str()
                        i = 0
                        for fault in faults_in_scenario :
                            
                            title += ' ' + str(fault)
                            i += 1
                            if i == 6 :
                                title += '\n'
                                i = 0
                        plt.annotate('Mmax : '+str(Mmax), xy=(0.1, 0.9), xycoords='axes fraction',size=6)
                        
                        plt.title(title)
                        plt.title(str(title))
                        plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set+'/'+str(index_scenario)+'.png',dpi=180)
                        file_names.append(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set+'/'+str(index_scenario)+'.png')
    
                        plt.close()
                        
                    
                    
                    
            # print the map for the model
            m = Basemap(projection='cyl',
                          llcrnrlon=llcrnrlon, 
                          llcrnrlat=llcrnrlat, 
                          urcrnrlon=urcrnrlon, 
                          urcrnrlat=urcrnrlat,resolution='h')

            Lon_bg = []
            Lat_bg = []
        
            # manually defined  in the file Background geometry
            geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)            
            column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
            index_model = np.where(np.array(column_model) == Model)[0]
            Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
            Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
            if len(Lon_bg) != 0 : #draw the background zone
                maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.1, 0.1, 'k')
                x, y = m(Lon_bg, Lat_bg) 
                m.plot(x,y,linewidth=0.2,color='k',linestyle = 'dashed')
            
                
            #draw the sub_areas
            
            #bbPath_sub_areas = []
            if os.path.exists(sub_area_file):
                read_sub_area_file = open(sub_area_file,'rU')
                lines_sub_area = read_sub_area_file.readlines()
                sub_area_names = []
                sub_area_coord = []
#                sub_area_lon = []
#                sub_area_lat = []
                for line in lines_sub_area:
                    model_sub_area = line.split('\t')[0]
                    if model_sub_area in Model_list:
                        sub_area_names.append(line.split('\t')[1])
                        sub_area_coord.append(line.split('\t')[2:])
                        sub_area_lon_i = []
                        sub_area_lat_i = []
                        for sub_area_coord_i in line.split('\t')[2:]:
                            if not '\n' in sub_area_coord_i.split(','):
                                if not '' in sub_area_coord_i.split(','):
                                    sub_area_lon_i.append(float(sub_area_coord_i.split(',')[1]))
                                    sub_area_lat_i.append(float(sub_area_coord_i.split(',')[0]))
                        maps_utils.draw_screen_poly(sub_area_lon_i, sub_area_lat_i,  m ,'k' , 0.01, 0.1, 'k')
                        x, y = m(sub_area_lon_i, sub_area_lat_i) 
                        m.plot(x,y,linewidth=0.2,color='k',linestyle = 'dotted')
                        

            #for each fault              
            for index_source in range(nb_sources):
                if fault_type[index_source] == 'sf':
                    dip_tan = math.tan(math.radians(Dip[index_source]))
                    
                    hdist_u = U_sism_depth[index_source] / dip_tan
                    hdist_d = L_sism_depth[index_source] / dip_tan
                    
                    azimuth = [0.,0.,0.,0.]
                    
                    lon_u = []
                    lat_u = []
                    lon_d = []
                    lat_d = []
                    if len(Lon[index_source]) < 9 :
                        azimuth = []
                        indexx = range(len(Lon[index_source]))
                        for i in indexx:
                            strike = abs(math.degrees(math.atan((
                            Lon[index_source][-1]-Lon[index_source][0])/(
                            Lat[index_source][-1]-Lat[index_source][0]))))
                            if strike > 180. :
                                strike = strike - 180.
                            azimuth.append((strike + 90.0) % 360)
                            
                        
                    else:
                        indexx = range(len(Lon[index_source]))
                        for i in indexx[4:-4]:
                            strike = abs(math.degrees(math.atan((
                            Lon[index_source][i+4]-Lon[index_source][i-4])/(
                            Lat[index_source][i+4]-Lat[index_source][i-4]))))
                            if strike > 180. :
                                strike = strike - 180.
                            azimuth.append((strike + 90.0) % 360)
                            
                        azimuth[0] = azimuth[4]
                        azimuth[1] = azimuth[4]
                        azimuth[2] = azimuth[4]
                        azimuth[3] = azimuth[4]
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                    
                    for i in indexx:
                        if (Lon[index_source][0])>(Lon[index_source][-1]):
                            x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                            y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
                            x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                            y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
                            lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
                            lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
                        if (Lon[index_source][0])<(Lon[index_source][-1]):
                            x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                            y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
                            x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                            y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
                            lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
                            lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
                        
                
                    
                    source_name_i = source_name[index_source].replace(Model+'_','')
                    
                    x, y = m(Lon[index_source], Lat[index_source]) 
                    #m.plot(x, y, 'D-', markersize=0.1, linewidth=0.2, color='k', markerfacecolor='b')
                    m.plot(x,y,linewidth=0.3,color='k')
                
                    poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                    poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                    
                    if len(Lon[index_source]) < 3 :
                        poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                        poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                    maps_utils.draw_screen_poly(poly_lons, poly_lats,  m ,'k' , 0.2, 0.05, 'k')
                if fault_type[index_source] == 'cf':
                    maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m ,'k' , 0.2, 0.2, 'k')
            
                
            m.drawcoastlines(linewidth=0.1)
            try:
                m.arcgisimage(service='World_Shaded_Relief', dpi = 400, alpha = 0.3, xpixels = 2000)
            except:
                m.fillcontinents(color='sienna',lake_color='w',alpha = 0.05)

            
            plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map.png',dpi=400)

            plt.close()

                     
            #'''        
            # map the NMS slip, the slip_rate, the seismic slip-rate and the Mmax
            if plot_sr_use == True:
                for MFD_type in MFD_type_list:
                    #print Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt'
                    slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt',
                                  dtype = [('U100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),
                                           ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t') 
                    fault_name_rep = list(map(lambda i : slip_rep_data[i][0], range(len(slip_rep_data))))
                    #print fault_name_rep
                    p_NMS = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))
                    
                    
                    m_nms = Basemap(projection='mill',
                              llcrnrlon=llcrnrlon, 
                              llcrnrlat=llcrnrlat, 
                              urcrnrlon=urcrnrlon, 
                              urcrnrlat=urcrnrlat,resolution='h')
                    
                    Lon_bg = []
                    Lat_bg = []
                
                    # manually defined  in the file Background geometry
                    geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
                    
                    column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
                    index_model = np.where(np.array(column_model) == Model)[0]
                    Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
                    Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
                    if len(Lon_bg) != 0 : #draw the background zone
                        maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_nms ,'g' , 0.02, 0.5, 'k')
                        
                    #for each fault  
                    
                    for index_source in range(nb_sources):
                        if fault_type[index_source] == 'sf':
                            dip_tan = math.tan(math.radians(Dip[index_source]))
                            
                            hdist_u = U_sism_depth[index_source] / dip_tan
                            hdist_d = L_sism_depth[index_source] / dip_tan
                            
                            azimuth = [0.,0.,0.,0.]
                            
                            lon_u = []
                            lat_u = []
                            lon_d = []
                            lat_d = []
                            if len(Lon[index_source]) < 9 :
                                azimuth = []
                                indexx = range(len(Lon[index_source]))
                                for i in indexx:
                                    strike = abs(math.degrees(math.atan((
                                    Lon[index_source][-1]-Lon[index_source][0])/(
                                    Lat[index_source][-1]-Lat[index_source][0]))))
                                    if strike > 180. :
                                        strike = strike - 180.
                                    azimuth.append((strike + 90.0) % 360)
                                    
                                
                            else:
                                indexx = range(len(Lon[index_source]))
                                for i in indexx[4:-4]:
                                    strike = abs(math.degrees(math.atan((
                                    Lon[index_source][i+4]-Lon[index_source][i-4])/(
                                    Lat[index_source][i+4]-Lat[index_source][i-4]))))
                                    if strike > 180. :
                                        strike = strike - 180.
                                    azimuth.append((strike + 90.0) % 360)
                                    
                                azimuth[0] = azimuth[4]
                                azimuth[1] = azimuth[4]
                                azimuth[2] = azimuth[4]
                                azimuth[3] = azimuth[4]
                                azimuth.append(azimuth[-4])
                                azimuth.append(azimuth[-4])
                                azimuth.append(azimuth[-4])
                                azimuth.append(azimuth[-4])
                            
                            for i in indexx:
                                if (Lon[index_source][0])>(Lon[index_source][-1]):
                                    x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                                    y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
                                    x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                                    y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
                                    lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
                                    lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
                                if (Lon[index_source][0])<(Lon[index_source][-1]):
                                    x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                                    y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
                                    x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                                    y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
                                    lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
                                    lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
                                
                        
                            
                            source_name_i = source_name[index_source].replace(Model+'_','')
                            #print source_name_i
                            index_fault = np.where(np.array(fault_name_rep)==source_name_i)[0][0]
                            NMS_i = p_NMS[index_fault]
    
                            x,y = m_nms(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                            if nb_on_maps == True :
                                plt.text(x,y,''+str(int(round(NMS_i)))+'',fontsize = 2)
                            
                            
                            
                            cmap = matplotlib.cm.get_cmap('rainbow')
                            if NMS_i >= 50.:
                                NMS_i = 50.
                            rgba = cmap(float(NMS_i*2.)/100.)
                            
                            x, y = m_nms(Lon[index_source], Lat[index_source]) 
                            
                            #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
                            
                            m_nms.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
                            
                        
                            poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                            poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                            
                            if len(Lon[index_source]) < 3 :
                                poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                                poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                            maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_nms ,rgba , 0.5, 1., rgba)
                        if fault_type[index_source] == 'cf':
                            source_name_i = source_name[index_source].replace(Model+'_','')
                            #print source_name_i
                            index_fault = np.where(np.array(fault_name_rep)==source_name_i)[0][0]
                            NMS_i = p_NMS[index_fault]
    
                            x,y = m_nms(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                            if nb_on_maps == True :
                                plt.text(x,y,''+str(int(round(NMS_i)))+'',fontsize = 2)
                            
                            
                            
                            cmap = matplotlib.cm.get_cmap('rainbow')
                            if NMS_i >= 50.:
                                NMS_i = 50.
                            rgba = cmap(float(NMS_i*2.)/100.)
                            maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_nms ,rgba , 0.5, 1., rgba)
                    
                    m_nms.drawcoastlines(linewidth=0.2)
                    m_nms.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
                    plt.title('NMS_'+MFD_type+'_'+str(scenario_set))
                    plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_NMS_'+MFD_type+'_'+str(scenario_set)+'.png',dpi=300,transparent=True)
        
                    plt.close()         
                        
                    
                    a = np.array([[0,1]])
                    pl.figure(figsize=(9, 1.5))
                    img = pl.imshow(a, cmap="rainbow")
                    pl.gca().set_visible(False)
                    cax = pl.axes([0.1, 0.2, 0.4, 0.3])
                    cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
                    cb.set_ticklabels(['0 %', '10 %', '20 %', '30 %', '40 %', '> 50 %'])
                    pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_NMS_ratio.png',dpi=180)
                    pl.close()
                
                
                
                     
            '''##################################################################        
            # map the slip_rate
            ###################################################################'''

            mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
                          dtype = [('U100'),('U100'),('U1000'),('f8'),('f8')], delimiter = '\t') 
            model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
            set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
            fault_name_mean_param =list( map(lambda i : mean_param[i][2], range(len(mean_param))))
#            print(set_mean_param[0],type(set_mean_param[0]))
#            print(fault_name_mean_param[0],type(fault_name_mean_param[0]))
            sr_mean =list( map(lambda i : mean_param[i][3], range(len(mean_param))))
            Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
            
            index_model = np.where(np.array(model_mean_param)==Model)[0]
            set_mean_param = np.take(set_mean_param,index_model)
            fault_name_mean_param = np.take(fault_name_mean_param,index_model)
            sr_mean = np.take(sr_mean,index_model)
            Mmax_mean = np.take(Mmax_mean,index_model)
            
            index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
            fault_name_mean_param = np.take(fault_name_mean_param,index_set)
            sr_mean = np.take(sr_mean,index_set)
            Mmax_mean = np.take(Mmax_mean,index_set)
            
            
            m_sr = Basemap(projection='mill',
                      llcrnrlon=llcrnrlon, 
                      llcrnrlat=llcrnrlat, 
                      urcrnrlon=urcrnrlon, 
                      urcrnrlat=urcrnrlat,resolution='h')
            
            
            
            
            Lon_bg = []
            Lat_bg = []
        
            # manually defined  in the file Background geometry
            geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
            
            column_model =list( map(lambda i : geom_bg[i][0],range(len(geom_bg))))
            index_model = np.where(np.array(column_model) == Model)[0]
            Lon_bg =list( map(lambda i : geom_bg[i][1],index_model))
            Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
            if len(Lon_bg) != 0 : #draw the background zone
                maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_sr ,'g' , 0.02, 0.05, 'k')
                
            #for each fault  
            
            for index_source in range(nb_sources):
                if fault_type[index_source] == 'sf':
                    dip_tan = math.tan(math.radians(Dip[index_source]))
                    
                    hdist_u = U_sism_depth[index_source] / dip_tan
                    hdist_d = L_sism_depth[index_source] / dip_tan
                    
                    azimuth = [0.,0.,0.,0.]
                    
                    lon_u = []
                    lat_u = []
                    lon_d = []
                    lat_d = []
                    if len(Lon[index_source]) < 9 :
                        azimuth = []
                        indexx = range(len(Lon[index_source]))
                        for i in indexx:
                            strike = abs(math.degrees(math.atan((
                            Lon[index_source][-1]-Lon[index_source][0])/(
                            Lat[index_source][-1]-Lat[index_source][0]))))
                            if strike > 180. :
                                strike = strike - 180.
                            azimuth.append((strike + 90.0) % 360)
                            
                        
                    else:
                        indexx = range(len(Lon[index_source]))
                        for i in indexx[4:-4]:
                            strike = abs(math.degrees(math.atan((
                            Lon[index_source][i+4]-Lon[index_source][i-4])/(
                            Lat[index_source][i+4]-Lat[index_source][i-4]))))
                            if strike > 180. :
                                strike = strike - 180.
                            azimuth.append((strike + 90.0) % 360)
                            
                        azimuth[0] = azimuth[4]
                        azimuth[1] = azimuth[4]
                        azimuth[2] = azimuth[4]
                        azimuth[3] = azimuth[4]
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                    
                    for i in indexx:
                        if (Lon[index_source][0])>(Lon[index_source][-1]):
                            x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                            y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
                            x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                            y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
                            lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
                            lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
                        if (Lon[index_source][0])<(Lon[index_source][-1]):
                            x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                            y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
                            x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                            y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
                            lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
                            lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
                        
                
                    
                    source_name_i = source_name[index_source].replace(Model+'_','')
#                    print(source_name_i)
#                    print(fault_name_mean_param)
                    index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
                    sr = sr_mean[index_fault]
                    
                    x,y = m_sr(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                    if nb_on_maps == True :
                        plt.text(x,y,''+str(round(sr,2))+'',fontsize = 2)
                    
                    cmap = matplotlib.cm.get_cmap('rainbow')
                    rgba = cmap(float(sr)/max(sr_mean))
                    
                    x, y = m_sr(Lon[index_source], Lat[index_source]) 
                    
                    #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
                    
                    m_sr.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
                    
                
                    poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                    poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                    
                    if len(Lon[index_source]) < 3 :
                        poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                        poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                    maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_sr ,rgba , 0.5, 1., rgba)
                if fault_type[index_source] == 'cf':
                    source_name_i = source_name[index_source].replace(Model+'_','')
                    index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
                    sr = sr_mean[index_fault]
                    
                    x,y = m_sr(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                    if nb_on_maps == True :
                        plt.text(x,y,''+str(round(sr,2))+'',fontsize = 2)
                    
                    cmap = matplotlib.cm.get_cmap('rainbow')
                    rgba = cmap(float(sr)/max(sr_mean))
                    maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_sr ,rgba , 0.5, 1., rgba)
            
            m_sr.drawcoastlines(linewidth=0.2)
            m_sr.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
            plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_sr_'+str(scenario_set)+'.png',dpi=300)

            plt.close()         
                
            
            a = np.array([[0,1]])
            pl.figure(figsize=(9, 1.5))
            img = pl.imshow(a, cmap="rainbow")
            pl.gca().set_visible(False)
            cax = pl.axes([0.1, 0.2, 0.4, 0.3])
            cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
            cb.set_ticklabels([str(round(min(sr_mean),2)),
                               str(round(min(sr_mean)+ 0.2*(max(sr_mean)-min(sr_mean)),2)),
                             str(round(min(sr_mean)+ 0.4*(max(sr_mean)-min(sr_mean)),2)),
                             str(round(min(sr_mean)+ 0.6*(max(sr_mean)-min(sr_mean)),2)),
                             str(round(min(sr_mean)+ 0.8*(max(sr_mean)-min(sr_mean)),2)),
                             str(round(max(sr_mean),1))])
            pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_sliprate.png',dpi=180)
            pl.close()
                
                
                
                     
            #'''        
            # map the seismic slip-rate
            if plot_sr_use == True:
                for MFD_type in MFD_type_list:
                    slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt',
                              dtype = [('U100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),
                                       ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t') 
                    fault_name_rep = list(map(lambda i : slip_rep_data[i][0], range(len(slip_rep_data))))
                    p_NMS = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))
                    
                    
                    mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
                                  dtype = [('U100'),('U100'),('U100'),('f8'),('f8')], delimiter = '\t') 
                    model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
                    set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
                    fault_name_mean_param = list(map(lambda i : mean_param[i][2], range(len(mean_param))))
                    sr_mean = list(map(lambda i : mean_param[i][3], range(len(mean_param))))
                    Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
                    
                    index_model = np.where(np.array(model_mean_param)==Model)[0]
                    set_mean_param = np.take(set_mean_param,index_model)
                    fault_name_mean_param = np.take(fault_name_mean_param,index_model)
                    sr_mean = np.take(sr_mean,index_model)
                    Mmax_mean = np.take(Mmax_mean,index_model)
                    
                    index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
                    fault_name_mean_param = np.take(fault_name_mean_param,index_set)
                    sr_mean = np.take(sr_mean,index_set)
                    Mmax_mean = np.take(Mmax_mean,index_set)
                    
                    
                    
                    m_sr_seismic = Basemap(projection='mill',
                              llcrnrlon=llcrnrlon, 
                              llcrnrlat=llcrnrlat, 
                              urcrnrlon=urcrnrlon, 
                              urcrnrlat=urcrnrlat,resolution='h')
                    
                    
                    
                    
                    Lon_bg = []
                    Lat_bg = []
                
                    # manually defined  in the file Background geometry
                    geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
                    
                    column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
                    index_model = np.where(np.array(column_model) == Model)[0]
                    Lon_bg =list( map(lambda i : geom_bg[i][1],index_model))
                    Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
                    if len(Lon_bg) != 0 : #draw the background zone
                        maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_sr_seismic ,'g' , 0.02, 0.05, 'k')
                        
                    #for each fault  
                    
                    for index_source in range(nb_sources):
                        if fault_type[index_source] == 'sf':
                            dip_tan = math.tan(math.radians(Dip[index_source]))
                            
                            hdist_u = U_sism_depth[index_source] / dip_tan
                            hdist_d = L_sism_depth[index_source] / dip_tan
                            
                            azimuth = [0.,0.,0.,0.]
                            
                            lon_u = []
                            lat_u = []
                            lon_d = []
                            lat_d = []
                            if len(Lon[index_source]) < 9 :
                                azimuth = []
                                indexx = range(len(Lon[index_source]))
                                for i in indexx:
                                    strike = abs(math.degrees(math.atan((
                                    Lon[index_source][-1]-Lon[index_source][0])/(
                                    Lat[index_source][-1]-Lat[index_source][0]))))
                                    if strike > 180. :
                                        strike = strike - 180.
                                    azimuth.append((strike + 90.0) % 360)
                                    
                                
                            else:
                                indexx = range(len(Lon[index_source]))
                                for i in indexx[4:-4]:
                                    strike = abs(math.degrees(math.atan((
                                    Lon[index_source][i+4]-Lon[index_source][i-4])/(
                                    Lat[index_source][i+4]-Lat[index_source][i-4]))))
                                    if strike > 180. :
                                        strike = strike - 180.
                                    azimuth.append((strike + 90.0) % 360)
                                    
                                azimuth[0] = azimuth[4]
                                azimuth[1] = azimuth[4]
                                azimuth[2] = azimuth[4]
                                azimuth[3] = azimuth[4]
                                azimuth.append(azimuth[-4])
                                azimuth.append(azimuth[-4])
                                azimuth.append(azimuth[-4])
                                azimuth.append(azimuth[-4])
                            
                            for i in indexx:
                                if (Lon[index_source][0])>(Lon[index_source][-1]):
                                    x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                                    y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
                                    x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                                    y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
                                    lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
                                    lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
                                if (Lon[index_source][0])<(Lon[index_source][-1]):
                                    x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                                    y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
                                    x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                                    y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
                                    lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
                                    lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                                    lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
                                
                        
                            
                            source_name_i = source_name[index_source].replace(Model+'_','')
                            #print source_name_i
                            index_fault = np.where(np.array(fault_name_rep)==source_name_i)[0][0]
                            NMS_i = p_NMS[index_fault]
                            
                            index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
                            sr = sr_mean[index_fault]
                            
                            sr_seismic = sr*(1. - float(NMS_i)/100.)
                            x,y = m_sr_seismic(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                            if nb_on_maps == True :
                                plt.text(x,y,''+str(round(sr_seismic,2))+'',fontsize = 2)
                            
                            
                            cmap = matplotlib.cm.get_cmap('rainbow')
                            rgba = cmap(sr_seismic/max(sr_mean))
                            
                            x, y = m_sr_seismic(Lon[index_source], Lat[index_source]) 
                            
                            #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
                            
                            m_sr_seismic.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
                            
                        
                            poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                            poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                            
                            if len(Lon[index_source]) < 3 :
                                poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                                poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                            maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_sr_seismic ,rgba , 0.5, 1., rgba)
                        if fault_type[index_source] == 'cf':
                            source_name_i = source_name[index_source].replace(Model+'_','')
                            #print source_name_i
                            index_fault = np.where(np.array(fault_name_rep)==source_name_i)[0][0]
                            NMS_i = p_NMS[index_fault]
                            
                            index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
                            sr = sr_mean[index_fault]
                            
                            sr_seismic = sr*(1. - float(NMS_i)/100.)
                            x,y = m_sr_seismic(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                            if nb_on_maps == True :
                                plt.text(x,y,''+str(round(sr_seismic,2))+'',fontsize = 2)
                            
                            
                            cmap = matplotlib.cm.get_cmap('rainbow')
                            rgba = cmap(sr_seismic/max(sr_mean))
                            maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_sr_seismic ,rgba , 0.5, 1., rgba)
                    
                    m_sr_seismic.drawcoastlines(linewidth=0.2)
                    m_sr_seismic.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
                    plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_sr_seismic_'+MFD_type+'_'+str(scenario_set)+'.png',dpi=300)
        
                    plt.close()         
                        
                    
                    a = np.array([[0,1]])
                    pl.figure(figsize=(9, 1.5))
                    img = pl.imshow(a, cmap="rainbow")
                    pl.gca().set_visible(False)
                    cax = pl.axes([0.1, 0.2, 0.4, 0.3])
                    cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
                    cb.set_ticklabels([str(round(min(sr_mean),2)),
                                       str(round(min(sr_mean)+ 0.2*(max(sr_mean)-min(sr_mean)),2)),
                                     str(round(min(sr_mean)+ 0.4*(max(sr_mean)-min(sr_mean)),2)),
                                     str(round(min(sr_mean)+ 0.6*(max(sr_mean)-min(sr_mean)),2)),
                                     str(round(min(sr_mean)+ 0.8*(max(sr_mean)-min(sr_mean)),2)),
                                     str(round(max(sr_mean),2))])
                    pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_sliprate_seismic.png',dpi=180)
                    pl.close()
                
                
                
                     
            #'''        
            # map the Mmax
            mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
                          dtype = [('U100'),('U100'),('U100'),('f8'),('f8')], delimiter = '\t') 
            model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
            set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
            fault_name_mean_param =list( map(lambda i : mean_param[i][2], range(len(mean_param))))
            sr_mean = list(map(lambda i : mean_param[i][3], range(len(mean_param))))
            Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
            
            index_model = np.where(np.array(model_mean_param)==Model)[0]
            set_mean_param = np.take(set_mean_param,index_model)
            fault_name_mean_param = np.take(fault_name_mean_param,index_model)
            sr_mean = np.take(sr_mean,index_model)
            Mmax_mean = np.take(Mmax_mean,index_model)
            
            index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
            fault_name_mean_param = np.take(fault_name_mean_param,index_set)
            sr_mean = np.take(sr_mean,index_set)
            Mmax_mean = np.take(Mmax_mean,index_set)
            
            
            m_mmax = Basemap(projection='mill',
                      llcrnrlon=llcrnrlon, 
                      llcrnrlat=llcrnrlat, 
                      urcrnrlon=urcrnrlon, 
                      urcrnrlat=urcrnrlat,resolution='h')
            
            
            
            
            Lon_bg = []
            Lat_bg = []
        
            # manually defined  in the file Background geometry
            geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
            
            column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
            index_model = np.where(np.array(column_model) == Model)[0]
            Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
            Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
            if len(Lon_bg) != 0 : #draw the background zone
                maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_mmax ,'g' , 0.02, 0.05, 'k')
                
            #for each fault  
            for index_source in range(nb_sources):
                if fault_type[index_source] == 'sf':
                    dip_tan = math.tan(math.radians(Dip[index_source]))
                    
                    hdist_u = U_sism_depth[index_source] / dip_tan
                    hdist_d = L_sism_depth[index_source] / dip_tan
                    
                    azimuth = [0.,0.,0.,0.]
                    
                    lon_u = []
                    lat_u = []
                    lon_d = []
                    lat_d = []
                    if len(Lon[index_source]) < 9 :
                        azimuth = []
                        indexx = range(len(Lon[index_source]))
                        for i in indexx:
                            strike = abs(math.degrees(math.atan((
                            Lon[index_source][-1]-Lon[index_source][0])/(
                            Lat[index_source][-1]-Lat[index_source][0]))))
                            if strike > 180. :
                                strike = strike - 180.
                            azimuth.append((strike + 90.0) % 360)
                            
                        
                    else:
                        indexx = range(len(Lon[index_source]))
                        for i in indexx[4:-4]:
                            strike = abs(math.degrees(math.atan((
                            Lon[index_source][i+4]-Lon[index_source][i-4])/(
                            Lat[index_source][i+4]-Lat[index_source][i-4]))))
                            if strike > 180. :
                                strike = strike - 180.
                            azimuth.append((strike + 90.0) % 360)
                            
                        azimuth[0] = azimuth[4]
                        azimuth[1] = azimuth[4]
                        azimuth[2] = azimuth[4]
                        azimuth[3] = azimuth[4]
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                        azimuth.append(azimuth[-4])
                    
                    for i in indexx:
                        if (Lon[index_source][0])>(Lon[index_source][-1]):
                            x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                            y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
                            x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                            y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
                            lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
                            lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
                        if (Lon[index_source][0])<(Lon[index_source][-1]):
                            x_u = hdist_u * math.sin(math.radians(azimuth[i]))
                            y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
                            x_d = hdist_d * math.sin(math.radians(azimuth[i]))
                            y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
                            lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
                            lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
                            lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
                        
                
                    
                    source_name_i = source_name[index_source].replace(Model+'_','')
                    index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
                    Mmax = Mmax_mean[index_fault]
                    
                    x,y = m_mmax(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                    if nb_on_maps == True :
                        plt.text(x,y,''+str(Mmax)+'',fontsize = 2)
                    
                    
                    cmap = matplotlib.cm.get_cmap('rainbow')
                    rgba = cmap((float(Mmax)-min(Mmax_mean))/(max(Mmax_mean)-min(Mmax_mean)))
                    
                    x, y = m_mmax(Lon[index_source], Lat[index_source]) 
                    
                    #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
                    
                    m_mmax.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
                    
                
                    poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                    poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                    
                    if len(Lon[index_source]) < 3 :
                        poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
                        poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
                    maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_mmax ,rgba , 0.5, 1., rgba)
                if fault_type[index_source] == 'cf':
                    source_name_i = source_name[index_source].replace(Model+'_','')
                    index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
                    Mmax = Mmax_mean[index_fault]
                    
                    x,y = m_mmax(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
                    if nb_on_maps == True :
                        plt.text(x,y,''+str(Mmax)+'',fontsize = 2)
                    
                    
                    cmap = matplotlib.cm.get_cmap('rainbow')
                    rgba = cmap((float(Mmax)-min(Mmax_mean))/(max(Mmax_mean)-min(Mmax_mean)))
                    maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_mmax ,rgba , 0.5, 1., rgba)
            
            m_mmax.drawcoastlines(linewidth=0.2)
            m_mmax.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
            plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_Mmax_'+str(scenario_set)+'.png',dpi=300)

            plt.close()         
                
            
            a = np.array([[0,1]])
            pl.figure(figsize=(9, 1.5))
            img = pl.imshow(a, cmap="rainbow")
            pl.gca().set_visible(False)
            cax = pl.axes([0.1, 0.2, 0.4, 0.3])
            cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
            cb.set_ticklabels([str(min(Mmax_mean)),
                               str(min(Mmax_mean)+ 0.2*(max(Mmax_mean)-min(Mmax_mean))),
                             str(min(Mmax_mean)+ 0.4*(max(Mmax_mean)-min(Mmax_mean))),
                             str(min(Mmax_mean)+ 0.6*(max(Mmax_mean)-min(Mmax_mean))),
                             str(min(Mmax_mean)+ 0.8*(max(Mmax_mean)-min(Mmax_mean))),
                             str(max(Mmax_mean))])
            pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_Mmax.png',dpi=180)
            pl.close()
                    

#def map_faults_2(Run_name,Model_name,list_FtF,file_fault_geom,Distance_criteria,
#               llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat): # without the XML
#
#    if not os.path.exists('../'+Run_name+'/prerun_analysis/'+Model_name+'/FtF_rutpures_'+str(Distance_criteria)):
#        os.makedirs('../'+Run_name+'/prerun_analysis/'+Model_name+'/FtF_rutpures_'+str(Distance_criteria))
#    file_names = []
#    InfosZonage = np.genfromtxt(file_fault_geom,dtype=[('U100'),('U100'),('f8'),('f8')],skip_header = 1)
#    Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
#    index_model = np.where(np.array(Column_model_name) == Model_name)
#    Column_Fault_name = list(map(lambda i : InfosZonage[i][1],index_model[0]))
#
#    Longitudes = list(map(lambda i : InfosZonage[i][2],index_model[0]))
#    Latitudes = list(map(lambda i : InfosZonage[i][3],index_model[0]))
#
#    print('creating the maps of each FtF')
#    print()
#    index_FtF = 1
#    for FtF in list_FtF:
#
#        plt.figure(figsize=(12,6))
#        m = Basemap(projection='mill',
#                      llcrnrlon=llcrnrlon,
#                      llcrnrlat=llcrnrlat,
#                      urcrnrlon=urcrnrlon,
#                      urcrnrlat=urcrnrlat,resolution='l')
#        m.drawcoastlines(linewidth=0.1)
#        #m.fillcontinents(color='grey',lake_color='w',alpha = 0.1, zorder=1)
#        x, y = m(Longitudes,Latitudes)
#        m.scatter(x ,y ,s=0.5,marker="o",color='k',alpha=0.2, zorder=1)
#        for fault in FtF:
#            index = np.where(np.array(Column_Fault_name) == fault)
#            Longitudes_fault = np.take(Longitudes,index)[0]
#            Latitudes_fault = np.take(Latitudes,index)[0]
#
#            x_fault, y_fault = m(Longitudes_fault, Latitudes_fault)
#
#            m.plot(x_fault, y_fault,  linewidth=1.5,color='r', zorder=2)
#
#        title = str()
#        i = 0
#        for fault in FtF :
#            title += ' ' + str(fault)
#            i += 1
#            if i == 6 :
#                title += '\n'
#                i = 0
#        plt.title(title)
#        plt.savefig('../'+Run_name+'/prerun_analysis/'+Model_name+'/FtF_rutpures_'+str(Distance_criteria)+'/'+str(index_FtF)+'.png',dpi=180)
#        file_names.append('../'+Run_name+'/prerun_analysis/'+Model_name+'/FtF_rutpures_'+str(Distance_criteria)+'/'+str(index_FtF)+'.png')
#        print(index_FtF, 'out of ', len(list_FtF))
#        index_FtF +=1
        

   
