#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""


import os,sys,inspect,copy
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
import utils.Geometry_scenario as Geometry_scenario
import maps.maps_utils as maps_utils
import maps.geom as geom
from sources.background import bg

import geojson
from geojson import LineString, Feature, FeatureCollection, dump


import utils.read_input as read_input


def map_faults(Run_name,Model_list,scenarios_names_list,
               ScL_complet_list, BG_hyp_list,
               sample_list,b_value_list,MFD_type_list,
               llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,File_bg,
               FileName_Prop,File_geom,plot_sr_use,visual_FtF,sub_area_file):
    nb_on_maps = False
    
    # use basemap option (create png figs)
    use_basemap = False
    
    #exctracting all the scenario set information
    available_sets = read_input.extract_sc_input('input/'+Run_name+'/ruptures.txt')
    
    
    for Model in Model_list :
    
        if not ".geojson" in FileName_Prop:
            # new method
            #Extraction of the faults and scenarios present in the model from the text file
            Prop = np.genfromtxt(FileName_Prop,
                                       dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
                                              ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
            Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
            Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
            index_model = np.where(np.array(Column_model_name) == Model)[0]
            Prop = np.take(Prop,index_model)
            faults_names = np.array(Column_fault_name[index_model[0]:index_model[-1]+1])
            faults_names = list(faults_names)
        
        else : #it's a geojson file
            with open(FileName_Prop) as f:
                gj = geojson.load(f)
            faults = gj['features']
            faults_names = []
            for fi in range(len(faults)):
                if faults[fi]['properties']['model'] == Model :
                    faults_names.append(str(faults[fi]['properties']['si']))
        
        geom_scenar = Geometry_scenario.Geom_scenar(faults_names,File_geom,Model)
        faults_lon = geom_scenar.faults_lon
        faults_lat = geom_scenar.faults_lat
        
        Column_Fault_name, Depths = geom.FaultGeometry(Model,File_geom)  #extract the geometries from the geometry file
        
        Lon_bg, Lat_bg  = bg.geom(Model,File_bg)
        
        # fault geom contains the geometry of the fault (trace and polygon)
        fault_geom = {}
        index_fault = 0
        for Fault_name in faults_names:
            lons = faults_lon[index_fault]
            lats = faults_lat[index_fault]
            # extract depth
            i_d = np.where(np.array(Column_Fault_name) == Fault_name)
            depth = list(map(lambda i : Depths[i],i_d[0]))
            
            dip, oriented, upper_sismo_depth, lower_sismo_depth = geom.FaultProperties(FileName_Prop,Fault_name,Model)
            
            if depth[0] == 'sf' :
                trace_lon, trace_lat, plot_trace, polygon = geom.get_sf_polygon(lons, lats,
                 upper_sismo_depth,
                 lower_sismo_depth,
                 dip,
                 oriented)
                 
            else : # it's a complexe fault
                trace_lon,trace_lat = [], []
                plot_trace = False
                lon_top, lon_bottom, lat_top, lat_bottom = [], [], [], []
                for lon_i, lat_i, d_i in zip(lons,lats,depth):
                    if d_i == min(depth):
                        lon_top.append(lon_i)
                        lat_top.append(lat_i)
                    if d_i == max(depth):
                        lon_bottom.append(lon_i)
                        lat_bottom.append(lat_i)
                        
                poly_lons = np.concatenate([lon_top,np.array(list(reversed(lon_bottom)))])
                poly_lats = np.concatenate([lat_top,np.array(list(reversed(lat_bottom)))])
                polygon = [poly_lons,poly_lats]
             
            fault_geom.update({index_fault:{'name':Fault_name,
            'trace_lon':trace_lon,
            'trace_lat':trace_lat,
            'plot_trace':plot_trace,
            'polygon':polygon}})
            
            index_fault += 1
            
        
        '''########################
        Print the faults activated for each scenario
        #########################'''
        if visual_FtF == True and use_basemap == True:
            for scenario_set in scenarios_names_list:
                file_names = []
                if not os.path.exists(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set):
                    os.makedirs(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set)
                    

                # extract the Mmax of the faults and the scenarios
#                log_Mmax_file = (str(Run_name)  + '/' + str(Model) + '/' + 'bg_' + str(BG_hyp_list[0]) + '/'
#                                + str(ScL_complet_list[0]) + '/sc_' +  str(scenario_set) + '/'
#                                + str(b_value_list[0]) + '/' + 'MFD_'+ str(MFD_type_list[0])
#                                +  '/Log/Mmax_sample_1.txt')
                log_Mmax_file = (str(Run_name)  + '/' + str(Model) +'/Log/Mmax_sample_'+str(ScL_complet_list[0])+'_'+scenario_set+'_1.txt')
                                                
                sources_names,sources_Mmax,sources_Lengths,sources_Areas = Read_file.read_file_Mmax_log(log_Mmax_file) #read the log of Mmax associated with the model
                rupture_set = available_sets['sc_'+scenario_set]

                m = Basemap(projection='mill',
                              llcrnrlon=llcrnrlon,
                              llcrnrlat=llcrnrlat,
                              urcrnrlon=urcrnrlon,
                              urcrnrlat=urcrnrlat,resolution='i')
                    
                
                index_rupt = 0
                for rupt in rupture_set:
                    if len(Lon_bg) != 0 : #draw the background zone
                        maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.2, 0.5, 'k')
                    m1 = copy.copy(m)
                    fault_colors = []
                    for fault in faults_names :
                        if fault in rupt:
                            fault_colors.append('r')
                        else :
                            fault_colors.append('k')
                            
                    title = str()
                    i = 0
                    for fault in rupt :
                        title += ' ' + str(fault)
                        i += 1
                        if i == 6 :
                            title += '\n'
                            i = 0
                    Mmax = sources_Mmax[len(faults_names)+index_rupt]
                    figpath = str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set+'/'+str(index_rupt)+'.png'
                    maps_utils.make_fault_map(m1,fault_geom,fault_colors,figpath,title,Mmax=Mmax)
                    index_rupt +=1
                    
                

        '''########################
        Print the map of the model
        #########################'''
        if use_basemap == True:
            m = Basemap(projection='mill',
                          llcrnrlon=llcrnrlon,
                          llcrnrlat=llcrnrlat,
                          urcrnrlon=urcrnrlon,
                          urcrnrlat=urcrnrlat,resolution='h')
                
            
            if len(Lon_bg) != 0 : #draw the background zone
                maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.2, 0.5, 'k')
            fault_colors = ['k' for fault in faults_names]
                    
            title = 'Map of the model : '+Model
            figpath = str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map.png'
                    

            m.drawcoastlines(linewidth=0.1)

            #draw the sub_areas
            if os.path.exists(sub_area_file):
                read_sub_area_file = open(sub_area_file,'rU')
                lines_sub_area = read_sub_area_file.readlines()
                sub_area_names = []
                sub_area_coord = []
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
                    
            maps_utils.make_fault_map(m,fault_geom,fault_colors,figpath,title,dpi=400,use_arcgis=True)
                
        
        '''########################
        Print the sliprate map
        #########################'''
        mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
                      dtype = [('U100'),('U100'),('U1000'),('f8'),('f8')], delimiter = '\t')
        model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
        set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
        fault_name_mean_param =list( map(lambda i : mean_param[i][2], range(len(mean_param))))
        sr_mean =list( map(lambda i : mean_param[i][3], range(len(mean_param))))
        Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))

        index_model = np.where(np.array(model_mean_param)==Model)[0]
        set_mean_param_model = np.take(set_mean_param,index_model)
        fault_name_mean_param_model = np.take(fault_name_mean_param,index_model)
        sr_mean_model = np.take(sr_mean,index_model)
        Mmax_mean_model = np.take(Mmax_mean,index_model)

#        index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
#        fault_name_mean_param = np.take(fault_name_mean_param,index_set)
#        sr_mean = np.take(sr_mean,index_set)
#        Mmax_mean = np.take(Mmax_mean,index_set)

        if use_basemap == True:
            m = Basemap(projection='mill',
                       llcrnrlon=llcrnrlon,
                       llcrnrlat=llcrnrlat,
                       urcrnrlon=urcrnrlon,
                       urcrnrlat=urcrnrlat,resolution='i')

            title = 'Slip rate : '+ Model
            figpath = str(Run_name) +'/analysis/figures/FtF/'+Model+'/'+'map_sliprate.png'

#        if len(Lon_bg) != 0 : #draw the background zone
#            maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.1, 0.5, 'k')
         
        fault_colors = []
        for fault in faults_names:
            index_fault = np.where(np.array(fault_name_mean_param_model)==fault)[0][0]
            sr = sr_mean_model[index_fault]

            cmap = matplotlib.cm.get_cmap('rainbow')
            rgba = cmap(float(sr)/max(sr_mean))
            fault_colors.append(rgba)

        if use_basemap == True:
            maps_utils.make_fault_map(m,fault_geom,fault_colors,figpath,title,dpi=180,use_arcgis=False)
                  
        '''########################
        Print the  Mmax map
        #########################'''
        if use_basemap == True:
            m = Basemap(projection='mill',
                          llcrnrlon=llcrnrlon,
                          llcrnrlat=llcrnrlat,
                          urcrnrlon=urcrnrlon,
                          urcrnrlat=urcrnrlat,resolution='i')
            
        for scenario_set in scenarios_names_list:
            if use_basemap == True:
                m1 = copy.copy(m)
                title = 'Slip rate : '+ Model + ' ' +scenario_set
                figpath = str(Run_name) +'/analysis/figures/FtF/'+Model+'/'+'map_Mmax_'+scenario_set+'.png'
                
            index_set = np.where(np.array(set_mean_param_model)==scenario_set)[0]
            fault_name_mean_param_set = np.take(fault_name_mean_param_model,index_set)
            sr_mean_set = np.take(sr_mean_model,index_set)
            Mmax_mean_set = np.take(Mmax_mean_model,index_set)
        

#            if len(Lon_bg) != 0 : #draw the background zone
#                maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.1, 0.5, 'k')
             
            fault_colors = []
            for fault in faults_names:
                index_fault = np.where(np.array(fault_name_mean_param_set)==fault)[0][0]
                Mmax = Mmax_mean_set[index_fault]

                cmap = matplotlib.cm.get_cmap('rainbow')
                rgba = cmap((float(Mmax)-min(Mmax_mean))/(max(Mmax_mean)-min(Mmax_mean)))
                fault_colors.append(rgba)
            if use_basemap == True:
                maps_utils.make_fault_map(m1,fault_geom
                ,fault_colors,figpath
                ,title,dpi=180,use_arcgis=False)
              
        '''########################
        Print the NMS
        #########################'''
        if use_basemap == True:
            m = Basemap(projection='mill',
                      llcrnrlon=llcrnrlon,
                      llcrnrlat=llcrnrlat,
                      urcrnrlon=urcrnrlon,
                      urcrnrlat=urcrnrlat,resolution='i')
            
        for scenario_set in scenarios_names_list:

            index_set = np.where(np.array(set_mean_param_model)==scenario_set)[0]
            fault_name_mean_param_set = np.take(fault_name_mean_param_model,index_set)
            Mmax_mean_set = np.take(Mmax_mean_model,index_set)
            
            for MFD_type in MFD_type_list:
                if use_basemap == True:
                    title = 'NMS : '+ Model +' '+ MFD_type +' '+ scenario_set
                    figpath = str(Run_name) +'/analysis/figures/FtF/'+Model+'/'+'map_NMS_'+MFD_type+'_'+str(scenario_set)+'.png'
                slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt',
                              dtype = [('U100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),
                                       ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t')
                fault_name_rep = list(map(lambda i : slip_rep_data[i][0], range(len(slip_rep_data))))
                p_NMS = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))
            
                if use_basemap == True:
                    m1 = copy.copy(m)

                if use_basemap == True:
                    fault_colors = []
                    for fault in faults_names:
                        index_fault = np.where(np.array(fault_name_rep)==fault)[0][0]
                        NMS_i = p_NMS[index_fault]

                        cmap = matplotlib.cm.get_cmap('rainbow')
                        if NMS_i >= 50.:
                            NMS_i = 50.
                        rgba = cmap(float(NMS_i*2.)/100.)
                        fault_colors.append(rgba)
                        maps_utils.make_fault_map(m,fault_geom,
                        fault_colors,figpath,
                        title,dpi=180,use_arcgis=False)

         

                '''########################
                Print the seismic sliprate
                #########################'''

#                slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt',
#                          dtype = [('U100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),
#                                   ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t')
#                fault_name_rep = list(map(lambda i : slip_rep_data[i][0], range(len(slip_rep_data))))
#                p_NMS = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))
                if use_basemap == True:
                    m2 = copy.copy(m)
                    title = 'Seismic slip rate : '+ Model +' '+ MFD_type +' '+ scenario_set
                    figpath = str(Run_name) +'/analysis/figures/FtF/'+Model+'/'+'map_seismicsliprate_'+MFD_type+'_'+str(scenario_set)+'.png'
#                slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt',
#                           dtype = [('U100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),
#                                    ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t')
#                fault_name_rep = list(map(lambda i : slip_rep_data[i][0], range(len(slip_rep_data))))
#                p_NMS = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))

#                if len(Lon_bg) != 0 : #draw the background zone
#                    maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m ,'g' , 0.1, 0.5, 'k')
                    
                    fault_colors = []
                    for fault in faults_names:
                        index_fault = np.where(np.array(fault_name_rep)==fault)[0][0]
                        NMS_i = p_NMS[index_fault]

                        index_fault = np.where(np.array(fault_name_mean_param)==fault)[0][0]
                        sr = sr_mean[index_fault]

                        sr_seismic_i = sr*(1. - float(NMS_i)/100.)
                        
                        cmap = matplotlib.cm.get_cmap('rainbow')
                        rgba = cmap(sr_seismic_i/max(sr_mean))
                        fault_colors.append(rgba)
                         
                        maps_utils.make_fault_map(m2,fault_geom,
                        fault_colors,figpath,
                        title,dpi=180,use_arcgis=False)
                  


                '''########################
                Build the geojson file
                #########################'''
                features = []
                id = 0
                for fault in faults_names:
                    index_fault = np.where(np.array(fault_name_rep)==fault)[0][0]
                    NMS_i = p_NMS[index_fault]

                    index_fault = np.where(np.array(fault_name_mean_param)==fault)[0][0]
                    sr = sr_mean[index_fault]

                    sr_seismic_i = sr*(1. - float(NMS_i)/100.)


                    index_fault = np.where(np.array(fault_name_mean_param_set)==fault)[0][0]
                    Mmax = Mmax_mean_set[index_fault]

                    trace = []
                    for lon_i,lat_i in zip(fault_geom[id]["trace_lon"],fault_geom[id]["trace_lat"]):
                        trace.append([lon_i,lat_i])
                    trace = LineString(trace)

                    features.append(Feature(geometry=trace, properties={"id": id,
                     "name": fault,
                     "Mmax": Mmax,
                     "sliprate": sr,
                     "NMS": NMS_i,
                     "sr_seismic": sr_seismic_i}))
                    
                    id +=1

                feature_collection = FeatureCollection(features)
                
                if not os.path.exists(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set):
                    os.makedirs(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+scenario_set)

                with open(str(Run_name) +'/analysis/figures/FtF/'+Model+'/'+MFD_type+'_'+str(scenario_set)+'.geojson', 'w') as f:
                   dump(feature_collection, f)





#
#
#
#
#
#
#        '''##################################################################
#        # map the slip_rate
#        ###################################################################'''
#
#        mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
#                      dtype = [('U100'),('U100'),('U1000'),('f8'),('f8')], delimiter = '\t')
#        model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
#        set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
#        fault_name_mean_param =list( map(lambda i : mean_param[i][2], range(len(mean_param))))
##            print(set_mean_param[0],type(set_mean_param[0]))
##            print(fault_name_mean_param[0],type(fault_name_mean_param[0]))
#        sr_mean =list( map(lambda i : mean_param[i][3], range(len(mean_param))))
#        Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
#
#        index_model = np.where(np.array(model_mean_param)==Model)[0]
#        set_mean_param = np.take(set_mean_param,index_model)
#        fault_name_mean_param = np.take(fault_name_mean_param,index_model)
#        sr_mean = np.take(sr_mean,index_model)
#        Mmax_mean = np.take(Mmax_mean,index_model)
#
#        index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
#        fault_name_mean_param = np.take(fault_name_mean_param,index_set)
#        sr_mean = np.take(sr_mean,index_set)
#        Mmax_mean = np.take(Mmax_mean,index_set)
#
#
#        m_sr = Basemap(projection='mill',
#                  llcrnrlon=llcrnrlon,
#                  llcrnrlat=llcrnrlat,
#                  urcrnrlon=urcrnrlon,
#                  urcrnrlat=urcrnrlat,resolution='i')
#
#
#        Lon_bg = []
#        Lat_bg = []
#
#        # manually defined  in the file Background geometry
#        geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
#
#        column_model =list( map(lambda i : geom_bg[i][0],range(len(geom_bg))))
#        index_model = np.where(np.array(column_model) == Model)[0]
#        Lon_bg =list( map(lambda i : geom_bg[i][1],index_model))
#        Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
#        if len(Lon_bg) != 0 : #draw the background zone
#            maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_sr ,'g' , 0.02, 0.05, 'k')
#
#        #for each fault
#
#        for index_source in range(nb_sources):
#            if fault_type[index_source] == 'sf':
#                dip_tan = math.tan(math.radians(Dip[index_source]))
#
#                hdist_u = U_sism_depth[index_source] / dip_tan
#                hdist_d = L_sism_depth[index_source] / dip_tan
#
#                azimuth = [0.,0.,0.,0.]
#
#                lon_u = []
#                lat_u = []
#                lon_d = []
#                lat_d = []
#                if len(Lon[index_source]) < 9 :
#                    azimuth = []
#                    indexx = range(len(Lon[index_source]))
#                    for i in indexx:
#                        strike = abs(math.degrees(math.atan((
#                        Lon[index_source][-1]-Lon[index_source][0])/(
#                        Lat[index_source][-1]-Lat[index_source][0]))))
#                        if strike > 180. :
#                            strike = strike - 180.
#                        azimuth.append((strike + 90.0) % 360)
#
#
#                else:
#                    indexx = range(len(Lon[index_source]))
#                    for i in indexx[4:-4]:
#                        strike = abs(math.degrees(math.atan((
#                        Lon[index_source][i+4]-Lon[index_source][i-4])/(
#                        Lat[index_source][i+4]-Lat[index_source][i-4]))))
#                        if strike > 180. :
#                            strike = strike - 180.
#                        azimuth.append((strike + 90.0) % 360)
#
#                    azimuth[0] = azimuth[4]
#                    azimuth[1] = azimuth[4]
#                    azimuth[2] = azimuth[4]
#                    azimuth[3] = azimuth[4]
#                    azimuth.append(azimuth[-4])
#                    azimuth.append(azimuth[-4])
#                    azimuth.append(azimuth[-4])
#                    azimuth.append(azimuth[-4])
#
#                for i in indexx:
#                    if (Lon[index_source][0])>(Lon[index_source][-1]):
#                        x_u = hdist_u * math.sin(math.radians(azimuth[i]))
#                        y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
#                        x_d = hdist_d * math.sin(math.radians(azimuth[i]))
#                        y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
#                        lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
#                        lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
#                    if (Lon[index_source][0])<(Lon[index_source][-1]):
#                        x_u = hdist_u * math.sin(math.radians(azimuth[i]))
#                        y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
#                        x_d = hdist_d * math.sin(math.radians(azimuth[i]))
#                        y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
#                        lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
#                        lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
#
#
#
#                source_name_i = source_name[index_source].replace(Model+'_','')
##                    print(source_name_i)
##                    print(fault_name_mean_param)
#                index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
#                sr = sr_mean[index_fault]
#
#                x,y = m_sr(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
#                if nb_on_maps == True :
#                    plt.text(x,y,''+str(round(sr,2))+'',fontsize = 2)
#
#                cmap = matplotlib.cm.get_cmap('rainbow')
#                rgba = cmap(float(sr)/max(sr_mean))
#
#                x, y = m_sr(Lon[index_source], Lat[index_source])
#
#                #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
#
#                m_sr.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
#
#
#                poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
#                poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
#
#                if len(Lon[index_source]) < 3 :
#                    poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
#                    poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
#                maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_sr ,rgba , 0.5, 1., rgba)
#            if fault_type[index_source] == 'cf':
#                source_name_i = source_name[index_source].replace(Model+'_','')
#                index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
#                sr = sr_mean[index_fault]
#
#                x,y = m_sr(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
#                if nb_on_maps == True :
#                    plt.text(x,y,''+str(round(sr,2))+'',fontsize = 2)
#
#                cmap = matplotlib.cm.get_cmap('rainbow')
#                rgba = cmap(float(sr)/max(sr_mean))
#                maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_sr ,rgba , 0.5, 1., rgba)
#
#        m_sr.drawcoastlines(linewidth=0.2)
#        m_sr.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
#        plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_sr_'+str(scenario_set)+'.png',dpi=300)
#
#        plt.close()
#
#
#        a = np.array([[0,1]])
#        pl.figure(figsize=(9, 1.5))
#        img = pl.imshow(a, cmap="rainbow")
#        pl.gca().set_visible(False)
#        cax = pl.axes([0.1, 0.2, 0.4, 0.3])
#        cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
#        cb.set_ticklabels([str(round(min(sr_mean),2)),
#                           str(round(min(sr_mean)+ 0.2*(max(sr_mean)-min(sr_mean)),2)),
#                         str(round(min(sr_mean)+ 0.4*(max(sr_mean)-min(sr_mean)),2)),
#                         str(round(min(sr_mean)+ 0.6*(max(sr_mean)-min(sr_mean)),2)),
#                         str(round(min(sr_mean)+ 0.8*(max(sr_mean)-min(sr_mean)),2)),
#                         str(round(max(sr_mean),1))])
#        pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_sliprate.png',dpi=180)
#        pl.close()
#
#
#
#
#        #'''
#        # map the seismic slip-rate
#        if plot_sr_use == True:
#            for MFD_type in MFD_type_list:
#                slip_rep_data = np.genfromtxt(Run_name + '/analysis/txt_files/slip_rep_on_faults_mean_'+str(Model)+'_'+ MFD_type +'_' +str(scenario_set)+'.txt',
#                          dtype = [('U100'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),
#                                   ('f8'),('f8'),('f8'),('f8'),('f8'),('f8'),('f8')], delimiter = '\t')
#                fault_name_rep = list(map(lambda i : slip_rep_data[i][0], range(len(slip_rep_data))))
#                p_NMS = list(map(lambda i : slip_rep_data[i][14], range(len(slip_rep_data))))
#
#
#                mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
#                              dtype = [('U100'),('U100'),('U100'),('f8'),('f8')], delimiter = '\t')
#                model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
#                set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
#                fault_name_mean_param = list(map(lambda i : mean_param[i][2], range(len(mean_param))))
#                sr_mean = list(map(lambda i : mean_param[i][3], range(len(mean_param))))
#                Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
#
#                index_model = np.where(np.array(model_mean_param)==Model)[0]
#                set_mean_param = np.take(set_mean_param,index_model)
#                fault_name_mean_param = np.take(fault_name_mean_param,index_model)
#                sr_mean = np.take(sr_mean,index_model)
#                Mmax_mean = np.take(Mmax_mean,index_model)
#
#                index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
#                fault_name_mean_param = np.take(fault_name_mean_param,index_set)
#                sr_mean = np.take(sr_mean,index_set)
#                Mmax_mean = np.take(Mmax_mean,index_set)
#
#
#
#                m_sr_seismic = Basemap(projection='mill',
#                          llcrnrlon=llcrnrlon,
#                          llcrnrlat=llcrnrlat,
#                          urcrnrlon=urcrnrlon,
#                          urcrnrlat=urcrnrlat,resolution='h')
#
#
#
#
#                Lon_bg = []
#                Lat_bg = []
#
#                # manually defined  in the file Background geometry
#                geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
#
#                column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
#                index_model = np.where(np.array(column_model) == Model)[0]
#                Lon_bg =list( map(lambda i : geom_bg[i][1],index_model))
#                Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
#                if len(Lon_bg) != 0 : #draw the background zone
#                    maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_sr_seismic ,'g' , 0.02, 0.05, 'k')
#
#                #for each fault
#
#                for index_source in range(nb_sources):
#                    if fault_type[index_source] == 'sf':
#                        dip_tan = math.tan(math.radians(Dip[index_source]))
#
#                        hdist_u = U_sism_depth[index_source] / dip_tan
#                        hdist_d = L_sism_depth[index_source] / dip_tan
#
#                        azimuth = [0.,0.,0.,0.]
#
#                        lon_u = []
#                        lat_u = []
#                        lon_d = []
#                        lat_d = []
#                        if len(Lon[index_source]) < 9 :
#                            azimuth = []
#                            indexx = range(len(Lon[index_source]))
#                            for i in indexx:
#                                strike = abs(math.degrees(math.atan((
#                                Lon[index_source][-1]-Lon[index_source][0])/(
#                                Lat[index_source][-1]-Lat[index_source][0]))))
#                                if strike > 180. :
#                                    strike = strike - 180.
#                                azimuth.append((strike + 90.0) % 360)
#
#
#                        else:
#                            indexx = range(len(Lon[index_source]))
#                            for i in indexx[4:-4]:
#                                strike = abs(math.degrees(math.atan((
#                                Lon[index_source][i+4]-Lon[index_source][i-4])/(
#                                Lat[index_source][i+4]-Lat[index_source][i-4]))))
#                                if strike > 180. :
#                                    strike = strike - 180.
#                                azimuth.append((strike + 90.0) % 360)
#
#                            azimuth[0] = azimuth[4]
#                            azimuth[1] = azimuth[4]
#                            azimuth[2] = azimuth[4]
#                            azimuth[3] = azimuth[4]
#                            azimuth.append(azimuth[-4])
#                            azimuth.append(azimuth[-4])
#                            azimuth.append(azimuth[-4])
#                            azimuth.append(azimuth[-4])
#
#                        for i in indexx:
#                            if (Lon[index_source][0])>(Lon[index_source][-1]):
#                                x_u = hdist_u * math.sin(math.radians(azimuth[i]))
#                                y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
#                                x_d = hdist_d * math.sin(math.radians(azimuth[i]))
#                                y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
#                                lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                                lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
#                                lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                                lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
#                            if (Lon[index_source][0])<(Lon[index_source][-1]):
#                                x_u = hdist_u * math.sin(math.radians(azimuth[i]))
#                                y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
#                                x_d = hdist_d * math.sin(math.radians(azimuth[i]))
#                                y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
#                                lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                                lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
#                                lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                                lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
#
#
#
#                        source_name_i = source_name[index_source].replace(Model+'_','')
#                        #print source_name_i
#                        index_fault = np.where(np.array(fault_name_rep)==source_name_i)[0][0]
#                        NMS_i = p_NMS[index_fault]
#
#                        index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
#                        sr = sr_mean[index_fault]
#
#                        sr_seismic = sr*(1. - float(NMS_i)/100.)
#                        x,y = m_sr_seismic(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
#                        if nb_on_maps == True :
#                            plt.text(x,y,''+str(round(sr_seismic,2))+'',fontsize = 2)
#
#
#                        cmap = matplotlib.cm.get_cmap('rainbow')
#                        rgba = cmap(sr_seismic/max(sr_mean))
#
#                        x, y = m_sr_seismic(Lon[index_source], Lat[index_source])
#
#                        #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
#
#                        m_sr_seismic.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
#
#
#                        poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
#                        poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
#
#                        if len(Lon[index_source]) < 3 :
#                            poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
#                            poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
#                        maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_sr_seismic ,rgba , 0.5, 1., rgba)
#                    if fault_type[index_source] == 'cf':
#                        source_name_i = source_name[index_source].replace(Model+'_','')
#                        #print source_name_i
#                        index_fault = np.where(np.array(fault_name_rep)==source_name_i)[0][0]
#                        NMS_i = p_NMS[index_fault]
#
#                        index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
#                        sr = sr_mean[index_fault]
#
#                        sr_seismic = sr*(1. - float(NMS_i)/100.)
#                        x,y = m_sr_seismic(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
#                        if nb_on_maps == True :
#                            plt.text(x,y,''+str(round(sr_seismic,2))+'',fontsize = 2)
#
#
#                        cmap = matplotlib.cm.get_cmap('rainbow')
#                        rgba = cmap(sr_seismic/max(sr_mean))
#                        maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_sr_seismic ,rgba , 0.5, 1., rgba)
#
#                m_sr_seismic.drawcoastlines(linewidth=0.2)
#                m_sr_seismic.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
#                plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_sr_seismic_'+MFD_type+'_'+str(scenario_set)+'.png',dpi=300)
#
#                plt.close()
#
#
#                a = np.array([[0,1]])
#                pl.figure(figsize=(9, 1.5))
#                img = pl.imshow(a, cmap="rainbow")
#                pl.gca().set_visible(False)
#                cax = pl.axes([0.1, 0.2, 0.4, 0.3])
#                cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
#                cb.set_ticklabels([str(round(min(sr_mean),2)),
#                                   str(round(min(sr_mean)+ 0.2*(max(sr_mean)-min(sr_mean)),2)),
#                                 str(round(min(sr_mean)+ 0.4*(max(sr_mean)-min(sr_mean)),2)),
#                                 str(round(min(sr_mean)+ 0.6*(max(sr_mean)-min(sr_mean)),2)),
#                                 str(round(min(sr_mean)+ 0.8*(max(sr_mean)-min(sr_mean)),2)),
#                                 str(round(max(sr_mean),2))])
#                pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_sliprate_seismic.png',dpi=180)
#                pl.close()
#
#
#
#
#        #'''
#        # map the Mmax
#        mean_param = np.genfromtxt(Run_name + '/analysis/txt_files/mean_parameters_faults.txt',
#                      dtype = [('U100'),('U100'),('U100'),('f8'),('f8')], delimiter = '\t')
#        model_mean_param = list(map(lambda i : mean_param[i][0], range(len(mean_param))))
#        set_mean_param = list(map(lambda i : mean_param[i][1], range(len(mean_param))))
#        fault_name_mean_param =list( map(lambda i : mean_param[i][2], range(len(mean_param))))
#        sr_mean = list(map(lambda i : mean_param[i][3], range(len(mean_param))))
#        Mmax_mean = list(map(lambda i : mean_param[i][4], range(len(mean_param))))
#
#        index_model = np.where(np.array(model_mean_param)==Model)[0]
#        set_mean_param = np.take(set_mean_param,index_model)
#        fault_name_mean_param = np.take(fault_name_mean_param,index_model)
#        sr_mean = np.take(sr_mean,index_model)
#        Mmax_mean = np.take(Mmax_mean,index_model)
#
#        index_set = np.where(np.array(set_mean_param)==scenario_set)[0]
#        fault_name_mean_param = np.take(fault_name_mean_param,index_set)
#        sr_mean = np.take(sr_mean,index_set)
#        Mmax_mean = np.take(Mmax_mean,index_set)
#
#
#        m_mmax = Basemap(projection='mill',
#                  llcrnrlon=llcrnrlon,
#                  llcrnrlat=llcrnrlat,
#                  urcrnrlon=urcrnrlon,
#                  urcrnrlat=urcrnrlat,resolution='h')
#
#
#
#
#        Lon_bg = []
#        Lat_bg = []
#
#        # manually defined  in the file Background geometry
#        geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
#
#        column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
#        index_model = np.where(np.array(column_model) == Model)[0]
#        Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
#        Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
#        if len(Lon_bg) != 0 : #draw the background zone
#            maps_utils.draw_screen_poly(Lon_bg, Lat_bg,  m_mmax ,'g' , 0.02, 0.05, 'k')
#
#        #for each fault
#        for index_source in range(nb_sources):
#            if fault_type[index_source] == 'sf':
#                dip_tan = math.tan(math.radians(Dip[index_source]))
#
#                hdist_u = U_sism_depth[index_source] / dip_tan
#                hdist_d = L_sism_depth[index_source] / dip_tan
#
#                azimuth = [0.,0.,0.,0.]
#
#                lon_u = []
#                lat_u = []
#                lon_d = []
#                lat_d = []
#                if len(Lon[index_source]) < 9 :
#                    azimuth = []
#                    indexx = range(len(Lon[index_source]))
#                    for i in indexx:
#                        strike = abs(math.degrees(math.atan((
#                        Lon[index_source][-1]-Lon[index_source][0])/(
#                        Lat[index_source][-1]-Lat[index_source][0]))))
#                        if strike > 180. :
#                            strike = strike - 180.
#                        azimuth.append((strike + 90.0) % 360)
#
#
#                else:
#                    indexx = range(len(Lon[index_source]))
#                    for i in indexx[4:-4]:
#                        strike = abs(math.degrees(math.atan((
#                        Lon[index_source][i+4]-Lon[index_source][i-4])/(
#                        Lat[index_source][i+4]-Lat[index_source][i-4]))))
#                        if strike > 180. :
#                            strike = strike - 180.
#                        azimuth.append((strike + 90.0) % 360)
#
#                    azimuth[0] = azimuth[4]
#                    azimuth[1] = azimuth[4]
#                    azimuth[2] = azimuth[4]
#                    azimuth[3] = azimuth[4]
#                    azimuth.append(azimuth[-4])
#                    azimuth.append(azimuth[-4])
#                    azimuth.append(azimuth[-4])
#                    azimuth.append(azimuth[-4])
#
#                for i in indexx:
#                    if (Lon[index_source][0])>(Lon[index_source][-1]):
#                        x_u = hdist_u * math.sin(math.radians(azimuth[i]))
#                        y_u = hdist_u * math.cos(math.radians(180. - azimuth[i]))
#                        x_d = hdist_d * math.sin(math.radians(azimuth[i]))
#                        y_d = hdist_d * math.cos(math.radians(180. - azimuth[i]))
#                        lon_u.append(Lon[index_source][i] + (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_u.append(Lat[index_source][i] + y_u/40007.*360.)
#                        lon_d.append(Lon[index_source][i] + (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_d.append(Lat[index_source][i] + y_d/40007.*360.)
#                    if (Lon[index_source][0])<(Lon[index_source][-1]):
#                        x_u = hdist_u * math.sin(math.radians(azimuth[i]))
#                        y_u = hdist_u * math.cos(math.radians(180. + azimuth[i]))
#                        x_d = hdist_d * math.sin(math.radians(azimuth[i]))
#                        y_d = hdist_d * math.cos(math.radians(180. + azimuth[i]))
#                        lon_u.append(Lon[index_source][i] - (x_u/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_u.append(Lat[index_source][i] - y_u/40007.*360.)
#                        lon_d.append(Lon[index_source][i] - (x_d/40075.*360.))#*((90.-abs(Lon[index_source][i]))/90.))
#                        lat_d.append(Lat[index_source][i] - y_d/40007.*360.)
#
#
#
#                source_name_i = source_name[index_source].replace(Model+'_','')
#                index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
#                Mmax = Mmax_mean[index_fault]
#
#                x,y = m_mmax(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
#                if nb_on_maps == True :
#                    plt.text(x,y,''+str(Mmax)+'',fontsize = 2)
#
#
#                cmap = matplotlib.cm.get_cmap('rainbow')
#                rgba = cmap((float(Mmax)-min(Mmax_mean))/(max(Mmax_mean)-min(Mmax_mean)))
#
#                x, y = m_mmax(Lon[index_source], Lat[index_source])
#
#                #go pick in the file for the NMS and do the colors with that. find a go color map (blue to green?)
#
#                m_mmax.plot(x, y, 'D-', markersize=1., linewidth=0.01, color=rgba, markerfacecolor=rgba, markeredgecolor = rgba)
#
#
#                poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
#                poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
#
#                if len(Lon[index_source]) < 3 :
#                    poly_lons = np.concatenate([lon_u,np.array(list(reversed(lon_d)))])
#                    poly_lats = np.concatenate([lat_u,np.array(list(reversed(lat_d)))])
#                maps_utils.draw_screen_poly(poly_lons, poly_lats,  m_mmax ,rgba , 0.5, 1., rgba)
#            if fault_type[index_source] == 'cf':
#                source_name_i = source_name[index_source].replace(Model+'_','')
#                index_fault = np.where(np.array(fault_name_mean_param)==source_name_i)[0][0]
#                Mmax = Mmax_mean[index_fault]
#
#                x,y = m_mmax(np.mean(Lon[index_source]),np.mean(Lat[index_source]))
#                if nb_on_maps == True :
#                    plt.text(x,y,''+str(Mmax)+'',fontsize = 2)
#
#
#                cmap = matplotlib.cm.get_cmap('rainbow')
#                rgba = cmap((float(Mmax)-min(Mmax_mean))/(max(Mmax_mean)-min(Mmax_mean)))
#                maps_utils.draw_screen_poly(Lon[index_source], Lat[index_source],  m_mmax ,rgba , 0.5, 1., rgba)
#
#        m_mmax.drawcoastlines(linewidth=0.2)
#        m_mmax.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
#        plt.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'map_Mmax_'+str(scenario_set)+'.png',dpi=300)
#
#        plt.close()
#
#
#        a = np.array([[0,1]])
#        pl.figure(figsize=(9, 1.5))
#        img = pl.imshow(a, cmap="rainbow")
#        pl.gca().set_visible(False)
#        cax = pl.axes([0.1, 0.2, 0.4, 0.3])
#        cb = pl.colorbar(orientation="horizontal", cax=cax,ticks=[0., 0.2, 0.4, 0.6, 0.8, 1])
#        cb.set_ticklabels([str(min(Mmax_mean)),
#                           str(min(Mmax_mean)+ 0.2*(max(Mmax_mean)-min(Mmax_mean))),
#                         str(min(Mmax_mean)+ 0.4*(max(Mmax_mean)-min(Mmax_mean))),
#                         str(min(Mmax_mean)+ 0.6*(max(Mmax_mean)-min(Mmax_mean))),
#                         str(min(Mmax_mean)+ 0.8*(max(Mmax_mean)-min(Mmax_mean))),
#                         str(max(Mmax_mean))])
#        pl.savefig(str(Run_name) + '/analysis/figures/FtF/'+Model+'/'+'colorbar_Mmax.png',dpi=180)
#        pl.close()
                

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
        

   
