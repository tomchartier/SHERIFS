#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""

import numpy as np
import math
import os
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.interpolate import interp1d
import random

mpl.interactive(False)
import sys


'''#############################
###############################
#     Extraction of 
#   the catalog 
###############################
##############################'''

def Geom_bg(Model_name,File_bg):
    Lon_bg = []
    Lat_bg = []

    # manually defined  in the file Background geometry
    geom_bg = np.genfromtxt(File_bg,dtype=[('U100'),('f'),('f')],skip_header = 1)
    
    column_model = list(map(lambda i : str(geom_bg[i][0]),range(len(geom_bg))))
    index_model = np.where(np.array(column_model) == Model_name)[0]
    Lon_bg = map(lambda i : geom_bg[i][1],index_model)
    Lat_bg = map(lambda i : geom_bg[i][2],index_model)
    return list(Lon_bg), list(Lat_bg)
    
    
    

def plt_catalog(do_catalog,Model_list,File_bg,catalog_file,Run_name,xmin,xmax,ymin,ymax,
                llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,completness_file,
                nb_inter,bining_in_mag,end_year_of_catalog,sub_area_file):
    if do_catalog == True:
        if not os.path.exists(str(Run_name) + '/analysis/figures/catalogue'):
            os.makedirs(str(Run_name) + '/analysis/figures/catalogue')
    
        catalog_cum_rate = []
        index_model = 0
        
        for model in Model_list : 
            file_catalog_rate_all = open(str(Run_name) + '/analysis/figures/catalogue/catalog_rates_all_'+model+'.txt','w')
            for mag in bining_in_mag[:-1]:
                file_catalog_rate_all.write(str(round(mag,1))+'\t')
            file_catalog_rate_all.write(str(round(bining_in_mag[-1],1))+'\n') 
                
            # extract the geometry of the zone ( geometry of the background)
            Lon_bg, Lat_bg = Geom_bg(model,File_bg)
            ColX = Lon_bg
            ColY = Lat_bg                       
            Poly = []   
            for x1,y1 in zip(ColX,ColY): # creation du polygon de la zone
                Poly.append((x1,y1))
            bbPath = mplPath.Path(Poly)
            
            bbPath_sub_areas = []
            if os.path.exists(sub_area_file):
                read_sub_area_file = open(sub_area_file,'rU')
                lines_sub_area = read_sub_area_file.readlines()
                sub_area_names = []
                sub_area_coord = []
                sub_area_lon = []
                sub_area_lat = []
                for line in lines_sub_area:
                    model_sub_area = line.split('\t')[0]
                    if model == model_sub_area:
                        sub_area_names.append(line.split('\t')[1])
                        sub_area_coord.append(line.split('\t')[2:])
                        sub_area_lon_i = []
                        sub_area_lat_i = []
                        for sub_area_coord_i in line.split('\t')[2:]:
                            if not '\n' in sub_area_coord_i.split(','):
                                if not '' in sub_area_coord_i.split(','):
                                    sub_area_lon_i.append(float(sub_area_coord_i.split(',')[1]))
                                    sub_area_lat_i.append(float(sub_area_coord_i.split(',')[0]))
                        sub_area_lon.append(sub_area_lon_i)
                        sub_area_lat.append(sub_area_lat_i)
                        if not os.path.exists(str(Run_name) + '/analysis/figures/catalogue/sub_area'):
                            os.makedirs(str(Run_name) + '/analysis/figures/catalogue/sub_area')  
                                     
                        Poly_sub = []   
                        for x1,y1 in zip(sub_area_lon_i,sub_area_lat_i): # creation du polygon de la zone
                            Poly_sub.append((x1,y1))    
                        bbPath_sub_areas.append(mplPath.Path(Poly_sub))
                            
            cat_data = np.genfromtxt(catalog_file,dtype=[('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),('S100'),
                                                                     ('S100'),('S100')],skip_header = 1)
            cat_lon =  list(map(lambda i : float(cat_data[i][5]), range(len(cat_data))))
            cat_lat = list( map(lambda i : float(cat_data[i][4]), range(len(cat_data))))
            cat_depth =  list(map(lambda i : float(cat_data[i][6]), range(len(cat_data))))
            cat_Mw =  list(map(lambda i : float(cat_data[i][7]), range(len(cat_data))))
            cat_sig_Mw = list(map(lambda i : float(cat_data[i][8]), range(len(cat_data))))
            cat_Yr =  list(map(lambda i : float(cat_data[i][0]), range(len(cat_data))))
            
            #selects only the earthquakes that are in the subarea
            indexes_in = []
            for index_cat in range(len(cat_lon)):
                if bbPath.contains_point((cat_lon[index_cat],cat_lat[index_cat])) == 1: #test pour savoir si la zone contient le point
                    indexes_in.append(index_cat)
            cat_lon = np.take(cat_lon,indexes_in)
            cat_lat = np.take(cat_lat,indexes_in)
            cat_depth = np.take(cat_depth,indexes_in)
            cat_Mw = np.take(cat_Mw,indexes_in)
            cat_sig_Mw = np.take(cat_sig_Mw,indexes_in)
            cat_Yr = np.take(cat_Yr,indexes_in)
    
    
            #extracting the completeness time of the catalogue
            completness = []
            weights_completness = []
            read_comp_file = open(completness_file,'rU')
            lines_of_the_file = read_comp_file.readlines()
            #lines_of_the_file.remove('\n')
            line_number=0
            for i in range(int(len(lines_of_the_file)/2)):
                if len(lines_of_the_file[line_number].split('\t'))!=0:
                    binning_comp_i = lines_of_the_file[line_number].split('\t')
                    if '\r\n' in binning_comp_i:
                        binning_comp_i.remove('\r\n')
                    if '\n' in binning_comp_i:
                        binning_comp_i.remove('\n')
                    binning_comp = []
                    for magnitudes in binning_comp_i[1:]:
                        mag = magnitudes.split(',')
                        binning_comp.append(float(mag[0]))
                        binning_comp.append(float(mag[1]))
                    #binning_comp = [float(x) for x in binning_comp]
                    binning_comp.append(bining_in_mag[-1])                
                    comp_value_i = lines_of_the_file[line_number+1].split('\t')
                    if '\r\n' in comp_value_i:
                        comp_value_i.remove('\r\n')
                    if '\n' in comp_value_i:
                        comp_value_i.remove('\n')
                    comp_value = []
                    for value in comp_value_i[1:]:
                        comp_value.append(float(value))
                        comp_value.append(float(value))
                    #comp_value = [float(x) for x in comp_value]
                    comp_value.append(comp_value[-1])
                    
                    completeness_interpolate = interp1d(binning_comp,comp_value)                
                    completness_i = []
                    for mag in bining_in_mag:
                        try :
                            completness_i.append(completeness_interpolate(mag))
                        except ValueError:
                            print('!!!!!!!!!!!!\n\n\nERROR in completeness file\n your minimum magnitude might not be low enough\n\n\n!!!!!!!!!!!!')
                    line_number += 2
                    completness.append(completness_i)
                    weights_completness.append(float(comp_value_i[0]))
            seismological_moment_rate = []
    
            number_of_earthquakes_for_rate = [] #nuMber of earthquake greater than this magnitude
            number_of_earthquakes_for_rate_sub_area = [] #nuMber of earthquake greater than this magnitude in the sub areas
            for poly_sub_area_i in bbPath_sub_areas:
                number_of_earthquakes_for_rate_sub_area.append([])
            
            #variables for the plot on the map
            file_catalog_for_map = open(str(Run_name) + '/analysis/figures/catalogue/catalog_for_map_'+str(model)+'.txt','w')
            earthquake_for_map = []
            nb_time_picked = [] #number of time this earthquake has been kept in the complete catalog
            
            #self.rate_cat_cum_kdeplot = []
            rate_cat_cum_density = []
            
            #for the sub_areas
            rate_in_sub_area_cum= []
            for poly_sub_area_i in bbPath_sub_areas:
                rate_in_sub_area_cum.append([])
                
            #loop on the interations to explore the uncertainties
            for i in range(nb_inter):
                cat_model_bin = np.zeros(len(bining_in_mag))
                rate_cat = np.zeros(len(bining_in_mag))        
                seismological_moment_rate_i = 0. #moment rate calculated using the catalog
                test_if_eq_here = [] #will be checked if the earthquake is not count twice
                
                #random picking of the completness used
                index_completness = np.random.choice(len(weights_completness), 1, p = weights_completness)[0] 
                completness_used = completness[index_completness]
    
                #for the sub areas
                cat_sub_area_bin = []
                rate_sub_area = []
                for poly_sub_area_i in bbPath_sub_areas:   
                    cat_sub_area_bin.append(np.zeros(len(bining_in_mag))  )
                    rate_sub_area.append(np.zeros(len(bining_in_mag))  )   
    
                    
                #random pick of the earthquakes magnitudes
                envents_magnitude = []
                index_cat=0
                for Yr in cat_Yr: #loop on all the earthquakes in the catalogue to pick the magnitude
                    if cat_sig_Mw[index_cat]==0:
                        event_magnitude = cat_Mw[index_cat]
                    elif cat_sig_Mw[index_cat] <= 3.:
                        event_magnitude = np.random.triangular(cat_Mw[index_cat]-cat_sig_Mw[index_cat]/2.,cat_Mw[index_cat],cat_Mw[index_cat]+cat_sig_Mw[index_cat]/2.)
                        event_magnitude = np.random.normal(loc=cat_Mw[index_cat],scale=cat_sig_Mw[index_cat])
    #                else :
    #                    event_magnitude = np.random.triangular(cat_sig_Mw[index_cat],cat_Mw[index_cat],cat_sig_Mw_plus[index_cat])
                    envents_magnitude.append(event_magnitude)
                    index_cat+=1
                
                
                index_mag=0
                #for mag_i,completness_i,sigma_completness_i in zip(bining_in_mag,completeness,sigma_completness) :
                for mag_i,completness_i in zip(bining_in_mag,completness_used) :
                    #loop on the magnitudes
                    index_cat = 0
                    
                    #picked_completness = completness_i + random_factor * sigma_completness_i
                    picked_completness = completness_i
                    
                    index_cat =0     
                    for Yr in cat_Yr: #loop on all the earthquakes in the catalogue
                        if cat_depth[index_cat] < 30.  or math.isnan(cat_depth[index_cat]): #check if the event is not on the subduction interface for Corinth
                            if Yr >= picked_completness and Yr <= end_year_of_catalog :
                                event_magnitude = envents_magnitude[index_cat]
                                if event_magnitude >= mag_i and event_magnitude < mag_i+0.099 :
                                    if bbPath.contains_point((cat_lon[index_cat],cat_lat[index_cat])) == 1: #test pour savoir si la zone contient le point
                                        if not (str(Yr)+str(cat_lon[index_cat])+str(cat_lat[index_cat])+str(cat_sig_Mw[index_cat])) in test_if_eq_here:
                                            #cat_SHARE_WCR_i.append(event_magnitude)
                                            test_if_eq_here.append(str(Yr)+str(cat_lon[index_cat])+str(cat_lat[index_cat]))
                                        
                                            cat_model_bin[index_mag] += 1  
                                            
                                            #for the map
                                            string_for_map = str(Yr) +'\t'+ str(cat_Mw[index_cat]) +'\t'+ str(cat_lon[index_cat]) +'\t'+  str(cat_lat[index_cat])
                                            if not string_for_map in earthquake_for_map : #the earthquake has not been pîcked in an earlier iteration
                                                earthquake_for_map.append(string_for_map)
                                                nb_time_picked.append(1)
                                            else :
                                                index_in_cat_for_map = np.where(np.array(earthquake_for_map)==string_for_map)[0][0]
                                                nb_time_picked[index_in_cat_for_map] += 1 #if earthquake arealy there, we count it one more time
                                                        
                                            #for the sub areas
                                            index_sub_area = 0
                                            for poly_sub_area_i in bbPath_sub_areas:
                                                if poly_sub_area_i.contains_point((cat_lon[index_cat],cat_lat[index_cat])) == 1: #test pour savoir si la zone contient le point
                                                    cat_sub_area_bin[index_sub_area][index_mag] += 1 
                                                index_sub_area+=1
                                                
                        index_cat += 1
                    rate_cat[index_mag] = cat_model_bin[index_mag] / (end_year_of_catalog - picked_completness)  
                    M0 = 10. ** (1.5 * mag_i + 9.1)  #calculation of the moment
                    rate_M0 = M0 * rate_cat[index_mag]
                    seismological_moment_rate_i += rate_M0
                    
                    
                    #for the sub areas
                    index_sub_area = 0
                    for poly_sub_area_i in bbPath_sub_areas:
                        rate_sub_area[index_sub_area][index_mag] = cat_sub_area_bin[index_sub_area][index_mag] / (end_year_of_catalog - picked_completness) 
                        index_sub_area+=1
                    
                    index_mag += 1
                    
                seismological_moment_rate.append(seismological_moment_rate_i)  
                
                number_of_earthquakes_for_rate_i=[]
                for i in range(len(cat_model_bin)): #calculate the cumulative number of earthquake
                    number_of_earthquakes_for_rate_i.append(int(np.sum(np.array(cat_model_bin)[-(len(cat_model_bin)-i):])))    
                number_of_earthquakes_for_rate.append(number_of_earthquakes_for_rate_i) 
            
                #for the sub areas
                index_sub_area = 0
                for poly_sub_area_i in bbPath_sub_areas:
                    number_of_earthquakes_for_rate_i=[]
                    for i in range(len(cat_model_bin)): #calculate the cumulative number of earthquake
                        number_of_earthquakes_for_rate_i.append(int(np.sum(np.array(cat_sub_area_bin[index_sub_area])[-(len(cat_sub_area_bin[index_sub_area])-i):])))  
                    number_of_earthquakes_for_rate_sub_area[index_sub_area].append(number_of_earthquakes_for_rate_i) 
                    index_sub_area+=1
                        
                rate_cat_cum = [] 
                for i in range(len(rate_cat)): #calculate the cumulative rate
                    rate_cat_cum.append(np.sum(np.array(rate_cat)[-(len(rate_cat)-i):]))
                rate_cat_cum_density.append(rate_cat_cum)
                
                #for the sub areas
                index_sub_area = 0
                for poly_sub_area_i in bbPath_sub_areas:
                    rate_sub_area_i_cum = [] 
                    for i in range(len(rate_sub_area[index_sub_area])): #calculate the cumulative rate
                        rate_sub_area_i_cum.append(np.sum(np.array(rate_sub_area[index_sub_area])[-(len(rate_sub_area[index_sub_area])-i):]))
                    rate_in_sub_area_cum[index_sub_area].append(rate_sub_area_i_cum)
                    index_sub_area+=1
               
                plt.scatter(bining_in_mag,rate_cat_cum,c = 'brown' , marker = '_', s = 50, alpha = 0.2)#, linewidth = 0
                
                for rate_i in rate_cat_cum[:-1]:
                    file_catalog_rate_all.write(str(rate_i)+'\t')
                file_catalog_rate_all.write(str(rate_cat_cum[-1])+'\n') 
            file_catalog_rate_all.close()
                
            
            axes = plt.gca()
            axes.set_xlim([xmin,xmax])
            axes.set_ylim([ymin,ymax])
            mean_nb_eq_in_cat = np.mean(number_of_earthquakes_for_rate, axis = 0)
                    
            nb_sample = 50
            time_simu = 1000000
            alphas=[]
            
            for index_mag in range(len(bining_in_mag)): 
                
                #does a simulation of the catalog in order to check the statistical validity
                
                nb_eq = mean_nb_eq_in_cat[index_mag]
                if nb_eq <= 10: #if there is less than 10 earthquakes, do the statistical représentativity
                    time_obs = int(round(end_year_of_catalog - completness[0][index_mag]))
                    
                    rate = float(nb_eq)/float(time_obs)
                    
                    cat = []
                    #build the catalog
                    for i in range(time_simu):
                        test = np.random.uniform(0.,1.)
                        if test<=rate:
                            cat.append(1)
                        else:
                            cat.append(0)
                            
                    rates_simu = []
                    for i in range(nb_sample):
                        yrs_start = random.choice(range(time_simu-time_obs))
                        cat_sampled = cat[yrs_start:yrs_start+time_obs]
                        #print sum(cat_sampled), sum(cat_sampled)/float(time_obs)
                        rates_simu.append(sum(cat_sampled)/float(time_obs))
                       
                    unc_rate = (np.std(rate_cat_cum_density,axis=0)[index_mag])/(np.square(mean_nb_eq_in_cat[index_mag]))
                    unc_rate = np.std(rates_simu)
                    
                    rate_plus = np.array(rate_cat_cum_density).mean(axis=0)[index_mag]+unc_rate
                    rate_minus = np.array(rate_cat_cum_density).mean(axis=0)[index_mag]-unc_rate
                    if rate_minus < 0. : # in order not to have negative values
                        rate_minus = np.array(rate_cat_cum_density).mean(axis=0)[index_mag]/100.
                        
                    #if this is smaller than the percentils, take percentiles
                    if np.percentile(rate_cat_cum_density,84,axis=0)[index_mag] >rate_plus:
                        rate_plus = np.percentile(rate_cat_cum_density,84,axis=0)[index_mag]
                    if np.percentile(rate_cat_cum_density,16,axis=0)[index_mag]<rate_minus:
                        rate_minus = np.percentile(rate_cat_cum_density,16,axis=0)[index_mag]
                        
                else:  #if there is more, use the percentiles of the distribution
                    rate_plus = np.percentile(rate_cat_cum_density,84,axis=0)[index_mag]
                    rate_minus = np.percentile(rate_cat_cum_density,16,axis=0)[index_mag]
                    
                mag = bining_in_mag[index_mag]
                mag_plus = mag+0.05
                mag_minus = mag-0.05
                verts = [(mag_minus, rate_minus ),
                         (mag_minus, rate_plus),
                         (mag_plus, rate_plus),
                         (mag_plus, rate_minus),
                         (mag_minus, rate_minus)]
                codes = [Path.MOVETO,
                         Path.LINETO,
                         Path.LINETO,
                         Path.LINETO,
                         Path.CLOSEPOLY]
                #alpha depending of the number of samples that are different to zero
                alpha = float(np.count_nonzero(rate_cat_cum_density,axis=0)[index_mag])/float(nb_sample) * 0.2
                alphas.append(alpha)
                         
                path_poly = Path(verts, codes)
                
                patch = patches.PathPatch(path_poly,facecolor = 'k', lw = 0., alpha = alpha)
                axes.add_patch(patch)
            
    #        plt.scatter(bining_in_mag,np.percentile(rate_cat_cum_density,50,axis=0),
    #                c='k', s=5, edgecolor='',marker = 'o',alpha = 0.8)
    #        plt.scatter(bining_in_mag,np.percentile(rate_cat_cum_density,16,axis=0),
    #            c='k', s=5, edgecolor='',marker = '+',alpha = 0.8)
    #        plt.scatter(bining_in_mag,np.percentile(rate_cat_cum_density,84,axis=0),
    #                    c='k', s=5, edgecolor='',marker = '+',alpha = 0.8)
    #        plt.scatter(bining_in_mag,np.percentile(rate_cat_cum_density,5,axis=0),
    #            c='k', s=4, edgecolor='',marker = 'o',alpha = 0.8)
    #        plt.scatter(bining_in_mag,np.percentile(rate_cat_cum_density,95,axis=0),
    #                    c='k', s=4, edgecolor='',marker = 'o',alpha = 0.8)
    
            rgba_colors = np.zeros((len(bining_in_mag),4))
            # for red the first column needs to be one
            rgba_colors[:,0] = 0.
            # the fourth column needs to be your alphas
            rgba_colors[:, 3] = alphas
            
            plt.scatter(bining_in_mag,np.array(rate_cat_cum_density).mean(axis=0),
                    c=rgba_colors, s=20, edgecolor='',marker = 's')
                
            plt.yscale('log')
            plt.title('earthquake catalog')
            plt.grid(alpha=0.3)
            plt.savefig(str(Run_name) + '/analysis/figures/catalogue/catalogue_'+str(model)+'.png' , dpi = 180, transparent=True)
            #plt.show()
            plt.close()
            
            #export the rates of the catalogue for future comparisons
            #line_cat_for model = [model, rate_cat_cum_density]
            catalog_cum_rate.append(np.array(rate_cat_cum_density).mean(axis=0))
            
            
            
            #for the sub areas
            index_sub_area = 0
            for poly_sub_area_i in bbPath_sub_areas:
                
                axes = plt.gca()
                axes.set_xlim([xmin,xmax])
                axes.set_ylim([ymin,ymax])
                mean_nb_eq_in_cat = np.mean(number_of_earthquakes_for_rate_sub_area[index_sub_area], axis = 0)
                alphas=[]
                for index_mag in range(len(bining_in_mag)): 
                    
                    #does a simulation of the catalog in order to check the statistical validity
                    
                    nb_eq = mean_nb_eq_in_cat[index_mag]
                    if nb_eq <= 10: #if there is less than 10 earthquakes, do the statistical représentativity
                        time_obs = int(round(end_year_of_catalog - completness[0][index_mag]))
                        
                        
                        rate = float(nb_eq)/float(time_obs)
                        
                        cat = []
                        #build the catalog
                        for i in range(time_simu):
                            test = np.random.uniform(0.,1.)
                            if test<=rate:
                                cat.append(1)
                            else:
                                cat.append(0)
                                
                        rates_simu = []
                        for i in range(nb_sample):
                            yrs_start = random.choice(range(time_simu-time_obs))
                            cat_sampled = cat[yrs_start:yrs_start+time_obs]
                            #print sum(cat_sampled), sum(cat_sampled)/float(time_obs)
                            rates_simu.append(sum(cat_sampled)/float(time_obs))
                           
                        
                        unc_rate = np.std(rates_simu)
                        rate_plus = np.array(rate_in_sub_area_cum[index_sub_area]).mean(axis=0)[index_mag]+unc_rate
                        rate_minus = np.array(rate_in_sub_area_cum[index_sub_area]).mean(axis=0)[index_mag]-unc_rate
                        if rate_minus < 0. : # in order not to have negative values
                            rate_minus = np.array(rate_in_sub_area_cum[index_sub_area]).mean(axis=0)[index_mag]/10.
                            
                        #if this is smaller than the percentils, take percentiles
                        if np.percentile(rate_in_sub_area_cum[index_sub_area],84,axis=0)[index_mag] >rate_plus:
                            rate_plus = np.percentile(rate_in_sub_area_cum[index_sub_area],84,axis=0)[index_mag]
                        if np.percentile(rate_in_sub_area_cum[index_sub_area],16,axis=0)[index_mag]<rate_minus:
                            rate_minus = np.percentile(rate_in_sub_area_cum[index_sub_area],16,axis=0)[index_mag]
                            
                    else:  #if there is more, use the percentiles of the distribution
                        rate_plus = np.percentile(rate_in_sub_area_cum[index_sub_area],84,axis=0)[index_mag]
                        rate_minus = np.percentile(rate_in_sub_area_cum[index_sub_area],16,axis=0)[index_mag]
                            
                    mag = bining_in_mag[index_mag]
                    mag_plus = mag+0.05
                    mag_minus = mag-0.05
                    verts = [(mag_minus, rate_minus ),
                             (mag_minus, rate_plus),
                             (mag_plus, rate_plus),
                             (mag_plus, rate_minus),
                             (mag_minus, rate_minus)]
                    codes = [Path.MOVETO,
                             Path.LINETO,
                             Path.LINETO,
                             Path.LINETO,
                             Path.CLOSEPOLY]
                             
                    path_poly = Path(verts, codes)
                        
                    #alpha depending of the number of samples that are different to zero
                    alpha = float(np.count_nonzero(rate_in_sub_area_cum[index_sub_area],axis=0)[index_mag])/float(nb_sample) * 0.2
                    alphas.append(alpha)
                
                    patch = patches.PathPatch(path_poly,facecolor = 'k', lw = 0., alpha = alpha)
                    axes.add_patch(patch)
                
                for i in range(nb_inter):
                    plt.scatter(bining_in_mag,rate_in_sub_area_cum[index_sub_area][i],c = 'brown' , marker = '_', s = 50, alpha = 0.2)#, linewidth = 0
                    
                rgba_colors = np.zeros((len(bining_in_mag),4))
                # for red the first column needs to be one
                rgba_colors[:,0] = 0.
                # the fourth column needs to be your alphas
                rgba_colors[:, 3] = alphas
                
                plt.scatter(bining_in_mag,np.array(rate_in_sub_area_cum[index_sub_area]).mean(axis=0),
                        c=rgba_colors, s=20, edgecolor='',marker = 's')
                
                
                    
                axes = plt.gca()
                axes.set_xlim([xmin,xmax])
                axes.set_ylim([ymin,ymax])
                plt.grid()
                plt.yscale('log')
                plt.title('catalogue sub area ' + sub_area_names[index_sub_area])
                plt.savefig(str(Run_name) + '/analysis/figures/catalogue/sub_area/catalogue_sub_area_'+sub_area_names[index_sub_area]+' Model '+str(model)+'.png' , dpi = 180, transparent=True)
                #plt.show()
                plt.close()
                
                file_cat_sub_area = open(str(Run_name) + '/analysis/figures/catalogue/sub_area/mean_rate_cat_'+sub_area_names[index_sub_area]+' Model '+str(model)+'.txt','w')
                for index_mag in range(len(bining_in_mag)): 
                    file_cat_sub_area.write(str(bining_in_mag[index_mag])+'\t'+str(np.array(rate_in_sub_area_cum[index_sub_area]).mean(axis=0)[index_mag])+'\n')
                file_cat_sub_area.close()
                
                index_sub_area+=1
            
            
            #print 'Number of earthquake for the rate in catalogue'
            file_nb_eq_in_cat = open(str(Run_name) + '/analysis/figures/catalogue/nb_eq_in_cat_'+str(model)+'.txt','w')
            file_nb_eq_in_cat.write('Number of earthquake for the rate in catalogue'+'\n')
            mean_nb_eq_in_cat = np.mean(number_of_earthquakes_for_rate, axis = 0)
            for i in range(len(bining_in_mag)) :
                #print bining_in_mag[i],int(round(mean_nb_eq_in_cat[i],0))
                file_nb_eq_in_cat.write(str(bining_in_mag[i]) + '\t' + str(round(mean_nb_eq_in_cat[i],0))+'\n')
            file_nb_eq_in_cat.close()
            
            plt.scatter(bining_in_mag,mean_nb_eq_in_cat)
            plt.axhline(y=1)
            plt.axhline(y=2)
            plt.axhline(y=5)
            plt.axhline(y=10)
            plt.title('nb of earthquakes (m>M) in the catalogue')
            plt.savefig(str(Run_name) + '/analysis/figures/catalogue/'+'nb_of_EQ_in_catalogue.png' , dpi = 180, transparent=True)
            plt.grid()
            plt.close()
            
            for string,nb in zip(earthquake_for_map,nb_time_picked) :
                #print string,nb
                #if nb >= nb_inter/2 : #we write in the catalog only if selected more than half the time
                file_catalog_for_map.write(string + '\t'+ str(nb)+'\n')
                
            file_catalog_for_map.close()
            
            
            #for the sub areas
            index_sub_area = 0
            for poly_sub_area_i in bbPath_sub_areas:
                plt.scatter(bining_in_mag,np.mean(number_of_earthquakes_for_rate_sub_area[index_sub_area], axis = 0))
                plt.axhline(y=1)
                plt.axhline(y=2)
                plt.axhline(y=5)
                plt.axhline(y=10)
                plt.title('nb of earthquakes (m>M) in the catalogue ' + sub_area_names[index_sub_area])
                plt.savefig(str(Run_name) + '/analysis/figures/catalogue/sub_area/'+'nb_of_EQ_in_ ' + sub_area_names[index_sub_area]+'.png' , dpi = 180, transparent=True)
                plt.grid()
                plt.close()
                index_sub_area+=1
                
                
            ###### map plot of the catalog #####
            if mean_nb_eq_in_cat[0] > 3:
                catalog_for_map = np.genfromtxt(str(Run_name) + '/analysis/figures/catalogue/catalog_for_map_'+str(model)+'.txt',
                                                dtype=[('f8'),('f8'),('f8'),('f8'),('f8')])
                yr_cat_for_map =  list(map(lambda i : round(float(catalog_for_map[i][0])), range(len(catalog_for_map))))
                M_cat_for_map = list( map(lambda i : float(catalog_for_map[i][1]), range(len(catalog_for_map))))
                lon_cat_for_map =  list(map(lambda i : float(catalog_for_map[i][2]), range(len(catalog_for_map))))
                lat_cat_for_map = list( map(lambda i : float(catalog_for_map[i][3]), range(len(catalog_for_map))))
                nb_time_picked = list( map(lambda i : float(catalog_for_map[i][4]), range(len(catalog_for_map))))
                
                m = Basemap(projection='mill',
                              llcrnrlon = llcrnrlon, 
                              llcrnrlat = llcrnrlat, 
                              urcrnrlon = urcrnrlon, 
                              urcrnrlat = urcrnrlat,resolution='i')
                
                for (Yr,M,lon,lat,nb) in zip(yr_cat_for_map,M_cat_for_map,lon_cat_for_map,lat_cat_for_map,nb_time_picked):
                    x, y = m(lon, lat)
                    # color scale
                    if M< 5. :
                        c='b'
                        s = 3
                    elif M< 6.0:
                        c='g'
                        s = 5
                    elif M< 7.:
                        c='orange'
                        s= 7
                    else:
                        c='r'
                        s = 10
                    
                    if nb >= 2*nb_inter/3 :
                        alpha = 1.
                    if nb <= 2*nb_inter/3 :
                        alpha = 0.8
                    if nb <= 1*nb_inter/2 :
                        alpha = 0.3
                    if nb < 1*nb_inter/3 :
                        alpha = 0.1
                        
                    c = plt.cm.jet((M-min(M_cat_for_map))/(max(M_cat_for_map)-min(M_cat_for_map)))
                    
                    m.plot(x, y, '.', markersize=s, markeredgewidth=0., linewidth=0.2, color = c, alpha = alpha)
                    
                    if M > max(M_cat_for_map) - 1.:
                        x_text, y_text = m(lon+0.005, lat+0.005)
                        plt.text(x_text,y_text,str(int(Yr))+','+str(int(nb)),fontsize = 4)
                                            
                                            
                m.drawcoastlines(linewidth=0.1)
                #m.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
                plt.gca().set_title('Earthquake catalog - Complete period')
                plt.savefig(str(Run_name) + '/analysis/figures/catalogue/'+'earthquake_map_'+str(model)+'.png' , dpi = 180, transparent=True)
                #plt.show()
                plt.close()
                
                
                ###### map plot of the catalog #####
                
                catalog_for_map = np.genfromtxt(str(Run_name) + '/analysis/figures/catalogue/catalog_for_map_'+str(model)+'.txt',dtype=[('f8'),('f8'),('f8'),('f8'),('f8')])
                yr_cat_for_map =  list(map(lambda i : round(float(catalog_for_map[i][0])), range(len(catalog_for_map))))
                M_cat_for_map = list( map(lambda i : float(catalog_for_map[i][1]), range(len(catalog_for_map))))
                lon_cat_for_map =  list(map(lambda i : float(catalog_for_map[i][2]), range(len(catalog_for_map))))
                lat_cat_for_map = list( map(lambda i : float(catalog_for_map[i][3]), range(len(catalog_for_map))))
                nb_time_picked = list( map(lambda i : float(catalog_for_map[i][4]), range(len(catalog_for_map))))
                
                m = Basemap(projection='mill',
                              llcrnrlon = llcrnrlon, 
                              llcrnrlat = llcrnrlat, 
                              urcrnrlon = urcrnrlon, 
                              urcrnrlat = urcrnrlat,resolution='l')
                
                for (Yr,M,lon,lat,nb) in zip(yr_cat_for_map,M_cat_for_map,lon_cat_for_map,lat_cat_for_map,nb_time_picked):
                    x, y = m(lon, lat)
                    # color scale
                    if M< 5. :
                        c='b'
                        s = 3
                    elif M< 6.0:
                        c='g'
                        s = 5
                    elif M< 7.:
                        c='orange'
                        s= 7
                    else:
                        c='r'
                        s = 10
                        
#                    if nb >= 2*nb_inter/3 :
#                        alpha = 1.
#                    if nb <= 2*nb_inter/3 :
#                        alpha = 0.8
#                    if nb <= 1*nb_inter/2 :
#                        alpha = 0.3
#                    if nb < 1*nb_inter/3 :
#                        alpha = 0.001
#                    
                    #c = plt.cm.jet((M-min(M_cat_for_map))/(max(M_cat_for_map)-min(M_cat_for_map)))
                    
                    if nb >= float(nb_inter)/2. : #plot the earthquake if it's piccked more than 50% of the iterations
                        m.plot(x, y, '.', markersize=s, markeredgewidth=0., linewidth=0., color = c, alpha = 0.6)
                    
                plt.savefig(str(Run_name) + '/analysis/figures/catalogue/'+'earthquake_map_just_EQ_'+str(model)+'.png' , dpi = 180, transparent=True)
                #plt.show()
                plt.close() 
            
        index_model += 1
    
        
        #write the mean rate of eq in the catalog
        file_eq_rate_cat = open(str(Run_name) + '/analysis/figures/catalogue/cat_cumulative_eq_rate.txt','w')
        file_eq_rate_cat.write('cumulative annual eq rate in the catalog'+'\n')
        for model_i, rate_i in zip(Model_list,catalog_cum_rate):
            file_eq_rate_cat.write(model_i)
            for r in rate_i :
                file_eq_rate_cat.write('\t'+str(r))
            file_eq_rate_cat.write('\n') 
        file_eq_rate_cat.close() 
        
    else: #read the files
        file_eq_rate_cat = str(Run_name) + '/analysis/figures/catalogue/cat_cumulative_eq_rate.txt'
        with open(file_eq_rate_cat) as f:#finds where to start reading
            lines_cat = f.readlines()
        Model_list=[]
        catalog_cum_rate=[]
        for line in lines_cat[1:]:
            cat_rates_i=line.split('\t')
            catalog_cum_rate.append([float(i) for i in cat_rates_i[1:]])
            
        seismological_moment_rate=[]
        for rates in catalog_cum_rate:
            rates_inc = []
            for i in range(len(rates)-1):
                rates_inc.append(rates[i]-rates[i+1])
            rates_inc.append(rates[-1])
            seismological_moment_rate_i=0
            for mag_i,rate_i in zip(bining_in_mag,rates_inc) :
                M0 = 10. ** (1.5 * mag_i + 9.1)  #calculation of the moment
                rate_M0 = M0 * rate_i
                seismological_moment_rate_i += rate_M0
        
            seismological_moment_rate.append(seismological_moment_rate_i)  
        
        
        
    return seismological_moment_rate, catalog_cum_rate#,yr_cat_for_map,M_cat_for_map,lon_cat_for_map,lat_cat_for_map
