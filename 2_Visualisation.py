# -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.1

Visualisation of the model you created using SHERIFS

@author: Thomas Chartier
"""

from __future__ import print_function
import numpy as np
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import time
#from glob import glob

mpl.interactive(False)
import sys
from lib.utils import sap


path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/lib'
sys.path.append(path_lib)
path_f = path_lib + '/file_writer'
sys.path.append(path_f)
path_f = path_lib + '/hm_visual'
sys.path.append(path_f)

from OQ_job_Creator import OQ_job_Creator
import Extract_data
import Plot_mfd
import Participation_rates
import plt_catalog
import plot_FtF_GIF
import Sampling_analysis
import moment_rate
import slip_rate_rep
import sherifs_report


def checking_the_input(input_file):
#    def __init__(self):
    debut = time.time()
    
    
    
    # Load the input file
    lines = open(input_file,'r').readlines()
    lines = [line.rstrip('\n') for line in lines]
    for line in lines:
        if "Run_Name" in line :
            Run_name = line.split(':')[1].replace(' ','')
        if "File_geom" in line :
            File_geom = line.split(':')[1].replace(' ','')
        if "File_prop" in line :
            FileName_Prop = line.split(':')[1].replace(' ','')
        if "File_bg" in line :
            File_bg = line.split(':')[1].replace(' ','')
        if "file_prop_bg" in line :
            file_prop_bg = line.split(':')[1].replace(' ','')
    #Example
    
    #Run_name = 'Example'
    #File_geom = 'data/Example/Faults_geometry.txt'
    #FileName_Prop = 'data/Example/Faults_properties.txt'
    #File_bg = 'data/Example/Background_geometry.txt'
    completness_file = 'data/Example/completness.txt'
    catalog_file = 'data/Example/catalog.txt'
    file_faults_data = 'data/Example/EQ_on_faults_data.txt'
#    sub_area_file = 'data/sub_areas.txt'
    sub_area_file = 'data/CHN/mfd_area.geojson'
    end_year_of_catalog = 2009

    llcrnrlon=21.6
    llcrnrlat=38.12
    urcrnrlon=22.5
    urcrnrlat=38.5
    
    xmax = 6.8
    ymin = 0.00001
    ymax = 5.
    
    
    

    print(Run_name)
    # booleans to know which figure to generat
    do_catalog = False
    #If it's the first time you run it, you have to do the catalog! Then you can save time and not run it again
    
    plot_mfd = True
    
    plot_mfd_detailled = True
    
    plot_Mmax = True
    
    plot_as_rep = True
    
    plot_rup_freq = True
    
    plot_sr_use = True
    
    plot_moment_rate = True
    
    visual_FtF = True   # !!! can take a very long time if there are a lot of FtF
    
    test_fit_to_data = False
            
    OQ_job = OQ_job_Creator(Run_name) # ask the info about the run and create the job.ini file
    Mmin = OQ_job.Mmin
    
    #plot mfd parameters
    xmin = float(Mmin)
    xmax = xmax
    ymin = ymin
    ymax = ymax
    
    '''########################
    # Extraction of the data
    ########################'''
    
    time_i = time.time()
    (mega_MFD, df_mega_MFD, scenarios_names_list, ScL_complet_list, ScL_list, Model_list,BG_hyp_list,
     dimension_used_list,faults_name_list,sample_list,b_value_list,MFD_type_list,m_Mmax,
     mega_bining_in_mag,a_s_model,b_sample,sm_sample,Mt_sample,sources_Lengths,sources_Areas,logictree) = Extract_data.extract(Run_name)
    print('\nTime to extract the data : ' + str(round(time.time() - time_i,2)) +' s.\n')
    
    
    '''########################
    # plot the Mmax distriution
    ########################'''
    if not os.path.exists(str(Run_name) + '/analysis/figures/Mmax'):
        os.makedirs(str(Run_name) + '/analysis/figures/Mmax')
        
    if plot_Mmax == True :
        plt.hist(m_Mmax,int(round(max(m_Mmax) - min(m_Mmax),1) * 10. + 1.))
        plt.savefig(str(Run_name) + '/analysis/figures/Mmax/Hist_Mmax.png',dpi = 180)
        plt.close()
        
        plt.hist(sources_Lengths,20)
        plt.xlabel('Length of the sources')
        plt.savefig(str(Run_name) + '/analysis/figures/Mmax/Hist_Lengths_all_models.png',dpi = 180)
        plt.close()
        
        plt.hist(sources_Areas,20)
        plt.xlabel('Area of the sources')
        plt.savefig(str(Run_name) + '/analysis/figures/Mmax/Hist_Areas_all_models.png',dpi = 180)
        plt.close()
     
    if xmax < max(m_Mmax):
        xmax = round(max(m_Mmax)+0.2,1)
    #bining_in_mag = np.linspace(xmin,xmax,int((xmax-xmin)*10)+2)
    bining_in_mag = [round(i,1) for i in np.arange(xmin,xmax+0.1,0.1)]
    
                
    if plot_Mmax == True :
        #plot the distribution of Mmax of each rupture for each fault for the mean branch
        i_model=0
        for Model in Model_list:
            for scenario in scenarios_names_list:
                log_Mmax_file = (str(Run_name)  + '/' + str(Model) + '/Log/' + 'Mmax_sample_'+ScL_complet_list[0]+'_sc_' +  str(scenario)+ '_1.txt')
                log_file = np.genfromtxt( log_Mmax_file,dtype=[('U10000'),('f8'),('U100'),('f8'),('f8'),('f8')],delimiter = '\t')
                sources_names =  list(map(lambda i : log_file[i][0], range(len(log_file))))
                sources_Mmax =  list(map(lambda i : log_file[i][5], range(len(log_file))))
                sources_names = [str(i) for i in sources_names]
                for fault in faults_name_list[0]:
                    Mmax_ruptures_fault=[]
                    for source,Mmax_i in zip(sources_names,sources_Mmax):
                        if fault == source or "['"+fault+"']" in source:
                            Mmax_ruptures_fault.append(Mmax_i)
                    n=[]
                    for mag in bining_in_mag:
                        n.append(sum(i > mag for i in Mmax_ruptures_fault))
                    plt.plot(bining_in_mag,n)
                    plt.title(fault+'   '+Model+'   '+scenario)
                    plt.xlabel('Mw')
                    plt.ylabel('Number of rupture with Mmax >= Mw')
                    
                    if not os.path.exists(str(Run_name) + '/analysis/figures/Mmax/'+Model):
                        os.makedirs(str(Run_name) + '/analysis/figures/Mmax/'+Model)
                    if not os.path.exists(str(Run_name) + '/analysis/figures/Mmax/'+Model+'/'+scenario):
                        os.makedirs(str(Run_name) + '/analysis/figures/Mmax/'+Model+'/'+scenario)
                    plt.savefig(str(Run_name) + '/analysis/figures/Mmax/'+Model+'/'+scenario+'/'+fault+'.png',dpi=80)
                    plt.close()
            i_model+=1
            
        
        
    '''########################
    # plot  Mmax vs NMS
    ########################'''
    plt.scatter(m_Mmax,a_s_model)
    plt.savefig(str(Run_name) + '/analysis/figures/Mmax/Mmax_vs_NMS.png',dpi = 180)
    plt.title('Mmax vs NMS')
    plt.close()
    

    '''#############################
    ###############################
    #     Extraction of
    #   the catalog
    ###############################
    ##############################'''
    
    time_i = time.time()
        
    if not os.path.exists(str(Run_name) + '/analysis/figures/mfd'):
        os.makedirs(str(Run_name) + '/analysis/figures/mfd')
        
        
    if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches'):
        os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches')
        
    nb_inter = 25
    try :
        (seismological_moment_rate, catalog_cum_rate,
         yr_cat_for_map,M_cat_for_map,lon_cat_for_map,lat_cat_for_map) = plt_catalog.plt_catalog(do_catalog,Model_list,File_bg,
                                                            catalog_file,Run_name,
                                                            xmin,xmax,ymin,ymax,
                                                            llcrnrlon,llcrnrlat,
                                                            urcrnrlon,urcrnrlat,
                                                            completness_file,nb_inter,bining_in_mag,end_year_of_catalog,
                                                            sub_area_file)
    except:
        seismological_moment_rate, catalog_cum_rate = plt_catalog.plt_catalog(do_catalog,Model_list,File_bg,
                                                            catalog_file,Run_name,
                                                            xmin,xmax,ymin,ymax,
                                                            llcrnrlon,llcrnrlat,
                                                            urcrnrlon,urcrnrlat,
                                                            completness_file,nb_inter,bining_in_mag,end_year_of_catalog,
                                                            sub_area_file)
    print('\nTime to plot the catalog : ' + str(round(time.time() - time_i,2)) +' s.\n')
    
    '''################################################
    ###################################################
    ### work on and plot the MFD of the sources in the model ###
    ###################################################
    ################################################'''
    
    
    time_i = time.time()
    (total_list_ScL,total_list_dimension_used,
     geologic_moment_rate,geologic_moment_rate_no_as,
     total_list_scenario_name,total_list_MFD_type,
     mega_mfd_cummulative,total_list_model,total_list_sample
     ,total_list_BG_hyp) = Plot_mfd.plt_mfd(Run_name,mega_MFD,scenarios_names_list,ScL_complet_list,
                                     ScL_list, Model_list,BG_hyp_list,dimension_used_list,faults_name_list,sample_list,b_value_list,
                                     MFD_type_list,m_Mmax,
                                     mega_bining_in_mag,a_s_model,b_sample,sm_sample,Mt_sample,plot_mfd,plot_as_rep,plot_Mmax,
                                     xmin,xmax,ymin,ymax, catalog_cum_rate,plot_mfd_detailled,bining_in_mag)
    print('\nTime to plot the MFDs : ' + str(round(time.time() - time_i,2)) +' s.\n')
       
#
#        print('size variables after mfd')
#        for var, obj in locals().items():
#            print(var, sys.getsizeof(obj))
    
    '''####################################
    # use of the slip rate per fault
    #######################################'''

    if plot_sr_use == True :
        time_i = time.time()
        slip_rate_rep.sr_rate(Run_name,scenarios_names_list,
        mega_MFD,Model_list,MFD_type_list,sub_area_file,File_geom,FileName_Prop)
    
        print('\nTime to see how the slip rate in distributed : ' + str(time.time() - time_i) +' s.\n')
        


    '''####################################
    #######################################
    #  FAULTS RUPTURES
    # plot return period of a fault of a mutli faults rupture to occure for
    # each set of rupture scenario used
    #
    #######################################
    #######################################'''

    
    '''####################################
    # rupture rate fault by fault
    #######################################'''
    time_i = time.time()
        
    if plot_rup_freq == True :
        Participation_rates.plt_EQ_rates(Run_name,mega_MFD,df_mega_MFD,scenarios_names_list,ScL_complet_list,
                                     ScL_list, Model_list,BG_hyp_list,dimension_used_list,faults_name_list,sample_list,b_value_list,MFD_type_list,m_Mmax,
                                     mega_bining_in_mag,a_s_model,b_sample,sm_sample,Mt_sample,plot_mfd,plot_as_rep,plot_Mmax,xmin,xmax,ymin,ymax,
                                     file_faults_data,File_bg,File_geom,sub_area_file,FileName_Prop )
#            #extract the faults data

    del mega_MFD,df_mega_MFD
    print('\nTime to plot the rupture rates of each faults : ' + str(time.time() - time_i) +' s.\n')
                

        
    '''####################################
    # impact of the sampling on the fit
    #######################################'''
    
    time_i = time.time()
    if test_fit_to_data == True :
        Sampling_analysis.sampling_analysis(Run_name,Model_list,m_Mmax,b_sample,a_s_model,
                                            mega_mfd_cummulative,catalog_cum_rate,xmin,xmax,ymin,ymax,
                                            total_list_model,bining_in_mag,total_list_MFD_type,
                                            total_list_scenario_name,file_faults_data,total_list_sample,
                                            total_list_BG_hyp)
    print('\nTime to do the sampling analysis : ' + str(round(time.time() - time_i,2)) +' s.\n')



    '''####################################
    # visualisation of the FtF ruptures
    #######################################'''
    
    time_i = time.time()
    plot_FtF_GIF.map_faults(Run_name,Model_list,scenarios_names_list,
       ScL_complet_list, BG_hyp_list,
       sample_list,b_value_list,MFD_type_list,
       llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,
       File_bg,FileName_Prop,File_geom ,plot_sr_use,visual_FtF,sub_area_file)
    print('\nTime to plot the different FtF ruptures : ' + str(round(time.time() - time_i,2)) +' s.\n')
       
    
    '''#############################
    ###############################
    #     comparison of the moment rate
    #   geologic, geodetic, sismologic
    ###############################
    ##############################'''
    moment_rate.moment_rate(Run_name,plot_moment_rate,geologic_moment_rate_no_as,geologic_moment_rate,
            seismological_moment_rate,scenarios_names_list,total_list_scenario_name,
            MFD_type_list,total_list_MFD_type)
            
    '''#############################
    ###############################
    #     CREATE PFD REPORT
    ###############################
    ##############################'''
    
    sherifs_report.create(Run_name,logictree)
    
    
    
    # Print time
    fin = time.time()-debut
    days = int(fin / 24. / 60. / 60.)
    hours = int((fin - days * 24. * 60. * 60.) / 60. / 60.)
    minutes = int((fin - days * 24. * 60. * 60. - hours* 60. * 60. ) / 60.)
    seconds = (fin - days * 24. * 60. * 60. - hours* 60. * 60.  - minutes * 60.)
    print("The calculation took: " + str(days) + ' days, ' + str(hours) + ' hours, ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds.')

        
#if __name__=="__main__":
#   app = checking_the_input()
   

def main(argv):
   """ Run SHERIFS analysis"""

   p = sap.Script(checking_the_input)
   p.arg(name='input_file', help='.txt file with the information concerning the run.')

   if len(argv) < 1:
       print(p.help())
   else:
       p.callfunc()


if __name__ == "__main__":
   main(sys.argv[1:])
