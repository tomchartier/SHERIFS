# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

This code is the main code for building OQ entries.
It creates the job.ini and launch the other codes creating logic trees (Sources and GMPEs) and sources models.

@author: Thomas Chartier
"""

import os
import numpy as np
#import tkinter as tk
#from tkinter import ttk


''''''''''''''''''''''''''''''''''''''''''
'''        Reading input data          '''
''''''''''''''''''''''''''''''''''''''''''
''' input data are located in the run.info and job.ini file'''

class OQ_job_Creator():    
    def __init__(self,Run_Name):  
        self.Run_Name = Run_Name
        
        self.input_data()
        self.write_ini()
        
    def input_data(self):
        
        info_run_file = 'input/'+str(self.Run_Name)+'/run.info'
        if os.path.exists(info_run_file): #the file exist, read it
            
            read_info_file = open(info_run_file,'r')
            lines_of_the_file = read_info_file.readlines()
            self.option_map = lines_of_the_file[1].split('\t')[1]
            lines_of_the_file = [x.strip('L\n') for x in lines_of_the_file] 
            lines_of_the_file = [x.strip('\r\n') for x in lines_of_the_file] 
            lines_of_the_file = [x.strip('\n') for x in lines_of_the_file] 
            if self.option_map == "yes":
                self.ll_lon = lines_of_the_file[2].split('\t')[1]
                self.ll_lat = lines_of_the_file[2].split('\t')[2]
                self.lr_lon = lines_of_the_file[2].split('\t')[3]
                self.lr_lat = lines_of_the_file[2].split('\t')[4]
                self.ur_lon = lines_of_the_file[2].split('\t')[5]
                self.ur_lat = lines_of_the_file[2].split('\t')[6]
                self.ul_lon = lines_of_the_file[2].split('\t')[7]
                self.ul_lat = lines_of_the_file[2].split('\t')[8]
                self.region_grid_spacing = lines_of_the_file[3].split('\t')[1]          
            else:    
                self.Longitude = lines_of_the_file[2].split('\t')[1]
                self.Latitude = lines_of_the_file[3].split('\t')[1]
            self.Vs30 = lines_of_the_file[4].split('\t')[1]
            self.Z1000 = lines_of_the_file[5].split('\t')[1]
            self.Z2500 = lines_of_the_file[6].split('\t')[1]
            self.nb_LT_samp = lines_of_the_file[7].split('\t')[1]
            self.rup_mesh = lines_of_the_file[8].split('\t')[1]
            self.source_discr = lines_of_the_file[9].split('\t')[1]
            self.investigation_time = lines_of_the_file[10].split('\t')[1]
            self.poes = lines_of_the_file[11].split('\t')[1]
            self.trunc_lvl = lines_of_the_file[12].split('\t')[1]
            self.max_dist = lines_of_the_file[13].split('\t')[1]
            self.nb_sample = int(lines_of_the_file[14].split('\t')[1])
            self.Mmin = float(lines_of_the_file[15].split('\t')[1])
            self.seed = int(lines_of_the_file[16].split('\t')[1])
            self.sr_correl = lines_of_the_file[17].split('\t')[1]
            if self.sr_correl == 'yes' or self.sr_correl == 'Yes' or self.sr_correl == 'True' or self.sr_correl == 'true'   :
                self.sr_correl = True
            else :
                self.sr_correl = False
            self.size_of_increment = float(lines_of_the_file[18].split('\t')[1]) 
            self.Mmax_range = [float(lines_of_the_file[19].split('\t')[1]),float(lines_of_the_file[19].split('\t')[2])]
            

            self.intensity_out = []
            self.min_int_out = []
            self.max_int_out = []
            self.nb_int_out = []            
            for line in lines_of_the_file[20::]:
                self.intensity_out.append(line.split('\t')[1])
                self.min_int_out.append(float(line.split('\t')[2]))
                self.max_int_out.append(float(line.split('\t')[3]))
                self.nb_int_out.append(int(line.split('\t')[4]))
                
            
            read_info_file.close()
            
            
        else : # the file doesn't exist, ask the information about the run
            print("ERROR : File run.info not found! \ Please refer to the user manual")
            exit()
#            self.w_ini = tk.Tk()
#            self.w_ini.grid()
#            self.w_ini.title('Entrer calculation information')
#
#            label = tk.Label(self.w_ini, text="Site lon")
#            label.grid(column=0,row=0)
#            self.Longitude = tk.Entry(self.w_ini)
#            self.Longitude.grid(column=0,row=1)
#            self.Longitude.insert(tk.INSERT,22.090575)
#
#            label = tk.Label(self.w_ini, text="Site lat")
#            label.grid(column=1,row=0)
#            self.Latitude = tk.Entry(self.w_ini)
#            self.Latitude.grid(column=1,row=1)
#            self.Latitude.insert(tk.INSERT,38.250372)
#
#            #self.option_map = 'no'
#
#            row_i = 3
#
#            row_i +=1
#
#
#            label = tk.Label(self.w_ini, text="\n SHERIFS Parameters")
#            label.grid(column=0,row=row_i)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Mmin")
#            label.grid(column=0,row=row_i)
#            self.Mmin = tk.Entry(self.w_ini)
#            self.Mmin.grid(column=1,row=row_i)
#            self.Mmin.insert(tk.INSERT,5.0)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Slip-rate increment size (mm/yr)")
#            label.grid(column=0,row=row_i)
#            self.size_of_increment = tk.Entry(self.w_ini)
#            self.size_of_increment.grid(column=1,row=row_i)
#            self.size_of_increment.insert(tk.INSERT,0.005)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Mmax range")
#            label.grid(column=0,row=row_i)
#            self.Mmaxmin = tk.Entry(self.w_ini)
#            self.Mmaxmin.grid(column=1,row=row_i)
#            self.Mmaxmin.insert(tk.INSERT,0.)
#            self.Mmaxmax = tk.Entry(self.w_ini)
#            self.Mmaxmax.grid(column=2,row=row_i)
#            self.Mmaxmax.insert(tk.INSERT,10.)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Nb sample (sr,b_value,Mmax)")
#            label.grid(column=0,row=row_i)
#            self.nb_sample = tk.Entry(self.w_ini)
#            self.nb_sample.grid(column=1,row=row_i)
#            self.nb_sample.insert(tk.INSERT,20)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Correlation of slip-rates")
#            label.grid(column=0,row=row_i)
#            self.sr_correl = tk.Entry(self.w_ini)
#            self.sr_correl.grid(column=1,row=row_i)
#            self.sr_correl.insert(tk.INSERT,'yes')
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Seed for random sampling")
#            label.grid(column=0,row=row_i)
#            self.seed = tk.Entry(self.w_ini)
#            self.seed.grid(column=1,row=row_i)
#            self.seed.insert(tk.INSERT,805)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="Option map")
#            label.grid(column=0,row=row_i)
#            self.option_map = tk.Entry(self.w_ini)
#            self.option_map.grid(column=1,row=row_i)
#            self.option_map.insert(tk.INSERT,'no')
#            row_i +=1
#
#
#            label = tk.Label(self.w_ini, text="\n Site Conditions")
#            label.grid(column=0,row=row_i)
#            row_i +=1
#            label = tk.Label(self.w_ini, text="VS30")
#            label.grid(column=0,row=row_i)
#            self.Vs30 = tk.Entry(self.w_ini)
#            self.Vs30.grid(column=1,row=row_i)
#            self.Vs30.insert(tk.INSERT,800.)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="depth_to_1pt0km_per_sec")
#            label.grid(column=0,row=row_i)
#            self.Z1000 = tk.Entry(self.w_ini)
#            self.Z1000.grid(column=1,row=row_i)
#            self.Z1000.insert(tk.INSERT,100.0)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="reference_depth_to_2pt5km_per_sec")
#            label.grid(column=0,row=row_i)
#            self.Z2500 = tk.Entry(self.w_ini)
#            self.Z2500.grid(column=1,row=row_i)
#            self.Z2500.insert(tk.INSERT,5.0)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="\n OpenQuake Parameters")
#            label.grid(column=0,row=row_i)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="rupture mesh spacing( km)")
#            label.grid(column=0,row=row_i)
#            self.rup_mesh = tk.Entry(self.w_ini)
#            self.rup_mesh.grid(column=1,row=row_i)
#            self.rup_mesh.insert(tk.INSERT,0.5)
#
#            label = tk.Label(self.w_ini, text="area source discretization (km)")
#            label.grid(column=2,row=row_i)
#            self.source_discr = tk.Entry(self.w_ini)
#            self.source_discr.grid(column=3,row=row_i)
#            self.source_discr.insert(tk.INSERT,5.0)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="investigation time")
#            label.grid(column=0,row=row_i)
#            self.investigation_time = tk.Entry(self.w_ini)
#            self.investigation_time.grid(column=1,row=row_i)
#            self.investigation_time.insert(tk.INSERT,50)
#
#            label = tk.Label(self.w_ini, text="max distance")
#            label.grid(column=2,row=row_i)
#            self.max_dist = tk.Entry(self.w_ini)
#            self.max_dist.grid(column=3,row=row_i)
#            self.max_dist.insert(tk.INSERT,200)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="truncation level (sigma)")
#            label.grid(column=0,row=row_i)
#            self.trunc_lvl = tk.Entry(self.w_ini)
#            self.trunc_lvl.grid(column=1,row=row_i)
#            self.trunc_lvl.insert(tk.INSERT,3)
#            row_i +=1
#
#            label = tk.Label(self.w_ini, text="LT sample for openquake")
#            label.grid(column=0,row=row_i)
#            self.nb_LT_samp =tk. Entry(self.w_ini)
#            self.nb_LT_samp.grid(column=1,row=row_i)
#            self.nb_LT_samp.insert(tk.INSERT,0)
#            row_i +=1
#
#
#            label = tk.Label(self.w_ini, text="Target PoE")
#            label.grid(column=0,row=row_i)
#            self.poes = tk.Entry(self.w_ini)
#            self.poes.grid(column=1,row=row_i)
#            self.poes.insert(tk.INSERT, 0.1)
#            row_i +=1
#
#
#
#            label = tk.Label(self.w_ini, text="\n Intensity levels")
#            label.grid(column=0,row=row_i)
#            label = tk.Label(self.w_ini, text="\n min (g)")
#            label.grid(column=1,row=row_i)
#            label = tk.Label(self.w_ini, text="\n max (g)")
#            label.grid(column=2,row=row_i)
#            label = tk.Label(self.w_ini, text="\n nb_points")
#            label.grid(column=3,row=row_i)
#            self.intensity = []
#            self.min_int = []
#            self.max_int = []
#            self.nb_int = []
#            for i in range(4):
#                if i == 0:
#                    self.intensity_i = tk.Entry(self.w_ini)
#                    self.intensity_i.grid(column=0,row=row_i+1)
#                    self.intensity_i.insert(tk.INSERT,'PGA')
#                    self.intensity.append(self.intensity_i)
#
#                    self.min_int_i = tk.Entry(self.w_ini)
#                    self.min_int_i.grid(column=1,row=row_i+1)
#                    self.min_int_i.insert(tk.INSERT,0.08)
#                    self.min_int.append(self.min_int_i)
#
#                    self.max_int_i = tk.Entry(self.w_ini)
#                    self.max_int_i.grid(column=2,row=row_i+1)
#                    self.max_int_i.insert(tk.INSERT,1.5)
#                    self.max_int.append(self.max_int_i)
#
#
#                    self.nb_int_i = tk.Entry(self.w_ini)
#                    self.nb_int_i.grid(column=3,row=row_i+1)
#                    self.nb_int_i.insert(tk.INSERT,10)
#                    self.nb_int.append(self.nb_int_i)
#
#                    row_i +=1
#
#                if i == 1:
#                    self.intensity_i = tk.Entry(self.w_ini)
#                    self.intensity_i.grid(column=0,row=row_i+1)
#                    self.intensity_i.insert(tk.INSERT,'SA(0.1)')
#                    self.intensity.append(self.intensity_i)
#
#                    self.min_int_i = tk.Entry(self.w_ini)
#                    self.min_int_i.grid(column=1,row=row_i+1)
#                    self.min_int_i.insert(tk.INSERT,0.08)
#                    self.min_int.append(self.min_int_i)
#
#                    self.max_int_i = tk.Entry(self.w_ini)
#                    self.max_int_i.grid(column=2,row=row_i+1)
#                    self.max_int_i.insert(tk.INSERT,1.5)
#                    self.max_int.append(self.max_int_i)
#
#
#                    self.nb_int_i = tk.Entry(self.w_ini)
#                    self.nb_int_i.grid(column=3,row=row_i+1)
#                    self.nb_int_i.insert(tk.INSERT,10)
#                    self.nb_int.append(self.nb_int_i)
#
#                    row_i +=1
#
#                else:
#                    self.intensity_i = tk.Entry(self.w_ini)
#                    self.intensity_i.grid(column=0,row=row_i+1)
#                    self.intensity_i.insert(tk.INSERT,'SA(0.1)')
#                    self.intensity.append(self.intensity_i)
#
#                    self.min_int_i = tk.Entry(self.w_ini)
#                    self.min_int_i.grid(column=1,row=row_i+1)
#                    self.min_int_i.insert(tk.INSERT,0.08)
#                    self.min_int.append(self.min_int_i)
#
#                    self.max_int_i = tk.Entry(self.w_ini)
#                    self.max_int_i.grid(column=2,row=row_i+1)
#                    self.max_int_i.insert(tk.INSERT,1.5)
#                    self.max_int.append(self.max_int_i)
#
#
#                    self.nb_int_i = tk.Entry(self.w_ini)
#                    self.nb_int_i.grid(column=3,row=row_i+1)
#                    self.nb_int_i.insert(tk.INSERT,10)
#                    self.nb_int.append(self.nb_int_i)
#
#                    row_i +=1
#
#            bou_ok = tk.Button(self.w_ini, text=u'OK', command = self.action)
#            bou_ok.grid(column=4,row= row_i+1)
#
#            self.w_ini.mainloop()
#
#            if self.option_map == 'yes': #if option map is yes, open new window
#                self.w_o_map = Tk()
#                self.w_o_map.grid()
#                self.w_o_map.title('Information about the map')
#
#                row_i = 1
#
#                label = tk.Label(self.w_o_map, text="Longitude")
#                label.grid(column=1,row=row_i)
#                label = tk.Label(self.w_o_map, text="Latitude")
#                label.grid(column=2,row=row_i)
#                row_i +=1
#
#                label = tk.Label(self.w_o_map, text="Lower Left corner")
#                label.grid(column=0,row=row_i)
#                self.ll_lon = tk.Entry(self.w_o_map)
#                self.ll_lon.grid(column=1,row=row_i)
#                self.ll_lon.insert(tk.INSERT,10.)
#                self.ll_lat = tk.Entry(self.w_o_map)
#                self.ll_lat.grid(column=2,row=row_i)
#                self.ll_lat.insert(tk.INSERT,43.)
#                row_i +=1
#
#                label = tk.Label(self.w_o_map, text="Lower Right corner")
#                label.grid(column=0,row=row_i)
#                self.lr_lon = tk.Entry(self.w_o_map)
#                self.lr_lon.grid(column=1,row=row_i)
#                self.lr_lon.insert(tk.INSERT,12.)
#                self.lr_lat = Entry(self.w_o_map)
#                self.lr_lat.grid(column=2,row=row_i)
#                self.lr_lat.insert(tk.INSERT,43.)
#                row_i +=1
#
#                label = tk.Label(self.w_o_map, text="Upper Right corner")
#                label.grid(column=0,row=row_i)
#                self.ur_lon = tk.Entry(self.w_o_map)
#                self.ur_lon.grid(column=1,row=row_i)
#                self.ur_lon.insert(tk.INSERT,12.)
#                self.ur_lat = tk.Entry(self.w_o_map)
#                self.ur_lat.grid(column=2,row=row_i)
#                self.ur_lat.insert(tk.INSERT,46.)
#                row_i +=1
#
#                label = tk.Label(self.w_o_map, text="Upper Left corner")
#                label.grid(column=0,row=row_i)
#                self.ul_lon = tk.Entry(self.w_o_map)
#                self.ul_lon.grid(column=1,row=row_i)
#                self.ul_lon.insert(tk.INSERT,10.)
#                self.ul_lat = tk.Entry(self.w_o_map)
#                self.ul_lat.grid(column=2,row=row_i)
#                self.ul_lat.insert(tk.INSERT,46.)
#                row_i +=1
#
#                label = tk.Label(self.w_o_map, text="Region grid spacing (km)")
#
#                label.grid(column=0,row=row_i)
#                self.region_grid_spacing = tk.Entry(self.w_o_map)
#                self.region_grid_spacing.grid(column=1,row=row_i)
#                self.region_grid_spacing.insert(tk.INSERT,5.)
#                row_i +=1
#
#                bou_ok_map = tk.Button(self.w_o_map, text=u'OK', command = self.action_map)
#                bou_ok_map.grid(column=4,row= row_i+1)
#
#                self.w_o_map.mainloop()
#
#
#            self.write_info()
            
        
    def write_ini(self):
        # write the job.ini for Openquake
        jobfile=open(str(self.Run_Name)+'/job.ini','w')
        Ligne='[general]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='description = \'' + str(self.Run_Name) + '\'\n'
        jobfile.write(Ligne)
        Ligne='calculation_mode = classical\n'
        jobfile.write(Ligne)
        Ligne='random_seed = '+ str(self.seed) + '\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='[geometry]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        if self.option_map == 'yes':        
            jobfile.write(Ligne)
            Ligne=('region = ' + str(self.ll_lon) + ' ' + str(self.ll_lat) 
            + ', ' + str(self.lr_lon) + ' '+ str(self.lr_lat) 
            + ', ' + str(self.ur_lon) + ' ' + str(self.ur_lat) 
            + ', ' + str(self.ul_lon) + ' ' + str(self.ul_lat) + '\n')
            jobfile.write(Ligne)
            Ligne='region_grid_spacing = ' + str(self.region_grid_spacing) + '\n'
        else:
            jobfile.write(Ligne)
            Ligne='sites = ' + str(self.Longitude)[:-1] + ' ' + str(self.Latitude) +'\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='[logic_tree]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='number_of_logic_tree_samples = '+str(self.nb_LT_samp)+'\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='[erf]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='rupture_mesh_spacing ='+str(self.rup_mesh)+'\n'
        jobfile.write(Ligne)
        Ligne='width_of_mfd_bin = 0.1\n'
        jobfile.write(Ligne)
        Ligne='area_source_discretization = '+str(self.source_discr)+'\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='[site_params]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='reference_vs30_type = measured\n'
        jobfile.write(Ligne)
        Ligne='reference_vs30_value = '+str(self.Vs30)+'\n'
        jobfile.write(Ligne)
        Ligne='reference_depth_to_2pt5km_per_sec = '+str(self.Z2500)+'\n'
        jobfile.write(Ligne)
        Ligne='reference_depth_to_1pt0km_per_sec = '+str(self.Z1000)+'\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='[calculation]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='source_model_logic_tree_file = Sources_Logic_tree.xml\n'
        jobfile.write(Ligne)
        Ligne='gsim_logic_tree_file = GMPE_Logic_tree.xml\n'
        jobfile.write(Ligne)
        Ligne='investigation_time = '+str(self.investigation_time)+'\n'
        jobfile.write(Ligne)
        
        Ligne='intensity_measure_types_and_levels = {'
        jobfile.write(Ligne)
        intensities_done=[]  #check to see if we don't put twice the same intensity
        for i in range(len(self.intensity_out)) :
        #'"PGA": [0.284, 0.397, 0.556, 0.778, 1.09, 1.52, 2.13], "SA(0.1)": [0.0527, 0.0738, 0.103, 0.145, 0.203, 0.284, 0.397, 0.556, 0.778, 1.09, 1.52, 2.13], "SA(1.0)": [0.0738, 0.103, 0.145, 0.203, 0.284, 0.397, 0.556, 0.778, 1.09, 1.52, 2.13]
            
            Ligne_1 = '"'+str(self.intensity_out[i])+'": [' 
            
            array_ints = np.logspace(np.log10(self.min_int_out[i]),np.log10(self.max_int_out[i]),self.nb_int_out[i])
            Ligne_2 = ''
            for values in array_ints:
                if len(Ligne_2) > 1 :
                    Ligne_2 += ', '
                Ligne_2 += str(values)
            
            if not self.intensity_out[i] in intensities_done :
                if len(intensities_done) != 0 :
                  jobfile.write(', ')  
                jobfile.write(Ligne_1 + Ligne_2+']')
            intensities_done.append(self.intensity_out[i])
        Ligne='}\n'
        jobfile.write(Ligne)
        
        Ligne='truncation_level = '+str(self.trunc_lvl)+'\n'
        jobfile.write(Ligne)
        Ligne='maximum_distance = '+str(self.max_dist)+'\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='[output]\n'
        jobfile.write(Ligne)
        Ligne='\n'
        jobfile.write(Ligne)
        Ligne='export_dir =  ./results\n'
        jobfile.write(Ligne)
        Ligne='mean_hazard_curves = true\n'
        jobfile.write(Ligne)
        Ligne='quantile_hazard_curves =\n'
        jobfile.write(Ligne)
        Ligne='uniform_hazard_spectra = false\n'
        jobfile.write(Ligne)
        
        Ligne='poes = '+str(self.poes)+'\n'
        jobfile.write(Ligne)
        
        
#    def write_info(self):
#        # write the input file containing all the information about the run
#        # this file will be re-used as an input if the program is runned again?
#        info_run_file = str(self.Run_Name)+'/run.info'
#        write_info_file = open(info_run_file,'w')
#        write_info_file.write('Information on run : '+ self.Run_Name +'\n')
#        write_info_file.write('Option map: \t'+ str(self.option_map) +'\n')
#        if self.option_map == 'yes':
#            write_info_file.write('region = \t' + str(self.ll_lon) + '\t' + str(self.ll_lat)
#            + '\t' + str(self.lr_lon) + '\t'+ str(self.lr_lat)
#            + '\t' + str(self.ur_lon) + '\t' + str(self.ur_lat)
#            + '\t' + str(self.ul_lon) + '\t' + str(self.ul_lat) + '\n')
#            write_info_file.write('region_grid_spacing = \t' + str(self.region_grid_spacing) + '\n')
#        else:
#            write_info_file.write('Site Longitude : \t'+ str(self.Longitude) +'\n')
#            write_info_file.write('Site Latitude : \t'+ str(self.Latitude) +'\n')
#        write_info_file.write('Vs30 : \t'+ str(self.Vs30) +'\n')
#        write_info_file.write('Site Z1000 : \t'+ str(self.Z1000) +'\n')
#        write_info_file.write('Site Z2500 : \t'+ str(self.Z2500) +'\n')
#        write_info_file.write('nb_LT_samp : \t'+ str(self.nb_LT_samp) +'\n')
#        write_info_file.write('rup_mesh : \t'+ str(self.rup_mesh) +'\n')
#        write_info_file.write('source_discr : \t'+ str(self.source_discr) +'\n')
#        write_info_file.write('investigation_time : \t'+ str(self.investigation_time) +'\n')
#        write_info_file.write('Probability of exceedance : \t'+ str(self.poes) +'\n')
#        write_info_file.write('trunc_lvl : \t'+ str(self.trunc_lvl) +'\n')
#        write_info_file.write('max_dist : \t'+ str(self.max_dist) +'\n')
#        write_info_file.write('nb sample (sr,b_value, Mmax) : \t'+ str(self.nb_sample) +'\n')
#        write_info_file.write('Mmin : \t'+ str(self.Mmin) +'\n')
#        write_info_file.write('Random seed : \t'+ str(self.seed) +'\n')
#        write_info_file.write('SR correl : \t'+ str(self.sr_correl) +'\n')
#        write_info_file.write('SR increment size (mm/yr) : \t'+ str(self.size_of_increment) +'\n')
#        write_info_file.write('Mmax range : \t'+ str(self.Mmaxmin) +'\t'+ str(self.Mmaxmax) +'\n')
#
#        for intensity_i,min_int_i,max_int_out,nb_int_i in zip(self.intensity_out,self.min_int_out,self.max_int_out,self.nb_int_out) :
#
#            write_info_file.write('intensity_i : \t'+ str(intensity_i) +'\t'+ str(min_int_i) +'\t'+ str(max_int_out) +'\t'+ str(nb_int_i) +'\n')
#
#
#        write_info_file.close()

#    def action(self): #saves the data put in  the window
#        self.Latitude = self.Latitude.get()
#        self.Longitude = self.Longitude.get()
#
#        self.Mmin = self.Mmin.get()
#        self.seed = int(self.seed.get())
#        self.sr_correl = self.sr_correl.get()
#        if self.sr_correl == 'yes' :
#            self.sr_correl = True
#        else :
#            self.sr_correl = False
#        self.size_of_increment = self.size_of_increment.get()
#
#        self.Mmaxmin = self.Mmaxmin.get()
#        self.Mmaxmax = self.Mmaxmax.get()
#        self.Mmax_range = [float(self.Mmaxmin),float(self.Mmaxmax)]
#
#        self.Vs30 = self.Vs30.get()
#        self.Z1000 = self.Z1000.get()
#        self.Z2500 = self.Z2500.get()
#
#        self.nb_LT_samp = self.nb_LT_samp.get()
#        self.rup_mesh = self.rup_mesh.get()
#        self.source_discr = self.source_discr.get()
#        self.investigation_time = self.investigation_time.get()
#        self.trunc_lvl = self.trunc_lvl.get()
#        self.max_dist = self.max_dist.get()
#        self.nb_sample = int(self.nb_sample.get())
#        self.option_map = self.option_map.get()
#        self.poes = self.poes.get()
#
#        self.intensity_out = []
#        for self.intensity_i in self.intensity :
#            self.intensity_out.append(self.intensity_i.get())
#        self.min_int_out = []
#        for self.min_int_i in self.min_int :
#            self.min_int_out.append(float(self.min_int_i.get()))
#        self.max_int_out = []
#        for self.max_int_i in self.max_int :
#            self.max_int_out.append(float(self.max_int_i.get()))
#        self.nb_int_out = []
#        for self.nb_int_i in self.nb_int :
#            self.nb_int_out.append(int(self.nb_int_i.get()))
#
#        self.w_ini.destroy()
#
#    def action_map(self): #saves the data put in  the window
#        self.ll_lon = self.ll_lon.get()
#        self.ll_lat = self.ll_lat.get()
#
#        self.lr_lon = self.lr_lon.get()
#        self.lr_lat = self.lr_lat.get()
#
#        self.ur_lon = self.ur_lon.get()
#        self.ur_lat = self.ur_lat.get()
#
#        self.ul_lon = self.ul_lon.get()
#        self.ul_lat = self.ul_lat.get()
#
#        self.region_grid_spacing = self.region_grid_spacing.get()
#
#        self.w_o_map.destroy()
#
#
#if __name__=="__main__":
#    app = S_LT()
