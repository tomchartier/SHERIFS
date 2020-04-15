# -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

This code is pretty much the core of SHERIFS. It converts the slip-rate into earthquake rates.

@author: Thomas Chartier
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy
from scipy.interpolate import interp1d
import scalling_laws
import populate_bins
import mfd_shape
import time
import warnings
import core_utils, target, rates
from seismic_moment import mag_to_M0

warnings.simplefilter("ignore",RuntimeWarning)


class EQ_on_faults_from_sr():
    def __init__(self,M_min,mfd_param,faults_names,faults_area,faults_length,faults_width,faults_slip_rates,
                 scenarios,faults_shear_mod,path,sample,selected_ScL,dimention_used,
                 use_all_ScL_data,faults_mecanism,bg_ratio,size_of_increment,mfd_hyp,count_reruns,
                 faults_lon,faults_lat,Mmax_range,calculation_log_file):
        self.M_min = M_min
        self.mfd_param = mfd_param
        self.faults_names = faults_names
        self.faults_area = faults_area
        self.faults_length = faults_length
        self.faults_width = faults_width
        self.faults_slip_rates = faults_slip_rates
        self.scenarios = scenarios
        self.faults_shear_mod = faults_shear_mod
        self.path = path
        self.sample = sample
        self.selected_ScL = selected_ScL
        self.dimention_used = dimention_used
        self.use_all_ScL_data = use_all_ScL_data
        self.faults_mecanism = faults_mecanism
        self.bg_ratio = bg_ratio
        self.size_of_increment = size_of_increment
        self.mfd_hyp =mfd_hyp
        self.count_reruns=count_reruns
        self.faults_lon=faults_lon
        self.faults_lat=faults_lat
        self.Mmax_range = Mmax_range
        self.calculation_log_file = calculation_log_file
       
        self.initialize()
    def initialize(self):
        #####################################################################
        # input fixed values
        #####################################################################        
        faults_shear_mod = self.faults_shear_mod
        mfd_param = self.mfd_param
        M_min = self.M_min
        Mmax_range = self.Mmax_range
        
        
        #####################################################################
        #faults data
        #####################################################################               
        faults_names = self.faults_names
        faults_areas = self.faults_area
        faults_length = self.faults_length
        faults_width = self.faults_width
        faults_slip_rates = self.faults_slip_rates
        
        # file containing a log of what happened during the calculation
        log_calculation_file = open(self.path +'/Log/calculation_sample_' + str(self.sample) + '.txt','w')
        
        log_sliprep_file = open(self.path +'/Log/sliprep_sample_' + str(self.sample) + '.txt','w')
        
        
        
        #####################################################################
        #possible scenarii in this model
        #####################################################################        
        scenarios_names = self.scenarios
              
        '''##################################################################
        # finds the faults belonging to each scenario
        ##################################################################'''
        index_faults_in_scenario = []  #contains the indexs of the faults in each scenario      
        for scenario in scenarios_names: 
            i_scenario = []
            for i in range(len(scenario)) :
                i_scenario = i_scenario + list(np.where(np.array(faults_names) == scenario["f_%s" % str(i+1)]))
            i_scenario = np.reshape(i_scenario,(1,len(scenario)))
            index_faults_in_scenario.append(i_scenario)
            
        self.index_faults_in_scenario = index_faults_in_scenario
        
        #check if some fault are isolated (not included in the large scenario)
        faults_alone = []
        faults_isolated = []
        len_max_section_scenario = 0
        for indexes in index_faults_in_scenario:
            if len(indexes[0])> len_max_section_scenario:
                len_max_section_scenario = len(indexes[0])
                
        index_fault = 0
        for fault_name in faults_names:
            fault_alone_bool = True
            fault_isolated_bool = True
            for indexes in index_faults_in_scenario:
                if index_fault in indexes[0] :
                    fault_alone_bool = False
                    if len(indexes[0]) > len_max_section_scenario/2. :
                        fault_isolated_bool = False
            if fault_isolated_bool == True and fault_alone_bool == False:
                faults_isolated.append(index_fault)
            if fault_alone_bool == True:
                faults_alone.append(index_fault)
            
            index_fault+=1
            
        
        '''##################################################################
        # area and length of each scenario
        ##################################################################'''        
        scenario_area = [] #area of each scenario
        
        #if two faults are overlaping (for example to describe a normal and strikeslip component), only take one area.
        index_faults_in_scenario_for_scl=[]
        for indexes_i in index_faults_in_scenario:
            list_of_points = []
            index_for_scenario = []
            for index_fault in indexes_i[0]:
                list_i = str(self.faults_lon[index_fault][0])+str(self.faults_lon[index_fault][-1])+str(self.faults_lat[index_fault][0])+str(self.faults_lat[index_fault][-1])
                if not list_i in list_of_points:
                    list_of_points.append(list_i)
                    index_for_scenario.append(index_fault)
            index_faults_in_scenario_for_scl.append(index_for_scenario)
            #print str(len(indexes_i)-len(index_for_scenario))
                
        
        for i in index_faults_in_scenario_for_scl:
            scenario_i_area = np.sum(np.take(faults_areas,i))
            scenario_area.append(scenario_i_area)
            
        scenario_length = [] #length of each scenario
        
        for i in index_faults_in_scenario_for_scl:
            scenario_i_length = np.sum(np.take(faults_length,i))
            scenario_length.append(scenario_i_length)
                
        scenario_width = [] #width of each scenario
        
        for i in index_faults_in_scenario_for_scl:
            scenario_i_width = np.mean(np.take(faults_width,i))
            scenario_width.append(scenario_i_width)
                
        
        
        '''##################################################################
        #####################################################################
        #
        # setting up the coefficients for the scalling law
        # and calculate the Mmax
        #
        #####################################################################
        ##################################################################'''
        # Mmax can be forced within a range by the user  Mmaxmin < Mmax < Mmaxmax
        
        Mmaxmin = Mmax_range[0]
        Mmaxmax = Mmax_range[1]
        Mmax = -1.
        loop_Mmax = 1
        while Mmax< Mmaxmin or Mmax > Mmaxmax :
            # file containing the log of the maximal magnitude of each fault and each scenario
            log_Mmax_file=open(self.path +'/Log/Mmax_sample_' + str(self.sample) + '.txt','w')
            if loop_Mmax == 1 :
                Mmaxs = scalling_laws.Calc_Mmax(self.faults_area,scenario_area,self.faults_length,scenario_length,self.faults_width,scenario_width,self.selected_ScL,
                                                self.dimention_used,self.use_all_ScL_data,self.faults_mecanism,index_faults_in_scenario, self.sample)
            else :
                Mmaxs = scalling_laws.Calc_Mmax(self.faults_area,scenario_area,self.faults_length,scenario_length,self.faults_width,scenario_width,self.selected_ScL,
                                                self.dimention_used,self.use_all_ScL_data,self.faults_mecanism,index_faults_in_scenario, 10)
                
            Mmax_faults = Mmaxs.Mmax_faults #Mmax of each fault
            Mmax_scenario = Mmaxs.Mmax_scenario #Mmax of each scenario    
            final_fault_length = Mmaxs.final_fault_length
            final_scenario_length = Mmaxs.final_scenario_length
            
            #writting in the log
            index_fault = 0
            for Mmax_i in Mmax_faults :
                line = (faults_names[index_fault] + '\t' + str(round(faults_length[index_fault]/1000.,1))   + '\t' + final_fault_length[index_fault] 
                + '\t' + str(round(self.faults_width[index_fault],1))  + '\t' + str(round(self.faults_area[index_fault]/1000000.,1)) 
                + '\t' + str(Mmax_i) + '\n')
                log_Mmax_file.write(line)
                index_fault += 1
            index_scenario = 0
            for Mmax_i in Mmax_scenario : 
                line = (str(scenarios_names[index_scenario]) + '\t' + str(round(scenario_length[index_scenario]/1000.,1))    + '\t' + final_scenario_length[index_scenario] 
                + '\t' + str(round(scenario_width[index_scenario],1)) + '\t' + str(round(scenario_area[index_scenario]/1000000.,1))  
                + '\t' + str(Mmax_i) + '\n')
                log_Mmax_file.write(line)
                index_scenario += 1
            log_Mmax_file.close()
                
            
            if np.size(scenarios_names) == 0 :
                Mmax = max(Mmax_faults)
            else :
                Mmax = max((max(Mmax_faults),max(Mmax_scenario)))
                
            loop_Mmax+=1
            
            if loop_Mmax == 30 :
                print('An Mmax incompatible with the ruptures is imposed!! Change it in run.info or change the rupture.txt file')
                
        if loop_Mmax >= 3 :
            print('Mmax imposed: '+str(Mmax))#+'   see EQ_on_faults.py for details.')
            self.calculation_log_file.write('\nMmax imposed: '+str(Mmax))#+'   see EQ_on_faults.py for details.')
        '''##################################################################
        #etablish the magnitude bin (0.1)
        ##################################################################'''
        bin_mag = np.linspace(M_min,Mmax,int(round((Mmax-M_min),1)*10. +1))
                
        '''#####################################################################
        # initializing incremental reccurence of the earthquakes for each rupture and each bin)
        #####################################################################'''
        rup_rates = {}
        index_rup = []
        for i in range(len(faults_names)):
            rup_rates.update({str(i):{'rup_id':i,'rup_name':faults_names[i],
            'involved_faults':[i],'Mmax':Mmax_faults[i] ,
            'area' : faults_areas[i],'rates':np.zeros(len(bin_mag))}})
            index_rup.append(i)
        i_end = i+1
        if np.size(scenarios_names) != 0 :
            for j in range(len(scenarios_names)):
                rup_rates.update({str(i_end+j):{'rup_id':i_end+j,'rup_name':scenarios_names[j],
                'involved_faults':index_faults_in_scenario[j][0],'Mmax':Mmax_scenario[j],
                'area' : scenario_area[j],'rates':np.zeros(len(bin_mag))}})
                index_rup.append(i_end+j)
                
        if str(self.sample) == '1':
            log_rup_file = open(self.path +'/Log/ruptures.txt','w')
            log_rup_file.write('rup_id\tinvolved_faults\n')
            for i in range(len(rup_rates)):
                log_rup_file.write(str(rup_rates.get(str(i)).get('rup_id'))+'\t')
                for j in rup_rates.get(str(i)).get('involved_faults'):
                    log_rup_file.write(str(j)+' ')
                log_rup_file.write('\n')
            log_rup_file.close()
        
        
#        OQ_entry_faults = np.zeros((len(faults_names),len(bin_mag)))
#
#        if np.size(scenarios_names) != 0 :
#            OQ_entry_scenarios = np.zeros((len(scenarios_names),len(bin_mag)))
        
        #####################################################################
        #delete thes scenarios that are too small and the flauts that are too small and are not used in any scenario
        #####################################################################
        
#        fault_or_scenario_deleted = False
#        #for the scenarios
#        index_scenario = 0
#        index_to_be_deleted = []
#        for scenario in scenarios_names :
#            if Mmax_scenario[index_scenario] < M_min:
#                print('scenario deleted : '+scenario+' Mmax : '+Mmax_scenario[index_scenario])
#                self.calculation_log_file.write('\nscenario deleted : '+scenario+' Mmax : '+Mmax_scenario[index_scenario])
#                index_to_be_deleted.append(index_scenario)
#                fault_or_scenario_deleted = True
#            index_scenario += 1
#
#        scenarios_names = [i for j, i in enumerate(scenarios_names) if j not in index_to_be_deleted]
#        Mmax_scenario = [i for j, i in enumerate(Mmax_scenario) if j not in index_to_be_deleted]
#        scenario_length = [i for j, i in enumerate(scenario_length) if j not in index_to_be_deleted]
#        scenario_area = [i for j, i in enumerate(scenario_area) if j not in index_to_be_deleted]
#        index_faults_in_scenario = [i for j, i in enumerate(index_faults_in_scenario) if j not in index_to_be_deleted]
#
#        index_fault = 0
#        index_to_be_deleted = []
#        for fault in faults_names :
#            is_fault_in_scenario = False
#            for index_faults_in_scenario_i in index_faults_in_scenario :
#                if index_fault in index_faults_in_scenario_i:
#                    is_fault_in_scenario = True
#            #print fault, is_fault_in_scenario
#            if Mmax_faults[index_fault] < M_min and is_fault_in_scenario == False :
#                log_calculation_file.write('fault deleted : ' + str(fault) + ' Mmax : '+str(Mmax_faults[index_fault]) + '\n')
#                print('fault too small to break on its own : '+fault+' Mmax : '+Mmax_faults[index_fault] )
#                self.calculation_log_file.write('\nfault too small to break on its own : '+fault+' Mmax : '+Mmax_faults[index_fault] )
#                index_to_be_deleted.append(index_fault)
#                fault_or_scenario_deleted = True
#            index_fault += 1
            
        self.index_faults_in_scenario = index_faults_in_scenario
        
        
        '''#####################################################################
        #for each bin, find which fault and which scenario populates it.
        #####################################################################'''
        rup_in_bin = populate_bins.pop(bin_mag,index_rup,rup_rates,M_min)
            
            
        '''##################################################################
        # Definition of the background - how much is not on the faults
        ##################################################################'''
            
        # proportion of the BG
        bin_mag_fault_prop = [ 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.]
        fault_prop_inc = self.bg_ratio
        
        bin_mag_fault_prop.append(10.)
        #print type(fault_prop_inc)
        fault_prop_inc = np.append(np.array(fault_prop_inc),1.)
        
        fault_prop = interp1d(bin_mag_fault_prop,fault_prop_inc) 
        
        EQ_rate_BG = np.zeros(len(bin_mag)) #rate of EQ in the BG for each magnitude
            
        
        '''##################################################################
        #Calculation of the total moment rate sum(slip-rate)*sum(area) + moment in the background
        ##################################################################'''
        
        Total_moment_faults_rate_init = 0. 
        
        index_fault = 0        
        for fault_area in faults_areas :   #adding the moment rate of each fault         
            Total_moment_faults_rate_init += faults_shear_mod[index_fault] * fault_area * faults_slip_rates[index_fault]
            index_fault += 1
        
        log_calculation_file.write('Moment_rate_fault_initial calculated using muAs' + '\t' + str(Total_moment_faults_rate_init) + '\t' + 'N.m' + '\n')
        
         
        Total_moment_rate_fault = 0. #will be use to compare with the initial moment rate
        
        '''##################################################################
        #etablish the target regional distribution
        ##################################################################'''
        
        #if GR ... create "p" which is taking into account for the random sample in the distribution later
        if self.mfd_hyp == 'GR' :
            p_MFD = mfd_shape.GR(mfd_param,bin_mag)
                
        if self.mfd_hyp == 'double_GR' :
            p_MFD = mfd_shape.double_GR(mfd_param,bin_mag)
                
        elif self.mfd_hyp == 'YC':
            p_MFD = mfd_shape.YC(mfd_param,bin_mag)
                
        elif self.mfd_hyp == 'YC_marmara':
            p_MFD = mfd_shape.YC_marmara(mfd_param,bin_mag)
                
        elif self.mfd_hyp == 'YC_modified':
            p_MFD = mfd_shape.YC_modified(mfd_param,bin_mag)
                
        elif self.mfd_hyp == 'UCERF_DV':
            # MFD calculated from the sum of the mfds of the faults of the UCERF
            # model for death valley
            p_MFD = mfd_shape.UCERF_DV(bin_mag)
            #can be used as a template for any shape numerically defined

        p_MFD = (p_MFD) / sum(p_MFD) #normalize the probability distribution

        
        '''##################################################################
        # Converting in term of moment
        ##################################################################'''
        p_MFD_MO = []  #distribution in terms of moment
        index_mag = 0
        for mag in bin_mag :
            M0 = mag_to_M0(mag)
            p_MFD_MO.append(fault_prop(mag) * p_MFD[index_mag] * M0) #considering the proportion that needs to be in the BG
            index_mag += 1
        p_MFD_MO = (p_MFD_MO) / sum(p_MFD_MO) #normalize the distribution
        
        '''##################################################################
        # Adding the moment rate of the background to the Total_moment_rate_init
        ##################################################################'''
        Total_moment_rate_init = Total_moment_faults_rate_init
        index_mag = 0
        for mag in bin_mag :
            moment_rate_fault_bin = p_MFD_MO[index_mag] * Total_moment_faults_rate_init # M0 of the faults for that bin
            moment_rate_BG_bin = ((1-fault_prop(mag)) * moment_rate_fault_bin) / fault_prop(mag) # M0 of the background for that bin
            Total_moment_rate_init += moment_rate_BG_bin
            index_mag += 1
        
        
        '''##################################################################
        # Etablish the target of moment rate per bin
        ##################################################################'''
        # distribution of probability multiply by the initial moment rate
        target_moment_per_bin = p_MFD_MO * Total_moment_rate_init
        # this is our target. This way we insure a GR distribution in the end
        
        
        '''##################################################################
        # Create the sliprate budget dictionary
        # This matrix contain a list of the fault name. each name is repeted 
        # by a number of time depending of the slip-rate of the fault.
        ##################################################################'''
        size_of_increment = float(self.size_of_increment) * 0.001  #in meter/year
        if self.count_reruns != 1:
            size_of_increment = size_of_increment/(float(self.count_reruns)*1.5-1.) #divide the size of the increment if there are re-runs
        
        faults_budget = {}
        for index_fault in range(len(faults_names)) :
            nb_f_name = int(round(faults_slip_rates[index_fault] / size_of_increment,0))  #nb of dsr to spend
            faults_budget.update({index_fault:nb_f_name})
            
        '''##################################################################
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        #   Populate the bin of magnitude of each fault and each scenario
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        ##################################################################'''
        
        M_slip_repartition = [] #how the slip rate is used
        for fault in faults_names :
            M_slip_repartition_i = []
            M_slip_repartition_i.append(fault)
            M_slip_repartition.append(M_slip_repartition_i)
            
        moment_rate_in_bin = np.zeros(len(bin_mag))   # moment aleardy present in the bin of magnitude
        model_MFD = []
        
        number_of_loops = 0   # variable to follow how the calculation is doing
        number_of_loops_before = 0 
        number_of_loops_for_nothing = 0
        number_of_loops_for_nothing_before = 0
        number_of_loops_last_checked = -1
        empty_bins = []
        bin_target_reached = []
        len_faults_budget = [] # will be use to check if the calculation is stuck and too slow (see under)
        aseismic_count = 0
        
        time_weighting_faults = 0
        time_checking_the_fit_2 = 0
        color_mag= []
        
        
        '''######################
        #####   MAIN LOOP   #####
        ######################'''
        nb_ss_to_spend = float(sum(faults_budget.values()))
        print('number of dsr to spend : '+ str(nb_ss_to_spend))
        self.calculation_log_file.write('\nnumber of dsr to spend : '+ str(nb_ss_to_spend)+'\n')
        print_percent = True
        print_target_set = True
        bool_target_set = False
        uniform_spending = True
        
        slip_rate_use_per_fault = np.zeros(len(faults_names))
        moment_rate_required = 0.
        moment_rate_left= Total_moment_faults_rate_init
        
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        
        test_mean_picked = []
        most_likely_pick = []
        
        while sum(faults_budget.values()) != 0 : # as long as there is some slip-rate to spend we keep going
            ratio_done = 1. - float(sum(faults_budget.values()))/nb_ss_to_spend
            if ratio_done > 0.01 :
                model_MFD, self.calculation_log_file,print_percent = core_utils.progress(rate_in_model,model_MFD,self.calculation_log_file,ratio_done,print_percent)
            
            number_of_loops += 1
            
            
            ''' Calculate the new target shape in each bin in terms of moment rate '''
            target_i = target.get_new_target(number_of_loops,moment_rate_in_bin,p_MFD_MO,target_moment_per_bin,bin_mag,empty_bins,bin_target_reached,rup_in_bin)
            #print(target_i)
            most_likely_pick.append(bin_mag[list(target_i).index(max(list(target_i)))])
            
            if len(empty_bins) != len(bin_mag):
                '''Pick the bin of magnitude to fill according to the current distribution '''
                # normalize the target mfd
                target_i = (target_i) / sum(target_i)
                picked_bin = np.random.choice(len(bin_mag), 1, p = target_i)[0]
                mag = bin_mag[picked_bin] #magnitude of that bin
                
#                test_mean_picked.append(picked_bin)
#                if str(round(ratio_done,2))[-1]=='9' and len(test_mean_picked) > 15:
#                    print(round(ratio_done,2),
#                    round(max(set(most_likely_pick[-10:]),key = most_likely_pick[-10:].count),1),
#                     round(np.mean(test_mean_picked[-10:]),1),
#                     len(bin_target_reached),len(bin_mag))
                
                if len(rup_in_bin[picked_bin]) != 0:
                    time_i = time.time()
                    '''Calculate the weight for sampling of the fault or scenario'''
                    weight_rup = core_utils.weight_fault_sampling(picked_bin,rup_in_bin,faults_names,faults_slip_rates,slip_rate_use_per_fault,faults_alone,scenarios_names,faults_isolated,index_faults_in_scenario,rup_rates)
                    
                    time_weighting_faults += time.time()-time_i
                    
                    '''#picking of the source'''
                    try :
                        picked_rup = np.random.choice(rup_in_bin[picked_bin],1,p=weight_rup)[0] #picked source
                    except ValueError:
                        picked_rup = np.random.choice(rup_in_bin[picked_bin])
                            
                    #index_fault = np.where(np.array(faults_names) == picked_rup)[0]
                    index_fault = rup_rates.get(str(picked_rup)).get('involved_faults')
                    if bool_target_set == False :
                        if (len(rup_in_bin[-3]) +  len(rup_in_bin[-2]) + len(rup_in_bin[-1])) == 0 : #if all the slip_rate has been spared for the Mmax - 0.3
                            rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
                            bool_target_set = True
                            #all the slip_rate has NOT been spared for the Mmax
                            #does other checks
                            # do we have enough moment to fit the shape?
                        if moment_rate_left <= (1. - 0.00001)* moment_rate_required:
                            #if not enough momen rate left, stop spending on the larger EQs
                            self.calculation_log_file.write('\n Not enough moment left ')
                            rup_in_bin[-3] = []
                            rup_in_bin[-2] = []
                            rup_in_bin[-1] = []
                            bool_target_set = True
                        # check if the rates of the antepenultimate bin is not too different of the two other ones
                        if len(rup_in_bin[-1])+ len(rup_in_bin[-2]) == 0 and bool_target_set == False:
                            if moment_rate_in_bin[-3] >= 20. * (moment_rate_in_bin[-2]+moment_rate_in_bin[-1]):
                                self.calculation_log_file.write('\n antepenultimate bin getting too high')
                                rup_in_bin[-3] = []
                                bool_target_set = True
   
                        #test if the moment rate left  to spend is enough to fit the target MFD
                        if number_of_loops > number_of_loops_last_checked + 10.:
                            time_ii = time.time()
                            number_of_loops_last_checked = number_of_loops
                            
                            moment_rate_left = Total_moment_faults_rate_init - Total_moment_rate_fault
                            rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
                            rate_Mmax_check = np.mean(rate_in_model[-3:])
                            
                            moment_rate_required = 0.
                            for index_mag in range(len(bin_mag)-3):#loop on all the magnitude bins except the three last ones
                                rate_Mi_check = rate_in_model[index_mag]
                                target_GR_i_check = rate_Mmax_check  * p_MFD[index_mag] / p_MFD[-2]
                                #difference between the target at this point and what is already filled in terms of moment rate
                                moment_rate_required += (mag_to_M0(bin_mag[index_mag]) * target_GR_i_check) - (mag_to_M0(bin_mag[index_mag]) * rate_Mi_check)
                                
                            if self.mfd_hyp == 'UCERF_DV_':
                                rate_Mmax = rate_f_in_model[-1]
                                moment_rate_required = 0.
                                for index_mag in range(len(bin_mag)-3):#loop on all the magnitude bins except the three last ones
                                    rate_Mi_check = rate_in_model[index_mag]
                                    target_GR_i_check = rate_Mmax_check  * p_MFD[index_mag] / p_MFD[-3]
                                    #difference between the target at this point and what is already filled in terms of moment rate
                                    moment_rate_required += (mag_to_M0(bin_mag[index_mag]) * target_GR_i_check) - (mag_to_M0(bin_mag[index_mag]) * rate_Mi_check)
        
                            time_checking_the_fit_2+=time.time()-time_ii
                            
                            
                        if print_target_set == True and bool_target_set ==True:
                            ####
                            # Setting the target
                            ####
                            print_target_set = False
                            print('- target set - ')
                            self.calculation_log_file.write('\n- target set - ')
                            rate_at_target_setting = rate_in_model
                            #record the MFD at the moment the target is set
                            
                            rate_Mmax = np.mean(rate_in_model[-3:])
                                
                            rate_Mi = rate_in_model[picked_bin]
                            target_GR_i = rate_Mmax  * p_MFD[picked_bin] / (p_MFD[-2]) #target of the picked bin in order to fit the original distribution
                            if self.mfd_hyp == 'UCERF_DV_':
                                rate_Mmax = rate_f_in_model[-1]
                                    
                                rate_Mi = rate_in_model[picked_bin]
                                target_GR_i = rate_Mmax  * p_MFD[picked_bin] / p_MFD[-3] #target of the picked bin in order to fit the original distribution
                                
                        
                    
                            
                    ''' spending the slip_rate increment '''
                    
                    index_fault = rup_rates.get(str(picked_rup)).get('involved_faults')
                    sr_to_spend = True #all of the faults involved still have slip rate to spend
                    shear_mod = 0
                    for index in index_fault :
                        shear_mod += faults_shear_mod[index]
                        if faults_budget[index] <= 0:
                            sr_to_spend = False #one of the faults don't have anymore budget
                            for index_mag in range(len(bin_mag)):
                                if picked_rup in rup_in_bin[index_mag]:
                                    rup_in_bin[index_mag].remove(picked_rup)
                    if sr_to_spend == True :
                        shear_mod = shear_mod / float(len(index_fault))
                        area = rup_rates.get(str(picked_rup)).get('area')
                        displacement = mag_to_M0(mag)/(shear_mod*area)
                        rate_i = size_of_increment/displacement
                        if bool_target_set == True:

                            rate_Mi = rate_in_model[picked_bin]
                            target_GR_i = rate_Mmax  * p_MFD[picked_bin] / (p_MFD[-2]) #target of the picked bin in order to fit the original distribution
                            if rate_Mi < target_GR_i : #for this bin, the GR hasn't been reached yet
                                #OQ_entry_scenarios[index_scenario[0]][picked_bin] = OQ_entry_scenarios[index_scenario[0]][picked_bin] + rate_i
                                
                                if uniform_spending == True or len(index_fault)==1:
                                    for index in index_fault :
                                        M_slip_repartition[index].append(picked_rup)
                                        faults_budget[index]+=-1
                                        slip_rate_use_per_fault[index] += size_of_increment
                                        rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i
                                        rate_in_model[picked_bin] += rate_i
                                        moment_rate_i = mag_to_M0(mag) * rate_i
                                else :
                                    M_slip_repartition,faults_budget,slip_rate_use_per_fault,nb_sdr_used = core_utils.variable_spending(index_fault,M_slip_repartition,faults_budget,slip_rate_use_per_fault,size_of_increment,faults_slip_rates,picked_rup)
                                        #adding to the rate
                                    rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i*(nb_sdr_used)
                                    rate_in_model[picked_bin] +=  rate_i*(nb_sdr_used)
                                    moment_rate_i = mag_to_M0(mag) * rate_i*(nb_sdr_used)
        
                                #substracting the moment used from the target
                                moment_rate_in_bin[picked_bin] += moment_rate_i
                                Total_moment_rate_fault += moment_rate_i
                            else : # for this bin, the GR has been reached, this slip rate needs to be aseismic
                                rup_in_bin[picked_bin] = []
                                if not picked_bin in bin_target_reached:
                                    bin_target_reached.append(picked_bin)
                                   
                        else : # the absolute target as not been reached yet

                            if uniform_spending == True or len(index_fault)==1:
                                for index in index_fault :
                                    M_slip_repartition[index].append(picked_rup)
                                    faults_budget[index]+=-1
                                    slip_rate_use_per_fault[index] += size_of_increment
                                    rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i
                                    rate_in_model[picked_bin] +=  rate_i
                                    moment_rate_i = mag_to_M0(mag) * rate_i
                            else :
                                M_slip_repartition,faults_budget,slip_rate_use_per_fault,nb_sdr_used = core_utils.variable_spending(index_fault,M_slip_repartition,faults_budget,slip_rate_use_per_fault,size_of_increment,faults_slip_rates,picked_rup)
                                    #adding to the rate
                                rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i*(nb_sdr_used)
                                rate_in_model[picked_bin] +=  rate_i*(nb_sdr_used)
                                moment_rate_i = mag_to_M0(mag) * rate_i
    
                            #substracting the moment used from the target
                            moment_rate_in_bin[picked_bin] += moment_rate_i
                            Total_moment_rate_fault += moment_rate_i
                    
                        #rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
                else :
                    if not picked_bin in empty_bins:
                        empty_bins.append(picked_bin)
                    
                    
                if number_of_loops > number_of_loops_before+200:
                    number_of_loops_before = number_of_loops
                    #if the fault has no more slip to spare, the scenarios and the fault are remove for the picking in order to have a faster picking
                    for index_mag in range(len(bin_mag)):
                        for i_rup in rup_in_bin[index_mag] :
                            remove_rup = False
                            for index_fault in rup_rates.get(str(i_rup)).get('involved_faults'):
                                if faults_budget[index_fault] <= 0:
                                    remove_rup = True
                            if remove_rup == True :
                                rup_in_bin[index_mag].remove(i_rup)
                
                                        
                len_faults_budget.append(sum(faults_budget.values()))
                #does a check to see if the remainingfaults that still have some slip rate fit in one of the bin of are too small
                if len(len_faults_budget) > 3 :
                    if len_faults_budget[-2] == len_faults_budget[-1] :
                        number_of_loops_for_nothing+=1
                    if number_of_loops_for_nothing_before<number_of_loops_for_nothing-100:
                        number_of_loops_for_nothing_before = number_of_loops_for_nothing
                        if len_faults_budget[-1] == len_faults_budget[-10] :
                            rup_still_used = []
                            for rup_i in rup_in_bin:
                                rup_still_used += rup_i
                            rup_still_used = set(rup_still_used)
                            fault_still_used = []
                            for rup_i in rup_still_used:
                                fault_still_used += list(rup_rates.get(str(rup_i)).get('involved_faults'))
                            fault_still_used = set(fault_still_used)
                            for fault,index_fault in zip(faults_names,range(len(faults_names))) :
                                if not (index_fault in fault_still_used) and (faults_budget[index_fault] > 0):
                                    while faults_budget[index_fault] > 0 :
                                        ratio_done = 1. - float(sum(faults_budget.values()))/nb_ss_to_spend
                                        if ratio_done > 0.01 :
                                            model_MFD, self.calculation_log_file,print_percent = core_utils.progress(rate_in_model,model_MFD,self.calculation_log_file,ratio_done,print_percent)
                                        faults_budget[index_fault]+=-1
                                        M_slip_repartition[index_fault].append('NMS')
                                        aseismic_count += 1
            else:
                print('-target filled-')       
                self.calculation_log_file.write('\n-target filled-')
                while sum(faults_budget.values())!=0: # as long as there is some slip-rate to spend we keep going
                
                    ratio_done = 1. - float(sum(faults_budget.values()))/nb_ss_to_spend
                    if ratio_done > 0.01 :
                        model_MFD, self.calculation_log_file,print_percent = core_utils.progress(rate_in_model,model_MFD,self.calculation_log_file,ratio_done,print_percent)
                    for index_fault in range(len(faults_names)):
                        if faults_budget[index_fault] > 0 :
                            faults_budget[index_fault]+=-1
                            M_slip_repartition[index_fault].append('NMS')
                            aseismic_count += 1
                    
        
        '''##################################################################  
        # 
        # Definition of the background seismicity
        # 
        ##################################################################'''
        rate_f_in_model, rate_bg_in_model = rates.get_rate_faults_n_bg(rup_rates,fault_prop,bin_mag)
        EQ_rate_BG = rate_bg_in_model
        '''##################################################################
        # printing
        ##################################################################''' 
        #printing the evolution of the model durring the incrementations
        rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        model_MFD.append(rate_in_model) 
        colors = ['gainsboro','darkgray','dimgray','dimgray','black','red']
        index_color = 0
        for MFD_i in model_MFD:
            plt.plot(bin_mag,MFD_i,color=colors[index_color])
            index_color+=1
        if print_target_set == True:
            rate_at_target_setting = model_MFD[-1]
        plt.plot(bin_mag,rate_at_target_setting,':g')
        target_i = []
        for index_mag in range(len(bin_mag)):
            rate_Mmax = (model_MFD[-1][-3]+model_MFD[-1][-2]+model_MFD[-1][-1])/3.
            target_GR_i = rate_Mmax  * p_MFD[index_mag] / p_MFD[-2]

            if self.mfd_hyp == 'UCERF_DV':
                rate_Mmax = model_MFD[-1][-3] #rate of the Mmax due to faults
                target_GR_i = rate_Mmax  * p_MFD[index_mag] / p_MFD[-3]
            target_i.append(target_GR_i)
        plt.plot(bin_mag,target_i,':b')
        plt.scatter(bin_mag[-3:],model_MFD[-1][-3:])
        plt.yscale('log')
        plt.savefig(self.path +'/Log/target_fit_' + str(self.sample) + '.png' ,dpi = 180, transparent=True)
        plt.close()
        index_5 = 0
        while bin_mag[index_mag] < 5.:
            index_5+=1
        self.ratio_test = np.mean([abs(target_i[index_5+0]/model_MFD[-1][index_5+0]),
                                   abs(target_i[index_5+1]/model_MFD[-1][index_5+1]),
                                    abs(target_i[index_5+2]/model_MFD[-1][index_5+2]),
                                    abs(target_i[index_5+3]/model_MFD[-1][index_5+3]),
                                    abs(target_i[index_5+4]/model_MFD[-1][index_5+4])])
        MFD_to_test = (model_MFD[-1])/sum(model_MFD[-1])
        target_to_test = (target_i)/sum(target_i)
        array_absolute_ratio = []
        for index_mag in range(len(bin_mag)):
            array_absolute_ratio.append(abs(MFD_to_test[index_mag]/target_to_test[index_mag]))
        self.ratio_test = np.mean(array_absolute_ratio)
        print('ratio between the target and the shape of the model : '+str(round(self.ratio_test,2)))
        self.calculation_log_file.write('\nratio between the target and the shape of the model : '+str(round(self.ratio_test,2)))
        '''##################################################################  
        # 
        # OUTPUT for Source_model_creator
        # 
        ##################################################################'''
        
        self.OQ_entry_faults = []
        self.OQ_entry_scenarios = []
        for i in range(len(faults_names)):
            self.OQ_entry_faults.append(rup_rates.get(str(i)).get('rates'))
        i_end = i+1
        if np.size(scenarios_names) != 0 :
            for j in range(len(scenarios_names)):
                self.OQ_entry_scenarios.append(rup_rates.get(str(i_end+j)).get('rates'))
        self.faults_names = faults_names
        self.scenarios_names = scenarios_names
        
        self.EQ_rate_BG = EQ_rate_BG
        
        self.bin_mag = bin_mag
        
        log_calculation_file.write('Moment_rate_faults_final calculated using 10^(1.5M+9.1)' + '\t' + str(Total_moment_rate_fault) + '\t' + 'N.m' + '\n')
        log_calculation_file.write('Ratio of NMS slip rate counted in the iterative process: ' + '\t' + str(round((100.) * (aseismic_count/nb_ss_to_spend))) + '\t %\n') 
        log_calculation_file.write('Moment_rate_model calculated using 10^(1.5M+9.1)' + '\t' + str(Total_moment_rate_fault) + '\t' + 'N.m' + '\n')
        log_calculation_file.close()
        self.ratio_NMS = round((100.) * (1. - Total_moment_rate_fault / Total_moment_faults_rate_init))
        print('ratio of NMS : '+str(round((100.) * (aseismic_count/nb_ss_to_spend))))
        self.calculation_log_file.write('\nratio of NMS : '+str(round((100.) * (aseismic_count/nb_ss_to_spend))))
        
        
        for rup in M_slip_repartition:
            for i in rup :
                log_sliprep_file.write(str(i)+' ')
            log_sliprep_file.write('\n')
        log_sliprep_file.close()
        
