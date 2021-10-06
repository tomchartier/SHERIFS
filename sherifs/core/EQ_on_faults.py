# -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

This code is pretty much the core of SHERIFS. It converts the slip-rate into
earthquake rates.

@author: Thomas Chartier
"""
import os
import time
import math
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from sherifs.core import (populate_bins, core_utils, mfd_shape, target, rates,
                          scaling_laws)
from sherifs.core.seismic_moment import mag_to_M0

warnings.simplefilter("ignore", RuntimeWarning)


class EQ_on_faults_from_sr():
    def __init__(self,Run_Name,M_min,mfd_param,faults_names,faults_area,faults_length,faults_width,faults_slip_rates,
                 scenarios,faults_shear_mod,path,pathlog,sample,selected_ScL,dimention_used,
                 use_all_ScL_data,faults_mecanism,bg_ratio,size_of_increment,mfd_hyp,count_reruns,
                 faults_lon,faults_lat,Mmax_range,calculation_log_file,branch,param):
        self.Run_Name = Run_Name
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
        self.pathlog = pathlog
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
        self.branch = branch
        self.param = param

        self.initialize()
    def initialize(self):
        #####################################################################
        # input fixed values
        #####################################################################
        faults_shear_mod = self.faults_shear_mod
        mfd_param = self.mfd_param
        M_min = self.M_min
        Mmax_range = self.Mmax_range


        inti_core_time = time.time()

        #####################################################################
        #faults data
        #####################################################################
        faults_names = self.faults_names
        faults_areas = self.faults_area
        faults_length = self.faults_length
        faults_width = self.faults_width
        faults_slip_rates = self.faults_slip_rates

        # file containing a log of what happened during the calculation
        log_calculation_file = open(self.pathlog +'/calculation_sample_' + str(self.sample) + '.txt','w')

        log_sliprep_file = self.pathlog +'/sliprep_sample_' + str(self.sample) + '.pkl'


        re_use = True
        #####################################################################
        #possible scenarii in this model
        #####################################################################
        scenarios_names = self.scenarios

        run_name = self.Run_Name
        model_name = self.branch["model"]

        # if not os.path.exists(run_name+'/'+model_name+'_Log'):
        #     os.makedirs(run_name+'/'+model_name +'_Log')
        scl_name = self.branch["scl"][0] + "_" + \
         self.branch["scl"][1] + "_" + self.branch["scl"][2]
        set_name = self.branch["set"]
        #f_bin_pop = self.path +'/Log/bin_pop_'+str(self.sample)+'.pkl'
        f_mmax = run_name+'/LOG/'+model_name +'_mmax_'+scl_name+'_'+set_name+'_'+str(self.sample)+'.pkl'
        if not os.path.isfile(f_mmax):
            re_use = False
        if re_use == False:
            print("Building scenarios and Mmax")
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

            '''
            check the max dimention for a single fault according to the
            aspect ratio
            '''
            f_len_ar, f_width_ar, f_area_ar = [], [], []
            ar = self.param["main"]["parameters"]["aspect_ratio"]
            for l, w in zip(self.faults_length,self.faults_width):
                if l < (w * ar) :
                    f_len_ar.append(l)
                    f_width_ar.append(l* ar)
                    f_area_ar.append(l*l* ar)
                else :
                    f_len_ar.append(l)
                    f_width_ar.append(w)
                    f_area_ar.append(l*w)


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
            while Mmax < Mmaxmin or Mmax > Mmaxmax:
                # file containing the log of the maximal magnitude of each fault and each scenario
#                log_Mmax_file=open(self.path +'/Log/Mmax_sample_' + str(self.sample) + '.txt','w')

                # TODO fix this using os.path.join
                log_Mmax_file=open(run_name+'/LOG/'+model_name +'_Log_Mmax_sample_'+scl_name+'_'+set_name+'_'+str(self.sample)+'.txt','w')

                if loop_Mmax == 1:
                    Mmaxs = scaling_laws.Calc_Mmax(f_area_ar,scenario_area,f_len_ar,scenario_length,f_width_ar,scenario_width,self.selected_ScL,
                                                    self.dimention_used,self.use_all_ScL_data,self.faults_mecanism,index_faults_in_scenario, self.sample)
                else :
                    Mmaxs = scaling_laws.Calc_Mmax(f_area_ar,scenario_area,f_len_ar,scenario_length,f_width_ar,scenario_width,self.selected_ScL,
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

            with open(f_mmax, 'wb') as f:

                dump_to_file = [Mmax_faults ,
                Mmax_scenario,
                final_fault_length,
                final_scenario_length,
                faults_alone,
                faults_isolated,
                index_faults_in_scenario,
                scenario_area]

                pickle.dump(dump_to_file, f)

            print("\t - scenario and max built")

        else :
            print('Reloading Mmax from data file')
            with open(f_mmax, 'rb') as f:
                load_from_file = pickle.load(f)
                Mmax_faults = load_from_file[0]
                Mmax_scenario = load_from_file[1]
                final_fault_length = load_from_file[2]
                final_scenario_length = load_from_file[3]
                faults_alone = load_from_file[4]
                faults_isolated = load_from_file[5]
                index_faults_in_scenario = load_from_file[6]
                scenario_area = load_from_file[7]

                if np.size(scenarios_names) == 0 :
                    Mmax = max(Mmax_faults)
                else :
                    Mmax = max((max(Mmax_faults),max(Mmax_scenario)))
                print("\t - scenario and max loaded")
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
            log_rup_file = open(self.pathlog +'/ruptures.txt','w')
            log_rup_file.write('rup_id\tinvolved_faults\n')
            for i in range(len(rup_rates)):
                log_rup_file.write(str(rup_rates.get(str(i)).get('rup_id'))+'\t')
                for j in rup_rates.get(str(i)).get('involved_faults'):
                    log_rup_file.write(str(j)+' ')
                log_rup_file.write('\n')
            log_rup_file.close()


        self.index_faults_in_scenario = index_faults_in_scenario


        '''#####################################################################
        # For each bin, find which fault and which scenario populates it.
        #####################################################################'''

        f_bin_pop = run_name+'/LOG/'+model_name +'_bin_pop_'+scl_name+'_'+set_name+'_'+str(self.sample)+'.pkl'
        if not os.path.isfile(f_bin_pop):
            re_use = False
        rup_in_bin = populate_bins.pop(bin_mag,index_rup,rup_rates,M_min,re_use,f_bin_pop)


        '''##################################################################
        # Definition of the background - how much is not on the faults
        ##################################################################'''

        # proportion of the BG
        bin_mag_fault_prop = [ 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.]
        fault_prop_inc = self.bg_ratio

        bin_mag_fault_prop.append(10.)
        fault_prop_inc = np.append(np.array(fault_prop_inc),1.)
        fault_prop = interp1d(bin_mag_fault_prop,fault_prop_inc)
        self.fault_prop = fault_prop



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

        elif self.mfd_hyp == 'tapered_GR' :
            # NOTE : The value of the corner magnitude could be computed in the future.
            # Kagan : Geophys. J. Int. (2002) 149, 731–754
            mfd_param['M_corner'] = Mmax-0.4
            p_MFD = mfd_shape.tapered_GR(mfd_param,bin_mag)

        elif self.mfd_hyp == 'double_GR' :
            p_MFD = mfd_shape.double_GR(mfd_param,bin_mag)

        elif self.mfd_hyp == 'YC':
            p_MFD = mfd_shape.YC(mfd_param,bin_mag,Mmax)

        elif self.mfd_hyp == 'YC_marmara':
            p_MFD = mfd_shape.YC_marmara(mfd_param,bin_mag,Mmax)

        elif self.mfd_hyp == 'YC_modified':
            p_MFD = mfd_shape.YC_modified(mfd_param,bin_mag,Mmax)

        elif self.mfd_hyp == 'UCERF_DV':
            # MFD calculated from the sum of the mfds of the faults of the UCERF
            # model for death valley
            p_MFD = mfd_shape.UCERF_DV(bin_mag)
            #can be used as a template for any shape numerically defined

        else:
            print("Error : Unknown MFD name!")

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
        target_moment_per_bin = p_MFD_MO * Total_moment_faults_rate_init
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
#        min_budget = 10**10
        for index_fault in range(len(faults_names)) :
            nb_dsr= int(round(faults_slip_rates[index_fault] / size_of_increment,0))  #nb of dsr to spend
            if nb_dsr == 0 :
                nb_dsr = 1
            faults_budget.update({index_fault:nb_dsr})
#            if nb_dsr < min_budget :
#                min_budget = nb_dsr
        min_budget = float(min(faults_budget.values())) + 1
        max_budget = float(max(faults_budget.values()))
        while max_budget/min_budget > 50.:
            min_budget *= 2.

        '''##################################################################
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        #   Populate the bin of magnitude of each fault and each scenario
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        ##################################################################'''

        M_slip_repartition = {} #how the slip rate is used
        # for each fault, contains all the rupture and the number of time each is picked.
        for fault,i in zip(faults_names,range(len(faults_names))) :

            dic_tmp = {}
            for rup_i in range(len(rup_rates)):
                if i in rup_rates.get(str(rup_i)).get('involved_faults'):
                    dic_tmp.update({str(rup_i):0})

            dic_tmp.update({'NMS':0})

            M_slip_repartition.update({str(fault):dic_tmp})

        moment_rate_in_bin = np.zeros(len(bin_mag))   # moment aleardy present in the bin of magnitude
        model_MFD = []

        number_of_loops = 0   # variable to follow how the calculation is doing
        number_of_loops_before = 0
        number_of_loops_for_nothing = 0
        number_of_loops_for_nothing_before = 0
        number_of_loops_last_checked = -1
        empty_bins = []
        empty_rups = []
        bin_target_reached = []
        len_faults_budget = [] # will be use to check if the calculation is stuck and too slow (see under)
        aseismic_count = 0

        color_mag= []


        '''######################
        #####   MAIN LOOP   #####
        ######################'''
        TARGET = []
        nb_ss_to_spend = float(sum(faults_budget.values()))
        sum_fault_budget = nb_ss_to_spend
        print("Number of dsr to spend : "+ str(nb_ss_to_spend))
        print("Min of sdr :",min(faults_budget.values()))
        print("Max of sdr :",max(faults_budget.values()))
        self.calculation_log_file.write("\nnumber of dsr to spend : "+ str(nb_ss_to_spend)+"\n")

        # !!!!!!!!!!!!!!!!!
        # Booleans used in the calculation, DO NOT MODIFY
        print_percent = True
        do_the_target = True
        bool_target_set = False
        # !!!!!!!!!!!!!!!!!

        # If uniform_spending is activated, the slip-rate budget
        # is not spent uniformly for the section participating to a rupture,
        # faster faults will slip-more
        uniform_spending = self.param["main"]["parameters"]["uniform_spending"]
        if uniform_spending in ["True", "true"] :
            uniform_spending = True
        else :
            uniform_spending = False

        # If deep analysis is turned on, display intermediate values in the prompt
        deep_analysis = self.param["main"]["parameters"]["deep_analysis"]
        if deep_analysis in ["True", "true"] :
            deep_analysis = True
        else :
            deep_analysis = False

        # Option to not to the weighting of the ruptures at every loop.
        faster_rup_weight = self.param["main"]["parameters"]["faster_rup_weight"]
        if faster_rup_weight in ["True", "true"] :
            faster_rup_weight = True
        else :
            faster_rup_weight = False


        # loop several time if the faults are much faster than the min ones
        # saving time on the random sampling and the calculation of the weights
        option_fast = self.param["main"]["parameters"]["option_fast"]
        if option_fast in ["True", "true"] :
            option_fast = True
        else :
            option_fast = False


        # if local MFD should also be respected
        # in this case, the
        local_MFD = self.param["main"]["parameters"]["local_MFD"]
        if local_MFD in ["True", "true"] :
            f_mfd_area = self.param["main"]["local_MFD_file"]
            local_mfds, associated_rup, associated_weight = core_utils.link_rup_mfd_area(rup_rates,f_mfd_area,self.faults_lon,self.faults_lat,bin_mag,self.bg_ratio)

        #log the time used for several parts
        time_weight_rupt = 0.
        time_clean_w_rupt = 0.
        time_target_building = 0.
        time_checking_target_reach = 0.
        time_spending_dsr = 0.
        time_checking_empty_bin = 0.
        time_checking_empty_faults = 0.


        # Does the temporary weighting
        budget_init = int(sum_fault_budget) #sum(faults_budget.values())

        if faster_rup_weight == True :
            nb_weigthings_rup_sampling = int(self.param["main"]["parameters"]["nb_weigthings_rup_sampling"])
            if "type_weigthings_rup_sampling" in self.param["main"]["parameters"].keys():
                type_weigthings_rup_sampling = self.param["main"]["parameters"]["type_weigthings_rup_sampling"]
            else :
                type_weigthings_rup_sampling = "lin"
            if type_weigthings_rup_sampling == "log":
                weigthing_built = [int(i) for i in np.logspace(0.,np.log10(budget_init),nb_weigthings_rup_sampling)]
            if type_weigthings_rup_sampling == "lin":
                weigthing_built = [int(i) for i in np.linspace(0.,budget_init,nb_weigthings_rup_sampling)]
        else :
            weigthing_built = [int(i) for i in range(budget_init)]

        weigthing_built.reverse()

        weigth_rup_sample = 0

        slip_rate_use_per_fault = np.zeros(len(faults_names))
        moment_rate_required = 0.
        moment_rate_left= Total_moment_faults_rate_init

        rate_tot_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        rate_in_model = np.zeros(len(bin_mag))

        test_mean_picked = []
        most_likely_pick = []

        n_w_work = 0
        n_w_crash = 0
        loop_last_rup_w = 0
        picked_empty_rup = 0
        old_percent = '0000'

        while sum_fault_budget != 0 : # sum(faults_budget.values()) as long as there is some slip-rate to spend we keep going
            ratio_done = 1. - float(sum_fault_budget)/nb_ss_to_spend
            if ratio_done > 0.01 :
                model_MFD, self.calculation_log_file,print_percent = core_utils.progress(model_MFD,self.calculation_log_file,ratio_done,print_percent,rup_rates,fault_prop,bin_mag)

            number_of_loops += 1

            if len(empty_bins) != len(bin_mag):
                ''' Calculate the new target shape in each bin in terms of moment rate '''
                tmp = time.time()
                target_i = target.get_new_target(number_of_loops,moment_rate_in_bin,p_MFD_MO,target_moment_per_bin,bin_mag,empty_bins,bin_target_reached,rup_in_bin)
                time_target_building += time.time() - tmp
                if sum(target_i) == 0.:
                    target_i = p_MFD_MO
                try :
                    most_likely_pick.append(bin_mag[list(target_i).index(max(list(target_i)))])
                except :
                    print(target_i)
                    most_likely_pick.append(bin_mag[list(target_i).index(max(list(target_i)))])

                ########
                # deep analysis mode : display intermediate values of variable here
                ########
                if deep_analysis == True:
                    percent = round((1.-float(sum_fault_budget)/float(nb_ss_to_spend)) * 100.)
                    percent = '{:04d}'.format(percent)
                    #if str(number_of_loops)[-4:] == "0000" or str(number_of_loops)[-4:] == "5000" :
                    if percent != old_percent :
                        old_percent = percent
                        print("\nnumber_of_loops",number_of_loops)
                        print("budget left : ",sum_fault_budget," | ",percent,"%")
                        time_str = core_utils.seconds_to_str(time_target_building)
                        print("time building target at time i : ",time_str)
                        time_str = core_utils.seconds_to_str(time_weight_rupt)
                        print("time weighting rupture pick : ",time_str)
                        time_str = core_utils.seconds_to_str(time_clean_w_rupt)
                        print("time cleaning weighting rupture weigth : ",time_str)
                        time_str = core_utils.seconds_to_str(time_checking_target_reach)
                        print("time checking target reach : ",time_str)
                        time_str = core_utils.seconds_to_str(time_checking_empty_bin)
                        print("time checking empty bins : ",time_str)
                        time_str = core_utils.seconds_to_str(time_checking_empty_faults)
                        print("time checking empty faults : ",time_str)
                        time_str = core_utils.seconds_to_str(time_target_building)
                        print("time spending dsr : ",time_str)
                        time_str = core_utils.seconds_to_str(time.time()-inti_core_time)
                        print("total core time : ",time_str)
                        tot_core_time = time.time()-inti_core_time
                        unaccounted_t = tot_core_time - (time_target_building+
                        time_weight_rupt+
                        time_clean_w_rupt+
                        time_checking_target_reach+
                        time_checking_empty_bin+
                        time_checking_empty_faults+
                        time_spending_dsr)
                        time_str = core_utils.seconds_to_str(unaccounted_t)
                        print("unaccounted time : ", time_str)
                        print("max target : ", round(max(target_i),4),"| last bin w : ",str(round(target_i[-1],6)))
                        print("Empty mag bins:")
                        print(empty_bins)
                        budget_last_bin = 0
                        for rup_i in rup_in_bin[-1]:
                            if not str(rup_i) in empty_rups:
                                for f_i in rup_rates.get(str(rup_i)).get('involved_faults'):
                                    budget_last_bin += faults_budget[f_i]
                        print("fault budget last bin :" , budget_last_bin)

                        fig, (ax0, ax1) = plt.subplots(ncols=2)
                        ax0.plot(bin_mag,target_i,label=("target_i"),marker='x')
                        ax0.plot(bin_mag,p_MFD_MO,label=("p_MFD_MO"),marker='x')
                        ax0.legend()

                        ax1.plot(bin_mag,target_moment_per_bin,label=("target_moment_per_bin"),marker='x')
                        ax1.plot(bin_mag,moment_rate_in_bin,label=("moment_rate_in_bin"),marker='x')
                        ax1.set_yscale('log')
                        ax1.legend()

                        plt.savefig(self.pathlog +'/Log_tmp_' + str(self.sample) + '.png' ,dpi = 80, transparent=True)
                        plt.close()


                '''Pick the bin of magnitude to fill according to the current distribution '''
                # normalize the target mfd
                target_i = (target_i) / sum(target_i)
                picked_bin = np.random.choice(len(bin_mag), 1, p = target_i)[0]
                mag = bin_mag[picked_bin] #magnitude of that bin

                tmp = time.time()
                # if len(empty_bins) != 0 :
                # testing if bin is not empty
                if not picked_bin in empty_bins:
                    empty_in_bin = set(rup_in_bin[picked_bin]) & set(empty_rups)
                    if len(empty_in_bin) == len(rup_in_bin[picked_bin]):
                        empty_bins.append(picked_bin)
                # not_empty_rup = []
                # for i_rup in rup_in_bin[picked_bin]:
                #     if not str(i_rup) in empty_rups :
                #         not_empty_rup.append(i_rup)
                # if len(not_empty_rup) == 0 :
                #     empty_bins.append(picked_bin)
                time_checking_empty_bin += time.time() - tmp


                if not picked_bin in empty_bins:
                    tmp = time.time()
                    '''Calculate the weight for sampling of the fault or scenario'''

                    do_rup_weight = False
                    if number_of_loops == 1 :
                        do_rup_weight = True
                    if sum_fault_budget < weigthing_built[weigth_rup_sample]:
                        do_rup_weight = True
                    # if number_of_loops - loop_last_rup_w > 50 :
                    #     do_rup_weight = True

                    if do_rup_weight == True :
                        loop_last_rup_w = number_of_loops
                        weigth_rup_sample += 1

                        if faster_rup_weight == True :
                            list_of_bins = range(len(bin_mag))
                        else :
                            list_of_bins = [picked_bin]

                        w_rup_binmag = []
                        for index_mag in list_of_bins:
                            if index_mag in empty_bins :
                                w_rup_binmag.append([])
                            else :
                                weight_rup_i = core_utils.weight_fault_sampling(index_mag,
                                rup_in_bin,faults_names,faults_slip_rates,slip_rate_use_per_fault,
                                faults_alone,scenarios_names,faults_isolated,index_faults_in_scenario,
                                rup_rates,empty_rups)

                                if local_MFD == True :
                                    # we check if local a local MFD is followed, if not, the ruptures that can help
                                    # to fit the local MFD are boosted or demotted.
                                    factor_on_weight = core_utils.check_local_mfd(rup_rates,
                                     rup_in_bin[index_mag], index_mag, bin_mag, local_mfds,
                                      associated_rup, associated_weight)
                                    weight_rup_i = np.array([i*w for i,w in zip(weight_rup_i,factor_on_weight)])
                                    weight_rup_i /= weight_rup_i.sum()

                                if sum(weight_rup_i) != 0. :
                                    weight_rup_i  = [float(i)/sum(weight_rup_i) for i in weight_rup_i]

                                w_rup_binmag.append(weight_rup_i)

                        time_weight_rupt += time.time() - tmp

                        # Cleaning up the weights
                        # to avoid nans and doesn't sum to 1
                        tmp = time.time()

                        for index_mag in list_of_bins:
                            if not index_mag in empty_bins :
                                weight_rup = w_rup_binmag[index_mag]
                                if (set(rup_in_bin[index_mag]) & set(empty_rups)) :
                                    i=0
                                    for i_rup in rup_in_bin[index_mag]:
                                        if str(i_rup) in empty_rups : #the rup is empty
                                            weight_rup[i]= 0.
                                        i+=1

                                weight_rup = list(weight_rup)
                                sum_weight_rup = sum(weight_rup)
                                if sum(weight_rup) == 0. :
                                    # print("Sum of the rupture weights in 0")
                                    # exit()
                                    empty_bins.append(index_mag)

                                if not index_mag in empty_bins :
                                    if sum_weight_rup != 1.0:
                                        weight_rup = [float(i)/sum_weight_rup for i in weight_rup]
                                    if math.isnan(sum(weight_rup)):
                                        print("WARNING : sum rup weight is nan")
                                        nb_nans = 0
                                        id = 0
                                        for i in weight_rup:
                                            if math.isnan(i):
                                                nb_nans += 1
                                            weight_rup[id] = 0.
                                            id += 1

                                w_rup_binmag[index_mag] = weight_rup

                        time_clean_w_rupt += time.time() - tmp





                    if faster_rup_weight == True :
                        weight_rup = w_rup_binmag[picked_bin]
                    else :
                        weight_rup = w_rup_binmag[0]
                    # i = 0
                    # for i_rup in rup_in_bin[picked_bin]:
                    #     if str(i_rup) in empty_rups : #the rup is empty
                    #         if weight_rup[i]!=0:
                    #             print("empty with weight",i_rup,weight_rup[i])
                    #     i+=1

                    #picking of the source'''
                    try :
                        #picked_rup = np.random.choice(rup_in_bin[picked_bin],1,p=weight_rup)[0] #picked source
                        i_picked = np.where(np.random.multinomial(1, weight_rup)==1)[0][0]
                        picked_rup = rup_in_bin[picked_bin][i_picked]
                        #print("index picked :",np.random.multinomial(1, weight_rup)[0],"  |  rup : ",picked_rup)
                        #n_w_work += 1

                    except ValueError:
                        print("rupt weights didn't work. sum:",sum(weight_rup))
                        picked_rup = np.random.choice(rup_in_bin[picked_bin])
                        n_w_crash += 1

                    #index_fault = np.where(np.array(faults_names) == picked_rup)[0]
                    index_fault = rup_rates.get(str(picked_rup)).get('involved_faults')
                    if bool_target_set == False :
                        tmp = time.time()
                        last_bins_empty = True
                        for bin_i in range(len(bin_mag))[-3:]:
                            if not bin_i in empty_bins :
                                last_bins_empty = False

                        if last_bins_empty == True : #if all the slip_rate has been spared for the Mmax - 0.3
                            rate_tot_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
                            bool_target_set = True
                            print("set target - limit on the largest magnitudes")
                            #all the slip_rate has NOT been spared for the Mmax
                            #does other checks
                            # do we have enough moment to fit the shape?
                        if moment_rate_left <= (1. - 0.00001)* moment_rate_required:
                            #if not enough momen rate left, stop spending on the larger EQs
                            self.calculation_log_file.write('\n Not enough moment left ')
                            for rup_i in rup_in_bin[-3]+rup_in_bin[-2]+rup_in_bin[-1]:
                                if not str(rup_i) in empty_rups :
                                    empty_rups.append(str(rup_i))
                            rup_in_bin[-3] = []
                            rup_in_bin[-2] = []
                            rup_in_bin[-1] = []
                            bool_target_set = True
                            print("set target - moment rate is limiting")
                        # check if the rates of the antepenultimate bin is not too different of the two other ones

                        antelast_bins_empty = True

                        for bin_i in range(len(bin_mag))[-2:]:
                            if not bin_i in empty_bins :
                                antelast_bins_empty = False

#                        for rup_i in rup_in_bin[-1]+rup_in_bin[-2]:
#                            if antelast_bins_empty == True :
#                                if not str(rup_i) in empty_rups :
#                                    antelast_bins_empty = False
                        if antelast_bins_empty == True and bool_target_set == False:
                            if moment_rate_in_bin[-3] >= 2. * (moment_rate_in_bin[-2]+moment_rate_in_bin[-1]):
                                self.calculation_log_file.write('\n antepenultimate bin getting too high')

                                for rup_i in rup_in_bin[-3]:
                                    if not str(rup_i) in empty_rups :
                                        empty_rups.append(str(rup_i))
                                empty_bins.append(range(len(bin_mag))[-3])
                                #rup_in_bin[-3] = []
                                bool_target_set = True

                                print("set target - antepenultimate bin getting too high")

                        #test if the moment rate left  to spend is enough to fit the target MFD
                        if number_of_loops > number_of_loops_last_checked + 50.:
#                            time_ii = time.time()
                            number_of_loops_last_checked = number_of_loops

                            moment_rate_left = Total_moment_faults_rate_init - Total_moment_rate_fault
                            rate_tot_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
                            rate_Mmax_check = np.mean(rate_tot_model[-3:])

                            moment_rate_required = 0.
                            for index_mag in range(len(bin_mag)-3):#loop on all the magnitude bins except the three last ones
                                rate_Mi_check = rate_tot_model[index_mag]
                                target_GR_i_check = rate_Mmax_check  * p_MFD[index_mag] / p_MFD[-2]
                                #difference between the target at this point and what is already filled in terms of moment rate
                                moment_rate_required += ((mag_to_M0(bin_mag[index_mag]) * target_GR_i_check) - (mag_to_M0(bin_mag[index_mag]) * rate_Mi_check)) * fault_prop(bin_mag[picked_bin])

                            if self.mfd_hyp == 'UCERF_DV_':
                                rate_Mmax = rate_tot_model[-1]
                                moment_rate_required = 0.
                                for index_mag in range(len(bin_mag)-3):#loop on all the magnitude bins except the three last ones
                                    rate_Mi_check = rate_tot_model[index_mag]
                                    target_GR_i_check = rate_Mmax_check  * p_MFD[index_mag] / p_MFD[-3]
                                    #difference between the target at this point and what is already filled in terms of moment rate
                                    moment_rate_required += ((mag_to_M0(bin_mag[index_mag]) * target_GR_i_check) - (mag_to_M0(bin_mag[index_mag]) * rate_Mi_check)) * fault_prop(bin_mag[picked_bin])

                        time_checking_target_reach += time.time() - tmp


                        if do_the_target == True and bool_target_set ==True:
                            ####
                            # Setting the target
                            ####

                            rate_tot_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)

                            do_the_target = False
                            print('- target set - ')
                            self.calculation_log_file.write('\n- target set - ')
                            rate_at_target_setting = rate_tot_model
                            #record the MFD at the moment the target is set

                            rate_Mmax = np.mean(rate_tot_model[-3:])

                            TARGET = []
                            for t_mag_bin in range(len(bin_mag)) :
                                TARGET.append(rate_Mmax  * p_MFD[t_mag_bin] / p_MFD[-2]) #target of the picked bin in order to fit the original distribution

                            if self.mfd_hyp == 'UCERF_DV_':
                                rate_Mmax = rate_tot_model[-1]

                                TARGET =[]
                                for t_mag_bin in range(len(bin_mag))  :
                                    rate_Mi = rate_tot_model[t_mag_bin]
                                    TARGET.append(rate_Mmax  * p_MFD[picked_bin] / p_MFD[-3]) #target of the picked bin in order to fit the original distribution


                    if picked_bin in bin_target_reached:
                        print("WHAT ARE YOU DOING HERE?",bin_mag[picked_bin])
                        exit()

                    ''' spending the slip_rate increment '''
                    tmp =time.time()
                    index_fault = rup_rates.get(str(picked_rup)).get('involved_faults')
                    sr_to_spend = True #all of the faults involved still have slip rate to spend
                    shear_mod = 0
                    for index in index_fault :
                        shear_mod += faults_shear_mod[index]
                        if faults_budget[index] <= 0.:
                            sr_to_spend = False #one of the faults don't have anymore budget

                    if sr_to_spend == False :
                        if not str(picked_rup) in empty_rups:
                            empty_rups.append(str(picked_rup))
                        else :
                            picked_empty_rup += 1

                    if sr_to_spend == True :
                        shear_mod = shear_mod / float(len(index_fault))
                        area = rup_rates.get(str(picked_rup)).get('area')
                        displacement = mag_to_M0(mag)/(shear_mod*area)
                        rate_i = size_of_increment/displacement

                        if option_fast == True :
                            min_budget_local = min([faults_budget[i] for i in index_fault])
                            nb_loop_spending = int(min_budget_local/min_budget)
                            if nb_loop_spending < 1:
                                nb_loop_spending = 1
                        else :
                            nb_loop_spending = 1


                        if bool_target_set == True:
                            rate_Mi = rate_in_model[picked_bin]
                            target_mfd_i = TARGET[picked_bin] * fault_prop(bin_mag[picked_bin])#target of the picked bin in order to fit the original distribution

                            if rate_Mi < target_mfd_i : #for this bin, the target hasn't been reached yet

                                if uniform_spending == True or len(index_fault)==1:
                                    moment_rate_i = 0.
                                    for index in index_fault :
                                        for loop_spending in range(nb_loop_spending):
                                            #M_slip_repartition[index].append(picked_rup)
                                            M_slip_repartition[str(faults_names[index])][str(picked_rup)] += 1
                                        faults_budget[index] += -1 * nb_loop_spending
                                        sum_fault_budget += -1 * nb_loop_spending
                                        slip_rate_use_per_fault[index] += size_of_increment * nb_loop_spending
                                    rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i * nb_loop_spending
                                    rate_in_model[picked_bin] += rate_i * nb_loop_spending
                                    moment_rate_i += mag_to_M0(mag) * rate_i * nb_loop_spending
                                else :
                                    moment_rate_i = 0.
                                    for loop_spending in range(nb_loop_spending):
                                        M_slip_repartition,faults_budget,slip_rate_use_per_fault,nb_sdr_used,sum_fault_budget = core_utils.variable_spending(index_fault,M_slip_repartition,faults_budget,slip_rate_use_per_fault,size_of_increment,faults_slip_rates,picked_rup,faults_names,sum_fault_budget)
                                            #adding to the rate
                                        rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i*(nb_sdr_used)
                                        rate_in_model[picked_bin] +=  rate_i*(nb_sdr_used)
                                        moment_rate_i += mag_to_M0(mag) * rate_i*(nb_sdr_used)

                                #substracting the moment used from the target
                                moment_rate_in_bin[picked_bin] += moment_rate_i #* ( 2.-fault_prop(mag))
                                Total_moment_rate_fault += moment_rate_i
                            else : # for this bin, the target has been reached, this slip rate needs to be aseismic
                                rup_in_bin[picked_bin] = []
                                if not picked_bin in bin_target_reached:
                                    bin_target_reached.append(picked_bin)


                        else : # the absolute target as not been reached yet

                            if uniform_spending == True or len(index_fault)==1:
                                moment_rate_i = 0.
                                for index in index_fault :
                                    for loop_spending in range(nb_loop_spending):
                                        #M_slip_repartition[index].append(picked_rup)
                                        M_slip_repartition[str(faults_names[index])][str(picked_rup)] += 1
                                    faults_budget[index]+=-1 * nb_loop_spending
                                    sum_fault_budget += -1 * nb_loop_spending
                                    slip_rate_use_per_fault[index] += size_of_increment * nb_loop_spending
                                rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i * nb_loop_spending
                                rate_in_model[picked_bin] +=  rate_i * nb_loop_spending
                                moment_rate_i += mag_to_M0(mag) * rate_i * nb_loop_spending
                            else :
                                moment_rate_i = 0.
                                for loop_spending in range(nb_loop_spending):
                                    M_slip_repartition,faults_budget,slip_rate_use_per_fault,nb_sdr_used,sum_fault_budget = core_utils.variable_spending(index_fault,M_slip_repartition,faults_budget,slip_rate_use_per_fault,size_of_increment,faults_slip_rates,picked_rup,faults_names,sum_fault_budget)
                                        #adding to the rate
                                    rup_rates[str(picked_rup)]['rates'][picked_bin] += rate_i*(nb_sdr_used)
                                    rate_in_model[picked_bin] +=  rate_i*(nb_sdr_used)
                                    moment_rate_i += mag_to_M0(mag) * rate_i*(nb_sdr_used)

                            #substracting the moment used from the target
                            moment_rate_in_bin[picked_bin] += moment_rate_i #* ( 2.-fault_prop(mag))
                            Total_moment_rate_fault += moment_rate_i

                        #rate_in_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
                    time_spending_dsr += time.time() - tmp

                else :
                    #print("picked empty bin")
                    if not picked_bin in empty_bins:
                        empty_bins.append(picked_bin)


#                if number_of_loops > number_of_loops_before+200:
#                    number_of_loops_before = number_of_loops
                    #if the fault has no more slip to spare, the scenarios and the fault are remove for the picking in order to have a faster picking


                tmp = time.time()
                # Checking is the bin is empty
                if number_of_loops > number_of_loops_before+500:
                    number_of_loops_before = number_of_loops
                    for index_mag in range(len(bin_mag)):
                        if not index_mag in empty_bins :
                            nb_rup_not_empty = 0
                            # empty_in_bin = set(rup_in_bin[index_mag]) & set(empty_rups)
                            # empty_in_bin = [i for i in rup_in_bin[index_mag] if i in empty_rups]
                            # if len(empty_in_bin) == len(rup_in_bin[index_mag]):
                            #     empty_bins.append(picked_bin)
                            for i_rup in rup_in_bin[index_mag] :
                                if nb_rup_not_empty == 0 :
                                    rup_is_empty = False
                                    for index_fault in rup_rates.get(str(i_rup)).get('involved_faults'):
                                        if faults_budget[index_fault] <= 0:
                                            rup_is_empty = True
                                    if rup_is_empty == False :
                                        nb_rup_not_empty += 1
                                        #rup_in_bin[index_mag].remove(i_rup)
                            if nb_rup_not_empty == 0 :
                                empty_bins.append(index_mag)
                time_checking_empty_bin += time.time() - tmp

                tmp = time.time()
                len_faults_budget.append(sum_fault_budget)
                #does a check to see if the remainingfaults that still have some slip rate fit in one of the bin of are too small
                if len(len_faults_budget) > 3 :
                    if len_faults_budget[-2] == len_faults_budget[-1] :
                        number_of_loops_for_nothing+=1
                    if number_of_loops_for_nothing_before<number_of_loops_for_nothing-100:
                        number_of_loops_for_nothing_before = number_of_loops_for_nothing
                        if len_faults_budget[-1] == len_faults_budget[-10] :
                            # rup_still_used = []
                            # for rup_i in rup_in_bin[picked_bin]:
                            #     rup_still_used += rup_i
                            # rup_still_used = set(rup_still_used)
                            rup_still_used = [i for i in rup_in_bin[picked_bin] if not i in empty_rups]
                            fault_still_used = []
                            for rup_i in rup_still_used:
                                fault_still_used += list(rup_rates.get(str(rup_i)).get('involved_faults'))
                            fault_still_used = set(fault_still_used)
                            for fault,index_fault in zip(faults_names,range(len(faults_names))) :
                                if not (index_fault in fault_still_used) and (faults_budget[index_fault] > 0):
                                    while faults_budget[index_fault] > 0 :
                                        ratio_done = 1. - float(sum_fault_budget)/nb_ss_to_spend
                                        if ratio_done > 0.01 :
                                            model_MFD, self.calculation_log_file,print_percent = core_utils.progress(model_MFD,self.calculation_log_file,ratio_done,print_percent,rup_rates,fault_prop,bin_mag)
                                        faults_budget[index_fault]+=-1
                                        sum_fault_budget+=-1
                                        #M_slip_repartition[index_fault].append('NMS')
                                        M_slip_repartition[str(faults_names[index_fault])]['NMS'] += 1
                                        aseismic_count += 1
                time_checking_empty_faults += time.time() - tmp
            else:
                print('-target filled-')
                self.calculation_log_file.write('\n-target filled-')
                while sum_fault_budget!=0: # as long as there is some slip-rate to spend we keep going

                    ratio_done = 1. - float(sum_fault_budget)/nb_ss_to_spend
                    if ratio_done > 0.01 :
                        model_MFD, self.calculation_log_file,print_percent = core_utils.progress(model_MFD,self.calculation_log_file,ratio_done,print_percent,rup_rates,fault_prop,bin_mag)
                    for index_fault in range(len(faults_names)):
                        if faults_budget[index_fault] > 0 :
                            faults_budget[index_fault]+=-1
                            sum_fault_budget+=-1
#                            M_slip_repartition[index_fault].append('NMS')
                            M_slip_repartition[str(faults_names[index_fault])]['NMS'] += 1
                            aseismic_count += 1
        ''' check if the TARGET as been set.
        if not, build it for comparing'''

        if TARGET == []:

            rate_tot_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)

            print('- target set at the end- ')
            self.calculation_log_file.write('\n- target set at the end - ')
            rate_at_target_setting = rate_tot_model
            #record the MFD at the moment the target is set

            rate_Mmax = np.mean(rate_tot_model[-3:])

            for t_mag_bin in range(len(bin_mag)) :
                TARGET.append(rate_Mmax  * p_MFD[t_mag_bin] / p_MFD[-2])

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
        rate_tot_model = rates.get_rate_model(rup_rates,fault_prop,bin_mag)
        model_MFD.append(rate_tot_model)
        colors = ['gainsboro','darkgray','dimgray','dimgray','black','red']
        index_color = 0
        for MFD_i in model_MFD:
            plt.plot(bin_mag,MFD_i,color=colors[index_color])
            index_color+=1
        if do_the_target == True:
            rate_at_target_setting = model_MFD[-1]
        plt.plot(bin_mag,rate_at_target_setting,':g')

#        target_i = []
#        for index_mag in range(len(bin_mag)):
##            rate_Mmax = (model_MFD[-1][-3]+model_MFD[-1][-2]+model_MFD[-1][-1])/3.
##            target_GR_i = rate_Mmax  * p_MFD[index_mag] / p_MFD[-3:]
#            target_GR_i = TARGET[index_mag]
#
##            if self.mfd_hyp == 'UCERF_DV':
##                rate_Mmax = model_MFD[-1][-3] #rate of the Mmax due to faults
##                target_GR_i = rate_Mmax  * p_MFD[index_mag] / p_MFD[-3]
#            target_i.append(target_GR_i)
        plt.plot(bin_mag,TARGET,':b')
        plt.scatter(bin_mag[-3:],rate_tot_model[-3:],c='k')
        plt.yscale('log')
        plt.savefig(self.pathlog +'/Log_target_fit_' + str(self.sample) + '.png' ,dpi = 180, transparent=True)
        plt.close()

        index_5 = 0
        while bin_mag[index_mag] < 5.:
            index_5+=1
        self.ratio_test = np.mean([abs(TARGET[index_5+0]/rate_tot_model[index_5+0]),
                                   abs(TARGET[index_5+1]/rate_tot_model[index_5+1]),
                                    abs(TARGET[index_5+2]/rate_tot_model[index_5+2]),
                                    abs(TARGET[index_5+3]/rate_tot_model[index_5+3]),
                                    abs(TARGET[index_5+4]/rate_tot_model[index_5+4])])
        MFD_to_test = (rate_tot_model)/sum(rate_tot_model)
        target_to_test = (TARGET)/sum(TARGET)
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
        self.rup_rates = rup_rates
        self.M_slip_repartition = M_slip_repartition


        log_calculation_file.write('Moment_rate_faults_final calculated using 10^(1.5M+9.1)' + '\t' + str(Total_moment_rate_fault) + '\t' + 'N.m' + '\n')
        log_calculation_file.write('Ratio of NMS slip rate counted in the iterative process: ' + '\t' + str(round((100.) * (aseismic_count/nb_ss_to_spend))) + '\t %\n')
        log_calculation_file.write('Moment_rate_model calculated using 10^(1.5M+9.1)' + '\t' + str(Total_moment_rate_fault) + '\t' + 'N.m' + '\n')
        log_calculation_file.close()
        self.ratio_NMS = round((100.) * (1. - Total_moment_rate_fault / Total_moment_faults_rate_init))
        print('ratio of NMS : '+str(round((100.) * (aseismic_count/nb_ss_to_spend))))
        self.calculation_log_file.write('\nratio of NMS : '+str(round((100.) * (aseismic_count/nb_ss_to_spend))))


#        for rup in M_slip_repartition:
#            for i in rup :
#                log_sliprep_file.write(str(i)+' ')
#            log_sliprep_file.write('\n')
#        log_sliprep_file.close()

        f = open(log_sliprep_file,"wb")
        pickle.dump(M_slip_repartition,f)
        f.close()
