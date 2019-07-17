#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches


def do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,
                 rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag): 
    for i in range(len(mfd_X)):
        plt.scatter(mega_bining_in_mag,mfd_X[i], c='darkcyan', s=50, edgecolor='',marker = '_',alpha = 0.5)
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    for index_mag in range(len(mega_bining_in_mag)): 
        rate_plus = np.percentile(mfd_X,84,axis=0)[index_mag]
        rate_minus = np.percentile(mfd_X,16,axis=0)[index_mag]
        mag = mega_bining_in_mag[index_mag]
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
        
        patch = patches.PathPatch(path_poly,facecolor = '#598556', lw = 0., alpha = 0.15)
        axes.add_patch(patch)
                    
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,50,axis=0),
                c='darkgreen', s=25, edgecolor='',marker = 'o',alpha = 0.8)
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,16,axis=0),
            c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,84,axis=0),
        c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
    plt.plot(mega_bining_in_mag,np.array(mfd_X).mean(axis=0),
                color='darkgreen', linewidth = 2)
    plt.grid()

    plt.yscale('log')
    plt.title('MFD ' + hyp_name)
    plt.savefig(path+ '/mdf_' + hyp_name + '_density.png' ,
                dpi = 180, transparent=True)
    #plt.show()
    plt.close()
    
    '''
    # plot linearely
    
    i_mag = 0

    for mag in [4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5]:
        #ploting the catalog
        for index_cat in range(len(rate_in_catalog)):
            mag_bin = [x + 0.01 for x in bining_in_mag]
#            print len(mag_bin)
#            print len(rate_in_catalog[index_cat])
            plt.scatter(mag_bin,rate_in_catalog[index_cat], c='k', s=50, edgecolor='',marker = '_', alpha = 0.25) 
           
        plt.scatter(mag_bin,np.percentile(rate_in_catalog,50,axis=0),
                    c='k', s=30, edgecolor='',marker = 'o',alpha = 0.8)
        plt.scatter(mag_bin,np.percentile(rate_in_catalog,16,axis=0),
                c='k', s=20, edgecolor='',marker = '+',alpha = 0.8)
        plt.scatter(mag_bin,np.percentile(rate_in_catalog,84,axis=0),
            c='k', s=20, edgecolor='',marker = '+',alpha = 0.8)
        plt.scatter(mag_bin,np.array(rate_in_catalog).mean(axis=0),
                    c='k', s=50, edgecolor='',marker = 's',alpha = 0.95)
        
        #plotting the modelled rates
        for i in range(len(mfd_X)):
            plt.scatter(mega_bining_in_mag,mfd_X[i], c='darkcyan', s=50, edgecolor='',marker = '_',alpha = 0.2)
        
        axes = plt.gca()
        axes.set_xlim([mag-0.05,mag+0.55])
        index_x0 = 0
        while bining_in_mag[index_x0] <= mag-0.05 and bining_in_mag[index_x0] != xmax :
            index_x0 +=1
            #print index_x0
        ymax = np.array(mfd_X).max(axis=0)[index_x0] + 0.2*np.array(mfd_X).max(axis=0)[index_x0]
        if ymax==0:
            ymax=0.00001
        axes.set_ylim([0.,ymax])
        mag_bin = [x - 0.01 for x in mega_bining_in_mag] 
        for index_mag in range(len(mag_bin)): 
            rate_plus = np.percentile(mfd_X,84,axis=0)[index_mag]
            rate_minus = np.percentile(mfd_X,16,axis=0)[index_mag]
            mag_i = mega_bining_in_mag[index_mag]
            mag_plus = mag_i+0.05
            mag_minus = mag_i-0.05
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
            
            patch = patches.PathPatch(path_poly,facecolor = '#598556', lw = 0., alpha = 0.15) 
            
               
            plt.scatter(mag_bin,np.percentile(mfd_X,50,axis=0),
                        c='darkgreen', s=25, edgecolor='',marker = 'o',alpha = 0.8)
            plt.scatter(mag_bin,np.percentile(mfd_X,16,axis=0),
                    c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
            plt.scatter(mag_bin,np.percentile(mfd_X,84,axis=0),
                c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
            plt.plot(mag_bin,np.array(mfd_X).mean(axis=0),
                        color='darkgreen', linewidth = 2)
                
            
    #        plt.scatter(mag_bin,np.percentile(mfd_X,50,axis=0),
    #                    c='darkgreen', s=30, edgecolor='',marker = 'o',alpha = 0.8)
    #        plt.scatter(mag_bin,np.percentile(mfd_X,16,axis=0),
    #                c='darkgreen', s=20, edgecolor='',marker = '+',alpha = 0.8)
    #        plt.scatter(mag_bin,np.percentile(mfd_X,84,axis=0),
    #            c='darkgreen', s=20, edgecolor='',marker = '+',alpha = 0.8)
    #        plt.scatter(mag_bin,np.array(mfd_X).mean(axis=0),
    #                    c='darkslateblue', s=50, edgecolor='',marker = 's',alpha = 0.95)
    #        axes = plt.gca()
        ymax = np.array(mfd_X).max(axis=0)[0] + 0.2 * np.array(mfd_X).max(axis=0)[0]
        axes.set_ylim([0.,ymax/float(12**i_mag)]) #depends of the b_value
        
        plt.grid() 
    
        #plt.yscale('log')
        plt.title('MFD ' + hyp_name)
        plt.savefig(path+ '/mdf_' + hyp_name + '_' + str(mag)+ '_' + str(mag+0.5) +'.png' ,
                    dpi = 100, transparent=True)
        #plt.show()
        plt.close()
        i_mag += 0.5
        '''
    if plot_as_rep == True :
        # plot the histogra of the aseismicity in this branch
        a_s_hyp = []
        for index in rows :
            a_s_hyp.append(a_s_model[index])
        
        plt.hist(a_s_hyp,20)
        plt.axis([0, 100, 0, max(plt.hist(a_s_hyp,20)[0] + 10)])
        plt.title('aseismic slip for bvalue ' +hyp_name)
        plt.savefig(path+ '/aseismic_slip.png' , dpi = 100, transparent=True)
        #plt.show()
        plt.close()

        
        
        
def plt_mfd(Run_name,mega_MFD, scenarios_names_list, ScL_complet_list, ScL_list, Model_list,BG_hyp_list,
             dimension_used_list,faults_name_list,sample_list,b_value_list,MFD_type_list,m_Mmax,
             mega_bining_in_mag,a_s_model,b_sample,sm_sample,Mt_sample,plot_mfd,plot_as_rep,plot_Mmax,xmin,xmax,ymin,ymax,
             catalog_cum_rate,plot_mfd_detailled,bining_in_mag):
    
    
    file_scenarios_MFD_name = str(Run_name) + '/analysis/txt_files/scenarios_MFD.txt'
    file_scenarios_MFD = open(file_scenarios_MFD_name,'w')
    
    if plot_mfd == True :
        for scenario in scenarios_names_list :
            
            mfds_scenario = []
            for mfd_i in mega_MFD:
                if mfd_i[8] == scenario:
                   mfds_scenario.append(mfd_i)
                   
            mfd_scenario_cumulative = []
            mfd_source_cummulative = []
            for mfd in mfds_scenario:
                mfd_i = mfd[11::].astype(np.float)
                
                mfd_source_cummulative_i = []
                for i in range(len(mfd_i)): #calculate the cumulative for each source
                    mfd_source_cummulative_i.append(np.sum(np.array(mfd_i)[-(len(mfd_i)-i):])) 
                mfd_source_cummulative.append(mfd_source_cummulative_i)
                
            for sample in sample_list:
                rows, cols = np.where(np.array(mfds_scenario) == sample) 
                
                mfds_scenario_sample = np.take(mfd_source_cummulative,rows,axis= 0)
                
                mfd_scenario_cumulative_sample = np.sum(mfds_scenario_sample,axis=0)
                mfd_scenario_cumulative.append(mfd_scenario_cumulative_sample)
                
                    
                file_scenarios_MFD.write(scenario + '\t' + str(mfd_scenario_cumulative_sample)+'\n')
        
    file_scenarios_MFD.close()
    
    
    #"#### plot for the whole tree
    file_branch_cumMFD_name = str(Run_name) + '/analysis/txt_files/branch_cumMFD.txt'        
    file_branch_cumMFD = open(file_branch_cumMFD_name,'w')
    
    mega_mfd_cummulative = [] #will contain the cummulative MFD for each model of the logic tree
    total_list_BG_hyp = []  #wil contain the list of the M_trunc for each model of the logic tree
    total_list_complet_ScL = []
    total_list_ScL = [] #wil contain the list of the ScL for each model of the logic tree
    total_list_dimension_used = [] #wil contain the list of the dimension used for each model of the logic tree
    total_list_b_value = []
    total_list_MFD_type = []
    total_list_scenario_name = []
    total_list_model = []
    total_list_sample = []
    
    geologic_moment_rate = [] # list of the moment rate of each model
    geologic_moment_rate_no_as = []  # list of the moment rate of each modelif no aseismic slip is considered
    
    selected_ScL = 'Init0'
    Dimention_used = 'Init0' 
    str_all_data = 'Init0' 
    Model = 'Init0' 
    BG_hyp = 'Init0'
    b_min = 'Init0'
    b_max = 'Init0' 
    MFD_type = 'Init0'
    scenario_name = 'Init0'
    sample = 'Init0'
    
    mfd_i = np.zeros(len(mega_MFD[0][11::]))
    
    index = 0
    for mega_mfd_i in mega_MFD :
        if (mega_mfd_i[0] == selected_ScL) and (mega_mfd_i[1] == Dimention_used) and (mega_mfd_i[2] == str_all_data) and (mega_mfd_i[3] == Model
        ) and (mega_mfd_i[4] == BG_hyp) and (mega_mfd_i[5] == b_min) and (mega_mfd_i[6] == b_max) and (mega_mfd_i[7] == MFD_type
        ) and (mega_mfd_i[8] == scenario_name)  and (mega_mfd_i[9] == sample): #same model, we add sources
            #print 'ok'
            mfd_i += mega_mfd_i[11::].astype(np.float)
            
        else : #it means it a new model
            if sum(mfd_i) != 0. : #we calculate the cumulative MFD
                mfd_cummulative_i = []
                geologic_moment_rate_i = 0.
                for i in range(len(mfd_i)): #calculate the cumulative for each source
                    mfd_cummulative_i.append(np.sum(np.array(mfd_i)[-(len(mfd_i)-i):]))
                    M0 = 10. ** (1.5 * mega_bining_in_mag[i] + 9.1)
                    rate_M0 = M0 * mfd_i[i]
                    geologic_moment_rate_i += rate_M0
                geologic_moment_rate.append(geologic_moment_rate_i)
                geologic_moment_rate_no_as.append(geologic_moment_rate_i * 100. / (100. - float(a_s_model[index])))
                
                mega_mfd_cummulative.append(mfd_cummulative_i)
                total_list_BG_hyp.append(BG_hyp)
                total_list_complet_ScL.append((str(selected_ScL) + '_' + str(Dimention_used) + '_' + str(str_all_data)))
                total_list_ScL.append(selected_ScL)
                total_list_dimension_used.append(Dimention_used)
                total_list_model.append(Model)
                total_list_b_value.append('bmin_'+str(b_min)+'_bmax_'+str(b_max))
                total_list_MFD_type.append(MFD_type)
                total_list_scenario_name.append(scenario_name)
                total_list_sample.append(sample)
                file_branch_cumMFD.write(str(Model) + '\t' + str(MFD_type) + '\t' + str(BG_hyp) + '\t' + str(scenario_name) + '\t' + str((str(selected_ScL) + '_' + str(Dimention_used) + '_' + str(str_all_data))) + '\t' + 'bmin_'+str(b_min)+'_bmax_'+str(b_max) + '\t' + str(sample) + '\t' + '\t'.join(map(str,mfd_cummulative_i)) + '\n')
                
                index += 1
                
            mfd_i = np.zeros(len(mega_mfd_i[11::]))
            selected_ScL = mega_mfd_i[0]
            Dimention_used = mega_mfd_i[1] 
            str_all_data = mega_mfd_i[2] 
            Model = mega_mfd_i[3] 
            BG_hyp = mega_mfd_i[4] 
            b_min = mega_mfd_i[5] 
            b_max = mega_mfd_i[6] 
            MFD_type = mega_mfd_i[7]
            #a_s = mega_mfd_i[8]
            scenario_name = mega_mfd_i[8] 
            sample = mega_mfd_i[9] 
            mfd_i += mega_mfd_i[11::].astype(np.float)

    #we write for the last model 
    mfd_cummulative_i = []
    geologic_moment_rate_i = 0.
    for i in range(len(mfd_i)): #calculate the cumulative for each source
        mfd_cummulative_i.append(np.sum(np.array(mfd_i)[-(len(mfd_i)-i):]))
        M0 = 10. ** (1.5 * mega_bining_in_mag[i] + 9.1)
        rate_M0 = M0 * mfd_i[i]
        geologic_moment_rate_i += rate_M0
    geologic_moment_rate.append(geologic_moment_rate_i)
    geologic_moment_rate_no_as.append(geologic_moment_rate_i * 100. / (100. - float(a_s_model[index])))       
    geologic_moment_rate.append(geologic_moment_rate_i)
    geologic_moment_rate_no_as.append(geologic_moment_rate_i * 100. / (100. - float(a_s_model[index])))
    
    mega_mfd_cummulative.append(mfd_cummulative_i)
    total_list_BG_hyp.append(BG_hyp)
    total_list_complet_ScL.append((str(selected_ScL) + '_' + str(Dimention_used) + '_' + str(str_all_data)))
    total_list_ScL.append(selected_ScL)
    total_list_dimension_used.append(Dimention_used)
    total_list_model.append(Model)
    total_list_b_value.append('bmin_'+str(b_min)+'_bmax_'+str(b_max))
    total_list_MFD_type.append(MFD_type)
    total_list_scenario_name.append(scenario_name)
    total_list_sample.append(sample)
    file_branch_cumMFD.write(str(Model) + '\t' + str(MFD_type) + '\t' + str(BG_hyp) + '\t' + str(scenario_name) + '\t' + str((str(selected_ScL) + '_' + str(Dimention_used) + '_' + str(str_all_data))) + '\t' + 'bmin_'+str(b_min)+'_bmax_'+str(b_max) + '\t' + str(sample) + '\t' + '\t'.join(map(str,mfd_cummulative_i)) + '\n')
    
    file_branch_cumMFD.close()
    
    
    if len(mega_mfd_cummulative) < 4 :
        plot_mfd = False       
    mfd_X = mega_mfd_cummulative
    for i in range(len(mfd_X)):
        plt.scatter(mega_bining_in_mag,mfd_X[i], c='darkcyan', s=50, edgecolor='',marker = '_',alpha = 0.5)
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    for index_mag in range(len(mega_bining_in_mag)): 
        rate_plus = np.percentile(mfd_X,84,axis=0)[index_mag]
        rate_minus = np.percentile(mfd_X,16,axis=0)[index_mag]
        mag = mega_bining_in_mag[index_mag]
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
        
        patch = patches.PathPatch(path_poly,facecolor = 'darkgreen', lw = 0., alpha = 0.15)
        axes.add_patch(patch)
                    
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,50,axis=0),
                c='darkgreen', s=25, edgecolor='',marker = 'o',alpha = 0.8)
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,16,axis=0),
            c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
    plt.scatter(mega_bining_in_mag,np.percentile(mfd_X,84,axis=0),
        c='darkgreen', s=60, edgecolor='',marker = '_',alpha = 0.8)
    plt.plot(mega_bining_in_mag,np.array(mfd_X).mean(axis=0),
                color='darkgreen', linewidth = 2)
    plt.grid()
    
    #plot the MFDs of the wholle tree with mean and percentiles    
#    for i in range(len(mega_mfd_cummulative)):
#        plt.scatter(mega_bining_in_mag,mega_mfd_cummulative[i], c='darkcyan', s=50, edgecolor='',marker = '_',alpha = 0.25)
#        
#    plt.scatter(mega_bining_in_mag,np.percentile(mega_mfd_cummulative,50,axis=0),
#                c='darkgreen', s=30, edgecolor='',marker = 'o',alpha = 0.8)
#    plt.scatter(mega_bining_in_mag,np.percentile(mega_mfd_cummulative,16,axis=0),
#            c='darkgreen', s=20, edgecolor='',marker = '+',alpha = 0.8)
#    plt.scatter(mega_bining_in_mag,np.percentile(mega_mfd_cummulative,84,axis=0),
#        c='darkgreen', s=20, edgecolor='',marker = '+',alpha = 0.8)
#    plt.scatter(mega_bining_in_mag,np.array(mega_mfd_cummulative).mean(axis=0),
#                c='darkslateblue', s=50, edgecolor='',marker = 's',alpha = 0.95) 
#    
#    
#    axes = plt.gca()
#    axes.set_xlim([xmin,xmax])
#    axes.set_ylim([ymin,ymax])
    plt.grid()
    plt.yscale('log')
    plt.title('MFD of the whole tree ')
    plt.savefig(str(Run_name) + '/analysis/figures/mfd/mdf_whole_tree.png' , dpi = 180, transparent=True)
    #plt.show()
    plt.close()
    
        
    rate_in_catalog = catalog_cum_rate

    #bining_in_mag = np.linspace(5.,7.5,26)
    '''##########################################
    #plot mfd for each scenario of the logic tree
    ############################################'''
    if len(scenarios_names_list)>1:
        index_model = 0
        for model in Model_list : 
            rate_in_catalog = catalog_cum_rate[index_model]
            for scenario in scenarios_names_list :
                
                if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/scenario_set/' + scenario):
                    os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/scenario_set/' + scenario)
                rows = np.where(np.array(total_list_scenario_name) == scenario)[0]  
                mfd_X = []
                for index in rows :  
                    mfd = mega_mfd_cummulative[index]
                    mfd_X.append(mfd)
                    
                #density plot
                if plot_mfd == True :   
                    
                    hyp_name = scenario
                    path = str(Run_name) + '/analysis/figures/analyze_branches/scenario_set/' + scenario
                    do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)

            index_model += 1
                
    '''##########################################
    #plot mfd for each model of the logic tree
    ############################################'''
    index_model = 0
    for model in Model_list : 
#        print catalog_cum_rate
#        print 
        rate_in_catalog = catalog_cum_rate[index_model]
        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model):
            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model)
        rows = np.where(np.array(total_list_model) == model)[0]  
        mfd_X = []
        for index in rows :  
            mfd = mega_mfd_cummulative[index]
            mfd_X.append(mfd)
            
        #density plot
        if plot_mfd == True :   
            hyp_name = model
            path = str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model
            do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)
        index_model +=1
    '''##########################################
    #plot mfd for each Background hypothesis  of the logic tree
    ############################################'''
    if len(BG_hyp_list) > 1:
        for BG_hyp in BG_hyp_list :        
            if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/BG/' + BG_hyp):
                os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/BG/' + BG_hyp)
            rows = np.where(np.array(total_list_BG_hyp) == BG_hyp)[0]  
            mfd_X = []
            index_check = 0
            for index in rows :  
                mfd = mega_mfd_cummulative[index]
                mfd_X.append(mfd)
                index_check += 1
                
            #density plot
            if plot_mfd == True :   
                hyp_name = BG_hyp
                path = str(Run_name) + '/analysis/figures/analyze_branches/BG/' + BG_hyp
                do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)

                
                                
    '''##########################################
    #plot mfd for each MFD  of the logic tree
    ############################################'''
    if len(MFD_type_list) > 1:
        for MFD_type in MFD_type_list :        
            if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/MFD_type/' + MFD_type):
                os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/MFD_type/' + MFD_type)
            rows = np.where(np.array(total_list_MFD_type) == MFD_type)[0]  
            mfd_X = []
            for index in rows :  
                mfd = mega_mfd_cummulative[index]
                mfd_X.append(mfd)
                
           #density plot
            if plot_mfd == True :  
                hyp_name = MFD_type
                path = str(Run_name) + '/analysis/figures/analyze_branches/MFD_type/' + MFD_type
                do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)

                
    '''##########################################
    #plot mfd for each bvalue  of the logic tree
    ############################################'''
    if len(b_value_list) > 1:
        for b in b_value_list :        
            if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/b_value/' + b):
                os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/b_value/' + b)
            rows = np.where(np.array(total_list_b_value) == b)[0]  
            mfd_X = []
            for index in rows :  
                mfd = mega_mfd_cummulative[index]
                mfd_X.append(mfd)
                
            #density plot
            if plot_mfd == True :   
                hyp_name = b
                path = str(Run_name) + '/analysis/figures/analyze_branches/b_value/' + b
                do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)

                
    '''##########################################
    #plot mfd for scalling law  of the logic tree
    ############################################'''
    if len(ScL_complet_list) > 1:
        for ScL in ScL_complet_list :        
            if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/ScL/' + ScL):
                os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/ScL/' + ScL)
            rows = np.where(np.array(total_list_complet_ScL) == ScL)[0]  
            mfd_X = []
            for index in rows :  
                mfd = mega_mfd_cummulative[index]
                mfd_X.append(mfd)
                
            #density plot
            if plot_mfd == True :   
                hyp_name = ScL
                path = str(Run_name) + '/analysis/figures/analyze_branches/ScL/' + ScL
                do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)

    
#            
    '''######################################
    #plot Mmax for each ScL of the logic tree
    ######################################'''
        
    for ScL in ScL_complet_list :
        rows = np.where(np.array(total_list_complet_ScL) == ScL)[0]
        #mfd_ScL_cumulative = []
        Mmax_m_ScL = []
        for index in rows :
            mfd = mega_mfd_cummulative[index]
            #mfd_ScL_cumulative.append(mfd)
            Mmax_m_ScL.append(m_Mmax[index])
            
        
        if not os.path.exists(str(Run_name) + '/analysis/figures/Mmax/for_each_ScL'):
            os.makedirs(str(Run_name) + '/analysis/figures/Mmax/for_each_ScL')
        
        if plot_Mmax == True :
            plt.hist(Mmax_m_ScL,int(round(max(m_Mmax) - min(m_Mmax),1) * 10. + 1.))
            plt.title(ScL)
            plt.savefig(str(Run_name) + '/analysis/figures/Mmax/for_each_ScL/Hist_Mmax_' + ScL +'.png',dpi = 100)
            #plt.show()
            plt.close()

#            
    '''######################################
    #plot Mmax for each scenario set of the logic tree
    ######################################'''
        
    for Sc_set in scenarios_names_list :
        rows = np.where(np.array(total_list_scenario_name) == Sc_set)[0]
        #mfd_Sc_set_cumulative = []
        Mmax_m_Sc_set = []
        for index in rows :
            mfd = mega_mfd_cummulative[index]
            #mfd_Sc_set_cumulative.append(mfd)
            Mmax_m_Sc_set.append(m_Mmax[index])
            
        
        if not os.path.exists(str(Run_name) + '/analysis/figures/Mmax/for_each_scenario_set'):
            os.makedirs(str(Run_name) + '/analysis/figures/Mmax/for_each_scenario_set')
        
        if plot_Mmax == True :
            plt.hist(Mmax_m_Sc_set,int(round(max(m_Mmax) - min(m_Mmax),1) * 10. + 1.))
            plt.title(Sc_set)
            plt.savefig(str(Run_name) + '/analysis/figures/Mmax/for_each_scenario_set/Hist_Mmax_' + Sc_set +'.png',dpi = 100)
            #plt.show()
            plt.close()
#
##            
#    '''######################################
#    # the magnitude of rupture in which each faults are involed, for each set of scenarios
#           work in kinda progress
#    ######################################'''
#
#    for fault in faults_name_list:
#        for Sc_set in scenarios_names_list :
#            rows = np.where(np.array(mega_MFD) == Sc_set)[0]
#            #mfd_Sc_set_cumulative = []
#            Mmax_m_Sc_set = []
#            for index in rows :
#                mfd = mega_mfd_cummulative[index]
#                Mmax_m_Sc_set.append(m_Mmax[index])
#                
#            
#            if not os.path.exists(str(Run_name) + '/analysis/figures/Mmax/for_each_scenario_set'):
#                os.makedirs(str(Run_name) + '/analysis/figures/Mmax/for_each_scenario_set')
#            
#            if plot_Mmax == True :
#                plt.hist(Mmax_m_Sc_set,int(round(max(m_Mmax) - min(m_Mmax),1) * 10. + 1.))
#                plt.title(Sc_set)
#                plt.savefig(str(Run_name) + '/analysis/figures/Mmax/for_each_scenario_set/Hist_Mmax_' + Sc_set +'.png',dpi = 100)
#                #plt.show()
#                plt.close()
#    
    '''######################################
    #########################################
    #    detailled plot for combinaison of
    #       hypothesis
    #########################################
    ######################################'''    

    '''##########################################
    # calculate the difference between the mean rate of the model and the mean rate of the catalog
    ############################################'''        
          
    if plot_mfd == True and plot_mfd_detailled == True:
        file_branch_to_catalog_name = str(Run_name) + '/analysis/txt_files/branch_vs_catalog.txt'
        file_branch_to_catalog = open(file_branch_to_catalog_name,'w')
        index_model = 0
        for model in Model_list : 
            rate_in_catalog = catalog_cum_rate[index_model]
            for MFD_type in MFD_type_list :  
                for scenario in scenarios_names_list :
                    for b_value in b_value_list :  
                        for BG_hyp in BG_hyp_list :  
                            for ScL in ScL_complet_list :  
                                rows_model = np.where(np.array(total_list_model) == model)[0]  
                                rows_mfd = np.where(np.array(total_list_MFD_type) == MFD_type)[0] 
                                rows_sc = np.where(np.array(total_list_scenario_name) == scenario)[0]  
                                rows_ScL = np.where(np.array(total_list_complet_ScL) == ScL)[0]  
                                rows_b = np.where(np.array(total_list_b_value) == b_value)[0]  
                                rows_bg = np.where(np.array(total_list_BG_hyp) == BG_hyp)[0]  
                                rows = list(set(rows_model).intersection(rows_mfd)) 
                                rows = list(set(rows).intersection(rows_sc)) 
                                rows = list(set(rows).intersection(rows_ScL)) 
                                rows = list(set(rows).intersection(rows_b)) 
                                rows = list(set(rows).intersection(rows_bg)) 
                                if len(rows) > 0:
                                    file_branch_to_catalog.write(str(model)+'\t')
                                    file_branch_to_catalog.write(str(MFD_type)+'\t')
                                    file_branch_to_catalog.write(str(scenario)+'\t')
                                    file_branch_to_catalog.write(str(b_value)+'\t')
                                    file_branch_to_catalog.write(str(BG_hyp)+'\t')
                                    file_branch_to_catalog.write(str(ScL)+'\t')
                                    mfd_X = []
                                    for index in rows :  
                                        mfd = mega_mfd_cummulative[index]
                                        mfd_X.append(mfd)
                                    mean_rate_model = np.array(mfd_X).mean(axis=0)
                                    mean_rate_catalog = np.array(rate_in_catalog)#.mean(axis=0)
                                    for i in range(len(mean_rate_catalog)):
                                        file_branch_to_catalog.write(str(mean_rate_model[i]/mean_rate_catalog[i]-1.)+'\t')
                                    file_branch_to_catalog.write('\n')
            index_model +=1
        file_branch_to_catalog.close()
        
    '''##########################################
    # plot the MFD for each branch of the logic tree (can be time consuming)
    ############################################'''        
          
    if plot_mfd == True and plot_mfd_detailled == True:
        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/detailled/'):
            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/detailled/')
        index_model = 0
        for model in Model_list : 
            rate_in_catalog = catalog_cum_rate[index_model]
            mean_rate_catalog = np.array(rate_in_catalog)#.mean(axis=0)
            for MFD_type in MFD_type_list :  
                for scenario in scenarios_names_list :
                    for b_value in b_value_list :  
                        for BG_hyp in BG_hyp_list :  
                            for ScL in ScL_complet_list :  
                                rows_model = np.where(np.array(total_list_model) == model)[0]  
                                rows_mfd = np.where(np.array(total_list_MFD_type) == MFD_type)[0] 
                                rows_sc = np.where(np.array(total_list_scenario_name) == scenario)[0]  
                                rows_ScL = np.where(np.array(total_list_complet_ScL) == ScL)[0]  
                                rows_b = np.where(np.array(total_list_b_value) == b_value)[0]  
                                rows_bg = np.where(np.array(total_list_BG_hyp) == BG_hyp)[0]  
                                rows = list(set(rows_model).intersection(rows_mfd)) 
                                rows = list(set(rows).intersection(rows_sc)) 
                                rows = list(set(rows).intersection(rows_ScL)) 
                                rows = list(set(rows).intersection(rows_b)) 
                                rows = list(set(rows).intersection(rows_bg)) 
                                if len(rows) > 0:
                                    mfd_X = []
                                    for index in rows :  
                                        mfd = mega_mfd_cummulative[index]
                                        mfd_X.append(mfd)
                                    hyp_name =  model + ' ' + MFD_type + ' ' + scenario + ' ' + ScL + ' ' + b_value + ' ' + BG_hyp
                                    path = str(Run_name) +'/analysis/figures/analyze_branches/detailled/' 
                                    do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)
            index_model +=1

    '''##########################################
    #plot mfd for each MFD shape hypothesis and scenario set
    ############################################'''        
          
    if plot_mfd == True and plot_mfd_detailled == True:
        if len(MFD_type_list) > 1 and len(scenarios_names_list)>1:
            index_model = 0
            for model in Model_list : 
                rate_in_catalog = catalog_cum_rate[index_model]
                for MFD_type in MFD_type_list :  
                    for scenario in scenarios_names_list :
                        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model):
                            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model)
                        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' + MFD_type):
                            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' + MFD_type)
                        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' + MFD_type+ '/' +scenario):
                            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' + MFD_type+ '/' +scenario)
                        rows_mfd = np.where(np.array(total_list_MFD_type) == MFD_type)[0] 
                        rows_sc = np.where(np.array(total_list_scenario_name) == scenario)[0]  
                        rows_i = list(set(rows_mfd).intersection(rows_sc)) 
                        rows_model = np.where(np.array(total_list_model) == model)[0]  
                        rows = list(set(rows_i).intersection(rows_model)) 
                        mfd_X = []
                        for index in rows :  
                            mfd = mega_mfd_cummulative[index]
                            mfd_X.append(mfd)
                            
                       #density plot
                        if plot_mfd == True :   
                            hyp_name =  model + ' ' + MFD_type + ' ' + scenario
                            path = str(Run_name) +'/analysis/figures/analyze_branches/Model/' + model+ '/' + MFD_type+ '/' +scenario
                            do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)
                index_model +=1
        
        
        
    '''##########################################
    #plot mfd for each background hypothesis and scenario set
    ############################################'''        
          
    if plot_mfd == True and plot_mfd_detailled == True:
        if len(BG_hyp_list) > 1 and len(scenarios_names_list)>1:
            index_model = 0
            for model in Model_list : 
                rate_in_catalog = catalog_cum_rate[index_model]
                for BG_hyp in BG_hyp_list :  
                    for scenario in scenarios_names_list :
                        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/'  + BG_hyp+ '/' +scenario):
                            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/'  + BG_hyp+ '/' +scenario)
                        rows_mfd = np.where(np.array(total_list_BG_hyp) == BG_hyp)[0] 
                        rows_sc = np.where(np.array(total_list_scenario_name) == scenario)[0]  
                        rows = list(set(rows_mfd).intersection(rows_sc)) 
                        rows_model = np.where(np.array(total_list_model) == model)[0]  
                        rows = list(set(rows).intersection(rows_model)) 
                        mfd_X = []
                        for index in rows :  
                            mfd = mega_mfd_cummulative[index]
                            mfd_X.append(mfd)
                            
                       #density plot
                        if plot_mfd == True :
                            hyp_name = BG_hyp + ' ' + scenario
                            path = str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/'  + BG_hyp+ '/' +scenario
                            #total_list_hyp = total_list_MFD_type
                            do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)
                index_model +=1
        
        
        
    '''##########################################
    #plot mfd for each model hypothesis and MFD
    ############################################'''        
          
          
    if plot_mfd == True and plot_mfd_detailled == True:
        if len(Model_list) > 1 and len(MFD_type_list)>1:
            index_model = 0
            for model in Model_list :  
                rate_in_catalog = catalog_cum_rate[index_model]
                for MFD_type in MFD_type_list :
                    if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' +MFD_type):
                        os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' +MFD_type)
                    rows_i = np.where(np.array(total_list_model) == model)[0] 
                    rows_j = np.where(np.array(total_list_MFD_type) == MFD_type)[0]  
                    rows = list(set(rows_i).intersection(rows_j)) 
                    mfd_X = []
                    for index in rows :  
                        mfd = mega_mfd_cummulative[index]
                        mfd_X.append(mfd)
                        
                   #density plot
                    if plot_mfd == True :
                        hyp_name = model + ' ' + MFD_type
                        path = str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/' +MFD_type
                        #total_list_hyp = total_list_MFD_type
                        do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)
                index_model +=1        
            
                
                
    '''##########################################
    #plot mfd for each background hypothesis and mfd
    ############################################'''        
          
    if plot_mfd == True and plot_mfd_detailled == True:
        if len(BG_hyp_list) > 1 and len(MFD_type_list)>1:
            index_model = 0
            for model in Model_list : 
                rate_in_catalog = catalog_cum_rate[index_model]
                for BG_hyp in BG_hyp_list :  
                    for MFD_type in MFD_type_list :
                        if not os.path.exists(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/'  + BG_hyp+ '/' +MFD_type):
                            os.makedirs(str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/'  + BG_hyp+ '/' +MFD_type)
                        rows_i = np.where(np.array(total_list_BG_hyp) == BG_hyp)[0] 
                        rows_j = np.where(np.array(total_list_MFD_type) == MFD_type)[0]  
                        rows = list(set(rows_i).intersection(rows_j)) 
                        rows_model = np.where(np.array(total_list_model) == model)[0]  
                        rows = list(set(rows).intersection(rows_model)) 
                        mfd_X = []
                        for index in rows :  
                            mfd = mega_mfd_cummulative[index]
                            mfd_X.append(mfd)
                            
                       #density plot
                        if plot_mfd == True :
                            hyp_name = BG_hyp + ' ' + MFD_type
                            path = str(Run_name) + '/analysis/figures/analyze_branches/Model/' + model+ '/'+ BG_hyp+ '/' +MFD_type
                            #total_list_hyp = total_list_MFD_type
                            do_the_plots(hyp_name,mfd_X,mega_bining_in_mag,xmin,xmax,ymin,ymax,Run_name,rate_in_catalog,plot_as_rep,a_s_model,rows,path,bining_in_mag)
                index_model +=1        
            
            
    return (total_list_ScL,total_list_dimension_used,geologic_moment_rate,
            geologic_moment_rate_no_as,total_list_scenario_name,total_list_MFD_type,
            mega_mfd_cummulative,total_list_model,total_list_sample,total_list_BG_hyp)
