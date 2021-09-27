#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: thomas
"""


import numpy as np
import os

from scipy.stats import chisquare
from scipy.stats import multivariate_normal

import matplotlib.pyplot as plt

def sampling_analysis(Run_name,Model_list,m_Mmax,b_sample,a_s_model,mega_mfd_cummulative,catalog_cum_rate,
                      xmin,xmax,ymin,ymax,total_list_model,bining_in_mag,total_list_MFD_type,
                      total_list_scenario_name,file_faults_data,total_list_sample,total_list_BG_hyp):
    
    if not os.path.exists(str(Run_name) + '/analysis/figures/sampling_analysis'):
        os.makedirs(str(Run_name) + '/analysis/figures/sampling_analysis')
        
        
    file_LT_metrics=open(str(Run_name) + '/analysis/txt_files/LT_metrics.txt','w')
    file_LT_metrics.write('ScL\tModel\tBG\tbvalue\tMFD\tSc\tsample\tmean_sr\tchi_score\tMmax_score\tNMS_score\tpaleo_score\n')
    

    
    

    '''
    #read the catalog cumulative rates for all sampling 
        and compare each branch to the catalog using the chi-squared test
    
    methodology:
        for each branch i of the logic tree
        30 random MFD calculated from the catalog are extracted for the comaprison
        the modeled rate of the branch i are compared the each on of these random samples
        the comparison is done using the following formula:
        we calculate the absolute value of the difference between the logs of the model minus the log of catalog rate
        we had 10 to this asolute value to make it accepatble or the chi-squarred test
        we run the chisquared test for an array of ten value corresponding to one unit of magnitude (ten bins of 0.1)
        bins of magnitude where one of the two rates are not defined are deleted
        the value are conpared to an array of value equal ten (the expected value is the model fits the data)
        we save the pvalue calculated 
        
        In order to get the p value for the whole MFD, we do the mean of the pvalue for each unit of magnitude 
        weighted by the number of filled bins in the range of magnitude.
        
        If the p value is close to 1, the two MFD are similar. 
        
        personal opinion: 
        p values superior to 0.9 seam like a good match
        p values superior to 0.8 seam like an acceptable match in most cases
        p values less than 0.7 make the match difficult to accept
        
        Warning! this method doesn't care if the two maximum magnitude are different,
        it will only take the bin where both MFDs are defined. 
        The fit in terms of Mmax need to rely on some other test. (hopefully I managed to provide one...)
        
    '''
    
    plot_fig=False
    
    
    index_Mmin=np.where(np.array(np.linspace(4.0,7.0,num=31).round(1))==xmin)[0][0]
    index_Mmax=np.where(np.array(np.linspace(4.0,10.0,num=61).round(1))==xmax)[0][0]+1
    
    file = open(str(Run_name) + '/analysis/txt_files/model_performance.txt','w')
    file.write('Model\tMFD type\tBG\tScenario Set\tsample\tFit to catalog\tFit to Paleo\tNMS score\n')
        
    index_model=0
    for model in Model_list:
    
        if not os.path.exists(str(Run_name) + '/analysis/figures/sampling_analysis/'+model):
            os.makedirs(str(Run_name) + '/analysis/figures/sampling_analysis/'+model)
        
        
        catfile_all=str(Run_name) + '/analysis/figures/catalogue/catalog_rates_all_'+model+'.txt'
        with open(catfile_all) as f:#finds where to start reading
            lines_cat = f.readlines()
        
        #ranges of magnitude where the test is made
        ranges_mag=[[4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9],
                [5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9],
                [6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9],
                [7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9],
                [8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9],
                [9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9]
                ]
        
        p_chi_branch=[]
        indexes_model=[]
        
        index_branch = 0
        for mfd,model_name_i in zip(mega_mfd_cummulative,total_list_model):
            if model_name_i ==model:
                indexes_model.append(index_branch)
                indexes_catalogs_to_test = np.random.choice(range(len(lines_cat))[1:],size=40)
                #indexes_catalogs_to_test = range(len(lines_cat)) #take them all, it doesn't take that long
                pvalues=[[],
                    [],
                    [],
                    [],
                    [],
                    []
                    ]
                weights_pvalues=[[],
                    [],
                    [],
                    [],
                    [],
                    []
                    ]
                
                if plot_fig==True:
                    f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
                    
                for i_cat in indexes_catalogs_to_test:
                    cat_rates_i=lines_cat[i_cat].split('\t')
                    cat_rates_i=[float(i) for i in cat_rates_i]
                    
                    if plot_fig==True:
                        ax1.scatter(bining_in_mag,cat_rates_i,c='k',alpha=0.1,s=0.5)
                        
                    index_range=0
                    for range_i in ranges_mag:
                        diff_rate=[]
                        target_value=[]
                        bining_i=[]
                        for model_rate_i,data_rate_i,mag_i in zip(mfd[index_Mmin:index_Mmax],cat_rates_i,bining_in_mag):
                            if model_rate_i!=0 and data_rate_i!=0 and mag_i in range_i:
        #                            model_rate.append(-np.log10(model_rate_i)*10.)
        #                            data_rate.append(-np.log10(data_rate_i)*10.)
                                diff_rate.append(abs(np.log10(model_rate_i)-np.log10(data_rate_i))*10. + 10.)
                                target_value.append(10.)
                                bining_i.append(mag_i)
                        if len(diff_rate)>=2:
                            pvalues[index_range].append(chisquare(diff_rate,f_exp=target_value)[1]) #pvalue for each range and sample
                            weights_pvalues[index_range].append(len(diff_rate)) #associated weight depending of number of bin in the range that are filled
                            
                            if plot_fig==True:
                                ax2.scatter(bining_i,diff_rate,c='r',alpha=0.2,s=2)
                                ax2.scatter(bining_i,target_value,c='k',alpha=0.1,s=2)
                        index_range+=1
                        
                if plot_fig==True:
                    ax1.scatter(bining_in_mag,mfd[index_Mmin:index_Mmax],c='r',s=0.5)
                    #ax1.set_title(str(round(np.mean(pvalues),3)))
                    ax1.set_yscale('log')
                    ax1.set_xlim([xmin,xmax])
                    ax1.set_ylim([ymin,ymax])
                p_total=[]
                weight_p=[]
                for range_i,p_i,w_i in zip(ranges_mag,pvalues,weights_pvalues):
                    if len(p_i)!=0 :
                        weight_p.append(np.mean(w_i))
                        p_total.append(round(np.mean(p_i),4))
                        
                        if plot_fig==True:
                            if round(np.mean(p_i),3) >= 0.9:
                                color='g'
                            elif round(np.mean(p_i),3) >= 0.8:
                                color='orange'
                            else:
                                color='r'
                            ax2.text(np.mean(range_i),25,str(round(np.mean(p_i),3)),fontsize=8,color=color)
                            
                p_chi_branch.append(round(np.average(p_total,weights=weight_p),3))
                if plot_fig==True:
                    if round(np.average(p_total,weights=weight_p),3) >= 0.9:
                        color='g'
                    elif round(np.average(p_total,weights=weight_p),3) >= 0.8:
                        color='orange'
                    else:
                        color='r'
                    ax1.set_title(str(round(np.average(p_total,weights=weight_p),3)),color=color)
                    ax2.set_xlim([xmin-0.1,xmax])
                    ax2.set_ylim([9,30])
                    plt.show()
                    plt.close()
            
            
            index_branch+=1
        
        
        
        '''
        #   Mmax fit to the Mmax in the catalog
        
        The rule is: The Mmax in the model should be at least the one in the catalog 
        but the catalog has some uncertainties on the magnitude of large historical EQs
        
        methodology:
            we calculate the cumulative density distribution of the Mmax in the catalog
            we associate the given density to each Mmax of the models
        
        '''
        Mmax_cat=[]
        bining_cat=lines_cat[0].split('\t')
        
        bining_cat=[float(i) for i in bining_cat]
        
        for i_cat in range(len(lines_cat)-1):
            cat_rates_i=lines_cat[i_cat+1].split('\t')
            cat_rates_i=[float(i) for i in cat_rates_i]
            
            i_test=0
            try :
                while cat_rates_i[i_test]!=0:
                    i_test+=1
            except:
                i_test=len(cat_rates_i)-1
            Mmax_cat.append(bining_cat[i_test])
            
        distribution_Mmax_cat=[]
        for mag_i in bining_cat:
            d_i = sum(i <= mag_i+0.1 for i in Mmax_cat)/len(Mmax_cat)
            distribution_Mmax_cat.append(d_i)
            
        plt.plot(bining_cat,distribution_Mmax_cat)
        plt.xlim([xmax-1.5,xmax])
        
        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/Mmax_distrib_in_the_cat.png',dpi = 180)
        plt.close()
        
        
        weight_model_Mmax=[]
        for Mmax_i,model_name_i in zip(m_Mmax,total_list_model):
            if model_name_i == model:
                index = np.where(np.array(bining_cat)==Mmax_i)[0][0]
                weight_model_Mmax.append(distribution_Mmax_cat[index])
            
        
        '''
        The NMS on a set of faults as a metric for judging the quality of a model
        '''
        fault_set=['F1','F2','F3']
        
        NMS_set=[]
        for fault in fault_set:
            NMS_set.append([])
        
        if len(NMS_set) != 0:

          sr_sample=[]
          for fault in fault_set:
              sr_sample.append([])

          score_nms=[]


          ####
          # extract the slip-rate of each target fault and does the mean for that branch. 
          # this can allow to see is some slip-rate values seam to work better
          ####
          srate_sample_file=str(Run_name) + '/analysis/txt_files/slip_rate_sampling.txt'
          with open(srate_sample_file) as f:#finds where to start reading
              lines_sr = f.readlines()



          srep_file=str(Run_name) + '/analysis/txt_files/slip_rep_on_faults_all_data.txt'
          try:
              with open(srep_file) as f:#finds where to start reading
                  lines = f.readlines()
              line_number=0    
              for line in lines:
                  #print(line)
                  if line.split('\t')[7] in fault_set and line.split('\t')[1]==model:
                      index_fault=np.where(np.array(fault_set)==line.split('\t')[7])[0][0]
                      NMS_set[index_fault].append(float(line.split('\t')[-1]))
                      sr_sample[index_fault].append(float(lines_sr[line_number].split('\t')[-1]))

                  line_number+=1
              if  np.sum(NMS_set) != 0. :       
                  #print('score NMS on target faults',np.mean(NMS_set,axis=0))
                  for i in range(len(p_chi_branch)):

                      '''
                      the score is 1 is MSN is less than 20%
                      the score is 0 if:
                          at least one of the NMS of the test faults if more than 50%
                          the mean is more the 40%
                      between 20 and 40 the score evolves linearily between 1 and 0
                      (this is very much open to discussion!)
                      '''
                      if np.mean(NMS_set,axis=0)[i] > 40.:
                          score_nms_i = 0.
                      elif np.mean(NMS_set,axis=0)[i] < 20.:
                          score_nms_i = 1.
                      else :
                          score_nms_i=2 - 1./20.*np.mean(NMS_set,axis=0)[i]
                      #print('score NMS on target faults',round(score_nms_i,2))
                      '''hard limit on acceptability'''    
                      for nms_row in NMS_set:
                          #print(nms_row[i])
                          if nms_row[i] > 50.:
                              score_nms_i = 0.

                      score_nms.append(score_nms_i)    
                      #print('score NMS on target faults',round(score_nms_i,2),'   NMS mean:',round(np.mean(NMS_set,axis=0)[i]))
          except FileNotFoundError:
              print('!!! you need to run the plot_sr_use if you want the NMS metric !!!')
              print('Default value = 1. ')
              for i in range(len(p_chi_branch)):
                  score_nms.append(1.)
        else:
          print('modify Sampling_analysis.py for the NMS metric')
          print('Default value = 1. ')
          for i in range(len(p_chi_branch)):
              score_nms.append(1.)
        
        #deos the mean sr of the faults for each branch
        mean_sr_branch = np.mean(sr_sample,axis=0)
#        plt.scatter(mean_sr_branch,p_chi_branch)
#        plt.show()
        
                
        '''#############################
        Weight based on the fit to the paleo rates
        
        and the RSQSim rates if they exist
        #######################################'''
        plot_paleo = False
        plot_rsqsim_pr =True
            
        #extract the faults data
        faults_data = np.genfromtxt(file_faults_data, dtype=[('model', 'U100000'),
                                                         ('fault_name', 'U100000'),
                                                         ('type', 'U100000'),
                                                         ('M', 'f8'),
                                                         ('sig_M', 'f8'),
                                                         ('rate', 'f8'),
                                                         ('sig_rate', 'f8')],
                                                 delimiter = '\t',skip_header = 1)
        # Dealing with one line files
        try:
            len_faults_data = len(faults_data)
        except TypeError:
            faults_data = faults_data.reshape((1,))
        
        
        rsqsim_pr=False
        RSQSim_pr_file = str(Run_name) + '/file_pr_rsqsim.txt'
        try:            
            with open(RSQSim_pr_file) as f:#finds where to start reading
                lines = f.readlines()
            bin_mag_rsqsim = [round(float(i),1) for i in lines[0].split('\t')[1:-1]]
            
            rqsim_pr_faults=[]
            faults_name_rsqsim = []
            for line in lines[1:]:
                faults_name_rsqsim.append(line.split('\t')[0])
                rqsim_pr_faults.append([float(i) for i in line.split('\t')[1:-1]]) #we don't take the last point of the MFD , too specific
            
            
            index_Mmin_rsqsim=np.where(np.array(bining_in_mag)==bin_mag_rsqsim[0])[0][0]
            index_Mmax_rsqsim=np.where(np.array(bining_in_mag)==bin_mag_rsqsim[-1])[0][0]+1
        except:
            rsqsim_pr=False
            
        
        #print faults_data
        data_model = list(map(lambda i : faults_data[i][0], range(len(faults_data))))
        data_fault_name =list( map(lambda i : faults_data[i][1], range(len(faults_data))))
        data_type =list( map(lambda i : faults_data[i][2], range(len(faults_data))))
        data_M =list( map(lambda i : float(faults_data[i][3]), range(len(faults_data))))
        data_sig_M =list( map(lambda i : float(faults_data[i][4]), range(len(faults_data))))
        data_rate = list(map(lambda i : float(faults_data[i][5]), range(len(faults_data))))
        data_sig_rate =list( map(lambda i : float(faults_data[i][6]), range(len(faults_data))))
        
        score_paleo = []
#        
#        score_paleo_per_fault = []
#        for fault in data_fault_name:
#            score_paleo_per_fault.append([])
        
        score_paleo_faults=[]
        faults_data=[]
        
        score_pr_rsqsim = []
        faults_rsqsim = []
        for fault,data_model_i in zip(data_fault_name,data_model):
            if data_model_i == model and fault not in faults_data:
                score_paleo_faults.append([])
                faults_data.append(fault)
                if rsqsim_pr == True:
                    if fault in faults_name_rsqsim and fault not in faults_rsqsim:
                        score_pr_rsqsim.append([])
                        faults_rsqsim.append(fault)
            
        participation_rate_file=str(Run_name) + '/analysis/figures/rupture_rate_for_each_fault_cum/' + model + '/file_for_comparison.txt'
        with open(participation_rate_file) as f:#finds where to start reading
            lines_pr = f.readlines()
            
        paleo_list_mfd = []
        paleo_list_bvalue = []
        paleo_list_bg = []
        paleo_list_scl = []
        paleo_list_scenario = []
        paleo_list_sample = []
        index_branch=0
        for line in lines_pr :
            index_fault=0
            for fault_name in faults_data:
                if line.split('\t')[0]==model and line.split('\t')[7]==fault_name:
                    if index_fault==0:
                        paleo_list_mfd.append(line.split('\t')[1])
                        paleo_list_scenario.append(line.split('\t')[2])
                        paleo_list_bg.append(line.split('\t')[3])
                        paleo_list_scl.append(line.split('\t')[4])
                        paleo_list_bvalue.append(line.split('\t')[5])
                        paleo_list_sample.append(line.split('\t')[6])
                    
                    mfd_i = [float(i) for i in list(line.split('\t')[(8+index_Mmin):(8+index_Mmax)])]
                    
                    #######
                    # COMPARE WITH THE PALEO
                    #######
                    self_data_M = []
                    self_data_sig_M = []
                    self_data_rate = []
                    self_data_sig_rate = []
        
                    index_fault_in_data = np.where(np.array(data_fault_name)==fault_name)[0]
                    for index_i in index_fault_in_data:
                        if data_model[index_i] == model and data_type[index_i] == 'pal':
                            self_data_M.append(data_M[index_i]) 
                            self_data_sig_M.append(data_sig_M[index_i])
                            self_data_rate.append(data_rate[index_i])
                            self_data_sig_rate.append(data_sig_rate[index_i]) 
                    
                    #calculating the paleoscore using a lognomral distribution for the paleouncertainties
                    paleo_score_i=[]
                    for m_i,sm_i,r_i,sr_i in zip(self_data_M,self_data_sig_M,self_data_rate,self_data_sig_rate):
                         
                        x, y = np.mgrid[4.5:7.5:.01, -5.:0.:.01]
                        pos = np.empty(x.shape + (2,))
                        pos[:, :, 0] = x; pos[:, :, 1] = y
                        #2D noraml * log normal law
                        rv = multivariate_normal([m_i, np.log10(r_i)], [sm_i+0.001, sr_i+0.0000001])
                        
                        #interpolates the MFD
                        detailed_bin_mag = np.linspace(bining_in_mag[0],bining_in_mag[-1],1000)
                        detailed_mfd_i = np.interp(detailed_bin_mag,bining_in_mag,np.log10(mfd_i))
                        
                        
                        if plot_paleo == True:
                            plt.contourf(x, y, rv.pdf(pos),alpha=0.5)
                            
                            plt.scatter(bining_in_mag,np.log10(mfd_i),c='k',marker='s',s=10,linewidths=0.01,alpha=0.7)
                            plt.scatter(detailed_bin_mag,detailed_mfd_i,c='k',marker='s',s=3,linewidths=0.01,alpha=0.7)
                            plt.xlim([5.,7.])
                            plt.ylim([-3,-1.])
                            plt.grid()
                            plt.show()
                        
                        paleo_score_i.append(max([rv.pdf([i,j])/rv.pdf([m_i,np.log10(r_i)]) for i,j in zip(detailed_bin_mag,detailed_mfd_i)]))
                        #print(max([rv.pdf([i,j])/rv.pdf([m_i,np.log10(r_i)]) for i,j in zip(detailed_bin_mag,detailed_mfd_i)]))
                    
                        
                    score_paleo_faults[index_fault].append(np.mean(paleo_score_i))
                    
                    #################
                    # Compare with RSQSim (if it exists)
                    # make the mean of the ration where both are defined expect for the last 2 bins 
                    # (the big drop in rate leads to very large ratios but actually it's small rates so it doesn't martter so much)
                    ################
                    
                    if rsqsim_pr == True and line.split('\t')[6] == '1':
                        pvalues = []
                        pshape = []
                        if fault_name in faults_rsqsim:
                            index_fault_rsqsim = np.where(np.array(faults_name_rsqsim)==fault_name)[0][0]
                            fault_pr_rsqsim = rqsim_pr_faults[index_fault_rsqsim]
                            
                            if plot_rsqsim_pr==True:
                                f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
                        
                            #for i_cat in indexes_catalogs_to_test:
                            #cat_rates_i=lines_cat[i_cat].split('\t')
                            #cat_rates_i=[float(i) for i in cat_rates_i]
                            
                            if plot_rsqsim_pr==True:
                                ax1.scatter(bin_mag_rsqsim[:-2],fault_pr_rsqsim[:-2],c='k',alpha=0.9,s=3)
                                ax1.scatter(bin_mag_rsqsim[-2:],fault_pr_rsqsim[-2:],c='k',alpha=0.5,s=3)
#                                
#                            index_range=0
#                            for range_i in ranges_mag:
                            diff_rate=[]
#                            target_value=[]
                            bining_i=[]
                            for model_rate_i,data_rate_i,mag_i in zip(mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2],fault_pr_rsqsim[:-2],bin_mag_rsqsim[:-2]):
                                if model_rate_i!=0 and data_rate_i!=0:
                                    if model_rate_i >= data_rate_i:
                                        diff_rate.append(model_rate_i/data_rate_i)
                                    else : 
                                        diff_rate.append(data_rate_i/model_rate_i)
                                    bining_i.append(mag_i)
                                    
                            pvalues.append(np.mean(diff_rate)) #pvalue for each range and sample
                            
                                
                            if plot_rsqsim_pr==True:
                                ax2.scatter(bining_i,diff_rate,c='b',alpha=0.8,s=2)
                                index_range+=1
                                    
                            if plot_rsqsim_pr==True:
                                ax1.scatter(bining_in_mag[index_Mmin_rsqsim:index_Mmax_rsqsim-2],mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2],c='r',s=0.5)
                                ax1.scatter(bining_in_mag[-2:],mfd_i[-2:],c='r',alpha=0.4,s=0.5)
                                ax1.set_yscale('log')
                                ax1.set_xlim([xmin+1.,xmax])
                                ax1.set_ylim([ymin,ymax/100.])
                            p_total=np.mean(diff_rate)
                            
                            #test on the shape (normalized mfd)
                            n_mfdi = [i/sum(mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2]) for i in mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2]]
                            n_mfd_rsqsim = [i/sum(fault_pr_rsqsim[:-2]) for i in fault_pr_rsqsim[:-2]]
                            diff_rate=[]
                            bining_i=[]
                            for model_rate_i,data_rate_i,mag_i in zip(n_mfdi,n_mfd_rsqsim,bin_mag_rsqsim[:-2]):
                                if model_rate_i!=0 and data_rate_i!=0:
                                    if model_rate_i >= data_rate_i:
                                        diff_rate.append(model_rate_i/data_rate_i)
                                    else : 
                                        diff_rate.append(data_rate_i/model_rate_i)
                                    bining_i.append(mag_i)
                            pshape.append(np.mean(diff_rate)) #pvalue for each range and sample
                            
                            if plot_rsqsim_pr==True:
                                ax2.scatter(bining_i,diff_rate,c='g',alpha=0.8,s=2)
                            
                            
                            if plot_rsqsim_pr==True:
                                if round(p_total,3) >= 1.3:
                                    color='r'
                                elif round(p_total,3) >= 1.2:
                                    color='orange'
                                else :
                                    color='g'
                                    
                                if round(np.mean(diff_rate),3) >= 1.3:
                                    color_shape='r'
                                elif round(np.mean(diff_rate),3) >= 1.2:
                                    color_shape='orange'
                                else :
                                    color_shape='g'
                                    
                                ax1.set_title(model +' '+ fault_name + ' '+str(round(p_total,2)),color=color)
                                ax2.set_title(str(round(np.mean(diff_rate),2)),color=color_shape)
                                ax1.set_xlim([xmin+1.,xmax])
                                ax2.set_xlim([xmin+1.,xmax])
                                ax2.set_ylim([0.9,3.])
                                plt.show()
                                plt.close()
                            
                
                index_fault +=1
        score_paleo = np.mean(score_paleo_faults,axis=0)
        
            
            
        '''###################""
        Compare with some other MFD at the system level (physics based for example)
        #####################"'''
        
        plot_fig_rsqsim=False
        RSQSim_MFD = str(Run_name) + '/mfd_RSQSim.txt'
        try:            
            with open(RSQSim_MFD) as f:#finds where to start reading
                lines = f.readlines()
            bin_mag_rsqsim = [round(float(i),1) for i in lines[0].split('\t')[1:-1]]
            mfd_rsqsim = [float(i) for i in lines[1].split('\t')[1:-1]]
            
            index_Mmin_rsqsim=np.where(np.array(bining_in_mag)==bin_mag_rsqsim[0])[0][0]
            index_Mmax_rsqsim=np.where(np.array(bining_in_mag)==bin_mag_rsqsim[-1])[0][0]+1
            
            index_branch = 0
            for mfd,model_name_i in zip(mega_mfd_cummulative,total_list_model):
                if model_name_i ==model:
                    
                    if plot_fig_rsqsim==True:
                        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
                        
                    if plot_fig_rsqsim==True:
                        ax1.scatter(bin_mag_rsqsim,mfd_rsqsim,c='k',alpha=0.9,s=3)
                        
                    pvalues = []
                    diff_rate=[]
                    bining_i=[]
                    mfd_i = mfd[index_Mmin:index_Mmax]
                    for model_rate_i,data_rate_i,mag_i in zip(mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2],mfd_rsqsim[:-2],bin_mag_rsqsim[:-2]):
                        if model_rate_i!=0 and data_rate_i!=0:
                            if model_rate_i >= data_rate_i:
                                diff_rate.append(model_rate_i/data_rate_i)
                            else : 
                                diff_rate.append(data_rate_i/model_rate_i)
                            bining_i.append(mag_i)
                    pvalues.append(np.mean(diff_rate))
                    p_total=np.mean(diff_rate)
                    
                    if plot_fig_rsqsim==True:
                        ax2.scatter(bining_i,diff_rate,c='r',alpha=0.9,s=3)
                        
                            
                    if plot_fig_rsqsim==True:
                        ax1.scatter(bining_in_mag[index_Mmin_rsqsim:index_Mmax_rsqsim],mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim],c='r',s=0.5)
                        
                                
                    #test on the shape (normalized mfd)
                    n_mfdi = [i/sum(mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2]) for i in mfd_i[index_Mmin_rsqsim:index_Mmax_rsqsim-2]]
                    n_mfd_rsqsim = [i/sum(mfd_rsqsim[:-2]) for i in mfd_rsqsim[:-2]]
                    diff_rate=[]
                    bining_i=[]
                    for model_rate_i,data_rate_i,mag_i in zip(n_mfdi,n_mfd_rsqsim,bin_mag_rsqsim[:-2]):
                        if model_rate_i!=0 and data_rate_i!=0:
                            if model_rate_i >= data_rate_i:
                                diff_rate.append(model_rate_i/data_rate_i)
                            else : 
                                diff_rate.append(data_rate_i/model_rate_i)
                            bining_i.append(mag_i)
                    pshape.append(np.mean(diff_rate)) #pvalue for each range and sample
                    
                    if plot_fig_rsqsim==True:
                        ax2.scatter(bining_i,diff_rate,c='g',alpha=0.8,s=3)
                                
                    if plot_fig_rsqsim==True:
                        if round(p_total,3) >= 1.3:
                            color='r'
                        elif round(p_total,3) >= 1.2:
                            color='orange'
                        else :
                            color='g'
                            
                        if round(np.mean(diff_rate),3) >= 1.4:
                            color_shape='r'
                        elif round(np.mean(diff_rate),3) >= 1.3:
                            color_shape='orange'
                        else :
                            color_shape='g'
                        ax1.set_title(model +' '+str(round(p_total,2)),color=color)
                        ax2.set_title(str(round(np.mean(diff_rate),2)),color=color_shape)
                        ax1.set_ylim([ymin/10.,ymax])
                        ax1.set_xlim([xmin+1.,xmax])
                        ax2.set_xlim([xmin+1.,xmax])
                        ax2.set_ylim([0.9,3.])
                        ax1.set_yscale('log')
                        plt.show()
                        plt.close()
                
                
                index_branch+=1
        
        except:
            pass
            #print('no rsqsim file')
        
        
        
        
                
        '''
        Setting the weight for each score
        '''
        #weight the different parameters
        #the sum must be one
        weight_chi=0.35
        weight_Mmax=0.05
        weight_NMS_faults_test=0.3
        weight_paleo = 0.3
        
        if len(score_nms)==0.:
            print('!!! no selected faults for the NMS metric !!!')
            print('Default value = 0.  Weight is set to 0.')
            weight_NMS_faults_test=0.
            weight_chi = weight_chi / (weight_chi+weight_Mmax+weight_NMS_faults_test+weight_paleo)
            weight_Mmax = weight_Mmax / (weight_chi+weight_Mmax+weight_NMS_faults_test+weight_paleo)
            weight_paleo = weight_paleo / (weight_chi+weight_Mmax+weight_NMS_faults_test+weight_paleo)
            for i in range(len(p_chi_branch)):
                score_nms.append(0.)
                
        if len(score_paleo)==0.:
            print('!!! no paleo data on the faults !!!')
            print('Default value = 0.  Weight is set to 0.')
            weight_paleo=0.
            weight_chi = weight_chi / (weight_chi+weight_Mmax+weight_NMS_faults_test+weight_paleo)
            weight_Mmax = weight_Mmax / (weight_chi+weight_Mmax+weight_NMS_faults_test+weight_paleo)
            weight_NMS_faults_test = weight_NMS_faults_test / (weight_chi+weight_Mmax+weight_NMS_faults_test+weight_paleo)
            for i in range(len(p_chi_branch)):
                score_paleo.append(0.)
        
        
        ''' 
        Builbing the text file
        '''
        lt_branch = []
        lt_i_before = 'truc'
        srep_file=str(Run_name) + '/analysis/txt_files/slip_rep_on_faults_all_data.txt'
        try:
            with open(srep_file) as f:#finds where to start reading
                lines = f.readlines()
            ordered_score_paleo = []
            i_lt=0    
            for line in lines:
                if line.split('\t')[1]==model:
                    lt_i=[]
                    for i in range(7): #add the branches parameters
                        lt_i.append(line.split('\t')[i])
                        
                    if str(lt_i) != lt_i_before:
                        lt_i_before = str(lt_i)
                        
                        lt_i.append(round(mean_sr_branch[i_lt],3))
                        lt_i.append(round(p_chi_branch[i_lt],3))
                        lt_i.append(round(weight_model_Mmax[i_lt],3))
                        lt_i.append(round(score_nms[i_lt],3))
                        
                        #oredering the score paleo
                        i1 = np.where(np.array(paleo_list_mfd)==line.split('\t')[4][4:])[0]
                        i2= np.where(np.array(paleo_list_scenario)==line.split('\t')[5])[0]
                        i3 = np.where(np.array(paleo_list_sample)==line.split('\t')[6].split('_')[1])[0]
                        i4= np.where(np.array(paleo_list_bvalue)==line.split('\t')[3])[0]
                        i5= np.where(np.array(paleo_list_bg)==line.split('\t')[2][3:])[0]
                        i6= np.where(np.array(paleo_list_scl)==line.split('\t')[0])[0]
                        i1 = np.intersect1d(i1,i2)
                        i1 = np.intersect1d(i1,i3)
                        i1 = np.intersect1d(i1,i4)
                        i1 = np.intersect1d(i1,i5)
                        i1 = np.intersect1d(i1,i6)
#                        print(line.split('\t')[3],line.split('\t')[2][3:],line.split('\t')[0])
#                        print(i1)
#                        index_score_paleo = np.where(np.logical_and(
#                                np.array(paleo_list_mfd)==line.split('\t')[4].split('_')[1],
#                                np.array(paleo_list_scenario)==line.split('\t')[5],
#                                np.array(paleo_list_sample)==line.split('\t')[6].split('_')[1]                                
#                                ))[0]
#                        print(index_score_paleo)
#                        print(paleo_list_mfd[0],line.split('\t')[4].split('_')[1])
#                        print(paleo_list_scenario[0],line.split('\t')[5])
#                        print(paleo_list_sample[0],line.split('\t')[6].split('_')[1])
#                        print(i1)
#                        print(len(paleo_list_scl),len(paleo_list_bg),len(paleo_list_bvalue),len(paleo_list_sample),len(paleo_list_scenario),len(paleo_list_mfd))
#                        print(score_paleo)
                        lt_i.append(round(np.mean(np.take(score_paleo,i1)),3))
                        ordered_score_paleo.append(round(np.mean(np.take(score_paleo,i1)),3))
                        
                        lt_branch.append(lt_i)
                        i_lt+=1
                        
            for lt_i in lt_branch:
                line=''
                for i in lt_i :
                    line+=str(i)+'\t'
                line=line[:-1]
                file_LT_metrics.write(line+'\n')
            
            
        except (FileNotFoundError, IndexError) as  e:
            print('!!! you need to run the plot_sr_use if you want the file with the metrics  and modify Sampling_analysis.py!!!')
    
        '''
        Calculataing the weighted score for each branch
        '''

        if ordered_score_paleo == []:
            ordered_score_paleo = [0 for i in range(len(p_chi_branch))]
        
        final_weigth = []
        for i in range(len(p_chi_branch)):
            final_weigth.append(p_chi_branch[i] * weight_chi+
                          weight_model_Mmax[i] * weight_Mmax+
                          score_nms[i] * weight_NMS_faults_test+
                          ordered_score_paleo[i] * weight_paleo
                          )
            
        '''
        Plotting section
        
        Weighted average of the different metric
        
        user defined weight for each metric. the figure give the weighted average as a final p value
        '''
        
        color_mfd=[]
        for MFD_type_i,model_name_i in zip(total_list_MFD_type,total_list_model):
            if model_name_i == model:
                if MFD_type_i =='GR':
                    color_mfd.append('darkblue')
                else:
                    color_mfd.append('darkgreen')
#          
#        
    
    
        
        f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey=True)
        
        ax1.axhspan(0.8, 1.1, facecolor='g', alpha=0.1)
        ax1.axhspan(0.6, 0.8, facecolor='orange', alpha=0.1)
        ax1.axhspan(-0.1, 0.6, facecolor='r', alpha=0.1)
        ax2.axhspan(0.8, 1.1, facecolor='g', alpha=0.1)
        ax2.axhspan(0.6, 0.8, facecolor='orange', alpha=0.1)
        ax2.axhspan(-0.1, 0.6, facecolor='r', alpha=0.1)
        ax3.axhspan(0.8, 1.1, facecolor='g', alpha=0.1)
        ax3.axhspan(0.6, 0.8, facecolor='orange', alpha=0.1)
        ax3.axhspan(-0.1, 0.6, facecolor='r', alpha=0.1)
        ax4.axhspan(0.8, 1.1, facecolor='g', alpha=0.1)
        ax4.axhspan(0.6, 0.8, facecolor='orange', alpha=0.1)
        ax4.axhspan(-0.1, 0.6, facecolor='r', alpha=0.1)
        ax5.axhspan(0.8, 1.1, facecolor='g', alpha=0.1)
        ax5.axhspan(0.6, 0.8, facecolor='orange', alpha=0.1)
        ax5.axhspan(-0.1, 0.6, facecolor='r', alpha=0.1)
        
        ax1.axhline(0.,linestyle=':',linewidth=0.2,color='k')
        ax1.axhline(1.,linestyle=':',linewidth=0.2,color='k')
        ax2.axhline(0.,linestyle=':',linewidth=0.2,color='k')
        ax2.axhline(1.,linestyle=':',linewidth=0.2,color='k')
        ax3.axhline(0.,linestyle=':',linewidth=0.2,color='k')
        ax3.axhline(1.,linestyle=':',linewidth=0.2,color='k')
        ax4.axhline(0.,linestyle=':',linewidth=0.2,color='k')
        ax4.axhline(1.,linestyle=':',linewidth=0.2,color='k')
        ax5.axhline(0.,linestyle=':',linewidth=0.2,color='k')
        ax5.axhline(1.,linestyle=':',linewidth=0.2,color='k')
        
        for i,j in zip(range(len(p_chi_branch)),indexes_model):
            if total_list_scenario_name[j] ==total_list_scenario_name[0]:
                if weight_model_Mmax[i]==0 or score_nms[i]==0 or ordered_score_paleo[i]< 0.25 or p_chi_branch[i]<0.3:
                    ax1.scatter(m_Mmax[j],p_chi_branch[i],c='darkred',marker='_',s=15,alpha=0.2,linewidth=1)
                    ax2.scatter(m_Mmax[j],weight_model_Mmax[i],c='darkred',marker='_',s=15,alpha=0.2,linewidth=1)
                    ax3.scatter(m_Mmax[j],score_nms[i],c='darkred',marker='_',s=15,alpha=0.2,linewidth=1)
                    ax4.scatter(m_Mmax[j],ordered_score_paleo[i],c='darkred',marker='_',s=15,alpha=0.2,linewidth=1)
                    ax5.scatter(m_Mmax[j],final_weigth[i],c='darkred',marker='_',s=15,alpha=0.2,linewidth=1)
                else:
                    ax1.scatter(m_Mmax[j],p_chi_branch[i],c=color_mfd[i],marker='_',s=15,alpha=0.9,linewidth=1)
                    ax2.scatter(m_Mmax[j],weight_model_Mmax[i],c=color_mfd[i],marker='_',s=15,alpha=0.9,linewidth=1)
                    ax3.scatter(m_Mmax[j],score_nms[i],c=color_mfd[i],marker='_',s=15,alpha=0.9,linewidth=1)
                    ax4.scatter(m_Mmax[j],ordered_score_paleo[i],c=color_mfd[i],marker='_',s=15,alpha=0.9,linewidth=1)
                    ax5.scatter(m_Mmax[j],final_weigth[i],c=color_mfd[i],marker='_',s=15,alpha=0.9,linewidth=1)
            else:
                if weight_model_Mmax[i]==0 or score_nms[i]==0 or ordered_score_paleo[i]< 0.25 or p_chi_branch[i]<0.3:
                    ax1.scatter(m_Mmax[j],p_chi_branch[i],c='darkred',marker='|',s=15,alpha=0.2,linewidth=1)
                    ax2.scatter(m_Mmax[j],weight_model_Mmax[i],c='darkred',marker='|',s=15,alpha=0.2,linewidth=1)
                    ax3.scatter(m_Mmax[j],score_nms[i],c='darkred',marker='|',s=15,alpha=0.2,linewidth=1)
                    ax4.scatter(m_Mmax[j],ordered_score_paleo[i],c='darkred',marker='|',s=15,alpha=0.2,linewidth=1)
                    ax5.scatter(m_Mmax[j],final_weigth[i],c='darkred',marker='|',s=15,alpha=0.2,linewidth=1)
                else:
                    ax1.scatter(m_Mmax[j],p_chi_branch[i],c=color_mfd[i],marker='|',s=15,alpha=0.9,linewidth=1)
                    ax2.scatter(m_Mmax[j],weight_model_Mmax[i],c=color_mfd[i],marker='|',s=15,alpha=0.9,linewidth=1)
                    ax3.scatter(m_Mmax[j],score_nms[i],c=color_mfd[i],marker='|',s=15,alpha=0.9,linewidth=1)
                    ax4.scatter(m_Mmax[j],ordered_score_paleo[i],c=color_mfd[i],marker='|',s=15,alpha=0.9,linewidth=1)
                    ax5.scatter(m_Mmax[j],final_weigth[i],c=color_mfd[i],marker='|',s=15,alpha=0.9,linewidth=1)
        ax1.set_xlabel('Mmax')
        ax1.set_ylabel('test value   '+str(model))
        ax1.set_ylim([-0.05,1.05])
        ax1.set_xlim([xmax-1.5,xmax])
        ax2.set_xlim([xmax-1.5,xmax])
        ax3.set_xlim([xmax-1.5,xmax])
        ax4.set_xlim([xmax-1.5,xmax])
        ax5.set_xlim([xmax-1.5,xmax])
        ax1.set_title('chi test')
        ax2.set_title('Mmax test')
        ax3.set_title('NMS test')
        ax4.set_title('Paleo test')
        ax5.set_title('weitghted total')
        
        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/model_performance.png',dpi = 180)
        plt.close()
    
        index_model+=1
    
        # records in a file
        for i,j in zip(range(len(p_chi_branch)),indexes_model):
            file.write(str(model)+'\t'+str(total_list_MFD_type[j])+'\t'
                       +str(total_list_BG_hyp[j])+'\t'+str(total_list_scenario_name[j])+'\t'+str(total_list_sample[j])
            +'\t'+str(round(p_chi_branch[i],2))+'\t'+str(round(ordered_score_paleo[i],2))+'\t'+str(round(score_nms[i],2))+'\n')
    
    
    
    
    
    file.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
#        
#        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#        
#        ax1.axhspan(0.8, 1., facecolor='g', alpha=0.1)
#        ax2.axhspan(0.8, 1.,xmax=0.3, facecolor='g', alpha=0.1)
#        ax1.axhspan(0.6, 0.8, facecolor='orange', alpha=0.1)
#        ax2.axhspan(0.6, 0.8,xmax=0.5, facecolor='orange', alpha=0.1)
#        ax1.axhspan(-0.1, 0.6, facecolor='r', alpha=0.1)
#        ax2.axhspan(-0.1, 0.6,xmax=0.5, facecolor='r', alpha=0.1)
#        #ax2.axvspan(0., 30.,ymin=0.5, facecolor='g', alpha=0.1)
#        ax2.axvspan(30., 50.,ymin=0.818, facecolor='orange', alpha=0.1)
#        ax2.axvspan(50., 100, facecolor='r', alpha=0.1)
#        for i,j in zip(range(len(p_chi_branch)),indexes_model):
#            final_pvalue = p_chi_branch[i]*weight_chi + weight_model_Mmax[i]*weight_Mmax
#            if total_list_scenario_name[j] ==total_list_scenario_name[0]:
#                ax1.scatter(m_Mmax[j],final_pvalue,c=color_mfd[i],marker='_',s=10,alpha=0.9,linewidth=1)
#                ax2.scatter(a_s_model[j],final_pvalue,c=color_mfd[i],marker='_',s=10,alpha=0.9,linewidth=1)
#            else:
#                ax1.scatter(m_Mmax[j],final_pvalue,c=color_mfd[i],marker='|',s=10,alpha=0.9,linewidth=1)
#                ax2.scatter(a_s_model[j],final_pvalue,c=color_mfd[i],marker='|',s=10,alpha=0.9,linewidth=1)
#        ax1.set_xlabel('Mmax')
#        ax1.set_ylabel('final p value (chi test, Mmax)')
#        ax2.set_xlabel('NMS')
#        ax2.set_ylim([-0.1,1.])
#        ax1.set_xlim([xmax-1.5,xmax])
#        ax2.set_xlim([0.,100.])
#        ax1.set_title(str(model))
#        ax2.set_title('w_chi : '+str(weight_chi)+'   w_Mmax : '+str(weight_Mmax))
#        
#        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/model_performance_small.png',dpi = 180)
#        plt.show()
#    
#        index_model+=1
          
#                
#        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
#        for i,j in zip(range(len(p_chi_branch)),indexes_model):
#            if total_list_scenario_name[j] ==total_list_scenario_name[0]:
#                ax1.scatter(a_s_model[j],weight_model_Mmax[i],c=color_mfd[i],marker='_',s=10,alpha=0.4,linewidth=1)
#                ax2.scatter(p_chi_branch[i],weight_model_Mmax[i],c=color_mfd[i],marker='_',s=10,alpha=0.4,linewidth=1)
#            else:
#                ax1.scatter(a_s_model[j],weight_model_Mmax[i],c=color_mfd[i],marker='|',s=10,alpha=0.4,linewidth=1)
#                ax2.scatter(p_chi_branch[i],weight_model_Mmax[i],c=color_mfd[i],marker='|',s=10,alpha=0.4,linewidth=1)
#        ax1.set_xlabel('NMS')
#        ax2.set_xlabel('pvalue')
#        ax1.set_ylabel('weight Mmax')
#        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/fig_1.png',dpi = 180)
#        plt.close()
#    
#        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
#        for i,j in zip(range(len(p_chi_branch)),indexes_model):
#            if total_list_scenario_name[j] ==total_list_scenario_name[0]:
#                ax1.scatter(b_sample[j],p_chi_branch[i],c=color_mfd[i],marker='_',s=10,alpha=0.6,linewidth=1)
#                ax2.scatter(m_Mmax[j],p_chi_branch[i],c=color_mfd[i],marker='_',s=10,alpha=0.6,linewidth=1)
#            else:
#                ax1.scatter(b_sample[j],p_chi_branch[i],c=color_mfd[i],marker='|',s=10,alpha=0.6,linewidth=1)
#                ax2.scatter(m_Mmax[j],p_chi_branch[i],c=color_mfd[i],marker='|',s=10,alpha=0.6,linewidth=1)
#        ax1.set_xlabel('b value')
#        ax2.set_xlabel('Mmax')
#        ax1.set_ylabel('p value')
#        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/fig_2.png',dpi = 180)
#        plt.close()
        
        
#        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
#        for i,j in zip(range(len(p_chi_branch)),indexes_model):
#            final_pvalue = p_chi_branch[i]*weight_chi + weight_model_Mmax[i]*weight_Mmax
#            if total_list_scenario_name[j] ==total_list_scenario_name[0]:
#                ax1.scatter(b_sample[j],final_pvalue,c=color_mfd[i],marker='_',s=10,alpha=0.6,linewidth=1)
#                ax2.scatter(m_Mmax[j],final_pvalue,c=color_mfd[i],marker='_',s=10,alpha=0.6,linewidth=1)
#            else:
#                ax1.scatter(b_sample[j],final_pvalue,c=color_mfd[i],marker='|',s=10,alpha=0.6,linewidth=1)
#                ax2.scatter(m_Mmax[j],final_pvalue,c=color_mfd[i],marker='|',s=10,alpha=0.6,linewidth=1)
#        ax1.set_xlabel('b value')
#        ax2.set_xlabel('Mmax')
#        ax1.set_ylabel('final p value (chi test, Mmax)')
#        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/fig_3.png',dpi = 180)
#        plt.close()
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#        
#    
#    index_model = 0
#    for model in Model_list:
#        indexes_model = np.where(np.array(total_list_model) == model)[0] 
#        mfd_X = []
#        for index in indexes_model :  
#            mfd = mega_mfd_cummulative[index]
#            mfd_X.append(mfd)
#        b_model = np.take(b_sample,indexes_model)
#        Mmax_model = np.take(m_Mmax,indexes_model)
#        as_model = np.take(a_s_model,indexes_model)
#        if not os.path.exists(str(Run_name) + '/analysis/figures/sampling_analysis/'+model):
#            os.makedirs(str(Run_name) + '/analysis/figures/sampling_analysis/'+model)
#        mean_rate_catalog = np.array(catalog_cum_rate[index_model])#.mean(axis=0)
#        #std_rate_catalog = np.std(np.array(catalog_cum_rate[index_model]),axis=0)
#        
#        err_rate = [] # line for models, column for magnitudes, ratio between the model and the mean rate of the catalog
#        for mfd in mfd_X:
#            err_rate_i = []
#            for i in range(len(mean_rate_catalog)):
#                err_rate_j = mfd[i]/mean_rate_catalog[i]-1.
#                err_rate_i.append(err_rate_j)
#            err_rate.append(err_rate_i)
#            
#        colors = ['royalblue','steelblue','powderblue','lightgreen','gold','darkorange','darkred']
#        labels = ['<-1','-1<...<-0.5','-0.5<...<-0.2','-0.2<...<0.2','0.2<...<0.5','0.5<...<1','...>1']
#            
#        for index_i in range(int(len(bining_in_mag)/10.)):
#            Mmax1, Mmax2, Mmax3, Mmax4, Mmax5, Mmax6, Mmax7 = [], [], [], [], [], [], []
#            b1, b2, b3, b4, b5, b6, b7 = [], [], [], [], [], [], []
#            index=0
#            for err in np.array(err_rate)[:,index_i*10] :
#                if err < -1. :
#                    Mmax1.append(Mmax_model[index])
#                    b1.append(b_model[index])
#                elif err < -0.5 :
#                    Mmax2.append(Mmax_model[index])
#                    b2.append(b_model[index])
#                elif err < - 0.2 :
#                    Mmax3.append(Mmax_model[index])
#                    b3.append(b_model[index])
#                elif err < 0.2 :
#                    Mmax4.append(Mmax_model[index])
#                    b4.append(b_model[index])
#                elif err < 0.5 :
#                    Mmax5.append(Mmax_model[index])
#                    b5.append(b_model[index])
#                elif err < 1. :
#                    Mmax6.append(Mmax_model[index])
#                    b6.append(b_model[index])
#                elif err > 1.  :
#                    Mmax7.append(Mmax_model[index])
#                    b7.append(b_model[index])
#                index+=1
#            for color, label, b, Mmax in zip(colors, labels, [b1, b2, b3, b4, b5, b6, b7], [Mmax1, Mmax2, Mmax3, Mmax4, Mmax5, Mmax6, Mmax7]):
#                plt.scatter(b,Mmax,c=color,s= 50 ,alpha=0.8,label=label)
#            plt.title('Modeled rate / Catalog rate')
#            plt.legend(loc=2,fontsize=6)
#            plt.xlabel('b value')
#            plt.ylabel('Mmax')
#            plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/'+str(model)+'_b_Mmax_vs_error_M'+str(bining_in_mag[index_i*10])+'.png',dpi = 100)
#            plt.close()
#            
#        colors = ['royalblue','steelblue','darkorange','red','darkred']
#        labels = ['<10%','10%<...<30%','30%<...<50%','50%<...<70%','...>70%']
#        Mmax1, Mmax2, Mmax3, Mmax4, Mmax5 = [], [], [], [], []
#        b1, b2, b3, b4, b5 = [], [], [], [], []
#        index=0
#        for NMS in as_model :
#            if NMS < 10. :
#                Mmax1.append(Mmax_model[index])
#                b1.append(b_model[index])
#            elif NMS < 30. :
#                Mmax2.append(Mmax_model[index])
#                b2.append(b_model[index])
#            elif NMS < 50. :
#                Mmax3.append(Mmax_model[index])
#                b3.append(b_model[index])
#            elif NMS < 70. :
#                Mmax4.append(Mmax_model[index])
#                b4.append(b_model[index])
#            elif NMS > 70. :
#                Mmax5.append(Mmax_model[index])
#                b5.append(b_model[index])
#            index+=1
#        for color, label, b, Mmax in zip(colors, labels, [b1, b2, b3, b4, b5], [Mmax1, Mmax2, Mmax3, Mmax4, Mmax5]):
#            plt.scatter(b,Mmax,c=color,s= 50 ,alpha=0.8,label=label)
#        plt.title('NMS is the model')
#        plt.legend(loc=2,fontsize=6)
#        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/'+str(model)+'_b_Mmax_vs_NMS.png',dpi = 100)
#        plt.close()
#    
#    
#    
    
    file_LT_metrics.close()
    
    
