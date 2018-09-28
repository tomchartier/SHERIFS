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

def sampling_analysis(Run_name,Model_list,m_Mmax,b_sample,a_s_model,mega_mfd_cummulative,catalog_cum_rate,
                      xmin,xmax,ymin,ymax,total_list_model,bining_in_mag):
    
    if not os.path.exists(str(Run_name) + '/analysis/figures/sampling_analysis'):
        os.makedirs(str(Run_name) + '/analysis/figures/sampling_analysis')
    index_model = 0
    for model in Model_list:
        indexes_model = np.where(np.array(total_list_model) == model)[0] 
        mfd_X = []
        for index in indexes_model :  
            mfd = mega_mfd_cummulative[index]
            mfd_X.append(mfd)
        b_model = np.take(b_sample,indexes_model)
        Mmax_model = np.take(m_Mmax,indexes_model)
        as_model = np.take(a_s_model,indexes_model)
        if not os.path.exists(str(Run_name) + '/analysis/figures/sampling_analysis/'+model):
            os.makedirs(str(Run_name) + '/analysis/figures/sampling_analysis/'+model)
        mean_rate_catalog = np.array(catalog_cum_rate[index_model]).mean(axis=0)
        std_rate_catalog = np.std(np.array(catalog_cum_rate[index_model]),axis=0)
        
        err_rate = [] # line for models, column for magnitudes, ratio between the model and the mean rate of the catalog
        for mfd in mfd_X:
            err_rate_i = []
            for i in range(len(mean_rate_catalog)):
                err_rate_j = mfd[i]/mean_rate_catalog[i]-1.
                err_rate_i.append(err_rate_j)
            err_rate.append(err_rate_i)
            
        colors = ['royalblue','steelblue','powderblue','lightgreen','gold','darkorange','darkred']
        labels = ['<-1','-1<...<-0.5','-0.5<...<-0.2','-0.2<...<0.2','0.2<...<0.5','0.5<...<1','...>1']
            
        for index_i in range(int(len(bining_in_mag)/10.)):
            Mmax1, Mmax2, Mmax3, Mmax4, Mmax5, Mmax6, Mmax7 = [], [], [], [], [], [], []
            b1, b2, b3, b4, b5, b6, b7 = [], [], [], [], [], [], []
            index=0
            for err in np.array(err_rate)[:,index_i*10] :
                if err < -1. :
                    Mmax1.append(Mmax_model[index])
                    b1.append(b_model[index])
                elif err < -0.5 :
                    Mmax2.append(Mmax_model[index])
                    b2.append(b_model[index])
                elif err < - 0.2 :
                    Mmax3.append(Mmax_model[index])
                    b3.append(b_model[index])
                elif err < 0.2 :
                    Mmax4.append(Mmax_model[index])
                    b4.append(b_model[index])
                elif err < 0.5 :
                    Mmax5.append(Mmax_model[index])
                    b5.append(b_model[index])
                elif err < 1. :
                    Mmax6.append(Mmax_model[index])
                    b6.append(b_model[index])
                elif err > 1.  :
                    Mmax7.append(Mmax_model[index])
                    b7.append(b_model[index])
                index+=1
            for color, label, b, Mmax in zip(colors, labels, [b1, b2, b3, b4, b5, b6, b7], [Mmax1, Mmax2, Mmax3, Mmax4, Mmax5, Mmax6, Mmax7]):
                plt.scatter(b,Mmax,c=color,s= 50 ,alpha=0.8,label=label)
            plt.title('Modeled rate / Catalog rate')
            plt.legend(loc=2,fontsize=6)
            plt.xlabel('b value')
            plt.ylabel('Mmax')
            plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/'+str(model)+'_b_Mmax_vs_error_M'+str(bining_in_mag[index_i*10])+'.png',dpi = 100)
            plt.close()
            
        colors = ['royalblue','steelblue','darkorange','red','darkred']
        labels = ['<10%','10%<...<30%','30%<...<50%','50%<...<70%','...>70%']
        Mmax1, Mmax2, Mmax3, Mmax4, Mmax5 = [], [], [], [], []
        b1, b2, b3, b4, b5 = [], [], [], [], []
        index=0
        for NMS in as_model :
            if NMS < 10. :
                Mmax1.append(Mmax_model[index])
                b1.append(b_model[index])
            elif NMS < 30. :
                Mmax2.append(Mmax_model[index])
                b2.append(b_model[index])
            elif NMS < 50. :
                Mmax3.append(Mmax_model[index])
                b3.append(b_model[index])
            elif NMS < 70. :
                Mmax4.append(Mmax_model[index])
                b4.append(b_model[index])
            elif NMS > 70. :
                Mmax5.append(Mmax_model[index])
                b5.append(b_model[index])
            index+=1
        for color, label, b, Mmax in zip(colors, labels, [b1, b2, b3, b4, b5], [Mmax1, Mmax2, Mmax3, Mmax4, Mmax5]):
            plt.scatter(b,Mmax,c=color,s= 50 ,alpha=0.8,label=label)
        plt.title('NMS is the model')
        plt.legend(loc=2,fontsize=6)
        plt.savefig(str(Run_name) + '/analysis/figures/sampling_analysis/'+model+'/'+str(model)+'_b_Mmax_vs_NMS.png',dpi = 100)
        plt.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    