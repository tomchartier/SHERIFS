# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np

def get_new_target(number_of_loops,moment_rate_in_bin,p_MFD_MO,target_moment_per_bin,bin_mag,empty_bins,bin_target_reached,rup_in_bin):
 if number_of_loops > 20:
     shape_mfd_i = (moment_rate_in_bin)/sum(moment_rate_in_bin)
     target_i = []
     for i in range(len(bin_mag)) :
         #build the target based on the original target modified according to the evolution of the calculation
         malus = 0.01
         ratio = (p_MFD_MO[i] - shape_mfd_i[i])/p_MFD_MO[i]
         
         if ratio > 0. :
             target_i.append(p_MFD_MO[i] * 1./malus)
         else :
             target_i.append(min(p_MFD_MO) * malus)
             
         if i in bin_target_reached:
             target_i[i] = 0.  # don't pick bins where target is
         elif i in empty_bins or len(rup_in_bin[i]) == 0:
             target_i[i] = 0.  #don't pick empty bins reached
         elif target_i[i] <= 0. :
             target_i[i] = 0. #can't have a negative number in the distribution so we put a very small one (10E-15 N.m)
         else :
            target_i[i] = p_MFD_MO[i]
 else:
     target_i = p_MFD_MO

 try :
    target_i = (target_i)/sum(target_i)# normalize the target to use it as a probability distribution
 except :
    target_i = p_MFD_MO
 return target_i
