#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

Return the maximum magnitude of the faults and scenarios using a scalling relationship

@author: thomas
"""
import numpy as np
import scipy.stats as stats

class Calc_Mmax():
    def __init__(self,faults_area,scenario_area,faults_length,scenario_length,faults_width,scenario_width,
                 selected_ScL,dimention_used,use_all_ScL_data,faults_mecanism,index_faults_in_scenario, sample):
        self.scenario_area = scenario_area
        self.faults_area = faults_area
        self.scenario_length = scenario_length
        self.faults_length = faults_length
        self.scenario_width = scenario_width
        self.faults_width = faults_width
        self.selected_ScL = selected_ScL
        self.dimention_used = dimention_used
        self.use_all_ScL_data = use_all_ScL_data
        self.faults_mecanism = faults_mecanism
        self.index_faults_in_scenario = index_faults_in_scenario
        self.sample = sample
        
       
        self.initialize()
    def initialize(self):
        self.Mmax_faults = [] #Mmax of each fault
        self.Mmax_scenario = [] #Mmax of each scenario
        
        # find the dominant fault kinematic of the FtF scenario 
        scenario_mecanism = []
        index_scenario = 0    
        for length in self.scenario_length:
            a= np.take(self.faults_mecanism,self.index_faults_in_scenario[index_scenario][0])
            unique,pos = np.unique(a,return_inverse=True) #Finds all unique elements and their positions
            counts = np.bincount(pos)                     #Count the number of each unique element
            maxpos = counts.argmax()                      #Finds the positions of the maximum count
            
            scenario_mecanism.append(unique[maxpos])
            
            index_scenario += 1
            
        # check if all the aspect ratio of the faults and scenario are okay. 
        # if some are not acceptable (according to TMG201- relationship)
        # need to ask the authors to provide the right equation for each fault kinematics
        # not accepatble is defined as two sigmas more than the median value of the regression
        #  if not, the dimension (area or length) is change to the maximum acceptable value
        # therefor, the Mmax will be reduced artificially even if the whole surface ruptures
        self.final_fault_length = []
        self.final_scenario_length = []
        print_ar = True
        i=0
        for l,w in zip([x / 1000. for x in self.faults_length],self.faults_width):
            string = 'aspect ratio acceptable'
            if self.faults_mecanism[i] == 'N':
                if np.log10(l) > np.log10(w)*(1.5+2.*0.21)-0.45:
                    string = 'aspect ratio too large'
#                    self.faults_length[i] = 10**(np.log10(w)*(1.5+2.*0.21)-0.45) *1000.
#                    self.faults_area[i] = 10**(np.log10(w)*(1.5+2.*0.21)-0.45) *1000. * w  *1000.
                    if print_ar == True :
#                        print w,l,10**(np.log10(w)*(1.5+2.*0.21)-0.45)
                        #print('Warning ! : Some faults or scenarios have a extremely large aspect ratio')
                        print_ar = False
                
            if self.faults_mecanism[i] == 'R':
                if np.log10(l) > np.log10(w)*(1.39+2.*0.09)-0.29:
                    string = 'aspect ratio too large'
#                    self.faults_length[i] = 10**(np.log10(w)*(1.39+2.*0.09)-0.29) *1000.
#                    self.faults_area[i] = 10**(np.log10(w)*(1.39+2.*0.09)-0.29) *1000. * w  *1000.
                    if print_ar == True :
                        #print('Warning ! : Some faults or scenarios have a extremely large aspect ratio')
                        print_ar = False
                
            if self.faults_mecanism[i] == 'S':
                if np.log10(l) > np.log10(w)*(2.62+2.*0.20)-1.52:
                    string = 'aspect ratio too large'
#                    self.faults_length[i] = 10**(np.log10(w)*(2.62+2.*0.20)-1.52) *1000.
#                    self.faults_area[i] = 10**(np.log10(w)*(2.62+2.*0.20)-1.52) *1000. * w  *1000.
                    if print_ar == True :
                        #print('Warning ! : Some faults or scenarios have a extremely large aspect ratio')
                        print_ar = False
            self.final_fault_length.append(string)
            i+=1
            
        i=0
        for l,w in zip([x / 1000. for x in self.scenario_length],self.scenario_width):
            string = 'aspect ratio acceptable'
            if scenario_mecanism[i] == 'N':
                if np.log10(l) > (np.log10(w)*(1.5+2.*0.21)-0.45):
                    string = 'aspect ratio too large'
#                    self.scenario_length[i] = 10**(np.log10(w)*(1.5+2.*0.21)-0.45) *1000.
#                    self.scenario_area[i] = 10**(np.log10(w)*(1.5+2.*0.21)-0.45)  *1000. * w  *1000.
                    if print_ar == True :
                        #print( w,l,10**(np.log10(w)*(1.5+2.*0.21)-0.45))
                        #print('Warning ! : Some faults or scenarios have a extremely large aspect ratio')
                        print_ar = False
                
            if scenario_mecanism[i] == 'R':
                if np.log10(l) > (np.log10(w)*(1.39+2.*0.09)-0.29):
                    string = 'aspect ratio too large'
#                    self.scenario_length[i] = 10**(np.log10(w)*(1.39+2.*0.09)-0.29) *1000.
#                    self.scenario_area[i] = 10**(np.log10(w)*(1.39+2.*0.09)-0.29)  *1000. * w  *1000.
                    if print_ar == True :
                        #print('Warning ! : Some faults or scenarios have a extremely large aspect ratio')
                        print_ar = False
                
            if scenario_mecanism[i] == 'S':
                if np.log10(l) > (np.log10(w)*(2.62+2.*0.20)-1.52):
                    string = 'aspect ratio too large'
#                    self.scenario_length[i] = 10**(np.log10(w)*(2.62+2.*0.20)-1.52) *1000.
#                    self.scenario_area[i] = 10**(np.log10(w)*(2.62+2.*0.20)-1.52)  *1000. * w  *1000.
                    if print_ar == True :
                        #print('Warning ! : Some faults or scenarios have a extremely large aspect ratio')
                        print_ar = False
            self.final_scenario_length.append(string)
                
            i+=1


        ######################################################################################""
        # in the case the scaling law Wells and Coppersmith 1994 is chosen
        if self.selected_ScL == 'WC1994' :
#            matrix, row 1 for normal; row 2 for reverse; row 3 for strikeslip; row 4 for ALL
#            each row contains aRup_Length, bRup_Length, sdRup_Length, aRup_Area, bRup_Area, sdRup_Area
            coeff_ref=np.array(([4.34, 1.54, 0.31, 3.93, 1.02, 0.25],[4.49, 1.49, 0.26, 4.33, 0.90, 0.25],[4.33, 1.49, 0.24, 3.98, 1.02, 0.23],[4.38, 1.49, 0.26, 4.07, 0.98, 0.24]))
            
            #random sampling of the uncertainty on the scaling law
            if self.sample == 1 : #for the first sample, the mean value is prefered
                
                coeff=np.array(([4.34, 1.54, 0., 3.93, 1.02, 0.],[4.49, 1.49, 0., 4.33, 0.90, 0.],[4.33, 1.49, 0., 3.98, 1.02, 0.],[4.38, 1.49, 0., 4.07, 0.98, 0.]))
            else:
                #triangular sampling
#                coeff=np.array(([4.34, 1.54, np.random.triangular(-coeff[0][2],0,coeff[0][2]), 3.93, 1.02, np.random.triangular(-coeff[0][5],0,coeff[0][5])]
#                ,[4.49, 1.49, np.random.triangular(-coeff[1][2],0,coeff[1][2]), 4.33, 0.90, np.random.triangular(-coeff[1][5],0,coeff[1][5])]
#                ,[4.33, 1.49, np.random.triangular(-coeff[2][2],0,coeff[2][2]), 3.98, 1.02, np.random.triangular(-coeff[2][5],0,coeff[2][5])]
#                ,[4.38, 1.49, np.random.triangular(-coeff[3][2],0,coeff[3][2]), 4.07, 0.98, np.random.triangular(-coeff[3][5],0,coeff[3][5])]))
                #normal sampling
                coeff=np.array(([4.34, 1.54, np.random.normal(0,coeff_ref[0][2]), 3.93, 1.02, np.random.normal(0,coeff_ref[0][5])]
                ,[4.49, 1.49, np.random.normal(0,coeff_ref[1][2]), 4.33, 0.90, np.random.normal(0,coeff_ref[1][5])]
                ,[4.33, 1.49, np.random.normal(0,coeff_ref[2][2]), 3.98, 1.02, np.random.normal(0,coeff_ref[2][5])]
                ,[4.38, 1.49, np.random.normal(0,coeff_ref[3][2]), 4.07, 0.98, np.random.normal(0,coeff_ref[3][5])]))
                #truncated normal sampling
                coeff=np.array(([4.34, 1.54, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[0][2]).rvs(1)[0], 3.93, 1.02, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[0][5]).rvs(1)[0]]
                ,[4.49, 1.49, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[1][2]).rvs(1)[0], 4.33, 0.90, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[1][5]).rvs(1)[0]]
                ,[4.33, 1.49, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[2][2]).rvs(1)[0], 3.98, 1.02, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[2][5]).rvs(1)[0]]
                ,[4.38, 1.49, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[3][2]).rvs(1)[0], 4.07, 0.98, stats.truncnorm(-2., 2., loc=0., scale=coeff_ref[3][5]).rvs(1)[0]]))
            
                
            if self.dimention_used == 'L' :  #the length is used 
                index_fault = 0
                for length in self.faults_length:
                    if self.use_all_ScL_data == True :
                        Mmax_i = coeff[3][0] + coeff[3][1] * np.log10(length/1000.) + coeff[3][2]
                    else:
                        if self.faults_mecanism[index_fault] == 'N':
                            Mmax_i = coeff[0][0] + coeff[0][1] * np.log10(length/1000.) + coeff[0][2]
                            
                        if self.faults_mecanism[index_fault] == 'R':
                            Mmax_i = coeff[1][0] + coeff[1][1] * np.log10(length/1000.) + coeff[1][2]
                            
                        if self.faults_mecanism[index_fault] == 'S':
                            Mmax_i = coeff[2][0] + coeff[2][1] * np.log10(length/1000.) + coeff[2][2]
                    self.Mmax_faults.append(round(Mmax_i,1))
                    
                    index_fault += 1
                    
                index_scenario = 0    
                for length in self.scenario_length:
                    if self.use_all_ScL_data == True :
                        Mmax_i = coeff[3][0] + coeff[3][1] * np.log10(length/1000.) + coeff[3][2]
                    else:
                        if scenario_mecanism[index_scenario] == 'N':
                            Mmax_i = coeff[0][0] + coeff[0][1] * np.log10(length/1000.) + coeff[0][2]
                            
                        if scenario_mecanism[index_scenario] == 'R':
                            Mmax_i = coeff[1][0] + coeff[1][1] * np.log10(length/1000.) + coeff[1][2]
                            
                        if scenario_mecanism[index_scenario] == 'S':
                            Mmax_i = coeff[2][0] + coeff[2][1] * np.log10(length/1000.) + coeff[2][2]
                    self.Mmax_scenario.append(round(Mmax_i,1))
                    
                    index_scenario += 1
                    
            
            if self.dimention_used == 'A' : #the area is used 
                index_fault = 0
                for area in self.faults_area:
                    if self.use_all_ScL_data == True :
                        Mmax_i = coeff[3][3] + coeff[3][4] * np.log10(area/1000000.) + coeff[3][5]
                    else:
                        if self.faults_mecanism[index_fault] == 'N':
                            Mmax_i = coeff[0][3] + coeff[0][4] * np.log10(area/1000000.) + coeff[0][5]
                            
                        if self.faults_mecanism[index_fault] == 'R':
                            Mmax_i = coeff[1][3] + coeff[1][4] * np.log10(area/1000000.) + coeff[1][5]
                            
                        if self.faults_mecanism[index_fault] == 'S':
                            Mmax_i = coeff[2][3] + coeff[2][4] * np.log10(area/1000000.) + coeff[2][5]
                    self.Mmax_faults.append(round(Mmax_i,1))
                    
                    index_fault += 1
                    
                index_scenario = 0 
                for area in self.scenario_area : 
                    if self.use_all_ScL_data == True :
                        Mmax_i = coeff[3][3] + coeff[3][4] * np.log10(area/1000000.) + coeff[3][5]
                    else : 
                        if scenario_mecanism[index_scenario] == 'N':
                            Mmax_i = coeff[0][3] + coeff[0][4] * np.log10(area/1000000.) + coeff[0][5]
                            
                        if scenario_mecanism[index_scenario] == 'R':
                            Mmax_i = coeff[1][3] + coeff[1][4] * np.log10(area/1000000.) + coeff[1][5]
                            
                        if scenario_mecanism[index_scenario] == 'S':
                            Mmax_i = coeff[2][3] + coeff[2][4] * np.log10(area/1000000.) + coeff[2][5]
                    self.Mmax_scenario.append(round(Mmax_i,1))
                    
                    index_scenario += 1
                    
           
                    
                    
                    
        #######################################################################
        # in the case the scaling law Leonard 2010 is chosen
        if self.selected_ScL == 'Le2010' :
#            matrix, row 1 for DS dip slip; row 2 for SS strike slip
#            each row contains aLength, baLengthmin, baLengthmax, aArea, bAreamin, bAreamax (Table 5 in Leonard, 2010)

            coeff=np.array(([2.5, 7.53, 8.51, 1.5, 5.69, 6.6],
                            [1.5, 12.01, 12.88, 1.5, 5.69, 6.47]))
            
            #random sampling of the uncertainty on the scaling law 
            #d_le10 is a number between 0 and 1, will locate the Mmax between the Mmax_min (d_le10 = 0) ad Mmax_max (d_le10 = 1)
            if self.sample == 1 : #for the first sample, the mean value is prefered
                d_le10 = 0.5
            else :
                d_le10 = np.random.triangular(0.,0.5,1.) #triangular picking
                         
            if self.dimention_used == 'L' : #the length is used 
                index_fault = 0
                for length in self.faults_length:
                    if self.faults_mecanism[index_fault] == 'N' or self.faults_mecanism[index_fault] == 'R':
                        if length > 5500. :
                            Mmax_min = 2. / 3. * (coeff[0][1] + coeff[0][0] * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (coeff[0][2] + coeff[0][0] * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        else:
                            Mmax_i = 2. / 3. * (6.1 + 3. * np.log10(length)) - 6.07
                    if self.faults_mecanism[index_fault] == 'S':
                        if length < 3400.:
                            Mmax_min = 2. / 3. * (5.65 + 3.0 * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (6.52 + 3.0 * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        elif length <45000.:
                            Mmax_min = 2. / 3. * (7.41 + 2.5 * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (8.28 + 2.5 * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        else :
                            Mmax_min = 2. / 3. (coeff[1][1] + coeff[1][0] * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. (coeff[1][2] + coeff[1][0] * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                    self.Mmax_faults.append(round(Mmax_i,1))
                    
                    index_fault += 1
                    
                index_scenario = 0    
                for length in self.scenario_length:
                    if scenario_mecanism[index_scenario] == 'N' or scenario_mecanism[index_scenario] == 'R':
                        if length > 5500. :
                            Mmax_min = 2. / 3. * (coeff[0][1] + coeff[0][0] * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (coeff[0][2] + coeff[0][0] * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        else:
                            Mmax_i = 2. / 3. * (6.1 + 3. * np.log10(length)) - 6.07
                        
                    if scenario_mecanism[index_scenario] == 'S':
                        if length < 3400.:
                            Mmax_min = 2. / 3. * (5.65 + 3.0 * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (6.52 + 3.0 * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        elif length <45000.:
                            Mmax_min = 2. / 3. * (7.41 + 2.5 * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (8.28 + 2.5 * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        else :
                            Mmax_min = 2. / 3. * (coeff[1][1] + coeff[1][0] * np.log10(length)) - 6.07
                            Mmax_max = 2. / 3. * (coeff[1][2] + coeff[1][0] * np.log10(length)) - 6.07
                            Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                    self.Mmax_scenario.append(round(Mmax_i,1))
                    
                    index_scenario += 1
                    
            
            if self.dimention_used == 'A' : #the area is used 
                index_fault = 0
                for area in self.faults_area:
                    if self.faults_mecanism[index_fault] == 'N' or self.faults_mecanism[index_fault] == 'R':
                        Mmax_min = 2. / 3. * (coeff[0][4] + coeff[0][3] * np.log10(area)) - 6.07
                        Mmax_max = 2. / 3. * (coeff[0][5] + coeff[0][3] * np.log10(area)) - 6.07
                        Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        
                    if self.faults_mecanism[index_fault] == 'S':
                        Mmax_min = 2. / 3. * (coeff[1][4] + coeff[1][3] * np.log10(area)) - 6.07
                        Mmax_max = 2. / 3. * (coeff[1][5] + coeff[1][3] * np.log10(area)) - 6.07
                        Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                    self.Mmax_faults.append(round(Mmax_i,1))
                    
                    index_fault += 1
                    
                index_scenario = 0 
                for area in self.scenario_area:
                    if scenario_mecanism[index_scenario] == 'N' or scenario_mecanism[index_scenario] == 'R':
                        Mmax_min = 2. / 3. * (coeff[0][4] + coeff[0][3] * np.log10(area)) - 6.07
                        Mmax_max = 2. / 3. * (coeff[0][5] + coeff[0][3] * np.log10(area)) - 6.07
                        Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                        
                    if scenario_mecanism[index_scenario] == 'S':
                        Mmax_min = 2. / 3. * (coeff[1][4] + coeff[1][3] * np.log10(area)) - 6.07
                        Mmax_max = 2. / 3. * (coeff[1][5] + coeff[1][3] * np.log10(area)) - 6.07
                        Mmax_i = Mmax_min + d_le10 * (Mmax_max - Mmax_min)
                    self.Mmax_scenario.append(round(Mmax_i,1))
                    
                    index_scenario += 1
        
                    
                    
        ####################################################################################            
        # in the case the scaling law Hanks and Bakun 2008 is chosen
        # only for Area
        # two cases depending if area is smaller or bigger than 537 km2
        if self.selected_ScL == 'HB08' :
            index_fault = 0
            if self.sample == 1 : #for the first sample, the mean value is prefered
                d_i = 0.
            else :
                d_i = np.random.normal(0,1.)
                d_i = stats.truncnorm(-2., 2., loc=0., scale=1.).rvs(1)[0]
            
            for area in self.faults_area:
                area = area/1000000.
                if area < 537.:
                  Mmax_i = np.log10(area) + 3.98 + d_i * 0.03
                else:
                  Mmax_i = 4./3.* np.log10(area) + 3.07 + d_i * 0.04
                self.Mmax_faults.append(round(Mmax_i,1))
                index_fault += 1
                
                
            index_scenario = 0 
            for area in self.scenario_area : 
                area = area/1000000.
                if area < 537.:
                  Mmax_i = np.log10(area) + 3.98 + d_i * 0.03
                else:
                  Mmax_i = 4./3.* np.log10(area) + 3.07 + d_i * 0.04
                self.Mmax_scenario.append(round(Mmax_i,1))
                index_scenario += 1
                
                
                
                
                
                
                
        ####################################################################################            
        # in the case the scaling law Thingbaijam Mai and Goda 2017
        #  mean value works -- need to find the good way to explore the uncertainties
        #
        if self.selected_ScL == 'TMG2017' :
#            matrix, row 1 for normal; row 2 for reverse; row 3 for strikeslip; row 4 for subduction interface
#            each row contains aRup_Length, bRup_Length, saRup_Length, sbRup_Length, aRup_Area, bRup_Area, saRup_Area, sbRup_Area
            coeff=np.array(([ -1.722, 0.485, 0.260, 0.036, -2.551, 0.808, 0.423, 0.059],
                            [ -2.693, 0.614, 0.292, 0.043, -4.362, 1.049, 0.445, 0.066],
                            [ -2.943, 0.681, 0.357, 0.052, -3.486, 0.942, 0.399, 0.058],
                            [ -2.412, 0.583, 0.288, 0.037, -3.292, 0.949, 0.377, 0.049]))
            
#            each row contains aRup_Length, bRup_Length, sigmaMRup_Length,  aRup_Area, bRup_Area, sigmaMRup_Area
            coeff=np.array(([ -1.722, 0.485, 0.128, -2.551, 0.808, 0.181],
                            [ -2.693, 0.614, 0.083, -4.362, 1.049, 0.121],
                            [ -2.943, 0.681, 0.151, -3.486, 0.942, 0.184],
                            [ -2.412, 0.583, 0.107, -3.292, 0.949, 0.150]))
            
            #random sampling of the uncertainty on the scaling law
            if self.sample == 1 : #for the first sample, the mean value is prefered
                coeff=np.array(([ -1.722, 0.485, 0., -2.551, 0.808, 0.],
                                [ -2.693, 0.614, 0., -4.362, 1.049, 0.],
                                [ -2.943, 0.681, 0., -3.486, 0.942, 0.],
                                [ -2.412, 0.583, 0., -3.292, 0.949, 0.]))
            else:
                coeff=np.array(([ -1.722, 0.485, np.random.normal(0,0.128), -2.551, 0.808, np.random.normal(0,0.181)],
                                [ -2.693, 0.614, np.random.normal(0,0.083), -4.362, 1.049, np.random.normal(0,0.121)],
                                [ -2.943, 0.681, np.random.normal(0,0.151), -3.486, 0.942, np.random.normal(0,0.184)],
                                [ -2.412, 0.583, np.random.normal(0,0.107), -3.292, 0.949, np.random.normal(0,0.150)]))
                #truncated normal distribution
                coeff=np.array(([ -1.722, 0.485, stats.truncnorm(-2., 2., loc=0., scale=0.128).rvs(1)[0], -2.551, 0.808, stats.truncnorm(-2., 2., loc=0., scale=0.181).rvs(1)[0]],
                                [ -2.693, 0.614, stats.truncnorm(-2., 2., loc=0., scale=0.083).rvs(1)[0], -4.362, 1.049, stats.truncnorm(-2., 2., loc=0., scale=0.121).rvs(1)[0]],
                                [ -2.943, 0.681, stats.truncnorm(-2., 2., loc=0., scale=0.151).rvs(1)[0], -3.486, 0.942, stats.truncnorm(-2., 2., loc=0., scale=0.184).rvs(1)[0]],
                                [ -2.412, 0.583, stats.truncnorm(-2., 2., loc=0., scale=0.107).rvs(1)[0], -3.292, 0.949, stats.truncnorm(-2., 2., loc=0., scale=0.150).rvs(1)[0]]))


            if self.dimention_used == 'L' :  #the length is used 
                index_fault = 0
                for length in self.faults_length:
                    if self.faults_mecanism[index_fault] == 'N':
                        Mmax_i = (np.log10(length/1000.) - (coeff[0][0]))/(coeff[0][1]) + coeff[0][2]
                        
                    if self.faults_mecanism[index_fault] == 'R':
                        Mmax_i = (np.log10(length/1000.) - (coeff[1][0]))/(coeff[1][1]) + coeff[1][2]
                        
                    if self.faults_mecanism[index_fault] == 'S':
                        Mmax_i = (np.log10(length/1000.) - (coeff[2][0]))/(coeff[2][1]) + coeff[2][2]

                    if self.faults_mecanism[index_fault] == 'Sub':
                        Mmax_i = (np.log10(length/1000.) - (coeff[3][0]))/(coeff[3][1]) + coeff[3][2]
                    self.Mmax_faults.append(round(Mmax_i,1))
                    
                    index_fault += 1
                    
                index_scenario = 0    
                for length in self.scenario_length:
                    if scenario_mecanism[index_scenario] == 'N':
                        Mmax_i = (np.log10(length/1000.) - (coeff[0][0]))/(coeff[0][1]) + coeff[0][2]
                        
                    if scenario_mecanism[index_scenario] == 'R':
                        Mmax_i = (np.log10(length/1000.) - (coeff[1][0]))/(coeff[1][1]) + coeff[1][2]
                        
                    if scenario_mecanism[index_scenario] == 'S':
                        Mmax_i = (np.log10(length/1000.) - (coeff[2][0]))/(coeff[2][1]) + coeff[2][2]

                    if scenario_mecanism[index_scenario] == 'Sub':
                        Mmax_i = (np.log10(length/1000.) - (coeff[3][0]))/(coeff[3][1]) + coeff[3][2]
                    self.Mmax_scenario.append(round(Mmax_i,1))
                    
                    index_scenario += 1
                    
            
            if self.dimention_used == 'A' : #the area is used 
                index_fault = 0
                for area in self.faults_area:
                    if self.faults_mecanism[index_fault] == 'N':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[0][3]))/(coeff[0][4]) + coeff[0][5]
                        
                    if self.faults_mecanism[index_fault] == 'R':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[1][3]))/(coeff[1][4]) + coeff[1][5]
                        
                    if self.faults_mecanism[index_fault] == 'S':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[2][3]))/(coeff[2][4]) + coeff[2][5]

                    if self.faults_mecanism[index_fault] == 'Sub':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[3][3]))/(coeff[3][4]) + coeff[3][5]
                    self.Mmax_faults.append(round(Mmax_i,1))
                    
                    index_fault += 1
                    
                index_scenario = 0 
                for area in self.scenario_area : 
                    
                    if scenario_mecanism[index_scenario] == 'N':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[0][3]))/(coeff[0][4]) + coeff[0][5]
                        
                    if scenario_mecanism[index_scenario] == 'R':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[1][3]))/(coeff[1][4]) + coeff[1][5]
                        
                    if scenario_mecanism[index_scenario] == 'S':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[2][3]))/(coeff[2][4]) + coeff[2][5]

                    if scenario_mecanism[index_scenario] == 'Sub':
                        Mmax_i = (np.log10(area/1000000.) - (coeff[3][3]))/(coeff[3][4]) + coeff[3][5]
                    self.Mmax_scenario.append(round(Mmax_i,1))
                    
                    index_scenario += 1
                
                
                
                
                
                
                
                
                
        ####################################################################################            
        # in the case the scaling law Thingbaijam Mai and Goda 2017
        #  mean value works -- need to find the good way to explore the uncertainties
        #
        if self.selected_ScL == 'Shaw2009mod' :
                
            if self.sample == 1 : #for the first sample, the mean value is prefered
                d_i = 0.
            else :
                d_i = np.random.normal(0,1.)
                d_i = stats.truncnorm(-2., 2., loc=0., scale=1.).rvs(1)[0]
            
            index_fault = 0
            for area,lenght in zip(self.faults_area,self.faults_length):
                area = area/1000000.
                length = length/1000.
                width = area/length
                Mmax_i = np.log10(area)+2./3.*np.log10(max(1.,np.sqrt(area/width**2.))/((1.+max(1.,area/(width**2.*7.4)))/2.))+3.98
                Mmax_i = Mmax_i + d_i * 0.2
                self.Mmax_faults.append(round(Mmax_i,1))
                index_fault += 1
                
                
            index_scenario = 0 
            for area,length in zip(self.scenario_area,self.scenario_length) : 
                area = area/1000000.
                Mmax_i = np.log10(area)+2./3.*np.log10(max(1.,np.sqrt(area/width**2.))/((1.+max(1.,area/(width**2.*7.4)))/2.))+3.98
                Mmax_i = Mmax_i + d_i * 0.2
                self.Mmax_scenario.append(round(Mmax_i,1))
                index_scenario += 1
                
                
                
                
                
                