#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems


Tool for weighting the Logic Tree.
To be used after each branch has been visualized.

@author: Thomas Chartier
"""


import numpy as np
#from Tkinter import Tk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton
#from ttk import Combobox
#from tkMessageBox import showerror
import tkinter as tk
from tkinter import ttk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton
from tkinter.ttk import Combobox
from tkinter import messagebox
import os
import sys
path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/lib'
sys.path.append(path_lib)
from file_writer.OQ_job_Creator import OQ_job_Creator


class Weight_LT():
    def __init__(self,):
        
        self.initialize()
    def initialize(self):
        Run_Name = 'Example'
        
        LT_log_name  =  str(Run_Name)+'/LT_log.txt' 

        file_log_LT = open(LT_log_name,'r')
        self.log_LT = file_log_LT.readlines()
        
        #print self.log_LT[1].split('\t')
        #print len(self.log_LT[1].split('\t'))
        
        OQ_job = OQ_job_Creator(Run_Name) # ask the info about the run and create the job.ini file
        
        nb_random_sampling = OQ_job.nb_sample               
        
        self.selected_Model = self.log_LT[1].split('\t')
        if '\n' in self.selected_Model:
            self.selected_Model.remove('\n')
        if '\r\n' in self.selected_Model:
            self.selected_Model.remove('\r\n')
            
        self.selected_ScL = self.log_LT[3].split('\t')
        if '\n' in self.selected_ScL:
            self.selected_ScL.remove('\n')
        if '\r\n' in self.selected_ScL:
            self.selected_ScL.remove('\r\n')
        
            
        index_advance = 0
        self.mfd_hyps = []
        self.b_values_hyps = []
        while self.log_LT[5+index_advance][0:3] == 'MFD':
            self.mfd_hyps.append(self.log_LT[5+index_advance].split('\t')[0])
            b_values_hyps_i = []
            for b_hyp in self.log_LT[5+index_advance].split('\t')[1:]:
                b_values_hyps_i.append(b_hyp)
            if '\n' in b_values_hyps_i:
                b_values_hyps_i.remove('\n')
            if '\r\n' in b_values_hyps_i:
                b_values_hyps_i.remove('\r\n')
            self.b_values_hyps.append(b_values_hyps_i)
            index_advance += 1
            
            
        self.bg_names = self.log_LT[6+index_advance].split('\t')
        if '\n' in self.bg_names:
            self.bg_names.remove('\n')
        if '\r\n' in self.bg_names:
            self.bg_names.remove('\r\n')
            
        self.sc_names = self.log_LT[8+index_advance].split('\t')
        if '\n' in self.sc_names:
            self.sc_names.remove('\n')
        if '\r\n' in self.sc_names:
            self.sc_names.remove('\r\n')
        if '' in self.sc_names:
            self.sc_names.remove('')
            
            
        self.window_weight()
            
        
        branches = []
        for model_i in self.selected_Model: 
            index_mfd = 0
            for mfd_i in self.mfd_hyps:               
                for bvalue in self.b_values_hyps[index_mfd]:            
                    for bg_hyp_i in self.bg_names:             
                        for sc_name in self.sc_names:  
                            for ScL_i in self.selected_ScL :  
                                ScL_i = ScL_i.split(' ')
                                ScL_name_i= ScL_i[0]
                                use_all_i = ScL_i[2][0]
                                dim_i = ScL_i[1]
                
                                branch_i = [model_i,ScL_name_i,use_all_i,dim_i,mfd_i,bvalue,bg_hyp_i,sc_name]
                                branches.append(branch_i)
                index_mfd += 1
                        

        LT_file = str(Run_Name)+'/Sources_Logic_tree.xml'
        XMLfile=open(LT_file,'w')
        Ligne='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
        XMLfile.write(Ligne)
        Ligne='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
        XMLfile.write(Ligne)
        Ligne='\txmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        XMLfile.write(Ligne)
        Ligne='\t<logicTree logicTreeID="lt1">\n'
        XMLfile.write(Ligne)
                 
        str_all_data = []
        id_number = 1
#        scenario_done = []
#        scenario_path = []
        
        Ligne='\t\t<logicTreeBranchingLevel branchingLevelID="bl' + str(id_number) + '">\n'
        XMLfile.write(Ligne)
        Ligne='\t\t\t<logicTreeBranchSet uncertaintyType="sourceModel"\n'
        XMLfile.write(Ligne)
        Ligne='\t\t\t\t\t\t\tbranchSetID="bs' + str(id_number) + '">\n'
        XMLfile.write(Ligne) 
        check_weights = []
        for branch in branches :
            Model = branch[0]
            selected_ScL = branch[1]
            dim_used = branch[3][0]

#            if branch[2] == True :
#                str_all_data = 'a' # ' a ' is for 'all data is used'
#            else :
#                str_all_data = 'm' # ' m ' is for 'mechanic specific data only'
            str_all_data = branch[2]
                
            bvalue = branch[5]
            mfd_hyp = str(branch[4])
            BG_hyp = branch[6]
            scenario_set = branch[7]
            for sample in range(1,nb_random_sampling+1):
                

                Ligne=('\t\t\t\t<logicTreeBranch branchID="' + str(Model) + '-' + str(BG_hyp) + '-' + str(selected_ScL) + '-' 
                       + str(dim_used) + '-' + str_all_data + '-' +  str(scenario_set) + '-' 
                        + str(bvalue)+ '-' + str(mfd_hyp)  + '-s_' + str(sample) + '">\n')
                XMLfile.write(Ligne)
                Ligne=('\t\t\t\t\t<uncertaintyModel>' + (str(Model) + '/' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                       + str(dim_used) + '_' + str_all_data + '/' +  str(scenario_set) + '/' 
                        + str(bvalue)+ '/' + str(mfd_hyp)) + '/Source_model_' 
                + str(sample) + '.xml</uncertaintyModel>\n')
                XMLfile.write(Ligne)
                
                index = np.where(np.array(self.selected_Model)==Model)[0][0]
                weight_model = float(self.weight_model[index])
                
                index = np.where(np.array(self.sc_names)==scenario_set)[0][0]
                weight_sc = float(self.weight_sc[index])
                
                index = np.where(np.array(self.selected_ScL)==branch[1]+' '+branch[3]+' '+branch[2])[0][0]
                weight_sl = float(self.weight_sl[index])
                
                index = np.where(np.array(self.bg_names)==BG_hyp)[0][0]
                weight_bg = float(self.weight_bg[index])
                
                index_mfd = np.where(np.array(self.mfd_hyps)==mfd_hyp)[0][0]
                weight_mfd = float(self.weight_mfd[index_mfd])
                
                index_bvalue = np.where(np.array(self.b_values_hyps[index_mfd])==bvalue)[0][0]
                weight_bv = float(self.weight_b_value[index_mfd][index_bvalue])
                
                weight_branch = weight_model * weight_bg * weight_sc * weight_sl * weight_mfd * weight_bv / float(nb_random_sampling)
                #print weight_branch
                Ligne='\t\t\t\t\t<uncertaintyWeight>' + str(round(float(weight_branch),6)) + '</uncertaintyWeight>\n'
                check_weights.append(round(float(weight_branch),6))

                XMLfile.write(Ligne)  
                print(str(Model) + '/' + str(BG_hyp) + '/' + str(selected_ScL) + '_' 
                       + str(dim_used) + '_' + str_all_data + '/' +  str(scenario_set) + '/' 
                        + str(bvalue)+ '/' + str(mfd_hyp)+ '   ' + str(round(float(weight_branch),7)))
                Ligne='\t\t\t\t</logicTreeBranch>\n'
                XMLfile.write(Ligne)   
                
                id_number += 1

        Ligne='\t\t\t</logicTreeBranchSet>\n'
        XMLfile.write(Ligne)            
        Ligne='\t\t</logicTreeBranchingLevel>\n'
        XMLfile.write(Ligne)            
        Ligne='\t</logicTree>\n'
        XMLfile.write(Ligne)
        Ligne='</nrml>\n'
        XMLfile.write(Ligne)
        
        XMLfile.close()
        
        if sum(check_weights) != 1.0 :
            print( '\n!!!!!!!!!!!!\nSum of weights is ', sum(check_weights) ,' \n!!!!!!!!!!!!\n' )

        
        
    def window_weight(self):
        self.weight = {}
        self.w = tk.Tk()
        self.w.grid()
        row_i = 1
        self.w.title('Please set the weights for the branches')
        
        label = Label(self.w, text='\nModels')
        label.grid(column=0,row=row_i)
        row_i +=1
        array_hypo = self.selected_Model
        for i in range(len(array_hypo)):
            label = Label(self.w, text="Model "+str(array_hypo[i]))
            label.grid(column=0,row=row_i)
            self.weight["w_model_{0}".format(i)] = Entry(self.w)
            self.weight["w_model_{0}".format(i)].grid(column=6,row= row_i) 
            self.weight["w_model_{0}".format(i)].insert(INSERT,str(round(1./float(len(array_hypo)),3)))
            row_i +=1
        
        label = Label(self.w, text='\nBackground hypotheses')
        label.grid(column=0,row=row_i)
        row_i +=1
        array_hypo = self.bg_names
        for i in range(len(array_hypo)):
            label = Label(self.w, text="Model "+str(array_hypo[i]))
            label.grid(column=0,row=row_i)
            self.weight["w_bg_{0}".format(i)] = Entry(self.w)
            self.weight["w_bg_{0}".format(i)].grid(column=6,row= row_i) 
            self.weight["w_bg_{0}".format(i)].insert(INSERT,str(round(1./float(len(array_hypo)),3)))
            row_i +=1
        
        label = Label(self.w, text='\nScalling Laws')
        label.grid(column=0,row=row_i)
        row_i +=1
        array_hypo = self.selected_ScL
        for i in range(len(array_hypo)):
            label = Label(self.w, text="Model "+str(array_hypo[i]))
            label.grid(column=0,row=row_i)
            self.weight["w_sl_{0}".format(i)] = Entry(self.w)
            self.weight["w_sl_{0}".format(i)].grid(column=6,row= row_i) 
            self.weight["w_sl_{0}".format(i)].insert(INSERT,str(round(1./float(len(array_hypo)),3)))
            row_i +=1
        
        label = Label(self.w, text='\nScenarios sets')
        label.grid(column=0,row=row_i)
        row_i +=1
        array_hypo = self.sc_names
        for i in range(len(array_hypo)):
            label = Label(self.w, text="Model "+str(array_hypo[i]))
            label.grid(column=0,row=row_i)
            self.weight["w_sc_{0}".format(i)] = Entry(self.w)
            self.weight["w_sc_{0}".format(i)].grid(column=6,row= row_i) 
            self.weight["w_sc_{0}".format(i)].insert(INSERT,str(round(1./float(len(array_hypo)),3)))
            row_i +=1
        
        label = Label(self.w, text='\nMDFs hypothesis')
        label.grid(column=0,row=row_i)
        row_i +=1
        array_hypo = self.mfd_hyps
        for i in range(len(array_hypo)):
            label = Label(self.w, text="Model "+str(array_hypo[i]))
            label.grid(column=0,row=row_i)
            self.weight["w_mfd_{0}".format(i)] = Entry(self.w)
            self.weight["w_mfd_{0}".format(i)].grid(column=6,row= row_i) 
            self.weight["w_mfd_{0}".format(i)].insert(INSERT,str(round(1./float(len(array_hypo)),3)))
            row_i +=1
            
        label = Label(self.w, text='\n b Values hypothesis')
        label.grid(column=0,row=row_i)
        row_i +=1
        array_hypo = self.b_values_hyps
        index_mfd = 0
        ii= 0
        for i in range(len(self.mfd_hyps)):
            for j in range(len(array_hypo[i])):
                label = Label(self.w, text="Model "+str(array_hypo[i][j]))
                label.grid(column=0,row=row_i)
                self.weight["w_bv_{0}".format(ii)] = Entry(self.w)
                self.weight["w_bv_{0}".format(ii)].grid(column=6,row= row_i) 
                self.weight["w_bv_{0}".format(ii)].insert(INSERT,str(round(1./float(len(array_hypo[i])),3)))
                row_i +=1
                ii+=1
            index_mfd += 1
            
        bou_ok = Button(self.w, text=u'OK', command = self.OK)
        bou_ok.grid(column=4,row=row_i+1)
                
        self.w.mainloop()
        
    def OK(self):
        
        self.weight_model = []
        array_hypo = self.selected_Model
        for i in range(len(array_hypo)):
            self.weight_model.append(self.weight["w_model_{0}".format(i)].get())     
            
        self.weight_bg = []
        array_hypo = self.bg_names
        for i in range(len(array_hypo)):
            self.weight_bg.append(self.weight["w_bg_{0}".format(i)].get())   
            
        self.weight_sl = []
        array_hypo = self.selected_ScL
        for i in range(len(array_hypo)):
            self.weight_sl.append(self.weight["w_sl_{0}".format(i)].get())    
            
        self.weight_sc = []
        array_hypo = self.sc_names
        for i in range(len(array_hypo)):
            self.weight_sc.append(self.weight["w_sc_{0}".format(i)].get())    
            
        self.weight_mfd = []
        array_hypo = self.mfd_hyps
        for i in range(len(array_hypo)):
            self.weight_mfd.append(self.weight["w_mfd_{0}".format(i)].get())    
            
        self.weight_b_value = []
        array_hypo = self.b_values_hyps
        ii=0
        for i in range(len(self.mfd_hyps)):
            weight_b_value_i = []
            for j in range(len(array_hypo[i])):
                weight_b_value_i.append(self.weight["w_bv_{0}".format(ii)].get())  
                ii+=1
            self.weight_b_value.append(weight_b_value_i)
        print( self.weight_b_value)
            
            
        self.w.destroy() 
        
if __name__=="__main__":
    app = Weight_LT()
