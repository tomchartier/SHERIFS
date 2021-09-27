# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: Thomas Chartier
"""

#!/usr/bin/python
#
#Selection et calcul des FAULT
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton
from tkinter.ttk import Combobox

class selecFAULT_tk():
    def __init__(self,Run_Name,Model_name,File_geom,File_faults_n_scenarios,scenario_set):
        self.Run_Name = Run_Name
        self.Model_name = Model_name
        self.File_geom = File_geom
        self.File_faults_n_scenarios = File_faults_n_scenarios
        self.scenario_set = scenario_set
        
        self.initialize()

    def initialize(self):
        self.IG = tk.Tk()
        self.IG.title('For ' + str(self.Model_name) + ' and ' + str(self.scenario_set))
        self.IG.grid()
        self.FaultGeometry()
        self.FAULTSelect = StringVar()
        self.choixFAULT = list(set(self.Column_Fault_name))
        self.listeFAULT = Combobox(self.IG, textvariable=self.FAULTSelect, values=self.choixFAULT,state = 'readonly')
        self.listeFAULT.grid(column=0,row=0,sticky='EW')
        self.listeFAULT.current(0)
        
        add_fault_button = Button(self.IG,text=u"Add fault to model",command=self.Add_fault_ButtonClick)
        add_fault_button.grid(column=1,row=0)
        
        add_all_faults_button = Button(self.IG,text=u"Add all faults to model",command=self.Add_all_faults_ButtonClick)
        add_all_faults_button.grid(column=4,row=0)
        
        suppr_fault_button = Button(self.IG,text=u"Delete fault from model",command=self.Delete_fault_ButtonClick)
        suppr_fault_button.grid(column=2,row=0)

        add_fault_to_scenario_button = Button(self.IG,text=u"Add fault to scenario",command=self.Add_fault_to_scenario_ButtonClick)
        add_fault_to_scenario_button.grid(column=2,row=1)
        
        suppr_fault_from_scenario_button = Button(self.IG,text=u"Delete fault from scenario",command=self.Delete_fault_from_scenario_ButtonClick)
        suppr_fault_from_scenario_button.grid(column=2,row=2)
        
        add_scenario_button = Button(self.IG,text=u"Add scenario to model",command=self.Add_scenario_ButtonClick)
        add_scenario_button.grid(column=6,row=1)
        suppr_scenario_button = Button(self.IG,text=u"Delete scenario from model",command=self.Delete_scenario_ButtonClick)
        suppr_scenario_button.grid(column=6,row=2)
  
   
        calcul_button = Button(self.IG,text=u"Build Hazard Models",command=self.CalculButtonClick)
        calcul_button.grid(column=14,row=10)
        self.ouverture_calcul = 0

        self.listechoix_fault = Listbox(self.IG)
        self.listechoix_fault.grid(column=0,row=1,columnspan=2,rowspan=3,sticky='EW')
        
        self.listechoix_scenario_tmp = Listbox(self.IG)
        self.listechoix_scenario_tmp.grid(column=4,row=1,columnspan=2,rowspan=3,sticky='EW')
        
        self.listechoix_scenario = Listbox(self.IG,width = 50)
        self.listechoix_scenario.grid(column=7,row=1,columnspan=4,rowspan=4,sticky='EW')
  
        self.IG.grid_columnconfigure(0,weight=1)
        self.IG.resizable(True,True)
        
        self.IG.mainloop()
        
    def Add_fault_ButtonClick(self):
        longueur_liste = self.listechoix_fault.size()
        compteur = 0
        for i in range(longueur_liste):
            if self.listeFAULT.get() == self.listechoix_fault.get(i):
               compteur = compteur + 1
        if compteur == 0:
            self.listechoix_fault.insert(END, self.listeFAULT.get())
        else:
            messagebox.showerror('Error','FAULT already selected') 
            
    def Add_all_faults_ButtonClick(self):
        longueur_liste = self.listechoix_fault.size()
        compteur = 0
        for i in range(longueur_liste):
            if self.listeFAULT.get() == self.listechoix_fault.get(i):
               compteur = compteur + 1
        if compteur == 0:
            for i in range(len(self.choixFAULT)):
                self.listechoix_fault.insert(END, self.choixFAULT[i])
        else:
            messagebox.showerror('Error','FAULT already selected') 
            
    def Add_fault_to_scenario_ButtonClick(self):
        items = self.listechoix_fault.curselection()
        longueur_liste = self.listechoix_scenario_tmp.size()
        compteur = 0
        for i in range(longueur_liste):
            if self.listechoix_fault.get(items) == self.listechoix_scenario_tmp.get(i):
               compteur = compteur + 1
        if compteur == 0:
            self.listechoix_scenario_tmp.insert(END, self.listechoix_fault.get(items))
        else:
            messagebox.showerror('Error','FAULT already selected') 
            
    def Add_scenario_ButtonClick(self):
        longueur_liste = self.listechoix_scenario.size()
        longueur_liste_tmp = self.listechoix_scenario_tmp.size()
        scenario = ''
        for i in range(longueur_liste_tmp):
            if scenario == '':
                scenario += self.listechoix_scenario_tmp.get(i)
            else:
                scenario += ' ' + self.listechoix_scenario_tmp.get(i)    
        compteur = 0
        for i in range(longueur_liste):
            if scenario == self.listechoix_scenario.get(i):
               compteur = compteur + 1
        if compteur == 0:
            self.listechoix_scenario.insert(END, scenario)
        else:
            messagebox.showerror('Error','Scenario already selected')               
        pos = 0
        for i in range(longueur_liste_tmp):    
            idx = int(i) - pos
            self.listechoix_scenario_tmp.delete(idx,idx)
            pos = pos + 1
    
    def Delete_fault_ButtonClick(self):
        items = self.listechoix_fault.curselection()
        pos=0
        for i in items:
            idx = int(i) - pos
            self.listechoix_fault.delete(idx,idx)
            pos = pos + 1
    def Delete_fault_from_scenario_ButtonClick(self):
        items = self.listechoix_scenario_tmp.curselection()
        pos=0
        for i in items:
            idx = int(i) - pos
            self.listechoix_scenario_tmp.delete(idx,idx)
            pos = pos + 1
    def Delete_scenario_ButtonClick(self):
        items = self.listechoix_scenario.curselection()
        pos=0
        for i in items:
            idx = int(i) - pos
            self.listechoix_scenario.delete(idx,idx)
            pos = pos + 1  
  
    def CalculButtonClick(self):
        #WaitWindow(self)
        if self.ouverture_calcul == 0:
            faults_n_scenar=open(self.File_faults_n_scenarios,'w')
            
            longueur_liste_faults = self.listechoix_fault.size()
            line_faults = ''
            for i in range(longueur_liste_faults):
                if line_faults == '':
                    line_faults += self.listechoix_fault.get(i)
                else:
                    line_faults += ' ' + self.listechoix_fault.get(i)
            faults_n_scenar.write(line_faults)
            
            longueur_liste_scenario = self.listechoix_scenario.size()
            for i in range(longueur_liste_scenario):
                line_scenario = '\n'+ self.listechoix_scenario.get(i)
                faults_n_scenar.write(line_scenario )
            faults_n_scenar.close()
            #WaitWindow(self) #Finn dances
            self.ouverture_calcul = 1 
            self.IG.destroy()
                
    def FaultGeometry(self):
        NomFichier_InfosZonage = self.File_geom
        InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8')],skip_header = 1)
        Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
        index_model = np.where(np.array(Column_model_name) == self.Model_name)
        self.Column_Fault_name = list(map(lambda i : InfosZonage[i][1],index_model[0]))


  
