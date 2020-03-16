# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: Thomas Chartier
"""
import tkinter as tk
from tkinter import ttk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton
from tkinter.ttk import Combobox
from tkinter import messagebox

class GMPE_lt():
    def __init__(self,available_GMPE,Domain_in_model,Run_Name):
        self.available_GMPE = available_GMPE
        self.Domain_in_model = Domain_in_model
        self.Run_Name = Run_Name
        
        self.initialize()
    def initialize(self):
        
        self.XMLfile=open(str(self.Run_Name)+'/GMPE_Logic_tree.xml','w')
        Ligne='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n\n'
        self.XMLfile.write(Ligne)
        Ligne='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
        self.XMLfile.write(Ligne)
        Ligne='\txmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        self.XMLfile.write(Ligne)
        Ligne='\t<logicTree logicTreeID="lt_gmpe">\n\n'
        self.XMLfile.write(Ligne)
        
        self.i_Domain = 0
        
        while self.i_Domain < len(self.Domain_in_model):
            self.Domain = self.Domain_in_model[self.i_Domain]
            self.weight = []
            self.fen1 = tk.Tk()
            self.GMPESelect = StringVar()
            Tex_box = Text(self.fen1, height = 2, width = 50)
            Tex_box.insert(INSERT,"Apply to  : " + str(self.Domain))
            Tex_box.pack()
            self.Box_GMPE = Combobox(self.fen1, textvariable=self.GMPESelect, values=self.available_GMPE,state = 'readonly', height = 5, width = 50)
            self.Box_GMPE.pack()
            bou_add = Button(self.fen1, text=u'Add GMPE', command = self.Add_gmpe)
            bou_add.pack()
            self.listechoix_gmpe = Listbox(self.fen1, width = 50) 
            self.listechoix_gmpe.pack()
            if not (self.i_Domain + 1) == len(self.Domain_in_model):             
                bou_quit = Button(self.fen1, text='Next Domain', command = self.Next_Domain) 
            else:
                bou_quit = Button(self.fen1, text='Build GMPE Logic Tree', command = self.Finalize_LT)
            bou_quit.pack()
            self.fen1.mainloop()
            self.i_Domain += 1
      
    def Add_gmpe(self):
        longueur_liste = self.listechoix_gmpe.size()
        compteur = 0
        for i in range(longueur_liste):
            if self.Box_GMPE.get() == self.listechoix_gmpe.get(i):
               compteur = compteur + 1
        if compteur == 0:
            self.listechoix_gmpe.insert(END, self.Box_GMPE.get())
            self.fen2 = tk.Tk()
            tex1 = Label(self.fen2, text='Associated weight:')
            tex1.pack()
            self.e = Entry(self.fen2)
            self.e.pack()
            bou_quit = Button(self.fen2, text='Apply weigth', command = self.Apply_weight)
            bou_quit.pack()
            self.fen2.mainloop()
        else:
            messagebox.showerror('Error','GMPE already selected') 
    
    def Apply_weight(self):
        self.weight.append(float(self.e.get()))
        self.fen2.destroy()
        
        
    def Next_Domain(self):
        if not sum(self.weight) == 1 :
            self.wind_err = tk.Tk()
            Text_box_err = Text(self.wind_err, height = 1, width = 80)
            Text_box_err.insert(INSERT,'Sum of weigths must be equal to one! Please enter GMPEs again.')
            Text_box_err.pack()
            bou_quit_box_err = Button(self.wind_err, text='OK', command = self.Delete_gmpe)
            bou_quit_box_err.pack()
            self.wind_err.mainloop()
        else :
            Ligne='\t\t<logicTreeBranchingLevel branchingLevelID="gmpe_bl' + str(self.i_Domain + 1 ) + '">\n'
            self.XMLfile.write(Ligne)        
            Ligne='\t\t\t<logicTreeBranchSet uncertaintyType="gmpeModel" branchSetID="' + str(self.Domain) + '"\n'
            self.XMLfile.write(Ligne)
            Ligne='\t\t\t\t\tapplyToTectonicRegionType="'+ str(self.Domain) + '">\n\n'
            self.XMLfile.write(Ligne)
            for i in range(self.listechoix_gmpe.size()):
                gmpe = self.listechoix_gmpe.get(i)
                Ligne='\t\t\t\t<logicTreeBranch branchID="'+ str(gmpe) + '">\n'
                self.XMLfile.write(Ligne)
                Ligne='\t\t\t\t\t<uncertaintyModel>'+ str(gmpe) + '</uncertaintyModel>\n'
                self.XMLfile.write(Ligne)
                Ligne='\t\t\t\t\t<uncertaintyWeight>'+ str(self.weight[i]) + '</uncertaintyWeight>\n'
                self.XMLfile.write(Ligne)
                Ligne='\t\t\t\t</logicTreeBranch>\n\n'
                self.XMLfile.write(Ligne)
                
            Ligne='\t\t\t</logicTreeBranchSet>\n'
            self.XMLfile.write(Ligne)
            Ligne='\t\t</logicTreeBranchingLevel>\n\n'
            self.XMLfile.write(Ligne)        
            self.fen1.destroy()
        
    def Finalize_LT(self):
        if not sum(self.weight) == 1 :
            self.wind_err = tk.Tk()
            Text_box_err = Text(self.wind_err, height = 1, width = 80)
            Text_box_err.insert(INSERT,'Sum of weigths must be equal to one! Please enter GMPEs again.')
            Text_box_err.pack()
            bou_quit_box_err = Button(self.wind_err, text='OK', command = self.Delete_gmpe)
            bou_quit_box_err.pack()
            self.wind_err.mainloop()
        else :        
            Ligne='\t\t<logicTreeBranchingLevel branchingLevelID="gmpe_bl' + str(self.i_Domain + 1) + '">\n'
            self.XMLfile.write(Ligne)        
            Ligne='\t\t\t<logicTreeBranchSet uncertaintyType="gmpeModel" branchSetID="' + str(self.Domain) + '"\n'
            self.XMLfile.write(Ligne)
            Ligne='\t\t\t\t\tapplyToTectonicRegionType="'+ str(self.Domain) + '">\n\n'
            self.XMLfile.write(Ligne)
            for i in range(self.listechoix_gmpe.size()):
                gmpe = self.listechoix_gmpe.get(i)
                Ligne='\t\t\t\t<logicTreeBranch branchID="'+ str(gmpe) + '">\n'
                self.XMLfile.write(Ligne)
                Ligne='\t\t\t\t\t<uncertaintyModel>'+ str(gmpe) + '</uncertaintyModel>\n'
                self.XMLfile.write(Ligne)
                Ligne='\t\t\t\t\t<uncertaintyWeight>'+ str(self.weight[i]) + '</uncertaintyWeight>\n'
                self.XMLfile.write(Ligne)
                Ligne='\t\t\t\t</logicTreeBranch>\n\n'
                self.XMLfile.write(Ligne)
                
            Ligne='\t\t\t</logicTreeBranchSet>\n'
            self.XMLfile.write(Ligne)
            Ligne='\t\t</logicTreeBranchingLevel>\n\n'
            self.XMLfile.write(Ligne)
            
            Ligne='\t</logicTree>\n'
            self.XMLfile.write(Ligne)
            Ligne='</nrml>\n'
            self.XMLfile.write(Ligne)                
            self.XMLfile.close()        
            self.fen1.destroy()
            
        
    def Delete_gmpe(self):
        for i in range(self.listechoix_gmpe.size()):
            self.listechoix_gmpe.delete(i)
        self.wind_err.destroy()
        self.i_Domain -= 1        
        self.fen1.destroy()
