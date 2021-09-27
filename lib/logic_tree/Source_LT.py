# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

This code open the interface to select the options explored in the logic tree

@author: Thomas Chartier
"""
import numpy as np

#import tkinter as tk
#from tkinter import ttk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton
#from tkinter.ttk import Combobox
#from tkinter import messagebox


class S_LT():
    def __init__(self,File_geom,File_prop,Run_Name):        
        self.File_geom = File_geom 
        self.File_prop = File_prop 
        self.Run_Name = Run_Name
        
        self.initialize()
    def initialize(self):
        
        self.get_available_scaling_laws()
        
        self.windows_nd_ScL()
        
        # structure of branch : [model,ScL,use_all_data,dimentio_used,µ,bmin_bmax,bg_hyp]
        self.branches = []
        
        
        for self.model_i in self.selected_Model:
            self.nb_mfd_hyp()
            
            self.bmin_hyp = []
            self.bmax_hyp = []
            for self.i in range(self.nb_of_mfd):
                
                self.nb_b_hyp()
                self.bmin_hyp.append(self.bmin_hyp_i)
                self.bmax_hyp.append(self.bmax_hyp_i)
                
            self.nb_sc_hyp()
            self.nb_bg_hyp()
            for ScL_i,use_all_i,dim_i in zip(self.selected_ScL,self.use_all_ScL_data,self.dimention_used) : 
                index_mfd = 0                 
                for mfd_i in self.mfd_hyp:     
                    for bmin_i,bmax_i in zip(self.bmin_hyp[index_mfd],self.bmax_hyp[index_mfd]):
                        for bg_hyp_i in self.bg_names: 
                            for sc_name in self.sc_names:
                                branch_i = [self.model_i,ScL_i,use_all_i,dim_i,mfd_i,str(bmin_i)+'_'+str(bmax_i),bg_hyp_i,sc_name]
                                self.branches.append(branch_i)
                    index_mfd += 1
         
        # file containing the branches names in the logic tree    
        LT_log_name  =  str(self.Run_Name)+'/LT_log.txt' 

        LT_log = open(LT_log_name, 'w')

        LT_log.write('Models\n')
        for self.model_i in self.selected_Model:
            LT_log.write(self.model_i+'\t')
        
        LT_log.write('\nSaling Laws\n')
        for ScL_i,use_all_i,dim_i in zip(self.selected_ScL,self.use_all_ScL_data,self.dimention_used) :         
            if use_all_i == True :
                str_all_data = 'a' # ' a ' is for 'all data is used'
            else :
                str_all_data = 'm' # ' m ' is for 'mechanic specific data only'
            LT_log.write(ScL_i+' '+dim_i+' '+str_all_data+'\t')
        
        LT_log.write('\nMFD\tb value\n')  
        index_mfd = 0
        for mfd_i in self.mfd_hyp:   
            LT_log.write('MFD_'+mfd_i+'\t')
            for bmin_i,bmax_i in zip(self.bmin_hyp[index_mfd],self.bmax_hyp[index_mfd]):
                LT_log.write('bmin_'+str(bmin_i)+'_bmax_'+bmax_i+'\t')   
            LT_log.write('\n')
            index_mfd += 1
            
            
        LT_log.write('Background\n') 
        for bg_hyp_i in self.bg_names: 
            LT_log.write('bg_'+bg_hyp_i+'\t')
            
        LT_log.write('\nscenario set\n') 
        for sc_name in self.sc_names:
            LT_log.write('sc_'+sc_name+'\t')

        LT_log.close()
            
    '''#####################'''        
    def windows_nd_ScL(self):
        self.w_ScL_nb = tk.Tk()
        self.w_ScL_nb.title('Number of scaling laws')
        label = Label(self.w_ScL_nb, text="\nHow many Scaling laws do you want to use?")
        label.pack()
        self.nb_of_scl = Entry(self.w_ScL_nb)
        self.nb_of_scl.pack()
        self.nb_of_scl.insert(INSERT,1)
        bou_ok_ScL = Button(self.w_ScL_nb, text=u'OK', command = self.OK_nb_Scl)
        bou_ok_ScL.pack()
        
        self.w_ScL_nb.mainloop()
        
    def OK_nb_Scl(self) :
        self.nb_of_scl = int(self.nb_of_scl.get())
        self.w_ScL_nb.destroy() 
        self.windows_ScL()
        
    def windows_ScL(self):        
        self.var = {}
        self.dimention_used_box = {}
        self.ScLSelect_dim_used = {}
        self.ScLname_used = {}
        self.ScL_name = {}
        self.w_ScL = tk.Tk()
        self.w_ScL.title('Scaling laws')
        self.w_ScL.grid()
        row_i = 1
        for i in range(self.nb_of_scl):   
            self.ScLname_used["ScLname{0}".format(i)] = StringVar()     
            self.ScL_name["ScLname{0}".format(i)] = Combobox(self.w_ScL, textvariable=self.ScLname_used["ScLname{0}".format(i)], values=self.available_scaling_laws,state = 'readonly', height = 5, width = 30)
            self.ScL_name["ScLname{0}".format(i)].current(0)
            self.ScL_name["ScLname{0}".format(i)].grid(column=0,row=row_i)
        
            self.var["use_all_ScL_data_{0}".format(i)] = StringVar()
            self.use_all_ScL_data_button = Checkbutton(self.w_ScL, text="Respect fault kinematic", variable=self.var["use_all_ScL_data_{0}".format(i)],onvalue="False", offvalue="True")
            self.use_all_ScL_data_button.grid(column=4,row= row_i)
            
            self.ScLSelect_dim_used["dimention_used_{0}".format(i)] = StringVar()
            self.dimention_used_box["dimention_used_{0}".format(i)] = Combobox(self.w_ScL, textvariable=self.ScLSelect_dim_used["dimention_used_{0}".format(i)], values=['Area','Length'],state = 'readonly', height = 5, width = 30)
            self.dimention_used_box["dimention_used_{0}".format(i)].current(0)
            self.dimention_used_box["dimention_used_{0}".format(i)].grid(column=8,row= row_i)
            row_i += 1
        
        bou_ok_ScL = Button(self.w_ScL, text=u'OK', command = self.OK_Scl)
        bou_ok_ScL.grid(column=8,row= row_i )
        
        self.w_ScL.mainloop()
        
    def OK_Scl(self) :
        self.dimention_used = []
        self.use_all_ScL_data = []
        self.selected_ScL = []
        for i in range(self.nb_of_scl):
            self.selected_ScL.append(self.ScL_name["ScLname{0}".format(i)].get())
            if self.var["use_all_ScL_data_{0}".format(i)].get() == 'False':                
                self.use_all_ScL_data.append(False)
            else :
                self.use_all_ScL_data.append(True)
            self.dimention_used.append(self.dimention_used_box["dimention_used_{0}".format(i)].get())
        
        self.w_ScL.destroy() 
        
        check_ScL = []
        for selected_ScL_i,use_all_ScL_data,dimention_used in zip(self.selected_ScL,self.use_all_ScL_data,self.dimention_used):
            line = str(selected_ScL_i)+str(use_all_ScL_data)+str(dimention_used)
            if not line in check_ScL:
                check_ScL.append(line)
            else:
                messagebox.showerror('Error','One scaling law has been selected twice\n please start again')   
        self.win_model()
    
    '''#####################''' 
    def win_model(self):
        self.get_available_models()
                    
        self.w_model = tk.Tk()
        self.ModelSelect = StringVar()
        
        label = Label(self.w_model, text="\nWhich model do you want to use?")
        label.pack()
        self.Box_Model = Combobox(self.w_model,  values=self.available_models,state = 'readonly', height = 5, width = 30)
        self.Box_Model.current(0)
        self.Box_Model.pack()
        bou_add_Model = Button(self.w_model, text=u'Add Model', command = self.Add_Model)
        bou_add_Model.pack()
        self.list_Model = Listbox(self.w_model, height = 5, width = 30) 
        self.list_Model.pack() 
        bou_del_Model = Button(self.w_model, text=u'Delete Selected Model', command = self.del_Model)
        bou_del_Model.pack() 
        label = Label(self.w_model, text="\n\n\nWhen ready click \"OK\"")
        label.pack()
        bou_ok_Model = Button(self.w_model, text=u'OK', command = self.ok_Model)
        bou_ok_Model.pack()
        self.w_model.mainloop()
        
        
    def Add_Model(self):
        len_list = self.list_Model.size()
        compteur = 0
        for i in range(len_list):
            if self.Box_Model.get() == self.list_Model.get(i):
               compteur = compteur + 1
        if compteur == 0:
            self.list_Model.insert(END, self.Box_Model.get())
        else:
            messagebox.showerror('Error','Model already selected') 
            
            
    def del_Model(self):
        items = self.list_Model.curselection()
        pos=0
        for i in items:
            idx = int(i) - pos
            self.list_Model.delete(idx,idx)
            pos = pos + 1
            
    def ok_Model(self):
        self.selected_Model = []
        longueur_liste = self.list_Model.size()
        for i in range(longueur_liste):
            if len(self.list_Model.get(i))!=0:
                self.selected_Model.append(self.list_Model.get(i))
        if len(self.selected_Model)==0:
            messagebox.showerror('Error','Select at least one model')
        self.w_model.destroy() 
        
    
    '''#####################''' 
    def nb_mfd_hyp(self):        
        self.w_mfd_nb = tk.Tk()
        self.w_mfd_nb.title('Number of MFD hypothesis for model : '+str(self.model_i))
        label = Label(self.w_mfd_nb, text='\nHow many MFD hypothesis for model : '+str(self.model_i)+' do you want to use?')
        label.pack()
        self.nb_of_mfd = Entry(self.w_mfd_nb)
        self.nb_of_mfd.pack()
        self.nb_of_mfd.insert(INSERT,1)
        bou_ok = Button(self.w_mfd_nb, text=u'OK', command = self.OK_nb_mfd_hyp)
        bou_ok.pack()
        
        self.w_mfd_nb.mainloop()
        
    def OK_nb_mfd_hyp(self) :
        self.nb_of_mfd = int(self.nb_of_mfd.get())
        self.w_mfd_nb.destroy() 
        self.mfd_hyp()
        
    def mfd_hyp(self):
        self.mfd = {}        
        self.w_mfd = tk.Tk()
        self.w_mfd.grid()
        self.w_mfd.title('Hypothesis on MFD for '+str(self.model_i))
        row_i = 1
        for i in range(self.nb_of_mfd):
            label = Label(self.w_mfd, text="Hypothesis "+str(i+1))
            label.grid(column=0,row=row_i)
            self.mfd["nb_mfd_{0}".format(i)] = Entry(self.w_mfd)
            self.mfd["nb_mfd_{0}".format(i)].grid(column=1,row= row_i) 
            if i == 0:
                self.mfd["nb_mfd_{0}".format(i)].insert(INSERT,'GR') 
            elif i == 1:
                self.mfd["nb_mfd_{0}".format(i)].insert(INSERT,'YC') 
            row_i +=1
        bou_ok = Button(self.w_mfd, text=u'OK', command = self.OK_mfd_hyp)
        bou_ok.grid(column=4,row=row_i+1)
        
        #self.w_shear_mod.mainloop()
        
        
    def OK_mfd_hyp(self) :
        self.mfd_hyp = []
        for i in range(self.nb_of_mfd):
            self.mfd_hyp.append(self.mfd["nb_mfd_{0}".format(i)].get())        
        self.w_mfd.destroy() 
        
    
    '''#####################''' 
    def nb_b_hyp(self):        
        self.w_b_nb = tk.Tk()
        self.w_b_nb.title('Number of b value distribution for model : '+str(self.model_i) +'\nand MFD : ' + str(self.mfd_hyp[self.i]))
        label = Label(self.w_b_nb, text='\nHow many b value distribution hypothesis for model : '+str(self.model_i)+'\nand MFD : ' + str(self.mfd_hyp[self.i])+' do you want to use?')
        label.pack()
        self.nb_of_b = Entry(self.w_b_nb)
        self.nb_of_b.pack()
        self.nb_of_b.insert(INSERT,1)
        bou_ok = Button(self.w_b_nb, text=u'OK', command = self.OK_nb_b_hyp)
        bou_ok.pack()
        
        self.w_b_nb.mainloop()
        
    def OK_nb_b_hyp(self) :
        self.nb_of_b = int(self.nb_of_b.get())
        self.w_b_nb.destroy() 
        self.b_hyp()
        
    def b_hyp(self):
        self.bmin = {}   
        self.bmax = {}       
        self.w_b = tk.Tk()
        self.w_b.grid()
        self.w_b.title('Hypothesis on b value')
        row_i = 0
        label = Label(self.w_b, text='\nHypothesis on b value for model '+str(self.model_i)+' and MFD : ' + str(self.mfd_hyp[self.i]))
        label.grid(column=0,row=row_i) 
        row_i +=1
        for i in range(self.nb_of_b):
            label = Label(self.w_b, text="Hypothesis "+str(i+1)+" for b min and b max")
            label.grid(column=0,row=row_i)
            self.bmin["nb_bmin_{0}".format(i)] = Entry(self.w_b)
            self.bmin["nb_bmin_{0}".format(i)].grid(column=4,row= row_i) 
            self.bmin["nb_bmin_{0}".format(i)].insert(INSERT,0.9)   
            self.bmax["nb_bmax_{0}".format(i)] = Entry(self.w_b)
            self.bmax["nb_bmax_{0}".format(i)].grid(column=6,row= row_i) 
            self.bmax["nb_bmax_{0}".format(i)].insert(INSERT,1.1)   
            row_i +=1
        bou_ok = Button(self.w_b, text=u'OK', command = self.OK_b_hyp)
        bou_ok.grid(column=4,row=row_i+1)
        
        self.w_b.mainloop()
        
    def OK_b_hyp(self) :
        self.bmin_hyp_i = []
        self.bmax_hyp_i = []
        for i in range(self.nb_of_b):
            self.bmin_hyp_i.append(self.bmin["nb_bmin_{0}".format(i)].get())  
            self.bmax_hyp_i.append(self.bmax["nb_bmax_{0}".format(i)].get())        
        self.w_b.destroy() 
        
    
    '''#####################''' 
    def nb_sc_hyp(self):        
        self.w_sc_nb = tk.Tk()
        self.w_sc_nb.title('Number of scenario sets for model : '+str(self.model_i))
        label = Label(self.w_sc_nb, text='\nHow many scenario sets for model : '+str(self.model_i)+' do you want to use?')
        label.pack()
        self.nb_of_sc = Entry(self.w_sc_nb)
        self.nb_of_sc.pack()
        self.nb_of_sc.insert(INSERT,1)
        bou_ok = Button(self.w_sc_nb, text=u'OK', command = self.OK_nb_sc_hyp)
        bou_ok.pack()
        
        self.w_sc_nb.mainloop()
        
    def OK_nb_sc_hyp(self) :
        self.nb_of_sc = int(self.nb_of_sc.get())
        self.w_sc_nb.destroy() 
        self.sc_hyp()
        
    def sc_hyp(self):
        self.sc = {}       
        self.w_sc = tk.Tk()
        self.w_sc.grid()
        self.w_sc.title('Hypothesis on scenario sets for '+str(self.model_i))
        row_i = 0
        for i in range(self.nb_of_sc):
            label = Label(self.w_sc, text="Hypothesis "+str(i+1))
            label.grid(column=0,row=row_i)
            self.sc["sc_{0}".format(i)] = Entry(self.w_sc)
            self.sc["sc_{0}".format(i)].grid(column=6,row= row_i) 
            self.sc["sc_{0}".format(i)].insert(INSERT,'Set_'+str(row_i+1))   
            row_i +=1
        bou_ok = Button(self.w_sc, text=u'OK', command = self.OK_sc_hyp)
        bou_ok.grid(column=4,row=row_i+1)
        
        self.w_sc.mainloop()
        
    def OK_sc_hyp(self) :
        self.sc_names = []
        for i in range(self.nb_of_sc):
            self.sc_names.append(self.sc["sc_{0}".format(i)].get())         
        self.w_sc.destroy() 
          
    
    '''#####################''' 
    def nb_bg_hyp(self):        
        self.w_bg_nb = tk.Tk()
        self.w_bg_nb.title('Number of background hypothesis for model '+str(self.model_i))
        label = Label(self.w_bg_nb, text='\nHow many background hypothesis for model '+str(self.model_i)+' do you want to use?')
        label.pack()
        self.nb_of_bg = Entry(self.w_bg_nb)
        self.nb_of_bg.pack()
        self.nb_of_bg.insert(INSERT,1)
        bou_ok = Button(self.w_bg_nb, text=u'OK', command = self.OK_nb_bg_hyp)
        bou_ok.pack()
        
        self.w_bg_nb.mainloop()
        
    def OK_nb_bg_hyp(self) :
        self.nb_of_bg = int(self.nb_of_bg.get())
        self.w_bg_nb.destroy() 
        self.bg_hyp()
        
    def bg_hyp(self):
        self.bg = {}       
        self.w_bg = tk.Tk()
        self.w_bg.grid()
        self.w_bg.title('Name of the background hypotheses for '+str(self.model_i))
        row_i = 0
        label = Label(self.w_bg, text='\nBackground hypotheses for '+str(self.model_i))
        label.grid(column=0,row=row_i) 
        row_i +=1
        for i in range(self.nb_of_bg):
            label = Label(self.w_bg, text="Hypothesis "+str(i+1))
            label.grid(column=0,row=row_i)
            self.bg["bg_{0}".format(i)] = Entry(self.w_bg)
            self.bg["bg_{0}".format(i)].grid(column=6,row= row_i) 
            self.bg["bg_{0}".format(i)].insert(INSERT,'BG_'+str(row_i))   
            row_i +=1
        bou_ok = Button(self.w_bg, text=u'OK', command = self.OK_bg_hyp)
        bou_ok.grid(column=4,row=row_i+1)
        
        self.w_bg.mainloop()
        
    def OK_bg_hyp(self) :
        self.bg_names = []
        for i in range(self.nb_of_bg):
            self.bg_names.append(self.bg["bg_{0}".format(i)].get())         
        self.w_bg.destroy()
          
     
    '''#####################'''    
    def get_available_models(self):
        NomFichier_InfosZonage = self.File_geom
        InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8')],skip_header = 1)
        Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
        self.available_models = list(np.unique(np.array(Column_model_name)))

    def get_available_scaling_laws(self):
        self.available_scaling_laws = ['WC1994','Le2010','HB08','TMG2017']
                

if __name__=="__main__":
    app = S_LT()

