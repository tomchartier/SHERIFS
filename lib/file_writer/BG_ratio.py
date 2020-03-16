# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

Open an interface fo the input of the background to fault ratio

@author: Thomas Chartier
"""

import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import ttk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton


class BG_ratio():
    def __init__(self,Run_Name,Model,BG_hyp,bg_ratio_file):  
        self. Run_Name = Run_Name
        self.Model = Model
        self.BG_hyp = BG_hyp
        self.bg_ratio_file = bg_ratio_file
        self.initialize()
    def initialize(self):
        
        self.Ratio_1 = {} #ratio of eq on faults
        self.Ratio_2 = {} #ratio of eq on faults
        self.Ratio_3 = {} #ratio of eq on faults
        self.Ratio_4 = {} #ratio of eq on faults
        self.Ratio_5 = {} #ratio of eq on faults
        self.Ratio_6 = {} #ratio of eq on faults
        self.Ratio_7 = {} #ratio of eq on faults
        self.Ratio_8 = {} #ratio of eq on faults
        self.Ratio_9 = {} #ratio of eq on faults
        
        self.w3 = tk.Tk()
        self.w3.title('Background hypotheis for : ' + str(self.Model) + ' and ' + str(self.BG_hyp))
        self.w3.grid()
        
        label = Label(self.w3, text="Ration of Earthquakes that happen on modelled faults")
        label.grid(column=1,row=1)
        
        label = Label(self.w3, text="Magnitude")
        label.grid(column=1,row=2)
        
        label = Label(self.w3, text="4.0")
        label.grid(column=2,row=2)
        self.Ratio_1["Ratio_1"] = Entry(self.w3)
        self.Ratio_1["Ratio_1"].grid(column=2,row= 3)
        self.Ratio_1["Ratio_1"].insert(INSERT,str(1.0)) 
        
        label = Label(self.w3, text="4.5")
        label.grid(column=3,row=2)
        self.Ratio_2["Ratio_2"] = Entry(self.w3)
        self.Ratio_2["Ratio_2"].grid(column=3,row= 3)
        self.Ratio_2["Ratio_2"].insert(INSERT,str(1.0)) 
        
        label = Label(self.w3, text="5.0")
        label.grid(column=4,row=2)
        self.Ratio_3["Ratio_3"] = Entry(self.w3)
        self.Ratio_3["Ratio_3"].grid(column=4,row= 3)
        self.Ratio_3["Ratio_3"].insert(INSERT,str(1.0)) 
        
        label = Label(self.w3, text="5.5")
        label.grid(column=5,row=2)
        self.Ratio_4["Ratio_4"] = Entry(self.w3)
        self.Ratio_4["Ratio_4"].grid(column=5,row= 3)
        self.Ratio_4["Ratio_4"].insert(INSERT,str(1.0))   
        
        label = Label(self.w3, text="6.0")
        label.grid(column=2,row=4)
        self.Ratio_5["Ratio_5"] = Entry(self.w3)
        self.Ratio_5["Ratio_5"].grid(column=2,row= 5)
        self.Ratio_5["Ratio_5"].insert(INSERT,str(1.0))  
        
        label = Label(self.w3, text="6.5")
        label.grid(column=3,row=4)
        self.Ratio_6["Ratio_6"] = Entry(self.w3)
        self.Ratio_6["Ratio_6"].grid(column=3,row= 5)
        self.Ratio_6["Ratio_6"].insert(INSERT,str(1.0))  
        
        label = Label(self.w3, text="7.0")
        label.grid(column=4,row=4)
        self.Ratio_7["Ratio_7"] = Entry(self.w3)
        self.Ratio_7["Ratio_7"].grid(column=4,row= 5)
        self.Ratio_7["Ratio_7"].insert(INSERT,str(1.0))  
        
        label = Label(self.w3, text="7.5")
        label.grid(column=5,row=4)
        self.Ratio_8["Ratio_8"] = Entry(self.w3)
        self.Ratio_8["Ratio_8"].grid(column=5,row= 5)
        self.Ratio_8["Ratio_8"].insert(INSERT,str(1.0))  
        
        label = Label(self.w3, text="8.0")
        label.grid(column=6,row=4)
        self.Ratio_9["Ratio_9"] = Entry(self.w3)
        self.Ratio_9["Ratio_9"].grid(column=6,row= 5)
        self.Ratio_9["Ratio_9"].insert(INSERT,str(1.0))   
        
        label = Label(self.w3, text="Ratio")
        label.grid(column=1,row=3)
        
        bou_ok = Button(self.w3, text=u'OK', command = self.ok)
        bou_ok.grid(column=3,row= 6)
            
        self.w3.mainloop()
          
    def ok(self):
        self.bg_ratio = []
        bg_ratio_file = open(self.bg_ratio_file,'w')

        self.bg_ratio.append(self.Ratio_1["Ratio_1"].get())
        bg_ratio_file.write(str(self.Ratio_1["Ratio_1"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_2["Ratio_2"].get())
        bg_ratio_file.write(str(self.Ratio_2["Ratio_2"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_3["Ratio_3"].get())
        bg_ratio_file.write(str(self.Ratio_3["Ratio_3"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_4["Ratio_4"].get())
        bg_ratio_file.write(str(self.Ratio_4["Ratio_4"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_5["Ratio_5"].get())
        bg_ratio_file.write(str(self.Ratio_5["Ratio_5"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_6["Ratio_6"].get())
        bg_ratio_file.write(str(self.Ratio_6["Ratio_6"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_7["Ratio_7"].get())
        bg_ratio_file.write(str(self.Ratio_7["Ratio_7"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_8["Ratio_8"].get())
        bg_ratio_file.write(str(self.Ratio_8["Ratio_8"].get())+ '\n')
        self.bg_ratio.append(self.Ratio_9["Ratio_9"].get())
        bg_ratio_file.write(str(self.Ratio_9["Ratio_9"].get())+ '\n')
            
        self.w3.destroy() 
        bg_ratio_file.close()

if __name__=="__main__":
    app = BG_ratio()

