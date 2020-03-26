# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
 
 
def progress(rate_in_model,model_MFD,calculation_log,ratio_done,print_percent):
    if ratio_done > 0.01 and ratio_done <= 0.25 and print_percent == True :
        print("1%")
        calculation_log.write("\n1%")
        model_MFD.append(rate_in_model)
        print_percent = False
    if ratio_done > 0.25 and ratio_done <= 0.5 and print_percent == False :
        print( "25%")
        calculation_log.write("\n25%")
        model_MFD.append(rate_in_model)
        print_percent = True
    if ratio_done > 0.5 and ratio_done <= 0.75 and print_percent == True :
        print( "50%")
        calculation_log.write("\n50%")
        model_MFD.append(rate_in_model)
        print_percent = False
    if ratio_done > 0.75 and ratio_done <= 0.9 and print_percent == False :
        print( "75%")
        calculation_log.write("\n75%")
        model_MFD.append(rate_in_model)
        print_percent = True
    if ratio_done > 0.9 and ratio_done <= 0.9999 and print_percent == True :
        print( "90%")
        calculation_log.write("\n90%")
        model_MFD.append(rate_in_model)
        print_percent = False
    return model_MFD,calculation_log,print_percent
