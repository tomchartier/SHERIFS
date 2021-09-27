# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
from fpdf import FPDF
#from PIL import Image


def create_title_page(Run_name,pdf):
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial', 'B', 18)
    pdf.cell(60)
    pdf.cell(90, 10, " ", 0, 2, 'C')
    pdf.cell(90, 10, " ", 0, 2, 'C')
    pdf.cell(90, 10, " ", 0, 2, 'C')
    pdf.cell(75, 10, Run_name, 0, 2, 'C')
    pdf.cell(90, 10, " ", 0, 2, 'C')
    pdf.cell(90, 10, " ", 0, 2, 'C')
    
    return pdf

def print_lt(logictree,pdf):
    # built a representation of the logic tree
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial', 'B', 8)
    pdf.cell(90, 10, " ", 0, 2, 'C')
    for key in logictree:
        for str_i in logictree[key]:
            pdf.cell(75, 10, str_i, 0, 2, 'C')
        pdf.cell(90, 10, " ", 0, 2, 'C')
        pdf.cell(90, 10, " ", 0, 2, 'C')
    
    
    return pdf
        
def compare_mfd_subareas(Run_name,pdf):
    
    
    return pdf
    
def create(Run_name,logictree):
    pdf = FPDF()
    
    # Create a nice title page
    pdf = create_title_page(Run_name,pdf)
    
    # Create a logic tree representation
    pdf = print_lt(logictree,pdf)
    
    # Compare the mfd for the subarea sources
    pdf = compare_mfd_subareas(Run_name,pdf)
    
    pdf.output(Run_name+'/report_'+Run_name+'.pdf', 'F')


#img = Image.open("data_mask_1354_2030.png")
#
#background = Image.open("background_1354_2030.png")
#
#background.paste(img, (0, 0), img)
#background.save('how_to_superimpose_two_images_01.png',"PNG")
