#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

Reads the different files containing all the information about the given run.

@author: thomas
"""
import numpy as np
import xml.etree.ElementTree as ET

def read_file_sr_log(log_sr_file):
    try:
        log_file = np.genfromtxt(log_sr_file,dtype=[('U100'),('f8')])
    except ValueError:
        log_file = np.genfromtxt(log_sr_file,dtype=[('U100'),('f8')],skip_footer = 1)
    faults_names =  list(map(lambda i : log_file[i][0], range(len(log_file))))
    faults_slip_rates =  list(map(lambda i : log_file[i][1], range(len(log_file))))
    faults_names = [str(i) for i in faults_names]
    return faults_names, faults_slip_rates
    
def read_file_as_log(log_as_file):
    log_file = np.genfromtxt( log_as_file,dtype=[('U100'),('f8'),('U100')],delimiter = '\t')
    try:
        a_s = log_file[2][1]
    except IndexError:
        print('\n\nThere is a problem with the following file : \n',log_as_file,'\n\n!!!!!')
    return a_s
    
def read_file_Mmax_log(log_Mmax_file):
    log_file = np.genfromtxt( log_Mmax_file,dtype=[('U10000'),('f8'),('U100'),('f8'),('f8'),('f8')],delimiter = '\t')
    sources_names =  list(map(lambda i : log_file[i][0], range(len(log_file))))
    sources_Areas =  list(map(lambda i : log_file[i][4], range(len(log_file))))
    sources_Lengths =  list(map(lambda i : log_file[i][1], range(len(log_file))))
    sources_Mmax =  list(map(lambda i : log_file[i][5], range(len(log_file))))
    sources_names = [str(i) for i in sources_names]
    return sources_names,sources_Mmax,sources_Lengths,sources_Areas
    
def read_file_general_param_log(log_general_param_file):
    log_file = np.genfromtxt( log_general_param_file,dtype=[('U100'),('f8')],delimiter = '\t')
    M_trunc_model =  log_file[0][1]
    b_value_model =  log_file[1][1]
    return M_trunc_model,b_value_model


def read_file_mfd_log(log_mfd_file):
    #print log_mfd_file
    log_file = np.genfromtxt( log_mfd_file,dtype = [('U1000'),('f8'),('U1000')],delimiter = '\t')
    sources_names =  list(map(lambda i : log_file[i][0], range(len(log_file))))
    Mmin =  list(map(lambda i : log_file[i][1], range(len(log_file)))        )
    MFD =  list(map(lambda i : log_file[i][2].split(' '), range(len(log_file))))
    sources_names = [str(i) for i in sources_names]
    return sources_names,Mmin,MFD
    

def read_hc_xml(hc_xml_file):
    tree = ET.parse( hc_xml_file)        
    nrml = tree.getroot()
   
    for hazardCurves in nrml:
        IMLs = str(hazardCurves[0].text)
        IMLs = np.array( IMLs.split(),float)
        poEs = str(hazardCurves[1][1].text)
        poEs = np.array( poEs.split(),float)     
        
         
    return IMLs,poEs



def read_LT_xml(Run_name):
    LT_file = str(Run_name) + '/Sources_Logic_tree.xml'
    tree = ET.parse(LT_file)        
    nrml = tree.getroot()
    
    Branch_names = []
    Branch_weight = []        
    
    for logicTree in nrml:
        for logicTreeBranchLevel in logicTree:
            for logicTreeBranchSet in logicTreeBranchLevel:
                for logicTreeBranch in logicTreeBranchSet:
                    logicTreeBranch_ID = logicTreeBranch.attrib['branchID']
                    Branch_names.append(logicTreeBranch_ID.split("-"))
                    Branch_weight.append(float(logicTreeBranch[1].text))
    return Branch_names,Branch_weight
    


def read_GMPE_LT_xml(Run_name):
    GMPE_LT_file = str( Run_name) + '/GMPE_Logic_tree.xml'
    tree = ET.parse(GMPE_LT_file)        
    nrml = tree.getroot()
    
    GMPE = []  
    Tectonic_region = []
    for logicTree in nrml:
        for logicTreeBranchLevel in logicTree:
            for logicTreeBranchSet in logicTreeBranchLevel:
                logicTreeBranchSet_ID = logicTreeBranchSet.attrib['branchSetID']
                for logicTreeBranch in logicTreeBranchSet:
                    logicTreeBranch_ID = logicTreeBranch.attrib['branchID']
                    GMPE.append(logicTreeBranch_ID) 
                    Tectonic_region.append(logicTreeBranchSet_ID)
    return GMPE,Tectonic_region

