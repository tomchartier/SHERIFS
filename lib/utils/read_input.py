# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""

def extract_bg_input(file):
    '''
    file : str, file containing the seismicity ratio bg/fault
    
    available_bg : dict, seismicity ratio fault/total
    '''
    available_bg = {}
    with open(file,'r') as f:
        lines = f.read().splitlines()
        for line,i in zip(lines,range(len(lines))):
            # check if the current line
            # starts with "#"
            if line.startswith("#"):
                pass
            elif line.startswith("bg"):
                hyp_name = line.split(' ')[1]
                ratios = lines[i+1].split(' ')
                ratios = [float(x) for x in ratios if x != '']
                available_bg.update({'bg_'+hyp_name:ratios})
            else:
                pass
    return available_bg

def extract_sc_input(file):
    '''
    file : str, file containing the complex fault ruptures
    
    available_set : dict, complexe ruptures for each hypothesis
    '''
    available_sets = {}
    with open(file,'r') as f:
        lines = f.read().splitlines()
        for line,i in zip(lines,range(len(lines))):
            # check if the current line
            # starts with "#"
            if line.startswith("#"):
                pass
            elif line.startswith("set "):
                hyp_name = line.split(' ')[1]
                set = []
                loop_in = True
                for j in range(len(lines))[i+1:]:
                    if lines[j].startswith("set ")==False and loop_in == True:
                        faults = lines[j].split(' ')
                        faults = [str(x) for x in faults if x != '']
                        if len(faults) > 1:
                            set.append(faults)
                    else :
                        loop_in = False
                available_sets.update({'sc_'+hyp_name:set})
            else:
                pass
    return available_sets
