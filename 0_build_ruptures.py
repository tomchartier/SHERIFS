# -*- coding: utf-8 -*-

"""
Prepare the input fils and ruptures for SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

using as input the geojson file from the global fault database
(oiler output)

@author: Thomas Chartier
contact : thomas.chartier@globalquakemodel.org
"""

import time
import os
import sys
path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/lib'
sys.path.append(path_lib)
path_f = path_lib + '/precom'
sys.path.append(path_f)
from lib.utils import sap
from read_precomp_files import read_oiler_file
from precomp_utils import *
from build_scenarios import *
from build_sections import *
from jumps import *

# If you are running SHERIFS with spyder define "input_file" here. Then run.

def build_rup(input_file):
    debut = time.time()

    
    print ('\nPreparting the rutpures for SHERIFS\n')

    # Load the input file
    lines = open(input_file,'r').readlines()
    lines = [line.rstrip('\n') for line in lines]
    apply_sr_reduction = 0.
    for line in lines:
        if "Run_Name" in line :
            Run_Name = line.split(':')[1].replace(' ','')
        if "File_Oiler" in line :
            File_Oiler = line.split(':')[1].replace(' ','')
        if "File_out" in line :
            File_out = line.split(':')[1].replace(' ','')
        if "File_Mmax_areas" in line :
            File_Mmax_areas = line.split(':')[1].replace(' ','')
        if "Model_name" in line :
            Model_name = line.split(':')[1].replace(' ','')
            
#        if "use_host_model" in line :
#            if "rue" in line :
#                use_host_model = True
#            elif "alse" in line :
#                use_host_model = False

                    
        if "File_Mu" in line :
            f_mu = line.split(':')[1].replace(' ','')
            
        if "jump_dist" in line :
            jump_dist = float(line.split(':')[1].replace(' ',''))
        if "apply_sr_reduction" in line :
            apply_sr_reduction = float(line.split(':')[1].replace(' ',''))
            
    path = "input/"+Run_Name
            
    # Reading things
    faults,nb_faults = read_oiler_file(File_Oiler)
    
    #bounding box for faults
    maxmin_pt_lon, maxmin_pt_lat = find_bounding_box(faults)
    
    #find possible fault asso
    assso_fault = find_possible_asso(maxmin_pt_lon,maxmin_pt_lat)
    
    # calc fault dimensions
    f_lengths, f_areas = calc_f_dims(faults)
    
    # cutting into smaller sections
    f_for_sherifs,id_sections_fault,sections_areas_tot,sections_lengths_tot = cut_faults(faults,
    f_lengths,
    f_areas,
    path)
    
    # force jumps
    force_jump_on_fault = force_jump_list()
    
    # find jumps
    section_jump = find_jumps(f_for_sherifs,
    assso_fault,
    id_sections_fault,
    jump_dist,
    path,
    force_jump_on_fault)
    
    # export section points
    export_sections_pts(f_for_sherifs,path)
    
    # find max rupture size for each section
    f_for_sherifs = find_sections_Mmax(f_for_sherifs,File_Mmax_areas)
    
    # create SHERIFS input dict for fault prperties
    f_for_sherifs = to_sherifs(f_for_sherifs,
    faults,
    Model_name,
    apply_sr_reduction,
    f_mu)
    
    # create binning in mag
    binning_in_mag, nb_rup_per_bin = mag_bin_distr(f_for_sherifs)
    
    # CORE : create the list of ruptures
    rup = build_scenarios(f_for_sherifs,
    id_sections_fault,
    sections_areas_tot,
    sections_lengths_tot,
    binning_in_mag,
    nb_rup_per_bin,
    section_jump)
    
    # write SHERIFS input file for faults
    write_section_json(f_for_sherifs,File_out)
    
    # write SHERIFS input file for ruptures
    write_rupt_file(rup,Run_Name)
    
    # Create visualization of ruptures
    #visu_rup(f_for_sherifs,rup,rups_length,rups_mag,path)

    fin = time.time()-debut
    days = int(fin / 24. / 60. / 60.)
    hours = int((fin - days * 24. * 60. * 60.) / 60. / 60.)
    minutes = int((fin - days * 24. * 60. * 60. - hours* 60. * 60. ) / 60.)
    seconds = (fin - days * 24. * 60. * 60. - hours* 60. * 60.  - minutes * 60.)
    print("The calculation took: " + str(days) + ' days, ' + str(hours) + ' hours, ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds.')



def main(argv):
    """ Run SHERIFS"""

    p = sap.Script(build_rup)
    p.arg(name='input_file', help='.txt file with the information concerning the run.')

    if len(argv) < 1:
        print(p.help())
    else:
        p.callfunc()


if __name__ == "__main__":
    main(sys.argv[1:])
