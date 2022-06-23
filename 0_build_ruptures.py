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
import toml

# If you are running SHERIFS with spyder define "input_file" here. Then run.

def build_rup(input_file):
    debut = time.time()


    print ('\nPreparting the rutpures for SHERIFS\n')

    # Load the input file
    param = toml.load(input_file)

    dirpath = param["dirpath"]
    if dirpath[-1] != "/":
        dirpath += "/"

    Run_Name = param["Run_Name"]
    Set_Name = param["pre"]["Set_Name"]
    File_Mmax_areas = dirpath + param["pre"]["File_Mmax_areas"]
    Model_name = param["pre"]["Model_name"]
    rupture_mesh_spacing = param["pre"]["rupture_mesh_spacing"]
    f_mu = dirpath + param["pre"]["File_Mu"]
    jump_dist = param["pre"]["jump_dist"]
    apply_sr_reduction = param["pre"]["apply_sr_reduction"]

    sectionning_param = {}
    if "max_section_length" in param["pre"].keys():
        sectionning_param.update({"max_section_length":param["pre"]["max_section_length"]})
    else :
        sectionning_param.update({"max_section_length":40.0})

    if "max_num_sections" in param["pre"].keys():
        sectionning_param.update({"max_num_sections":param["pre"]["max_num_sections"]})
    else :
        sectionning_param.update({"max_num_sections":6})

    if "distance_resolution" in param["pre"].keys():
        sectionning_param.update({"distance_resolution":param["pre"]["distance_resolution"]})
    else :
        sectionning_param.update({"distance_resolution":0.5})


    path = dirpath + "input/"+Run_Name

    if  "File_Oiler" in param["pre"].keys() :
        do_sectionning = True
        File_faults = dirpath + param["pre"]["File_Oiler"]
        File_out = dirpath + param["pre"]["File_out"]

    else :
        do_sectionning = False
        File_faults = dirpath + param["pre"]["File_sections"]
        File_out = File_faults


    # Reading things
    faults,nb_faults = read_oiler_file(File_faults)

    #bounding box for faults
    maxmin_pt_lon, maxmin_pt_lat = find_bounding_box(faults)

    #find possible fault asso
    assso_fault = find_possible_asso(maxmin_pt_lon,maxmin_pt_lat)

    # calc fault dimensions
    f_lengths, f_areas = calc_f_dims(faults,)

    if do_sectionning == True :
        # cutting into smaller sections
        f_for_sherifs,id_sections_fault,sections_areas_tot,sections_lengths_tot = cut_faults(faults,
        f_lengths,
        f_areas,
        path,
        rupture_mesh_spacing,
        sectionning_param)
    else :
        # converts faults to sections
        f_for_sherifs,id_sections_fault,sections_areas_tot,sections_lengths_tot = converts_to_sections(faults,
        f_lengths,
        f_areas,
        path,
        rupture_mesh_spacing,
        sectionning_param)

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
    rup,rup_param = build_scenarios(f_for_sherifs,
    id_sections_fault,
    sections_areas_tot,
    sections_lengths_tot,
    binning_in_mag,
    nb_rup_per_bin,
    section_jump)

    # write SHERIFS input file for faults
    write_section_json(f_for_sherifs,File_out)

    # write SHERIFS input file for ruptures
    write_rupt_file(rup,Run_Name,Set_Name,dirpath)

    # Create visualization of ruptures
    visu_rup(f_for_sherifs,rup,rup_param[0],rup_param[1],path)

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
