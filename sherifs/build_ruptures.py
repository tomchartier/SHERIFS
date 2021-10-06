# -*- coding: utf-8 -*-

"""
Prepare the input fils and ruptures for SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

using as input the geojson file from the global fault database
(oiler output)

@author: Thomas Chartier
contact : thomas.chartier@globalquakemodel.org
"""

import os
import sys
import time
import toml
from pathlib import Path

from sherifs.utils import sap
from sherifs.precom.read_precomp_files import read_oiler_file
from sherifs.precom.precomp_utils import *
from sherifs.precom.build_scenarios import *
from sherifs.precom.build_sections import *
from sherifs.precom.jumps import *

# If you are running SHERIFS with spyder define "input_file" here. Then run.


def build_rup(input_file):
    debut = time.time()

    print('\nPreparing the ruptures for SHERIFS\n')

    # Load the input file
    param = toml.load(input_file)

    # This is the path to the .toml file. All the paths in this file are
    # relative to its position
    root = os.path.dirname(input_file)

    # This is the path where the output from build_ruptures will be created
    dirpath = os.path.join(root, param["dirpath"])
    Path(dirpath).mkdir(parents=True, exist_ok=True)

    Run_Name = param["Run_Name"]
    Set_Name = param["pre"]["Set_Name"]

    File_Oiler = os.path.join(root, param["pre"]["File_Oiler"])
    File_out = os.path.join(root, param["pre"]["File_out"])
    File_Mmax_areas = os.path.join(root, param["pre"]["File_Mmax_areas"])
    Model_name = param["pre"]["Model_name"]
    rupture_mesh_spacing = param["pre"]["rupture_mesh_spacing"]
    f_mu = os.path.join(root, param["pre"]["File_Mu"])
    jump_dist = param["pre"]["jump_dist"]
    apply_sr_reduction = param["pre"]["apply_sr_reduction"]

    path = os.path.join(dirpath, Run_Name)

    # Reading things
    faults, nb_faults = read_oiler_file(File_Oiler)

    #bounding box for faults
    maxmin_pt_lon, maxmin_pt_lat = find_bounding_box(faults)

    #find possible fault asso
    assso_fault = find_possible_asso(maxmin_pt_lon,maxmin_pt_lat)

    # calc fault dimensions
    f_lengths, f_areas = calc_f_dims(faults,)

    # cutting into smaller sections
    f_for_sherifs,id_sections_fault,sections_areas_tot,sections_lengths_tot = cut_faults(faults,
    f_lengths,
    f_areas,
    path,
    rupture_mesh_spacing)

    # force jumps
    force_jump_on_fault = force_jump_list()

    # find jumps
    section_jump = find_jumps(f_for_sherifs, assso_fault, id_sections_fault,
                              jump_dist, path, force_jump_on_fault)

    # export section points
    export_sections_pts(f_for_sherifs, path)

    # find max rupture size for each section
    f_for_sherifs = find_sections_Mmax(f_for_sherifs,File_Mmax_areas)

    # create SHERIFS input dict for fault prperties
    f_for_sherifs = to_sherifs(f_for_sherifs, faults, Model_name,
                               apply_sr_reduction, f_mu)

    # create binning in mag
    binning_in_mag, nb_rup_per_bin = mag_bin_distr(f_for_sherifs)

    # CORE : create the list of ruptures
    rup, rup_param = build_scenarios(f_for_sherifs, id_sections_fault,
                                     sections_areas_tot, sections_lengths_tot,
                                     binning_in_mag, nb_rup_per_bin,
                                     section_jump)

    # write SHERIFS input file for faults
    write_section_json(f_for_sherifs, File_out)

    # write SHERIFS input file for ruptures
    write_rupt_file(path, rup, Run_Name, Set_Name)

    # Create visualization of ruptures
    visu_rup(f_for_sherifs, rup, rup_param[0], rup_param[1], path)

    fin = time.time()-debut
    days = int(fin / 24. / 60. / 60.)
    hours = int((fin - days * 24. * 60. * 60.) / 60. / 60.)
    minutes = int((fin - days * 24. * 60. * 60. - hours* 60. * 60. ) / 60.)
    seconds = (fin - days * 24. * 60. * 60. - hours* 60. * 60.  - minutes * 60.)
    print("The calculation took: " + str(days) + ' days, ' + str(hours) +
          ' hours, ' + str(minutes) + ' minutes and ' + str(seconds) +
          ' seconds.')


def main(argv):
    """ Run SHERIFS"""

    p = sap.Script(build_rup)
    msg = '.txt file with the information concerning the run.'
    p.arg(name='input_file', help=msg)

    if len(argv) < 1:
        print(p.help())
    else:
        p.callfunc()


if __name__ == "__main__":
    main(sys.argv[1:])
