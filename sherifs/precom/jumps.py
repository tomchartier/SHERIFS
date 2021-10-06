# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

@author: Thomas Chartier
"""
from sherifs.precom.precomp_utils import *


def force_jump_list():
    # find the faults where a specific distance
    # is imposed a priori (to force a jump for example)
    # indicate the fault name from the oiler file
    force_jump_on_fault = ['Main Himalaya Thrust']
    return force_jump_on_fault


def find_jumps(f_for_sherifs, assso_fault, id_sections_fault, jump_dist,
               path, force_jump_on_fault):
    # section_jump contains the list of each sections that can rupture with the section
    section_jump = []
    d_jump = jump_dist

    file_jumps_1 = open(path+"/qgis/jumps_1.csv", 'w')
    file_jumps_1.write('lon,lat'+'\n')
    file_jumps_2 = open(path+"/qgis/jumps_2.csv", 'w')
    file_jumps_2.write('lon,lat'+'\n')
    file_jumps_line = open(path+"/qgis/jumps_lines.csv", 'w')
    file_jumps_line.write('section;wkt'+'\n')

    for si in range(len(f_for_sherifs)):
        si_jump = []
        fi = f_for_sherifs[si]["oiler_id"]
        lons_i = f_for_sherifs[si]["lons"]
        lats_i = f_for_sherifs[si]["lats"]


        for fj in assso_fault[fi]:
            if fi != fj :
                for sj in id_sections_fault[fj]:

                    shortest_pt = 100000
                    jump = False
                    maybe = False
                    lons_j = f_for_sherifs[sj]["lons"]
                    lats_j = f_for_sherifs[sj]["lats"]

                    for lon_i, lat_i in zip(lons_i,lats_i):
                        for lon_j, lat_j in zip(lons_j,lats_j):
                            dist = distance(lon_i, lat_i, lon_j, lat_j)
                            if dist<shortest_pt:
                                shortest_pt = dist
                                save = [lon_i, lat_i, lon_j, lat_j]
                    if shortest_pt < d_jump :
                        jump=True

                    if jump == True :
                        azi = calculate_initial_compass_bearing([lats_i[0],lons_i[0]],
                                                                [lats_i[-1],lons_i[-1]])
                        azj = calculate_initial_compass_bearing([lats_j[0],lons_j[0]],
                                                                [lats_j[-1],lons_j[-1]])
                        if abs(azi-azj) > 75. and abs(azi-azj) < 115. :
                            jump = False
                        if abs(azi-azj) > 255. and abs(azi-azj) < 185. :
                            jump = False

                    # check if the force need to be forced
                    if (jump == False
                    and f_for_sherifs[si]["oiler_name"] == f_for_sherifs[sj]["oiler_name"] and
                    f_for_sherifs[si]["oiler_name"] != None and
                       f_for_sherifs[si]["oiler_name"] in force_jump_on_fault):
                        shortest_pt_between_faults = 100000
                        for sk in id_sections_fault[fj]:
                            if (f_for_sherifs[si]["oiler_name"] == f_for_sherifs[sk]["oiler_name"]):
                                lons_k = f_for_sherifs[sk]["lons"]
                                lats_k = f_for_sherifs[sk]["lats"]
                                for sl in id_sections_fault[fi]:
                                    if (f_for_sherifs[si]["oiler_name"] == f_for_sherifs[sk]["oiler_name"]):
                                        lons_l = f_for_sherifs[sl]["lons"]
                                        lats_l = f_for_sherifs[sl]["lats"]
                                        for lon_l, lat_l in zip(lons_l,lats_l):
                                            for lon_k, lat_k in zip(lons_k,lats_k):
                                                dist = distance(lon_l, lat_l, lon_k, lat_k)
                                                if dist<shortest_pt_between_faults:
                                                    shortest_pt_between_faults = dist
                        if shortest_pt_between_faults == shortest_pt :
                            jump = True
    #                         if f_for_sherifs[si]["oiler_fid"] == 11:
#                            print("force jump",si,sj)
    #                         print(shortest_pt_between_faults,shortest_pt)
    #                         print(f_for_sherifs[si]["oiler_fid"],f_for_sherifs[sj]["oiler_fid"])
    #                         print()
    #                     else :
    #                         if f_for_sherifs[si]["oiler_fid"] == 11:
    #                             print("dont",si,sj)
    #                             print(shortest_pt_k,shortest_pt)
    #                             print(f_for_sherifs[si]["oiler_fid"],f_for_sherifs[sj]["oiler_fid"])
    #                             print()

                    if jump == True :
                        si_jump.append(sj)
                        file_jumps_1.write(str(save[0])+','+str(save[1])+'\n')
                        file_jumps_2.write(str(save[2])+','+str(save[3])+'\n')
                        file_jumps_line.write(str(si)+";LINESTRING("+
                                              str(save[0])+' '+
                                              str(save[1])+','+
                                              str(save[2])+' '+
                                              str(save[3])+')\n')
        section_jump.append(si_jump)
    file_jumps_1.close()
    file_jumps_2.close()
    file_jumps_line.close()

    print("There are ",sum([len(i) for i in section_jump]),"jumps.")

    return section_jump

