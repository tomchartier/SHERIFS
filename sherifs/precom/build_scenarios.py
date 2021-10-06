# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

@author: Thomas Chartier
"""
import os
import time
import numpy as np
from sherifs.precom.precomp_utils import *
from shapely.geometry import MultiPoint
from geojson import Feature, FeatureCollection, dump


def mag_bin_distr(f_for_sherifs):
    # prepare the 2d list for counting the number of
    # ruptures with Mmax in each bin of magnitude for each section
    binning_in_mag = np.array([round(i,1) for i in list(np.linspace(5,10,51))])
    nb_rup_per_bin = []
    for si in range(len(f_for_sherifs)):
        nb_for_i = []
        for j in binning_in_mag:
            nb_for_i.append(0)
        nb_rup_per_bin.append(nb_for_i)
    return binning_in_mag, nb_rup_per_bin


def build_scenarios(f_for_sherifs,id_sections_fault,sections_areas_tot,sections_lengths_tot,binning_in_mag,nb_rup_per_bin,section_jump):
    t0 = time.time()

    # parameters
    #max_length = 350
    max_fault_jumps = 100
    max_fault_invo = 100
    max_zig_zag = 2

    nb_sections = len(f_for_sherifs)

    rup = []
    rup_paired = []
    full_fault_rup = []
    rup_f_jump = []
    #rup_id_for_s = [[] for _ in range(len(f_for_sherifs))]
    for si in range(nb_sections):
        f_for_sherifs[si]["rup_id"] = []
        #f_for_sherifs[si]["rup_to_pair"] = []
    rups_length = []
    rups_area = []
    rups_rake = []
    rups_mag = []

    target_section =False
    if target_section == True:
        # this is to pring details about one given section
        look_at = 917

    too_long = 0
    too_big = 0

    rup_id = 0
    #solo ruptures
    for si in range(len(f_for_sherifs)):
        new_rup = [si]
        rup.append(new_rup)
        rup_paired.append([])
        rup_area = sections_areas_tot[si]
        rake = f_for_sherifs[si]["rake"]
        mag = wc1994_median_mag( rup_area, rake)
        index_mag = np.where(binning_in_mag==round(mag,1))[0][0]

        for sk in new_rup:
            f_for_sherifs[sk]["rup_id"].append(rup_id)
            nb_rup_per_bin[sk][index_mag]+=1
            f_for_sherifs[sk]["max_length"] = sections_lengths_tot[si]

        rups_length.append(sections_lengths_tot[si])
        rups_area.append(rup_area)
        rups_rake.append(f_for_sherifs[si]["rake"])
        rups_mag.append(mag)


        full_fault = False
        len_2 = len(id_sections_fault[f_for_sherifs[new_rup[0]]["oiler_id"]])
        if len_2 == len(new_rup):
            full_fault = True
        full_fault_rup.append(full_fault)
        rup_f_jump.append([])

        rup_id+=1

    #ruptures in the faults
    for si in range(nb_sections):
        fi = f_for_sherifs[si]["oiler_id"]
        sec_fi = id_sections_fault[fi] #sections of the fault fi
        if len(sec_fi)> 1:
            id_si_in_fi = sec_fi.index(si)
            for i in range(len(sec_fi)-id_si_in_fi):
                new_rup = sec_fi[id_si_in_fi:id_si_in_fi+i+1]
                new_rup = sorted(list(set(new_rup)))
                rup_length = sum([sections_lengths_tot[sk] for sk in new_rup])
                max_length = min([f_for_sherifs[i]["max_possible_length"] for sk in new_rup])

                rup_area = sum([sections_areas_tot[sk] for sk in new_rup])
                max_mmax = min([f_for_sherifs[sk]["max_possible_Mmax"] for sk in new_rup])
                rake = sum([sections_areas_tot[sk]*f_for_sherifs[sk]["rake"] for sk in new_rup])/rup_area
                mag = wc1994_median_mag( rup_area, rake)
                index_mag = np.where(binning_in_mag==round(mag,1))[0][0]


                if len(new_rup)>1 and rup_length < max_length:
                    if mag <= max_mmax :
                        if not new_rup in rup :
                            # number of different faults in the rupture
                            nb_fault_invo = len(set([f_for_sherifs[i]["oiler_id"] for sk in new_rup]))
                            if nb_fault_invo <= max_fault_invo:
                                rup.append(new_rup)
                                rup_paired.append([])
                                for sk in new_rup:
                                    nb_rup_per_bin[sk][index_mag]+=1
                                    f_for_sherifs[sk]["rup_id"].append(rup_id)
                                    #f_for_sherifs[sk]["rup_to_pair"].append(rup_id)
                                    if rup_length > f_for_sherifs[sk]["max_length"]:
                                        f_for_sherifs[sk]["max_length"] = rup_length

                                rups_length.append(rup_length)
                                rups_area.append(rup_area)
                                rups_rake.append(rake)
                                rups_mag.append(mag)


                                full_fault = False
                                if nb_fault_invo == 1:
                                    # number of sections in the fault
                                    len_2 = len(id_sections_fault[f_for_sherifs[new_rup[0]]["oiler_id"]])
                                    if len_2 == len(new_rup):
                                        full_fault = True
                                full_fault_rup.append(full_fault)
                                rup_f_jump.append([])

                                rup_id+=1
                    else :
                        too_big += 1
                else :
                    too_long += 1

    if target_section == True:
        print("look at section :", look_at,
              "\nother sections in the same fault",id_sections_fault[f_for_sherifs[look_at]["oiler_id"]],
              "\njumps with",section_jump[look_at],
              "\nruptures",f_for_sherifs[look_at]["rup_id"])
        print("ruptures :")
        for rup_i in [rup[i] for i in f_for_sherifs[look_at]["rup_id"]]:
            print(rup_i)


    print("\nNumber of rup without jumps: ",len(rup),"\n")

    last_l = -1
    loop = 1
    already_created = 0
    already_paired = 0
    too_many_faults = 0
    too_many_jumps = 0
    too_many_zig_zag = 0
    not_diverse_enough = 0

    while last_l != len(rup):
        print("\n##\n-> Loop",loop)
        last_l = len(rup)
        for si in range(nb_sections):
            sec_fi = id_sections_fault[fi] #sections of the fault fi
            jumps = section_jump[si]
            si_rup = f_for_sherifs[si]["rup_id"]
            for sj in jumps :
                sj_rup = f_for_sherifs[sj]["rup_id"]

                # removing rupture id if the other fault is already in it.
                sj_rup = [i for i in sj_rup if not si in rup[i]]
                sii_rup = [i for i in si_rup if not sj in rup[i]]

                for id_rup_i in sii_rup :

                    # Filter if the rupture to pair is alredy quite large
                    # If it is large, it can only pair with full fault ruptures
                    if rups_length[id_rup_i] > max_length/2.:
                        sjj_rup = [j for j in sj_rup if full_fault_rup[j] == True]
                    else:
                        sjj_rup = [j for j in sj_rup if rups_length[j] < max_length/2.]

                    for id_rup_j in sjj_rup :
                        # boolean to know if we need to build to rupture
                        build = True

                        # Looking if the number of jumps doesn't exceed the max
                        nb_jumps = len(rup_f_jump[id_rup_j]+rup_f_jump[id_rup_i])+1
                        if nb_jumps <= max_fault_jumps:
                            build = True
                        else:
                            build = False
                            too_many_jumps += 1

                        # Looking if the number of zigzag doesn't exceed the max
                        if build == True:
                            jump_i = sorted([f_for_sherifs[si]["oiler_id"],f_for_sherifs[sj]["oiler_id"]])
                            nb_zig_zag = rup_f_jump[id_rup_j].count(jump_i)+rup_f_jump[id_rup_i].count(jump_i)
                            if not nb_zig_zag > max_zig_zag :
                                build = True
                            else:
                                build = False
                                too_many_zig_zag += 1

                        # Looking if these ruptures have been paired before already
                        if build == True:
                            if id_rup_j in rup_paired[id_rup_i]:
                                build = False
                                already_paired += 1
                        if build == True:
                            if id_rup_i in rup_paired[id_rup_j] :
                                build = False
                                already_paired += 1

                        # Looking if the length isn't larger than the max
                        if build == True:
                            new_rup = rup[id_rup_i]+rup[id_rup_j]
                            new_rup = sorted(list(set(new_rup)))
                            rup_length = rups_length[id_rup_i] + rups_length[id_rup_j]
                            max_length = min([f_for_sherifs[sk]["max_possible_length"] for sk in new_rup])


                            #rup_area = sum([sections_areas_tot[i] for i in new_rup])
                            rup_area = rups_area[id_rup_i] + rups_area[id_rup_j]
                            max_mmax = min([f_for_sherifs[sk]["max_possible_Mmax"] for sk in new_rup])
                            rake =  (rups_rake[id_rup_i]*rups_area[id_rup_i]
                                     + rups_rake[id_rup_j]*rups_area[id_rup_j]
                                    )/rup_area
                            mag = wc1994_median_mag( rup_area, rake)
                            index_mag = np.where(binning_in_mag==round(mag,1))[0][0]


                            if rup_length > max_length :
                                build = False
                                too_long += 1
                                rup_paired[id_rup_i].append(id_rup_j)
                                rup_paired[id_rup_j].append(id_rup_i)
                            if build == True and mag > max_mmax :
                                build = False
                                too_big += 1
                                rup_paired[id_rup_i].append(id_rup_j)
                                rup_paired[id_rup_j].append(id_rup_i)


                        if build == True:
                            if not new_rup in rup:
                                build = True
                                # number of different faults in the rupture
                                nb_fault_invo = len(set([f_for_sherifs[i]["oiler_id"] for i in new_rup]))
                            else :
                                build = False
                                already_created += 1
                                rup_paired[id_rup_i].append(id_rup_j)
                                rup_paired[id_rup_j].append(id_rup_i)

                        if build == True:
                            if nb_fault_invo <= max_fault_invo:
                                build = True
                            else :
                                build = False
                                too_many_faults += 1
                                rup_paired[id_rup_i].append(id_rup_j)
                                rup_paired[id_rup_j].append(id_rup_i)

                        # Test is ruptures are getting too diverse
                        if build == True:
                            test_using_length = False
                            test_using_mag = True

                            if test_using_length == True :
                                nb_rup_sec = [len(f_for_sherifs[i]["rup_id"]) for i in new_rup]
                                max_nb_rup_i = max(nb_rup_sec)
                                mean_nb_rup = len(rup)/nb_sections

                                # Test if the new rupture is longer for at least one of the faults
                                min_max_length_i = min([f_for_sherifs[i]["max_length"] for i in new_rup])
                                mean_max_length_i = np.mean([f_for_sherifs[i]["max_length"] for i in new_rup])

                                hogger = False
        #                         if max_nb_rup_i > 40\
        #                         and rup_length < min_max_length_i * 1.05:
        #                             hogger = True
                                if max_nb_rup_i > 40\
                                and rup_length < mean_max_length_i * 1.5:
                                    hogger = True

                                if hogger == True:
                                    nb_i_to_remove = 0
                                    if rup_length < min_max_length_i/2. :
                                        if len(rup[id_rup_j]) < 2 or len(rup[id_rup_i]) < 2:
                                            build = False
                                        else :
                                            nb_i_to_remove = 2
                                    elif rup_length < min_max_length_i*3./4. :
                                        if len(rup[id_rup_j]) < 3 or len(rup[id_rup_i]) < 3:
                                            build = False
                                        else :
                                            nb_i_to_remove = 3
                                    else :
                                        if len(rup[id_rup_j]) < 4 or len(rup[id_rup_i]) < 4:
                                            build = False
                                        else :
                                            nb_i_to_remove = 4


                                    if build == True and nb_i_to_remove != 0:
                                        r = list(combinations(new_rup, nb_i_to_remove))
                                        grp_i = 0
                                        while build == True and grp_i<len(r):
                                            if len(r[grp_i]) == nb_i_to_remove :
                                                test_rup = [i for i in new_rup if not i in r[grp_i]]
                                                if test_rup in rup:
                                                    build = False
                                                    if target_section == True:
                                                        if si ==look_at :
                                                            print(id_rup_j+id_rup_i,test_rup,
                                                                  "already there|nb_i_to_remove:",nb_i_to_remove)
                                            grp_i+=1

                            elif test_using_mag == True :
                                #mean number of rup already in this bin for the sections involved
                                nb_i = np.mean([nb_rup_per_bin[sk][index_mag] for sk in new_rup])
                                if nb_i > 7 :
                                    # nb_j is the number of sections in the new rup for
                                    # which the number of rupture in the bin is part of the
                                    # bins with already the most ruptures (percentile 80)
                                    nb_j = 0
                                    for sk in new_rup :
                                        tmp_list = np.trim_zeros(nb_rup_per_bin[sk])
                                        p=np.percentile(tmp_list,66)
                                        if nb_rup_per_bin[sk][index_mag]>=p:
                                            nb_j+=1
                                    if nb_j >= len(new_rup)/2.:
                                        build = False

                            if build == False :
                                not_diverse_enough += 1
                                rup_paired[id_rup_i].append(id_rup_j)
                                rup_paired[id_rup_j].append(id_rup_i)


                        ##################
                        # Add the rupture
                        ##################
                        if build == True:
                            rup.append(new_rup)
                            rup_paired.append([])
                            rup_paired[id_rup_i].append(id_rup_j)
                            rup_paired[id_rup_j].append(id_rup_i)
                            l = len(rup)
                            for sk in new_rup:
                                f_for_sherifs[sk]["rup_id"].append(rup_id)
                                nb_rup_per_bin[sk][index_mag]+=1
                                if rup_length > f_for_sherifs[sk]["max_length"]:
                                    f_for_sherifs[sk]["max_length"] = rup_length
                            rups_length.append(rup_length)
                            rups_area.append(rup_area)
                            rups_rake.append(rake)
                            rups_mag.append(mag)
                            full_fault_rup.append(False)
                            rup_f_jump.append(rup_f_jump[id_rup_j]+rup_f_jump[id_rup_i]+[jump_i])



                            if "000" == str(l)[-3:]:# or "500" == str(l)[-3:]:
                                print(l, "ruptures   |   active section :",si,
                                      " |    not_diverse_enough :",not_diverse_enough)

                            rup_id+=1

            if target_section == True:
                if si ==look_at :
                    print(si,len(f_for_sherifs[si]["rup_id"]))
        print("\n -> Number of rup with ",loop," loops: ",len(rup))
        loop += 1


    print("\n########\n")
    print("Total ruptures : ", len(rup))
    print("The longest rupture is : ", round(max(rups_length))," km.")
    print("The largest magnitude is : ", round(max(rups_mag),1),".")
    print("Max number of ruptures for a section :",
          max([len(f_for_sherifs[i]["rup_id"]) for i in range(len(f_for_sherifs))]))
    print("Number of ruptures created for nothing because already there : ", already_created)
    print("Number of ruptures tried for nothing because already paired : ", already_paired)
    print("Number of ruptures not used because too many faults : ", too_many_faults)
    print("Number of ruptures not used because too many jumps : ", too_many_jumps)
    print("Number of ruptures not used because too many zigzag : ", too_many_zig_zag)
    print("Number of ruptures not used because not diverse enough : ", not_diverse_enough)
    print("Number of ruptures not used because too long : ", too_long)
    print("Number of ruptures not used because magnitude was too large : ", too_big)

    if target_section == True:
        print("\n########\n")
        print("ruptures for the look_at :","    look_at id is ",look_at)
        print("number of ruptures : ",len(f_for_sherifs[look_at]["rup_id"]))
        for i in f_for_sherifs[look_at]["rup_id"]:
            if rup.count(rup[i]) != 1:
                print("error")
            print(i,rup[i], round(rups_length[i]),'km | M :', rups_mag[i])

    duration = time.time() - t0
    print("\n\nIt took ",round(duration,3),"seconds to find all the ruptures.")

    return rup, [rups_length, rups_mag]


def write_rupt_file(root, rup, Run_Name, Set_Name):

    f_name = os.path.join(root, 'ruptures.txt')
    f = open(f_name, 'w')
    f.write("set "+Set_Name+"\n")
    for rup_i in rup:
        if len(rup_i) > 1:
            line = ""
            for si in rup_i:
                line += str(si)+" "
            line = line[:-1]
            f.write(line+"\n")
    f.close()

    print("Rupture file built")

def visu_rup(f_for_sherifs, rup, rups_length, rups_mag, path):

    north_shift = 0.08
    for s in range(len(f_for_sherifs)):
        if len(f_for_sherifs[s]["rup_id"]) > 3. :
            features = []
            z = 0
            id_rup = 0
            for rup_i in rup :
                add = False
                for si in rup_i :
                    if str(si) == str(s):
                        add = True
                if add == True :
                    geom = []
                    for si in rup_i:
                        for lon_i,lat_i in zip(f_for_sherifs[si]["lons"],f_for_sherifs[si]["lats"]):
                            geom.append([lon_i,lat_i+z*north_shift])
                    geom = MultiPoint(geom)
                    features.append(Feature(geometry=geom, properties={"id_rup": id_rup,
                                                                      "length":rups_length[id_rup],
                                                                      "mag":rups_mag[id_rup]}))
                    z+=1
                id_rup+=1
            feature_collection = FeatureCollection(features)

            if not os.path.exists(path+'/qgis/rups'):
                os.makedirs(path+'/qgis/rups')
            with open(path+'/qgis/rups/rup_'+str(s)+'.geojson', 'w') as f:
               dump(feature_collection, f)
