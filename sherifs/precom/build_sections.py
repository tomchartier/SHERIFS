# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

@author: Thomas Chartier
"""
import os
import geojson
import numpy as np
import matplotlib.pyplot as plt
from geojson import LineString, Feature, FeatureCollection, dump
from math import radians, sin, degrees, acos
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from sherifs.precom.precomp_utils import (distance, fault_length,
                                          calculate_initial_compass_bearing)


def cut_faults(faults, f_lengths, f_areas, path, rupture_mesh_spacing):

    # for all faults longer than 50 km cut them either in roughly 25 km long sections
    # or in 5 sections if the faults is longer than 25 x 5 = 125 km
    f_for_sherifs = {}
    f_lons = []
    f_lats = []
    f_id = -1
    sections_lengths_tot = []
    sections_areas_tot = []
    id_sections_fault = []

    if not os.path.exists(path+"/qgis"):
        os.makedirs(path+"/qgis")
    file_section_tips = open(path+"/qgis/sect_tips.csv",'w')
    file_section_tips.write("lon,lat\n")

    # check if the smallest fault is large enough for discretization
    error_msg = "A fault is too small for the discretization. Min length :"
    error_msg += str(round(min(f_lengths),3))
    error_msg += " km."
    assert (min(f_lengths)> 2. * rupture_mesh_spacing), error_msg

    nb_faults = len(faults)
    for fi in range(nb_faults):
        lons_tot = [i[0] for i in faults[fi]['geometry']["coordinates"]]
        lats_tot = [i[1] for i in faults[fi]['geometry']["coordinates"]]
        id_sections_fault_i = []

        # if the number of points is too small for the length, we interpolate new points
        if len(lons_tot) < f_lengths[fi]/3.:
            inter_lon_i, inter_lat_i = [], []

            for i in range(len(lons_tot)-1):
                di = distance(lons_tot[i],lats_tot[i],lons_tot[i+1],lats_tot[i+1])
                nb_extra_pts = int(di/3.)
                if i == len(lons_tot)-2 :
                    if nb_extra_pts > 1:
                        inter_lon_i += list(np.linspace(lons_tot[i],
                                                        lons_tot[i+1],
                                                        nb_extra_pts))
                        inter_lat_i += list(np.linspace(lats_tot[i],
                                                        lats_tot[i+1],
                                                        nb_extra_pts))
                    else :
                        inter_lon_i += [lons_tot[i],lons_tot[i+1]]
                        inter_lat_i += [lats_tot[i],lats_tot[i+1]]
                else :
                    if nb_extra_pts > 1:
                        inter_lon_i += list(np.linspace(lons_tot[i],
                                                        lons_tot[i+1],
                                                        nb_extra_pts))[:-1]
                        inter_lat_i += list(np.linspace(lats_tot[i],
                                                        lats_tot[i+1],
                                                        nb_extra_pts))[:-1]
                    else :
                        inter_lon_i += [lons_tot[i],lons_tot[i+1]][:-1]
                        inter_lat_i += [lats_tot[i],lats_tot[i+1]][:-1]
            lons_tot, lats_tot = inter_lon_i, inter_lat_i


        if f_lengths[fi] < 35. : #TODO clean up hard variables
            f_name = str(fi)
            f_id += 1
            f_for_sherifs.update({f_id:{'f_name' : f_name,
                                 'lons' : lons_tot,
                                 'lats' : lats_tot,
                                 'oiler_id' : fi,
                                'oiler_name' : faults[fi]['properties']["name"],
                                    'oiler_fid' : faults[fi]['properties']["fid"],
                                 'length' : f_lengths[fi]}})

            file_section_tips.write(str(lons_tot[0])+","+str(lats_tot[0])+"\n")
            file_section_tips.write(str(lons_tot[-1])+","+str(lats_tot[-1])+"\n")

            sections_lengths_tot.append(f_lengths[fi])
            sections_areas_tot.append(f_areas[fi])
            id_sections_fault_i.append(f_id)
        else :
            # Find maximal distance between two points in the fault
            dists = []
            for i in range(len(lons_tot)-1):
                dists.append(distance(lons_tot[i],lats_tot[i],
                                     lons_tot[i+1],lats_tot[i+1]))
            if f_lengths[fi] < 125. :
                if max(dists) > 25. :
                    section_lenght = max(dists)
                    nb_sections = int(round(f_lengths[fi]/section_lenght))
                else :
                    nb_sections = 1
                    section_lenght = f_lengths[fi]/float(nb_sections)
                    while section_lenght > 25.  :
                        nb_sections += 1
                        section_lenght = f_lengths[fi]/float(nb_sections)

            else :
                if max(dists) > f_lengths[fi]/5. :
                    nb_sections = 5
                else :
                    nb_sections = 6

            if nb_sections != 1 :
                # Find cutting points
                av_section_len = f_lengths[fi]/float(nb_sections)
                index_0 = 0
                for cut in range(nb_sections-1):
                    dist = 0 #if dist is large, force cut
                    index_cut = 0+index_0
                    while dist < av_section_len :
                        try :
                            dist += distance(lons_tot[index_cut],lats_tot[index_cut],
                                             lons_tot[index_cut+1],lats_tot[index_cut+1])
                            index_cut += 1
                            do_last = True

                            # check for major azimut change
                            if dist > av_section_len * 0.3 and dist > 2. * rupture_mesh_spacing:
                                az1 = calculate_initial_compass_bearing([lons_tot[index_cut-1],lats_tot[index_cut-1]], [lons_tot[index_cut],lats_tot[index_cut]])
                                az2 = calculate_initial_compass_bearing([lons_tot[index_cut],lats_tot[index_cut]], [lons_tot[index_cut+1],lats_tot[index_cut+1]])
                                if abs(az1-az2) % 360 > 60. :
                                    dist = 100000.
                        except IndexError:
                            index_cut = len(lons_tot)-1
                            do_last = False
                            dist = 100000.

                    # feed the dico
                    lons = lons_tot[index_0:index_cut+1]
                    lats = lats_tot[index_0:index_cut+1]
                    length = 0.
                    for i in range(len(lons)-1):
                        length += distance(lons[i],lats[i],lons[i+1],lats[i+1])

                    width =((faults[fi]['properties']['lsd']-
                              faults[fi]['properties']['usd'])/
                              sin(radians(faults[fi]['properties']['dip'])))
                    area = length*width
                    f_name = str(fi)+"_"+str(cut)
                    f_id += 1
                    f_for_sherifs.update({f_id:{'f_name' : f_name,
                                         'lons' : lons,
                                         'lats' : lats,
                                         'oiler_id' : fi,
                                        'oiler_name' : faults[fi]['properties']["name"],
                                    'oiler_fid' : faults[fi]['properties']["fid"],
                                         'length' : length}})

                    file_section_tips.write(str(lons[0])+","+str(lats[0])+"\n")
                    file_section_tips.write(str(lons[-1])+","+str(lats[-1])+"\n")

                    index_0 = index_cut
                    sections_lengths_tot.append(length)
                    sections_areas_tot.append(area)
                    id_sections_fault_i.append(f_id)
                #for the last one :
                if do_last == True:
                    lons = lons_tot[index_cut:]
                    lats = lats_tot[index_cut:]
                    f_name = str(fi)+"_"+str(cut)
                    f_id += 1
                    length = fault_length(lons,lats)
                    width =((faults[fi]['properties']['lsd']-
                              faults[fi]['properties']['usd'])/
                              sin(radians(faults[fi]['properties']['dip'])))

                    # check if the width is large enought for the discretization
                    error_msg = "A fault is too thin for the discretization. Fault name :"
                    error_msg += str(faults[fi]['properties']["fid"])
                    error_msg += " width : "
                    error_msg += str(round(width,2))
                    error_msg += " km."
                    assert (width > 2. * rupture_mesh_spacing), error_msg

                    area = length*width
        #             length = 0.
        #             for i in range(len(lons)-1):
        #                 length += distance(lons[i],lats[i],lons[i+1],lats[i+1])
                    f_for_sherifs.update({f_id:{'f_name' : f_name,
                                         'lons' : lons,
                                         'lats' : lats,
                                         'oiler_id' : fi,
                                        'oiler_name' : faults[fi]['properties']["name"],
                                    'oiler_fid' : faults[fi]['properties']["fid"],
                                         'length' : length}})

                    file_section_tips.write(str(lons[0])+","+str(lats[0])+"\n")
                    file_section_tips.write(str(lons[-1])+","+str(lats[-1])+"\n")

                    sections_lengths_tot.append(length)
                    sections_areas_tot.append(area)
                    id_sections_fault_i.append(f_id)

            else : # in this case, there is only one section
                f_name = str(fi)
                f_id += 1
                f_for_sherifs.update({f_id:{'f_name' : f_name,
                                     'lons' : lons_tot,
                                     'lats' : lats_tot,
                                     'oiler_id' : fi,
                                    'oiler_name' : faults[fi]['properties']["name"],
                                    'oiler_fid' : faults[fi]['properties']["fid"],
                                     'length' : f_lengths[fi]}})

                file_section_tips.write(str(lons_tot[0])+","+str(lats_tot[0])+"\n")
                file_section_tips.write(str(lons_tot[-1])+","+str(lats_tot[-1])+"\n")

                sections_lengths_tot.append(f_lengths[fi])
                sections_areas_tot.append(f_areas[fi])
                id_sections_fault_i.append(f_id)

            #print("#####\n fault",fi)
            if round(sum(sections_lengths_tot[-nb_sections:])) != round(f_lengths[fi]):
                print("#####\n fault",fi)
                print("Failing!")
                print(round(sum(sections_lengths_tot[-nb_sections:])), round(f_lengths[fi]), "missing", round(-sum(sections_lengths_tot[-nb_sections:]))+ round(f_lengths[fi]))
                #print("len",len([i[0] for i in faults[fi]['geometry']["coordinates"]]),len(lons_tot))
                print("sections:", nb_sections)
    #         else :
    #             print("is okay. sections:", nb_sections)
        id_sections_fault.append(id_sections_fault_i)

    plt_fig = False
    if plt_fig == True :
        plt.hist(sections_lengths_tot)
        plt.ylabel("nb of sections")
        plt.xlabel("Length of the section (km)")
        plt.show()

    file_section_tips.close()

    print("And a total of ",len(f_for_sherifs)," sections.")

    if round(sum(sections_lengths_tot)) != round(sum(f_lengths)):
        print("!!!!!\n ERROR! There are ",round(sum(f_lengths)-sum(sections_lengths_tot)), " km of faults missing.")

    return f_for_sherifs,id_sections_fault,sections_areas_tot,sections_lengths_tot

def export_sections_pts(f_for_sherifs,path):
    # Export the section points
    sections_points = open(path+"/qgis/sections_pt.csv",'w')
    sections_points.write('lon,lat'+'\n')
    for si in range(len(f_for_sherifs)):
        lons_i = f_for_sherifs[si]["lons"]
        lats_i = f_for_sherifs[si]["lats"]
        for lon_i, lat_i in zip(lons_i,lats_i):
            sections_points.write(str(lon_i)+','+str(lat_i)+'\n')
    sections_points.close()

    return

def find_sections_Mmax(f_for_sherifs,zones_json):
    # Find the max length for each section
    nb_sections = len(f_for_sherifs)
    for si in range(nb_sections):
        f_for_sherifs[si]["max_possible_length"] = 10000.
        f_for_sherifs[si]["max_possible_Mmax"] = 11.

    # Load the zones
    with open(zones_json) as f:
        gj = geojson.load(f)
    zones_length = gj['features']

    # Loop on the zone to find the sections within
    for zone_i in zones_length:
        poly = []
        for pt in zone_i["geometry"]["coordinates"][0][0]:
            poly.append((pt[0],pt[1]))
        polygon = Polygon(poly)
        max_i = zone_i["properties"]["max_length"]
        Mmax_i = zone_i["properties"]["max_possible_Mmax"]
        for si in range(nb_sections):
            if f_for_sherifs[si]["max_possible_length"]>max_i:
                lons_si = f_for_sherifs[si]['lons']
                lats_si = f_for_sherifs[si]['lats']
                for lon_i,lat_i in zip(lons_si,lats_si):
                    if polygon.contains(Point(lon_i, lat_i)):
                        f_for_sherifs[si]["max_possible_length"]=max_i

        for si in range(nb_sections):
            if f_for_sherifs[si]["max_possible_Mmax"]>Mmax_i:
                lons_si = f_for_sherifs[si]['lons']
                lats_si = f_for_sherifs[si]['lats']
                for lon_i,lat_i in zip(lons_si,lats_si):
                    if polygon.contains(Point(lon_i, lat_i)):
                        f_for_sherifs[si]["max_possible_Mmax"]=Mmax_i

    return f_for_sherifs

def to_sherifs(f_for_sherifs,faults,Model_name,apply_sr_reduction,f_mu):
    # Attribute the parameters usable by SHERIFS for each section
    nb_sections = len(f_for_sherifs)

    # if the shear modulus is modified for some faults
    modify_mu = []
    mu_value = []
    if f_mu != None:
        # Load the zones
        with open(f_mu) as f:
            gj = geojson.load(f)
        zones_mu = gj['features']
        # Loop on the zone to find the sections within
        for zone_i in zones_mu:
            poly = []
            for pt in zone_i["geometry"]["coordinates"][0][0]:
                poly.append((pt[0],pt[1]))
            polygon = Polygon(poly)
            mu_i = zone_i["properties"]["shear_mod"]
            for si in range(nb_sections):
                lons_si = f_for_sherifs[si]['lons']
                lats_si = f_for_sherifs[si]['lats']
                isin = False
                for lon_i,lat_i in zip(lons_si,lats_si):
                    if isin == False and polygon.contains(Point(lon_i, lat_i)):
                        isin = True
                        modify_mu.append(si)
                        mu_value.append(mu_i)

    i_mu = 0
    # top - bottom - dip - orientation - defaults
    for si in range(nb_sections):
        fi = f_for_sherifs[si]["oiler_id"]
        f_for_sherifs[si]["up_s_d"] = faults[fi]['properties']["usd"]
        f_for_sherifs[si]["lo_s_d"] = faults[fi]['properties']["lsd"]
        f_for_sherifs[si]["dip"] = faults[fi]['properties']["dip"]
        f_for_sherifs[si]["oriented"] = faults[fi]['properties']["dip_dir"]

        # default for now but can be changed later
        f_for_sherifs[si]["Domain"] = "Active Shallow Crust"

        if si in modify_mu:
            f_for_sherifs[si]["shear_modulus"] = mu_value[i_mu]
            i_mu += 1
        else :
            f_for_sherifs[si]["shear_modulus"] = "30"

    # slip rate and mechanism
    for si in range(nb_sections):
        fi = f_for_sherifs[si]["oiler_id"]
        if "v_rl" in faults[fi]['properties'].keys():
            #using the format of the GEM GFDB
            v_rl = abs(faults[fi]['properties']["v_rl"])
            v_ex = faults[fi]['properties']["v_ex"]
            slip_rate_moy = (v_rl**2 + v_ex**2)**0.5
            err = (faults[fi]['properties']["e_rl"]**2 + faults[fi]['properties']["e_ex"]**2)**0.5
            slip_rate_min = slip_rate_moy - err
            slip_rate_max = slip_rate_moy + err
            if slip_rate_min < 0. :
                print("ERROR - slip rate min is negative !")

            if v_ex != 0. :
                if abs(v_ex) < v_rl :
                    rake = degrees(acos(abs(v_ex) / v_rl))
                else :
                    rake = degrees(acos(v_rl / abs(v_ex)))
                if v_ex < 0 :
                    rake = -rake
            else :
                rake = 0.
        elif "sr" in faults[fi]['properties'].keys():
            #using slip-rate direclty
            slip_rate_moy = abs(faults[fi]['properties']["sr"])
            err = abs(faults[fi]['properties']["e_sr"])
            slip_rate_min = slip_rate_moy - err
            slip_rate_max = slip_rate_moy + err
            rake = abs(faults[fi]['properties']["rake"])

        else :
            print("Please check the format of the input geojson file")



        f_for_sherifs[si]["rake"] = rake

        f_for_sherifs[si]["slip_rate_min"] = slip_rate_min * (1.-apply_sr_reduction)
        f_for_sherifs[si]["slip_rate_moy"] = slip_rate_moy * (1.-apply_sr_reduction)
        f_for_sherifs[si]["slip_rate_max"] = slip_rate_max * (1.-apply_sr_reduction)


        f_for_sherifs[si]["model"] = Model_name

    return f_for_sherifs

def write_section_json(f_for_sherifs,f_name):
    # Write in a Geojson

    features = []

    for si in range(len(f_for_sherifs)):
        geom = []
        for lon_i,lat_i in zip(f_for_sherifs[si]["lons"],f_for_sherifs[si]["lats"]):
            geom.append([lon_i,lat_i])
        geom = LineString(geom)
        nb_rup = len(f_for_sherifs[si]["rup_id"])
        max_rup_length = f_for_sherifs[si]["max_length"]
        up_s_d = f_for_sherifs[si]["up_s_d"]
        lo_s_d = f_for_sherifs[si]["lo_s_d"]
        dip = f_for_sherifs[si]["dip"]
        oriented = f_for_sherifs[si]["oriented"]
        Domain = f_for_sherifs[si]["Domain"]
        shear_modulus = f_for_sherifs[si]["shear_modulus"]
        rake = f_for_sherifs[si]["rake"]
        sr_min = f_for_sherifs[si]["slip_rate_min"]
        sr_mean = f_for_sherifs[si]["slip_rate_moy"]
        sr_max = f_for_sherifs[si]["slip_rate_max"]

        model = f_for_sherifs[si]["model"]

        features.append(Feature(geometry=geom, properties={"si": si,
                                                           "nb_rup": nb_rup,
                                                           "max_rup_length":max_rup_length,
                                                          "up_s_d":up_s_d,
                                                          "lo_s_d":lo_s_d,
                                                          "dip":dip,
                                                          "oriented":oriented,
                                                          "Domain":Domain,
                                                          "shear_modulus":shear_modulus,
                                                          "rake":rake,
                                                          "sr_min":sr_min,
                                                          "sr_mean":sr_mean,
                                                          "sr_max":sr_max,
                                                          "model":model}))


    feature_collection = FeatureCollection(features)

    with open(f_name, 'w') as f:
       dump(feature_collection, f)
