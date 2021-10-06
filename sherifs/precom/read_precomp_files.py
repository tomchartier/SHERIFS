# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
import geojson


def read_oiler_file(oiler_json):
    with open(oiler_json) as f:
        gj = geojson.load(f)
    faults = gj['features']
    nb_faults = len(faults)
    print("There are ",nb_faults," faults in this model.")

    # TODO option that can be useful in the future
    sub_selec = False
    app = ""
    if sub_selec == True :
        faults = gj['features'][333:]
        nb_faults = len(faults)
        app = "_subsample"
        print("But only ",nb_faults,"are selected.")
    
    return faults,nb_faults
