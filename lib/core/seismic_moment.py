# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np

def mag_to_M0(mag):
    #returns Mo for a given mag
    M0 = 10. ** (1.5 * mag + 9.1)
    return M0
