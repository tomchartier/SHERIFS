# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap


def draw_screen_poly(lons, lats,  m , color , op, linewidth, edgecolor):
    x, y = m( lons, lats )
    xy = list(zip(x,y))
    poly = Polygon( xy, facecolor=color, alpha= op , linewidth=linewidth, edgecolor = edgecolor)
    plt.gca().add_patch(poly)



def make_fault_map(m,fault_geom,fault_colors,figpath,title,Mmax=0.,dpi=80,use_arcgis=False):
    for index_fault in range(len(fault_geom.items())):
        if fault_geom[index_fault]['plot_trace'] == True :
            x, y = m(fault_geom[index_fault]['trace_lon'], fault_geom[index_fault]['trace_lat'])
            m.plot(x, y, '-', markersize=1, linewidth=1, color=fault_colors[index_fault], markerfacecolor=fault_colors[index_fault], markeredgewidth = 1)
        draw_screen_poly(fault_geom[index_fault]['polygon'][0], fault_geom[index_fault]['polygon'][1],  m ,fault_colors[index_fault] , 0.3, 0.5, fault_colors[index_fault])
        
    m.drawcoastlines(linewidth=0.1)
    if use_arcgis == True :
        try:
            m.arcgisimage(service='World_Shaded_Relief', dpi = 400, alpha = 0.3, xpixels = 2000)
        except:
            m.fillcontinents(color='sienna',lake_color='w',alpha = 0.05)
    else :
        m.fillcontinents(color='grey',lake_color='w',alpha = 0.2)
         
    if Mmax != 0. :
        plt.annotate('Mmax : '+str(Mmax), xy=(0.1, 0.9), xycoords='axes fraction',size=6)

    plt.title(str(title))
    plt.savefig(figpath,dpi=dpi)

    plt.close()
