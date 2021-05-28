# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.4

@author: Thomas Chartier
"""

import numpy as np
import matplotlib.pyplot as plt

def add_data(x,y):
    y_cum = list(np.cumsum(y.reversed()))
    y_cum = y_cum.reversed()
    plt.scatter(x,y_cum, c='grey', s=50, edgecolor='',marker = 'o',alpha = 0.7)
    plt.plot(x,y,'grey', linewidth = 1 ,alpha = 0.8)


def plot(x,y,lim,axis,data,path,title):
    '''
    x : list, bining in magnitude
    y : list, mfd values (same length as x)
    lim : list, 2D of len 2, [[xmin,xmax],[ymin,ymax]]
    axis : list, 2 strings for the axis title x and y
    data : bool or list, False or mfd value of the catalog
    path : bool, destination path of the figure
    title : str, title of the figure

    '''

    y_cum = list(np.cumsum(y.reversed()))
    y_cum = y_cum.reversed()
    plt.scatter(x,y_cum, c='darkcyan', s=50, edgecolor='',marker = 's',alpha = 0.7)
    plt.plot(x,y,'darkgreen', linewidth = 2 ,alpha = 0.8)
    if not data == False :
        add_data(x,data)
    plt.yscale('log')
    plt.title(title)
    plt.grid()
    plt.savefig(path,dpi = 180, transparent=True)
    plt.close()
