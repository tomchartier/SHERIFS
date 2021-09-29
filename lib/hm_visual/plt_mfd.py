# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.4

@author: Thomas Chartier
"""

import numpy as np
import matplotlib.pyplot as plt

def add_data(data):
    x = data[0]
    y = data[1]
    #plt.scatter(x,y, c='grey', s=50, edgecolor='',marker = 'o',alpha = 0.7)
    plt.plot(x,y,'k', linewidth = 3 ,alpha = 0.8, label='cumulative catalogue MFD')


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
    y = list(y)
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.scatter(x,y, c='darkcyan', s=50, marker = 's',alpha = 0.7, label='SHERIFS incremental MFD')
    plt.plot(x,y_cum,'darkgreen', linewidth = 2 ,alpha = 0.8, label='SHERIFS cumulative MFD')
    if not data == False :
        add_data(data)
    plt.yscale('log')
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.savefig(path,dpi = 180, transparent=False)
    plt.close()

def plot_bg_ft(x,ys,lim,axis,path,title):
    '''
    x : list, bining in magnitude
    y : list, mfd values for whole model, faults, bg (3 x same length as x)
    lim : list, 2D of len 2, [[xmin,xmax],[ymin,ymax]]
    axis : list, 2 strings for the axis title x and y
    data : bool or list, False or mfd value of the catalog
    path : bool, destination path of the figure
    title : str, title of the figure

    '''
    y = list(ys[0])
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.plot(x,y_cum,'darkgreen', linewidth = 2 ,
    alpha = 0.8, label='SHERIFS cumulative MFD')


    y = list(ys[1])
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.plot(x,y_cum,'red', linewidth = 2 ,
    alpha = 0.8, label='Faults cumulative MFD')


    y = list(ys[2])
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.plot(x,y_cum,'blue', linewidth = 2 ,
    alpha = 0.8, label='Background cumulative MFD')

    plt.yscale('log')
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.savefig(path,dpi = 180, transparent=False)
    plt.close()


def local(x,ys,data,lim,axis,path,title):
    '''
    x : list, bining in magnitude
    y : list, mfd values for whole model, faults, bg
                and smooth before scalling (4 x same length as x)
    lim : list, 2D of len 2, [[xmin,xmax],[ymin,ymax]]
    axis : list, 2 strings for the axis title x and y
    data : bool or list, False or mfd value of the catalog
    path : bool, destination path of the figure
    title : str, title of the figure

    '''
    y = list(ys[0])
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.plot(x,y_cum,'darkgreen', linewidth = 2 ,
    alpha = 0.8, label='SHERIFS cumulative MFD')


    y = list(ys[1])
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.plot(x,y_cum,'red', linewidth = 2 ,
    alpha = 0.8, label='Faults cumulative MFD')


    y = list(ys[2])
    y_cum = list(np.cumsum(np.array(y[::-1])))
    y_cum = y_cum[::-1]
    plt.plot(x,y_cum,'blue', linewidth = 2 ,
    alpha = 0.8, label='Background cumulative MFD')

    # y = list(ys[3])
    # y_cum = list(np.cumsum(np.array(y[::-1])))
    # y_cum = y_cum[::-1]
    # plt.plot(x,y_cum,'blue', linewidth = 0.5 ,
    # alpha = 0.8, label='smoothing MFD (preSHERIFS)')

    if not data == False :
        add_data(data)

    plt.yscale('log')
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.savefig(path,dpi = 180, transparent=False)
    plt.close()
