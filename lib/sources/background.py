# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np

class bg():
    """
    Extract the geometry and properities of the background.
    """
    def geom(model_name,file_geom):
        Lon_bg = []
        Lat_bg = []

        # manually defined  in the file Background geometry
        geom_bg = np.genfromtxt(file_geom,dtype=[('U100'),('f8'),('f8')],skip_header = 1)
        
        column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
        index_model = np.where(np.array(column_model) == model_name)[0]
        Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
        Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))
        
        
        if Lon_bg == 0:
            print('Error!! Check your input background geometry')
        
        return Lon_bg, Lat_bg
            

    def prop(model_name,file_prop):
        prop_bg = open(file_prop,'r').readlines()
        # background general parameters read from the input file
        nodalPlanes = []
        hypoDepths = []
        for line in prop_bg :
            if line.split('\t')[0] == model_name:
                if line.split('\t')[1] == 'upperSeismoDepth':
                    upperSeismoDepth = float(line.split('\t')[2])
                if line.split('\t')[1] == 'lowerSeismoDepth':
                    lowerSeismoDepth = float(line.split('\t')[2])
                if line.split('\t')[1] == 'ruptAspectRatio':
                    ruptAspectRatio = float(line.split('\t')[2])
                if line.split('\t')[1] == 'nodalPlane':
                    nodalPlanes.append([float(line.split('\t')[2]),
                                                   float(line.split('\t')[3]),
                                                         float(line.split('\t')[4]),
                                                               float(line.split('\t')[5])])
                if line.split('\t')[1] == 'hypoDepth':
                    hypoDepths.append([float(line.split('\t')[2]),
                                                   float(line.split('\t')[3])])
        if len(str(nodalPlanes))==0:
            print('Error!! Verify your Background parameters file')

        return upperSeismoDepth, lowerSeismoDepth, ruptAspectRatio, nodalPlanes, hypoDepths
