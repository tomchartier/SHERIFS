# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.0 

@author: Thomas Chartier
"""

import numpy as np
import math
import sys
from math import radians, cos, sin, asin, sqrt

class Geom_scenar:
    def __init__(self,faults_names,File_geom,Model_name):
        self.faults_names = faults_names
        #self.scenarios = scenarios
        self.File_geom = File_geom
        self.Model_name = Model_name
        
        self.initialize()
    def initialize(self):
        faults_names = self.faults_names
        self.FaultGeometry() # extract the faults geometries
        length = []
        names = []        
        #####################################################################
        #find de geometry of each fault in the scenario
        #####################################################################
        faults_lon_ini = np.zeros((len(self.faults_names),1000))
        faults_lat_ini = np.zeros((len(self.faults_names),1000))        

        index_fault = 0
        for fault in faults_names:
            index = np.where(np.array(self.Column_Fault_name) == fault)
            lon_fault = np.array(list(map(lambda i : self.Longitudes[i],index[0])))
            lon_fault.resize(faults_lon_ini[index_fault][:].shape)
            lat_fault = np.array(list(map(lambda i : self.Latitudes[i],index[0])))
            lat_fault.resize(faults_lat_ini[index_fault][:].shape)
            faults_lon_ini[index_fault][:] = faults_lon_ini[index_fault][:] + lon_fault
            faults_lat_ini[index_fault][:] = faults_lat_ini[index_fault][:] + lat_fault
                    
            
            index_fault = index_fault +1
        #rearranging the arrays to delete the zeros
        faults_lon = []
        faults_lat = []
        index_fault = 0
        tmp_f_names = []
        for fault in faults_names:
            index_zeros = np.where((faults_lat_ini[index_fault] == 0.)&(faults_lon_ini[index_fault] == 0.))  
            tmp = np.array(faults_lon_ini[index_fault])
            tmp = np.delete(tmp,index_zeros)
            faults_lon.append(tmp)
            tmp = np.array(faults_lat_ini[index_fault])
            tmp = np.delete(tmp,index_zeros)
            faults_lat.append(tmp)
            index = np.where(np.array(self.Column_Fault_name) == fault)[0]
            depths_fault = np.take(self.Depths,index)
            #calculating the length of the fault
            dist = 0
            try :
                depths_fault[0]
            except:
                sys.exit("Error with the geometry of fault :"+fault+"\n"
                +"Please define a geometry for all faults in the model")
            if depths_fault[0] ==  'sf':
                for i in range(len(faults_lon[index_fault])-1):
                    dist += self.distance(faults_lon[index_fault][i],faults_lat[index_fault][i],faults_lon[index_fault][i+1],faults_lat[index_fault][i+1])
            else: 
                for i in range(len(faults_lon[index_fault])-1):
                    if depths_fault[i+1] == depths_fault[0]:
                        dist += self.distance(faults_lon[index_fault][i],faults_lat[index_fault][i],faults_lon[index_fault][i+1],faults_lat[index_fault][i+1])
             
            length.append(dist)
            names.append(fault)
            
            index_fault += 1
            
        self.length = length
        self.faults_lon = faults_lon
        self.faults_lat = faults_lat

                
    def FaultGeometry(self): 
        #CritereDistance = 3. 
        NomFichier_InfosZonage = self.File_geom
        InfosZonage = np.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8'),('U100')],skip_header = 1)
        Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
        index_model = np.where(np.array(Column_model_name) == self.Model_name)
        self.Column_Fault_name = list(map(lambda i : InfosZonage[i][1],index_model[0]))
        self.Longitudes = list(map(lambda i : InfosZonage[i][2],index_model[0]))
        self.Latitudes = list(map(lambda i : InfosZonage[i][3],index_model[0]))
        self.Depths = list(map(lambda i : InfosZonage[i][4],index_model[0]))
        
        ZoneSelec = self.Column_Fault_name
        DicoZone = dict([(k,ZoneSelec.count(k)) for k in set(ZoneSelec)])
        Longitudes = []
        Latitudes = []
        Depths = []
        Column_Fault_name = []
        for cle in DicoZone.keys():
            indices_ZonesSelec = np.where(np.array(self.Column_Fault_name) == cle)
            ColonneNomZone_inter = np.take(self.Column_Fault_name,indices_ZonesSelec)
            Longitudes_inter = np.take(self.Longitudes,indices_ZonesSelec)
            Latitudes_inter = np.take(self.Latitudes,indices_ZonesSelec)
            depth_inter = np.take(self.Depths,indices_ZonesSelec)
        
            Longitudes_inter = Longitudes_inter[0].tolist()
            Latitudes_inter = Latitudes_inter[0].tolist()
            depth_inter = depth_inter[0].tolist()
            ColonneNomZone_inter = ColonneNomZone_inter[0].tolist()
            compt = 0
            for xx,yy,nn,dd in zip(Longitudes_inter,Latitudes_inter,ColonneNomZone_inter,depth_inter):
                compt+=1
                Longitudes.append(xx)
                Latitudes.append(yy)
                Depths.append(dd)
                Column_Fault_name.append(nn)
                
        self.Longitudes =Longitudes
        self.Latitudes =Latitudes
        self.Depths =Depths
        self.Column_Fault_name = Column_Fault_name
        self.Nb_data_per_zone = dict([(k,self.Column_Fault_name.count(k)) for k in set(self.Column_Fault_name)])

    def distance(self,lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points 
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians 
        lon1, lat1, lon2, lat2 = list(map(radians, [lon1, lat1, lon2, lat2]))
        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 
        km = 6367 * c
        return km        

