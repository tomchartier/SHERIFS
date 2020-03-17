# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import numpy as np
from math import pi, cos, radians , sin, asin, sqrt, atan2, degrees

def reproject(latitude, longitude):
    """Returns the x & y coordinates in meters using a sinusoidal projection"""
    earth_radius = 6371009 # in meters
    lat_dist = pi * earth_radius / 180.0
    y = [lat * lat_dist for lat in latitude]
    x = [long * lat_dist * cos(radians(lat))
                for lat, long in zip(latitude, longitude)]
    return x, y



def points_aligned(a, b, c):
    crossproduct = (c[1] - a[1]) * (b[0] -a[0]) - (c[0] -a[0]) * (b[1] - a[1])
    epsilon = 10000000.
    if abs(crossproduct) > epsilon:
        return False
    dotproduct = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1])*(b[1] - a[1])
    if dotproduct < 0:
        return False
    
    squaredlengthba = (b[0]-a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1]-a[1])
    if dotproduct > squaredlengthba:
        return False
    
    return True



    
  
def area_of_polygon(x, y):
    """Calculates the area of an arbitrary polygon given its verticies"""
    # first, does the list of vertices
    x_vertices = []
    y_vertices = []
    inn = []
    for i in range(len(x)):
        if i == 0:
            if points_aligned([x[-1],y[-1]],[x[0],y[0]], [x[1],y[1]]) == False:
                 x_vertices.append(x[i])
                 y_vertices.append(y[i])
                 inn.append(1)
            else:
                inn.append(0)
        elif i == len(x)-1:
            if points_aligned([x[-2],y[-2]], [x[-1],y[-1]],[x[0],y[0]]) == False:
                 x_vertices.append(x[i])
                 y_vertices.append(y[i])
                 inn.append(1)
            else:
                 inn.append(0)
        else:
             if points_aligned([x[i-1],y[i-1]], [x[i],y[i]],[x[i+1],y[i+1]]) == False:
                 x_vertices.append(x[i])
                 y_vertices.append(y[i])
                 inn.append(1)
             else:
                 inn.append(0)
    #print len(x),len(x_vertices),inn
    
    area = 0.0
    for i in range(-1, len(x_vertices)-1):
        area += x_vertices[i] * (y_vertices[i+1] - y_vertices[i-1])
    return abs(area) / 2.0

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       #raise Exception('lines do not intersect')
       x = 'no_intesection'
       y = 'no_intesection'
    else:
        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
    return x, y


                          
def calculate_initial_compass_bearing(pointA, pointB):
    """
    Calculates the bearing between two points.

    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))

    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees

    :Returns:
      The bearing in degrees

    :Returns Type:
      float
    """
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = radians(pointA[0])
    lat2 = radians(pointB[0])

    diffLong = radians(pointB[1] - pointA[1])

    x = sin(diffLong) * cos(lat2)
    y = cos(lat1) * sin(lat2) - (sin(lat1)
            * cos(lat2) * cos(diffLong))

    initial_bearing = atan2(x, y)

    # Now we have the initial bearing but atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def distance(lon1, lat1, lon2, lat2):
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

