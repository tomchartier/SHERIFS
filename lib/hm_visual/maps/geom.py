# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

These functions are all from the OpenQuake engine
they have been modified for usage in SHERIFS

Version 1.2

@author: Thomas Chartier
"""
import numpy, math

EARTH_RADIUS = 6371.0

def _prepare_coords(lons1, lats1, lons2, lats2):
    """
    Convert two pairs of spherical coordinates in decimal degrees
    to numpy arrays of radians. Makes sure that respective coordinates
    in pairs have the same shape.
    """
    lons1 = numpy.radians(lons1)
    lats1 = numpy.radians(lats1)
    assert lons1.shape == lats1.shape
    lons2 = numpy.radians(lons2)
    lats2 = numpy.radians(lats2)
    assert lons2.shape == lats2.shape
    return lons1, lats1, lons2, lats2
    
def geodetic_distance(lons1, lats1, lons2, lats2, diameter=2*EARTH_RADIUS):
    """
    Calculate the geodetic distance between two points or two collections
    of points.

    Parameters are coordinates in decimal degrees. They could be scalar
    float numbers or numpy arrays, in which case they should "broadcast
    together".

    Implements http://williams.best.vwh.net/avform.htm#Dist

    :returns:
        Distance in km, floating point scalar or numpy array of such.
    """
    lons1, lats1, lons2, lats2 = _prepare_coords(lons1, lats1, lons2, lats2)
    distance = numpy.arcsin(numpy.sqrt(
        numpy.sin((lats1 - lats2) / 2.0) ** 2.0
        + numpy.cos(lats1) * numpy.cos(lats2)
        * numpy.sin((lons1 - lons2) / 2.0) ** 2.0
    ))
    return diameter * distance
    
    
def find_azimuth(lons1, lats1, lons2, lats2):
    """
    Calculate the azimuth between two points or two collections of points.

    Parameters are the same as for :func:`geodetic_distance`.

    Implements an "alternative formula" from
    http://williams.best.vwh.net/avform.htm#Crs

    :returns:
        Azimuth as an angle between direction to north from first point and
        direction to the second point measured clockwise in decimal degrees.
    """
    lons1, lats1, lons2, lats2 = _prepare_coords(lons1, lats1, lons2, lats2)
    cos_lat2 = numpy.cos(lats2)
    true_course = numpy.degrees(numpy.arctan2(
        numpy.sin(lons1 - lons2) * cos_lat2,
        numpy.cos(lats1) * numpy.sin(lats2)
        - numpy.sin(lats1) * cos_lat2 * numpy.cos(lons1 - lons2)
    ))
    return (360 - true_course) % 360

def average_azimuth(lons,lats):
    """
    Calculate and return weighted average azimuth of all line's segments
    in decimal degrees.

    Uses formula from
    http://en.wikipedia.org/wiki/Mean_of_circular_quantities

    >>> from openquake.hazardlib.geo.point import Point as P
    >>> '%.1f' % Line([P(0, 0), P(1e-5, 1e-5)]).average_azimuth()
    '45.0'
    >>> '%.1f' % Line([P(0, 0), P(0, 1e-5), P(1e-5, 1e-5)]).average_azimuth()
    '45.0'
    >>> line = Line([P(0, 0), P(-2e-5, 0), P(-2e-5, 1.154e-5)])
    >>> '%.1f' % line.average_azimuth()
    '300.0'
    """
#    if len(self.points) == 2:
#        return self.points[0].azimuth(self.points[1])
#    lons = numpy.array([point.longitude for point in self.points])
#    lats = numpy.array([point.latitude for point in self.points])
    lons = numpy.array(lons)
    lats = numpy.array(lats)
    
    azimuths = find_azimuth(lons[:-1], lats[:-1], lons[1:], lats[1:])
    distances = geodetic_distance(lons[:-1], lats[:-1],
                                           lons[1:], lats[1:])
    azimuths = numpy.radians(azimuths)
    # convert polar coordinates to Cartesian ones and calculate
    # the average coordinate of each component
    avg_x = numpy.mean(distances * numpy.sin(azimuths))
    avg_y = numpy.mean(distances * numpy.cos(azimuths))
    # find the mean azimuth from that mean vector
    azimuth = numpy.degrees(numpy.arctan2(avg_x, avg_y))
    if azimuth < 0:
        azimuth += 360
    return azimuth


def point_at(lon, lat, azimuth, distance):
    """
    Perform a forward geodetic transformation: find a point lying at a given
    distance from a given one on a great circle arc defined by azimuth.

    :param float lon, lat:
        Coordinates of a reference point, in decimal degrees.
    :param azimuth:
        An azimuth of a great circle arc of interest measured in a reference
        point in decimal degrees.
    :param distance:
        Distance to target point in km.
    :returns:
        Tuple of two float numbers: longitude and latitude of a target point
        in decimal degrees respectively.

    Implements the same approach as :func:`npoints_towards`.
    """
    # this is a simplified version of npoints_towards().
    # code duplication is justified by performance reasons.
    lon, lat = numpy.radians(lon), numpy.radians(lat)
    tc = numpy.radians(360 - azimuth)
    sin_dists = numpy.sin(distance / EARTH_RADIUS)
    cos_dists = numpy.cos(distance / EARTH_RADIUS)
    sin_lat = numpy.sin(lat)
    cos_lat = numpy.cos(lat)

    sin_lats = sin_lat * cos_dists + cos_lat * sin_dists * numpy.cos(tc)
    lats = numpy.degrees(numpy.arcsin(sin_lats))

    dlon = numpy.arctan2(numpy.sin(tc) * sin_dists * cos_lat,
                         cos_dists - sin_lat * sin_lats)
    lons = numpy.mod(lon - dlon + numpy.pi, 2 * numpy.pi) - numpy.pi
    lons = numpy.degrees(lons)

    return lons, lats

def orient_fault(lons,lats,orientated):

    # orienting the arrays in order to respect OQ right hand rule
    compass_bearing = calculate_initial_compass_bearing((lats[0],lons[0]),(lats[-1],lons[-1]))
    
    if str('N') in str(orientated):
        if compass_bearing < 180. :
            lons = reversed(lons)
            lats = reversed(lats)
    if str('S') in str(orientated):
        if compass_bearing > 180. :
            lons = reversed(lons)
            lats = reversed(lats)
    if str('E') in str(orientated):
        if compass_bearing > 90. and compass_bearing < 270. :
            lons = reversed(lons)
            lats = reversed(lats)
    if str('W') in str(orientated):
        if compass_bearing < 90. or compass_bearing > 270. :
            lons = reversed(lons)
            lats = reversed(lats)
    return list(lons), list(lats)

def get_sf_polygon(lons, lats, upper_seismogenic_depth,
                    lower_seismogenic_depth, dip, oriented):
    """
    Create and return a fault surface using fault source data.

    :param lons:
        list of the longitude of each point of the fault.
    :param lats:
        list of the latitude of each point of the fault.
    :param upper_seismo_depth:
        Minimum depth ruptures can reach, in km (i.e. depth
        to fault's top edge).
    :param lower_seismo_depth:
        Maximum depth ruptures can reach, in km (i.e. depth
        to fault's bottom edge).
    :param dip:
        Dip angle (i.e. angle between fault surface
        and earth surface), in degrees.
    :param oriented:
        Orientation of the dip. N S E or W.
    :returns:
        Pylygon.
    """
#    cls.check_fault_data(fault_trace, upper_seismogenic_depth,
#                         lower_seismogenic_depth, dip, mesh_spacing)
    # Loops over points in the top edge, for each point
    # on the top edge compute corresponding point on the bottom edge, then
    # computes equally spaced points between top and bottom points.
    
    # organize the points to follow the right hand rule
    lons, lats = orient_fault(lons,lats,oriented)

    vdist_top = upper_seismogenic_depth
    vdist_bottom = lower_seismogenic_depth

    hdist_top = vdist_top / math.tan(math.radians(dip))
    hdist_bottom = vdist_bottom / math.tan(math.radians(dip))
    
    # if the trace is too far from the rupture plan, don't plot it
    if hdist_top > (hdist_bottom-hdist_top)/2.:
        plot_trace = False
    else :
        plot_trace = True
        
    trace_lon = lons
    trace_lat = lats
    
    strike = average_azimuth(lons,lats)
    azimuth = (strike + 90.0) % 360

    lon_top, lon_bottom, lat_top, lat_bottom = [], [], [], []
    for lon,lat in zip(lons,lats):
        lon_t, lat_t = point_at(lon,lat, azimuth, hdist_top)
        lon_b, lat_b = point_at(lon,lat, azimuth, hdist_bottom)
            
        lon_top.append(lon_t)
        lon_bottom.append(lon_b)
        lat_top.append(lat_t)
        lat_bottom.append(lat_b)
        
    poly_lons = numpy.concatenate([lon_top,numpy.array(list(reversed(lon_bottom)))])
    poly_lats = numpy.concatenate([lat_top,numpy.array(list(reversed(lat_bottom)))])
    polygon = [poly_lons,poly_lats]
        
    return trace_lon, trace_lat, plot_trace, polygon

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

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
            * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def FaultProperties(File_prop,Name_of_fault,Model):
    FileName_Prop = File_prop
    Prop = numpy.genfromtxt(FileName_Prop,
                               dtype=[('U100'),('U100'),('f8'),('U100'),('U100'),('f8'),('f8'),('f8'),
                                      ('f8'),('f8'),('U100'),('f8')],skip_header = 1)
    Column_model_name = list(map(lambda i : Prop[i][0],range(len(Prop))))
    Column_fault_name = list(map(lambda i : Prop[i][1],range(len(Prop))))
    index_model = numpy.where(numpy.array(Column_model_name) == Model)[0]
    
    Prop = numpy.take(Prop,index_model)
    index_fault = numpy.where(numpy.array(Column_fault_name[index_model[0]:index_model[-1]+1]) == Name_of_fault)
    Indexfault_final = index_fault[0]

    dip = Prop[Indexfault_final][0][2]
    oriented = Prop[Indexfault_final][0][3]
    upper_sismo_depth = Prop[Indexfault_final][0][5]
    lower_sismo_depth = Prop[Indexfault_final][0][6]
    
    return dip, oriented, upper_sismo_depth, lower_sismo_depth

def FaultGeometry(Model,File_geom):
    NomFichier_InfosZonage = File_geom
    InfosZonage = numpy.genfromtxt(NomFichier_InfosZonage,dtype=[('U100'),('U100'),('f8'),('f8'),('U100')],skip_header = 1)
    Column_model_name = list(map(lambda i : InfosZonage[i][0],range(len(InfosZonage))))
    index_model = numpy.where(numpy.array(Column_model_name) == Model)
    Column_Fault_name_tmp = list(map(lambda i : InfosZonage[i][1],index_model[0]))
    Longitudes_tmp = list(map(lambda i : InfosZonage[i][2],index_model[0]))
    Latitudes_tmp = list(map(lambda i : InfosZonage[i][3],index_model[0]))
    Depths_tmp = list(map(lambda i : InfosZonage[i][4],index_model[0]))
    
    ZoneSelec = Column_Fault_name_tmp
    DicoZone = dict([(k,ZoneSelec.count(k)) for k in set(ZoneSelec)])
    Longitudes = []
    Latitudes = []
    Depths = []
    Column_Fault_name = []
    for cle in DicoZone.keys():
        indices_ZonesSelec = numpy.where(numpy.array(Column_Fault_name_tmp) == cle)
        ColonneNomZone_inter = numpy.take(Column_Fault_name_tmp,indices_ZonesSelec)
        Longitudes_inter = numpy.take(Longitudes_tmp,indices_ZonesSelec)
        Latitudes_inter = numpy.take(Latitudes_tmp,indices_ZonesSelec)
        depth_inter = numpy.take(Depths_tmp,indices_ZonesSelec)

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
    Depths =Depths
    Column_Fault_name = Column_Fault_name
    
    return Column_Fault_name, Depths
