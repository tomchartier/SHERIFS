# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

@author: Thomas Chartier
"""
import os, csv
import numpy as np
import geojson
import xml.etree.ElementTree as ET
import mfd_shape

from geometry_tools import distance

class bg():
    """
    Extract the geometry and properities of the background.
    """
    def geom(model_name,file_geom):
        Lon_bg = []
        Lat_bg = []
        if not ".geojson" in file_geom:

            # manually defined  in the file Background geometry
            geom_bg = np.genfromtxt(file_geom,dtype=[('U100'),('f8'),('f8')],skip_header = 1)

            column_model = list(map(lambda i : geom_bg[i][0],range(len(geom_bg))))
            index_model = np.where(np.array(column_model) == model_name)[0]
            Lon_bg = list(map(lambda i : geom_bg[i][1],index_model))
            Lat_bg = list(map(lambda i : geom_bg[i][2],index_model))


            if Lon_bg == 0:
                print('Error!! Check your input background geometry')

        else : #it's a geojson file
            with open(file_geom) as f:
                gj = geojson.load(f)
            bgs = gj['features']
            for bg_i in range(len(bgs)):
                if bgs[bg_i]['properties']['model'] == model_name :
                    Lon_bg = [i[0] for i in bgs[bg_i]['geometry']["coordinates"][0][0]]
                    Lat_bg = [i[1] for i in bgs[bg_i]['geometry']["coordinates"][0][0]]
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

    def get_multipoints(EQ_rate_BG,
                        M_min,
                        bbPath_BG,
                        list_bg_xml,
                        include_all_faults,
                        outside_faults,
                        faults_data,
                        OQ_entry_faults,
                        pathssm,
                        pathlog,
                        list_src_files):
        """
        Import the multipoint sources from a list of xml and create a list
        of PointSources with the modified rates to fit the background rates
        calculated by SHERIFS. Then creates the xml bg files with the modified
        rates.

        input:
        EQ_rate_BG : list, MFD for the background seismicity
        M_min : float, Mmin for the model
        bbPath_BG: mplPath.Path, ploygon of the background
        list_bg_xml : str, path to the smooth seismicity background xml files
        include_all_faults : bool, if faults outside of the bg are included
        outside_faults : list, list of the faults outside of the bg
        faults_data : dict, fault information
        OQ_entry_faults : list, mfd for single fault ruptres
        pathssm : str, path to the the folder containing a list of xml files
        pathlog : str, log folder for the branch
        list_src_files : list, list of path of the xml source files for oq


        """

        sharp_dist_criteria = 15. # km (TODO move out)

        Mmin_checked = False

        Mmax = M_min+len(EQ_rate_BG)*0.1-1
        mags = np.linspace(M_min,Mmax,len(EQ_rate_BG))

        log_data = [['lon','lat','M5_init','M5_final','M6_init','M6_final','M7_init','M7_final']]

        # read the xml and stores the list of aValues
        pts_list = {} # points in the background area
        pts_out_list = {} # points outside the background area
        list_pt_loc_in, list_pt_loc_out = [], []
        faults_out_pts = []
        #sum_rates = 0.
        sum_rates_in = [0. for _ in mags]
        id_bg = 0
        bg_info = {}
        print()
        for fbg in list_bg_xml:
            tree = ET.parse(fbg)
            ET.register_namespace('', "http://openquake.org/xmlns/nrml/0.5")
            nrml = tree.getroot()
            i_mltptsr =0
            point_data = {}
            tectonicRegion = nrml[0][0].get('tectonicRegion')
            for multiPointSource in nrml[0][0]:
                i_child = 0
                for child in multiPointSource:
                    if "multiPointGeometry" in str(child) :
                        i_geom = i_child
                        i_grandchild = 0
                        for grandchild in nrml[0][0][i_mltptsr][i_geom]:

                            if "posList" in str(grandchild) :
                                s_tmp = multiPointSource[i_geom][i_grandchild].text
                                s_tmp=s_tmp.replace('\n','')
                                list_pts=[float(i) for i in s_tmp.split(' ') if i != '']
                                for id_pt in range(int(len(list_pts)/2)):
                                    pt_in_BG = False
                                    close_to_out_fault = False
                                    lon, lat = list_pts[id_pt*2], list_pts[id_pt*2+1]

                                    str_loc = str(lon)+'_'+str(lat)
                                    if bbPath_BG.contains_point((lon,lat)) == 1:
                                        pt_in_BG = True
                                        list_pt_loc_in.append(str_loc)
                                    else :
                                        list_pt_loc_out.append(id_pt)


                                        # if point is close to an outside fault
                                        # first, check roughly
                                        pt_very_far = True
                                        rough_dist_crateria = 80 #km
                                        for i_fault in outside_faults :
                                            if pt_very_far == True :
                                                mean_lon = np.mean(faults_data[i_fault]['lon'])
                                                mean_lat = np.mean(faults_data[i_fault]['lat'])
                                                dist_i = distance(mean_lon, mean_lat, lon, lat)
                                                if dist_i < rough_dist_crateria :
                                                    pt_very_far = False

                                        # find distance to fault mean location
                                        if pt_very_far == False :
                                            dist = 1000000.
                                            closest_fault = "nope"
                                            for i_fault in outside_faults :
                                                for lon_f,lat_f in zip(faults_data[i_fault]['lon'],faults_data[i_fault]['lat']):
                                                    dist_i = distance(lon_f, lat_f, lon, lat)
                                                    if dist_i < dist :
                                                        dist = dist_i
                                                        closest_fault = i_fault
                                            closest_dist = dist
                                            if closest_dist < sharp_dist_criteria :
                                                close_to_out_fault = True
                                                faults_out_pts.append(closest_fault)

                                    point_data.update({id_pt:{'lon':lon,
                                                              'lat':lat,
                                                            'pt_in_BG':pt_in_BG,
                                                            'id_bg':id_bg,
                                                            'str_loc':str_loc,
                                                            'close_to_out_fault':close_to_out_fault}})
                                    if close_to_out_fault :
                                        point_data[id_pt].update({'closest_fault':closest_fault})

                            if "upperSeismoDepth" in str(grandchild) :
                                upperSeismoDepth = multiPointSource[i_geom][i_grandchild].text
                            if "lowerSeismoDepth" in str(grandchild) :
                                lowerSeismoDepth = multiPointSource[i_geom][i_grandchild].text
                            i_grandchild += 1

                    if "multiMFD" in str(child) :
                        i_mfd = i_child
                        i_grandchild = 0
                        for grandchild in nrml[0][0][i_mltptsr][i_child]:
                            if "a_val" in str(grandchild) :
                                s_tmp = nrml[0][0][i_mltptsr][i_child][i_grandchild].text
                                s_tmp=s_tmp.replace('\n','')
                                list_a=[float(i) for i in s_tmp.split(' ') if i != '']
                            if "b_val" in str(grandchild) :
                                b_val = nrml[0][0][i_mltptsr][i_child][i_grandchild].text
                                b_val = float(b_val)
                            if "bin_width" in str(grandchild) :
                                bin_width = nrml[0][0][i_mltptsr][i_child][i_grandchild].text
                                bin_width = float(bin_width)
                            if "max_mag" in str(grandchild) :
                                max_mag = nrml[0][0][i_mltptsr][i_child][i_grandchild].text
                                max_mag = float(max_mag)
                            if "min_mag" in str(grandchild) :
                                min_mag = nrml[0][0][i_mltptsr][i_child][i_grandchild].text
                                min_mag = float(min_mag)


                            i_grandchild += 1

                        if Mmin_checked == False :
                            if float(min_mag) < M_min :
                                print("!!!!")
                                print("WARNING : BG has a smaller Mmin than the SHERIFS input")
                                print("!!!!")
                                Mmin_checked = True
                        mfd_kind = nrml[0][0][i_mltptsr][i_child].get('kind')
                    if "magScaleRel" in str(child) :
                        magScaleRel = multiPointSource[i_child].text
                    if "ruptAspectRatio" in str(child) :
                        ruptAspectRatio = multiPointSource[i_child].text
                    if "nodalPlaneDist" in str(child) :
                        nodalPlaneDist = []
                        for nodalPlane in multiPointSource[i_child]:
                            nodalPlaneDist.append([nodalPlane.get('dip'),
                                                   nodalPlane.get('probability'),
                                                   nodalPlane.get('rake'),
                                                   nodalPlane.get('strike')])
                    if "hypoDepthDist" in str(child) :
                        hypoDepthDist = []
                        for hypoDepth in multiPointSource[i_child]:
                            hypoDepthDist.append([hypoDepth.get('depth'),
                                                  hypoDepth.get('probability')])
                    i_child+=1

            # defining the mfd for each point
            mfd_param = {'b_value':b_val}
            bin_mag = [M_min+float(i)/10. for i in range(int((max_mag-M_min)*10.)+1)]

            for id_pt in range(len(point_data)):
                p_mfd = mfd_shape.GR(mfd_param,bin_mag)
                r_Min = 10.**(list_a[id_pt]-b_val*(M_min-0.05)) - \
                10.**(list_a[id_pt]-b_val*(M_min+0.05))
                mfd = [i * r_Min for i in p_mfd]
                point_data[id_pt].update({'mfd':mfd})
                if point_data[id_pt]['pt_in_BG'] :
                    sum_rates_in = [i+j for i,j in zip(sum_rates_in,mfd)]
            bg_info.update({id_bg:{'tectonicRegion':tectonicRegion,
                                    'magScaleRel':magScaleRel,
                                    'ruptAspectRatio':ruptAspectRatio,
                                    'nodalPlaneDist':nodalPlaneDist,
                                    'hypoDepthDist':hypoDepthDist,
                                    'upperSeismoDepth':upperSeismoDepth,
                                    'lowerSeismoDepth':lowerSeismoDepth,
                                    'point_data':point_data}})
            id_bg += 1



        # Normalize to fit the BG MFD and writting
        sum_bg_min = 0.
        i_bg = 0
        for fbg in list_bg_xml:
            bg_file = pathssm +'/bg_' + str(i_bg) + '.xml'
            bg_file = open(bg_file,'w')
            # Initiate the xml file
            line='<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
            line+='<nrml xmlns:gml="http://www.opengis.net/gml"\n'
            line+='\txmlns="http://openquake.org/xmlns/nrml/0.5">\n'
            line+='\t<sourceModel name="background '+str(i_bg)+'">\n'
            line+='\t\t<sourceGroup name="bg_'+str(i_bg)+'" rup_interdep="indep" src_interdep="indep" tectonicRegion="' + str(bg_info[i_bg]['tectonicRegion']) + '" >\n'
            point_data = bg_info[i_bg]['point_data']
            for id_pt in range(len(point_data)):
                line+='\t\t\t<pointSource id="'+str(i_bg)+'_'+str(id_pt)+'" name="">\n'
                line+='\t\t\t\t<pointGeometry>\n'
                line+='\t\t\t\t\t<gml:Point>\n'
                line+='\t\t\t\t\t\t<gml:pos>\n'
                line+='\t\t\t\t\t\t\t '+str(point_data[id_pt]['lon'])+' '+str(point_data[id_pt]['lat'])+'\n'
                line+='\t\t\t\t\t\t</gml:pos>\n'
                line+='\t\t\t\t\t</gml:Point>\n'
                line+='\t\t\t\t<upperSeismoDepth>'
                line+='\t\t\t\t'+str(bg_info[i_bg]['upperSeismoDepth'])
                line+='</upperSeismoDepth>\n'
                line+='\t\t\t\t<lowerSeismoDepth>'
                line+='\t\t\t\t'+str(bg_info[i_bg]['lowerSeismoDepth'])
                line+='</lowerSeismoDepth>\n'
                line+='\t\t\t\t</pointGeometry>\n'
                line+='\t\t\t\t<magScaleRel>'
                line+='\t\t\t\t'+str(bg_info[i_bg]['magScaleRel'])
                line+='</magScaleRel>\n'
                line+='\t\t\t\t<ruptAspectRatio>'
                line+='\t\t\t\t'+str(bg_info[i_bg]['ruptAspectRatio'])
                line+='</ruptAspectRatio>\n'

                # WORKING OF THE MFD
                line+='\t\t\t<incrementalMFD binWidth=\"0.10\" minMag="'+ str(M_min)+'">\n'
                mfd = point_data[id_pt]['mfd']
                mfd_init = mfd
                if point_data[id_pt]['pt_in_BG'] :
                    mod_mfd = []
                    for i_m in range(min(len(mfd),len(sum_rates_in))):
                        weight_pt = mfd[i_m]/sum_rates_in[i_m]
                        mod_mfd.append(EQ_rate_BG[i_m] * weight_pt)
                    mfd = mod_mfd
                if include_all_faults : # if there are faults outside the bg
                    if point_data[id_pt]['close_to_out_fault'] :
                        ### ONLY DOING THE REDUCTION CONSIDERING SINGLE FAULT RUPTURES
                        # get closest fault
                        i_fault = point_data[id_pt]["closest_fault"]
                        # find ruptures of the closest faults
                        mfd_single_clst_f = OQ_entry_faults[i_fault]
                        mod_mfd = []

                        for i_m in range(min(len(mfd),len(mfd_single_clst_f))):
                            nb_pt_close = faults_out_pts.count(i_fault)
                            reduction = mfd_single_clst_f[i_m]/float(nb_pt_close)
                            mod_ri = mfd[i_m]-reduction
                            # the reduction is up to 50 %, not more
                            if mod_ri >= mfd[i_m] * 0.5 :
                                mod_mfd.append(mod_ri)
                            else :
                                mod_mfd.append(mfd[i_m]* 0.5)
                        mfd = mod_mfd



                line+='\t\t\t<occurRates> ' + ' '.join(list(map(str, mfd))) + '</occurRates>\n'

                # log the mfd values for csv
                log_mags = [5.0,6.0,7.0]
                r_log = []
                for mag in log_mags:
                    if mag < M_min :
                        r_log.append(0.)
                        r_log.append(0.)
                    elif mag > M_min + 0.1 * len(mfd):
                        r_log.append(0.)
                        r_log.append(0.)
                    else :
                        i_mag  = int((mag-M_min)*10.)
                        r_log.append(sum(mfd_init[i_mag:]))
                        r_log.append(sum(mfd[i_mag:]))
                log_data.append([point_data[id_pt]['lon'],point_data[id_pt]['lat']]+r_log)



                line+='\t\t\t</incrementalMFD>\n'
                line+='\t\t\t\t<nodalPlaneDist>\n'
                for nodalPlane in bg_info[i_bg]['nodalPlaneDist']:
                    line+='\t\t\t\t\t<nodalPlane dip="'+nodalPlane[0]+\
                    '" probability="'+nodalPlane[1]+\
                    '" rake="'+nodalPlane[2]+\
                    '" strike="'+nodalPlane[3]+\
                    '"/>\n'
                line+='\t\t\t\t</nodalPlaneDist>\n'
                line+='\t\t\t\t<hypoDepthDist>\n'
                for hypoDepth in bg_info[i_bg]['hypoDepthDist']:
                    line+='\t\t\t\t\t<hypoDepth depth="'+hypoDepth[0]+\
                    '" probability="'+hypoDepth[1]+\
                    '"/>\n'
                line+='\t\t\t\t</hypoDepthDist>\n'
                line+='\t\t\t</pointSource>\n'

            line+='\t\t</sourceGroup>\n'
            line+='\t</sourceModel>\n'
            line+='</nrml>\n'

            bg_file.write(line)
            bg_file.close()

            if not pathssm +'/bg_' + str(i_bg) + '.xml' in list_src_files :
                list_src_files.append(pathssm +'/bg_' + str(i_bg) + '.xml')

            i_bg += 1

        if not os.path.exists(pathlog+'/bg'):
            os.makedirs(pathlog+'/bg')
        with open(pathlog + "/bg/bg_rates.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(log_data)
