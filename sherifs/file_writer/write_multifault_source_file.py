# -*- coding: utf-8 -*-
"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

@author: Thomas Chartier
"""
import numpy as np

from sherifs.utils.geometry_tools import *
from decimal import Decimal, getcontext
getcontext().prec = 10


def build(f, txt):
    '''
    f : path to file name
    txt : str containing the file info
    '''
    f = open(f,'w')
    f.write(txt)
    f.close()

def start(explo_time, trt):
    '''
    txt : str containing the file info
    explo_time : float or int, exploration time of the model
    trt : str, tectonic region type
    '''
    txt = ''
    # Initiate the xml file
    txt += '<?xml version=\'1.0\' encoding=\'utf-8\'?>\n'
    txt += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
    txt += '\txmlns="http://openquake.org/xmlns/nrml/0.5">\n'
    explo_time = str(round(explo_time,1))
    txt += '<sourceModel name="Hazard Model"'
    txt += ' investigation_time="'+explo_time+'">\n'
    txt += '	<sourceGroup name="group 1" rup_interdep="indep"\n'
    txt += '        src_interdep="indep" '
    txt += ' tectonicRegion="'+trt+'">\n'

    return txt

def wrt_rupture(txt,mag,l,explo_time,rake,sections_id):
    '''
    txt : str containing the file info
    '''

    t = explo_time

    # Poisson probability
    p_occur_1 = np.float32(
    (l*t)**1*np.exp(-l*t)/(np.math.factorial(1))
    )

    xxx1 = Decimal('{:.8f}'.format(np.float32(p_occur_1)))
    p_not_occur = (Decimal('1')-xxx1)
    p_not_occur = '{:.8f}'.format(p_not_occur)

    txt += '\t\t\t\t<multiPlanesRupture probs_occur="'
    txt += str(p_not_occur) +' '+ str(xxx1) +'">\n'
    txt += '\t\t\t\t\t<magnitude>'+ str(mag) +'</magnitude>\n'

    list_sections = ','.join(str(i) for i in sections_id)
    txt += '\t\t\t\t\t<sectionIndexes indexes="'
    txt += str(list_sections)
    txt += '"/>\n'

    txt += '\t\t\t\t\t<rake>'+str(rake)+'</rake>\n'
    txt += '\t\t\t\t</multiPlanesRupture>\n'


    return txt

def start_multifault_source(txt,name,trt,sec_f,ID_number):

    txt += '        <multiFaultSource id="'+str(ID_number)+'" name="'+name+'">\n'
    return txt

def end_multifault_source(txt):
    '''
    txt : str containing the file info
    '''
    txt += '	    </multiFaultSource>\n'
    return txt

def wrt_multifault_source(txt,MFD,Mmin,explo_time,rake,sections_id):
    '''
    txt : str containing the file info
    name : str, rupture names
    trt : str, tectonic region type

    '''

    bin_mag = np.linspace(Mmin,Mmin+0.1*len(MFD)+0.1,num = 2+len(MFD))

    for mag,i_mag in zip(bin_mag,range(len(bin_mag))):
        if i_mag <= len(MFD)-1:
            # annual_rate
            l = MFD[i_mag]
            if l!=0. :
                txt = wrt_rupture(txt,mag,l,explo_time,rake,
                sections_id)

    return txt


def end(txt):
    '''
    txt : str containing the file info
    '''
    txt += '	</sourceGroup>\n'
    txt += '    </sourceModel>\n'
    txt += '</nrml>\n'
    return txt


    return txt
