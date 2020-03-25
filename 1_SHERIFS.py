# -*- coding: utf-8 -*-

"""
SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.2

The Seismic Hazard and Earthquake Rates In Fault Systems (SHERIFS) program, is an open source collection of 
tools for calculating the rates of earthquakes on each fault of a fault system 
allowing complex Fault to Fault ruptures following the methodology presented 
in Chartier et al 2017. It includes a set of tools for checking and visualizing input and outputs
and producing JPEG illustrations.It is released under the GNU Lesser General Public License.
The SHERIFS program is a code developed in the framework of the PhD thesis of Thomas Chartier 
under the supervision of Oona Scotti (IRSN) and Hélène Lyon-Caen (ENS).


@author: Thomas Chartier
contact : chartier@geologie.ens.fr
"""

import time
import os
import sys
debut = time.time()

path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/lib'
sys.path.append(path_lib)
path_f = path_lib + '/logic_tree'
sys.path.append(path_f)
path_f = path_lib + '/file_writer'
sys.path.append(path_f)

from GMPE_Logic_Tree_Creator import GMPE_Logic_Tree_Creator
from Sources_Logic_Tree_Creator import Sources_Logic_Tree_Creator
from OQ_job_Creator import OQ_job_Creator

print ('\nRunning SHERIFS version 1.2\n')




'''###########################'''
'''       Input files         '''
'''###########################'''

Run_Name = 'Example'
File_geom = 'data/Example/Faults_geometry.txt'
File_prop = 'data/Example/Faults_properties.txt'
File_bg = 'data/Example/Background_geometry.txt'
file_prop_bg = 'data/Example/Background_properties.txt'
host_model_file = 'delete_share.xml'




# if rerunning the calcuation do you want to overwrite the existing source model files
overwrite_files = True

use_host_model = False

#if you want to build a GMPE logic tree for OpenQuake, can can always do it later before running the hazard calculation
build_GMPE_LT = False

if not os.path.exists(str(Run_Name)):
    os.makedirs(str(Run_Name))
if not os.path.exists(str(Run_Name) + '/results'):
     os.makedirs(str(Run_Name) + '/results')

Domain_in_model = []    

OQ_job_Creator = OQ_job_Creator(Run_Name) # ask the info about the run and create the job.ini file

seed = OQ_job_Creator.seed
nb_random_sampling = OQ_job_Creator.nb_sample
Mmin = OQ_job_Creator.Mmin
sr_correl = OQ_job_Creator.sr_correl
size_of_increment = OQ_job_Creator.size_of_increment
Mmax_range = OQ_job_Creator.Mmax_range
fit_quality = 5 #maximum misfit between the model and the target (in %)

calculation_log_file = open(Run_Name+'/calculation_log.txt','w')
Sources_Logic_Tree_Creator = Sources_Logic_Tree_Creator(Run_Name,File_geom,
                                                        File_prop,File_bg,file_prop_bg,Domain_in_model,
                                                        nb_random_sampling,seed,Mmin,sr_correl,
                                                        size_of_increment,Mmax_range,overwrite_files,float(fit_quality)/100.,
                                                        calculation_log_file,use_host_model,host_model_file) 
                                                        #create the source models logic tree
                                                        
calculation_log_file.close()

Domain_in_model = Sources_Logic_Tree_Creator.Domain_in_model

if build_GMPE_LT == True:
    GMPE_Logic_Tree_Creator(Run_Name,Domain_in_model) #create the logic tree of GMPEs


fin = time.time()-debut
days = int(fin / 24. / 60. / 60.)
hours = int((fin - days * 24. * 60. * 60.) / 60. / 60.)
minutes = int((fin - days * 24. * 60. * 60. - hours* 60. * 60. ) / 60.)
seconds = (fin - days * 24. * 60. * 60. - hours* 60. * 60.  - minutes * 60.)
print("The calculation took: " + str(days) + ' days, ' + str(hours) + ' hours, ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds.')


