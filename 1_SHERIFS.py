# -*- coding: utf-8 -*-

"""
SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

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
from lib.utils import sap
import toml
import shutil

# If you are running SHERIFS with spyder define "input_file" here. Then run.

def SHERIFS(input_file):
    debut = time.time()

    path_actuel=os.path.dirname(os.path.abspath(__file__))
    path_lib = path_actuel + '/lib'
    sys.path.append(path_lib)
    path_f = path_lib + '/logic_tree'
    sys.path.append(path_f)
    path_f = path_lib + '/file_writer'
    sys.path.append(path_f)

    #from GMPE_Logic_Tree_Creator import GMPE_Logic_Tree_Creator
    from Sources_Logic_Tree_Creator import Sources_Logic_Tree_Creator
    from OQ_job_Creator import OQ_job_Creator

    print ('\nRunning SHERIFS version 1.3\n')

    '''###########################'''
    '''       Input files         '''
    '''###########################'''
    # Load the input file
    param = toml.load(input_file)
    '''
    # Load the input file (old)
    lines = open(input_file,'r').readlines()
    lines = [line.rstrip('\n') for line in lines]
    list_fbg = []
    for line in lines:
        if "Run_Name" in line :
            Run_Name = line.split(':')[1].replace(' ','')
        if "File_geom" in line :
            File_geom = line.split(':')[1].replace(' ','')
        if "File_prop" in line :
            File_prop = line.split(':')[1].replace(' ','')
        if "File_bg" in line :
            File_bg = line.split(':')[1].replace(' ','')
        if "file_prop_bg" in line :
            file_prop_bg = line.split(':')[1].replace(' ','')
        if "host_model_file" in line :
            host_model_file = line.split(':')[1].replace(' ','')
        if "overwrite_files" in line :
            if "rue" in line :
                overwrite_files = True
            elif "alse" in line :
                overwrite_files = False
        if "use_host_model" in line :
            if "rue" in line :
                use_host_model = True
            elif "alse" in line :
                use_host_model = False
        #maximum misfit between the model and the target (in %)
        if "fit_quality" in line :
            fit_quality = float(line.split(':')[1].replace(' ',''))

        # read the list of bg files to use
        if "bgf" in line :
            tmp = line.split(':')[1]
            fbgpath = tmp.replace(" ", "")
            if os.path.isdir(fbgpath):
                list_fbg = [f for f in listdir(fbgpath) if isfile(join(fbgpath, f))]
            else :
                list_fbg = tmp.split(' ')
                while '' in list_fbg:
                    list_fbg.remove('')

        else :
            fbgpath = None
    '''
    Run_Name = param["Run_Name"]

    path = param["dirpath"]
    if len(path)!= 0:
        if not path[-1]=="/":
            path += "/"
            param["dirpath"] = path

    # if rerun is forced, delete everything
    if param['main']['parameters']['force_rerun'] in ['true','True']:
        if os.path.exists(path+str(Run_Name)):
            shutil.rmtree(path+str(Run_Name))
    #create folder structure
    if not os.path.exists(path+str(Run_Name)):
        os.makedirs(path+str(Run_Name))
    if not os.path.exists(path+str(Run_Name) + '/results'):
         os.makedirs(path+str(Run_Name) + '/results')
    if not os.path.exists(path+str(Run_Name) + '/LOG'):
         os.makedirs(path+str(Run_Name) + '/LOG')
    if not os.path.exists(path+str(Run_Name) + '/ssm'):
         os.makedirs(path+str(Run_Name) + '/ssm')

    #Domain_in_model = []

    '''
    OQ_job_Creator = OQ_job_Creator(Run_Name) # ask the info about the run and create the job.ini file

    seed = OQ_job_Creator.seed
    nb_random_sampling = OQ_job_Creator.nb_sample
    Mmin = OQ_job_Creator.Mmin
    sr_correl = OQ_job_Creator.sr_correl
    size_of_increment = OQ_job_Creator.size_of_increment
    Mmax_range = OQ_job_Creator.Mmax_range
    #fit_quality = 5 #maximum misfit between the model and the target (in %)
    '''
    calculation_log_file = open(path+Run_Name+'/calculation_log.txt','w')

    '''
    sltc = Sources_Logic_Tree_Creator(Run_Name,File_geom,
                                    File_prop,File_bg,file_prop_bg,Domain_in_model,
                                    nb_random_sampling,seed,Mmin,sr_correl,
                                    size_of_increment,Mmax_range,overwrite_files,float(fit_quality)/100.,
                                    calculation_log_file,use_host_model,host_model_file
                                    ,list_fbg,fbgpath)
                                    #create the source models logic tree
    '''

    Sources_Logic_Tree_Creator(param,calculation_log_file)
                                #create the source models logic tree

    calculation_log_file.close()

    #Domain_in_model = sltc.Domain_in_model

    #if you want to build a GMPE logic tree for OpenQuake, you can always do it later before running the hazard calculation
    #build_GMPE_LT = False
    #if build_GMPE_LT == True:
    #    GMPE_Logic_Tree_Creator(Run_Name,Domain_in_model) #create the logic tree of GMPEs


    fin = time.time()-debut
    days = int(fin / 24. / 60. / 60.)
    hours = int((fin - days * 24. * 60. * 60.) / 60. / 60.)
    minutes = int((fin - days * 24. * 60. * 60. - hours* 60. * 60. ) / 60.)
    seconds = (fin - days * 24. * 60. * 60. - hours* 60. * 60.  - minutes * 60.)
    print("The calculation took: " + str(days) + ' days, ' + str(hours) + ' hours, ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds.')



def main(argv):
    """ Run SHERIFS"""

    p = sap.Script(SHERIFS)
    p.arg(name='input_file', help='.txt file with the information concerning the run.')

    if len(argv) < 1:
        print(p.help())
    else:
        p.callfunc()


if __name__ == "__main__":
    main(sys.argv[1:])
