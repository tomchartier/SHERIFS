# -*- coding: utf-8 -*-

"""
SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Version 1.3

The Seismic Hazard and Earthquake Rates In Fault Systems (SHERIFS) program, is
an open source collection of tools for calculating the rates of earthquakes on
each fault of a fault system allowing complex Fault to Fault ruptures following
the methodology presented in Chartier et al 2017. It includes a set of tools
for checking and visualizing input and outputs and producing JPEG
illustrations. It is released under the GNU Lesser General Public License.  The
SHERIFS program is a code developed in the framework of the PhD thesis of
Thomas Chartier under the supervision of Oona Scotti (IRSN) and Hélène
Lyon-Caen (ENS).

@author: Thomas Chartier
contact : chartier@geologie.ens.fr
"""

import os
import sys
import time
import toml
from pathlib import Path
from sherifs.utils import sap


def SHERIFS(input_file):
    debut = time.time()

    from sherifs.logic_tree.Sources_Logic_Tree_Creator import (
        Sources_Logic_Tree_Creator)

    print('\nRunning SHERIFS version 1.3\n')

    # Load the input file
    param = toml.load(input_file)
    Run_Name = param["Run_Name"]

    root = os.path.dirname(input_file)
    path = os.path.join(root, param["dirpath"])
    param['path'] = path

    # Create folder structure
    folder = os.path.join(path, Run_Name, 'results')
    Path(folder).mkdir(parents=True, exist_ok=True)
    folder = os.path.join(path, Run_Name, 'LOG')
    Path(folder).mkdir(parents=True, exist_ok=True)
    folder = os.path.join(path, Run_Name, 'ssm')
    Path(folder).mkdir(parents=True, exist_ok=True)

    tmp = os.path.join(path, Run_Name, 'calculation_log.txt')
    calculation_log_file = open(tmp, 'w')

    Sources_Logic_Tree_Creator(param, calculation_log_file)
    calculation_log_file.close()

    fin = time.time()-debut
    print(type(fin))
    days = int(fin / 24. / 60. / 60.)
    hours = int((fin - days * 24. * 60. * 60.) / 60. / 60.)
    minutes = int((fin - days * 24. * 60. * 60. - hours* 60. * 60. ) / 60.)
    seconds = (fin - days * 24. * 60. * 60. - hours* 60. * 60.  - minutes * 60.)
    print("The calculation took: " + str(days) + ' days, ' + str(hours) + ' hours, ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds.')


def main(argv):
    """ Run SHERIFS"""

    p = sap.Script(SHERIFS)
    msg = '.txt or .toml file with the information concerning the run.'
    p.arg(name='input_file', help=msg)

    if len(argv) < 1:
        print(p.help())
    else:
        p.callfunc()


if __name__ == "__main__":
    main(sys.argv[1:])
