#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for SHERIFS
If this test works, SHERIFS should run fine as well.
@author: thomas
"""
import time
import os
import sys
path_actuel=os.path.dirname(os.path.abspath(__file__))
path_lib = path_actuel + '/lib'
sys.path.append(path_lib)
path_f = path_lib + '/logic_tree'
sys.path.append(path_f)
path_f = path_lib + '/file_writer'
sys.path.append(path_f)
path_f = path_lib + '/core'
sys.path.append(path_f)
path_f = path_lib + '/hm_visual'
sys.path.append(path_f)
path_f = path_lib + '/utils'
sys.path.append(path_f)

from logic_tree.GMPE_Logic_Tree_Creator import GMPE_Logic_Tree_Creator
from logic_tree.Sources_Logic_Tree_Creator import Sources_Logic_Tree_Creator
from OQ_job_Creator import OQ_job_Creator
import xml.etree.ElementTree as ET
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Polygon
import matplotlib.patches as patches
import matplotlib.path as mplPath
from scipy.stats import gaussian_kde
mpl.interactive(False)
import Extract_data
import Plot_mfd
import Participation_rates
import plt_catalog
import plot_FtF_GIF
import pandas as pd
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from glob import glob
import xml.etree.ElementTree as ET
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import toml
import pickle

print('\n\n\n\nEverything is OK\n\n\n\n')
print('')
print('You can use SHERIFS!\n\n')
