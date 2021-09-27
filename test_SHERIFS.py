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
from mpl_toolkits.basemap import Basemap
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
import tkinter as tk
from tkinter import ttk, Label, Text, INSERT,END, StringVar,Listbox,Button,Entry,Checkbutton
from tkinter.ttk import Combobox
from tkinter import messagebox
import pandas as pd
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from glob import glob
import xml.etree.ElementTree as ET
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection


w_sc_nb = tk.Tk()
w_sc_nb.title('OK?')
label = Label(w_sc_nb, text='\n\n\n\n\n\n        Everything is working well           \n\nYou can close this window\nand start with SHERIFS\n\n\n\n')
label.pack()
w_sc_nb.mainloop()

print('\n\n\n\nEverything is OK\n\n\n\n')
print('')
print('You can use SHERIFS!\n\n')
