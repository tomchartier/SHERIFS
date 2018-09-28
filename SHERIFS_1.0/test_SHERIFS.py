#!/usr/bin/env python2
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
path_dossier = path_actuel + '/lib'
sys.path.append(path_dossier)
path_actuel=os.path.dirname(os.path.abspath(__file__))
path_dossier = path_actuel + '/lib/input_checks'
sys.path.append(path_dossier)
from GMPE_Logic_Tree_Creator import GMPE_Logic_Tree_Creator
from Sources_Logic_Tree_Creator import Sources_Logic_Tree_Creator
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
import EQ_rates_for_each_fault
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
label = Label(w_sc_nb, text='\n\n\n\n\n\n        Everything is working            \n\n\n\n\n\n')
label.pack()
w_sc_nb.mainloop()

print('\n\n\n\n      Everything is working\n\n\n\n')
print('')
print('      You can use SHERIFS!\n\n')

        