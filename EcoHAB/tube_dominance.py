from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utility_functions as utils
import numpy as np
import matplotlib.pyplot as plt
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot
from numba import jit
from collections import OrderedDict
nbins = 10
homepath = os.path.expanduser("~/")

opposite_antenna = { 1:2,
                     2:1,
                     3:4,
                     4:3,
                     5:6,
                     6:5,
                     7:8,
                     8:7}


#1 check the mouse 2 before the readout of mouse 1
def check_mouse2_before_mouse1_reading(antenna_m1, t_m1, mouse2_antennas, mouse2_times):
    idx = np.where(np.array(times2) < t_m1)[0]
    antenna_m2 = mouse2_antennas[idx]
    if 
