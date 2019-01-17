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


def get_idx_pre(t0, times):
    return np.where(np.array(times) < t0)[0]

#1 check the mouse 2 before the readout of mouse 1
def check_mouse2_antenna_pre_time1(antenna_m1, t_m1, mouse2_antennas, mouse2_times):
    idx = get_idx_pre(t_m1, mouse2_times)
    antenna_m2 = mouse2_antennas[idx]
    if antenna_m2 == opposite_antenna[antenna_m1]:
        return True
    return False

def check_mouse2_before_mouse1_reading(idx_m1_t1, mouse1_antennas, mouse1_times, mouse2_antennas, mouse2_times):
    a_m1 = mouse1_antennas[idx_m1_t1]
    t1_m1 = mouse1_times[idx_m1_t1]
    if check_mouse2_antenna_pre_time1(a_m1, t1_m1, mouse2_antennas, mouse2_times):
        idx_m2 = get_idx_pre(t1_m1, mouse2_times)
        nexta_m1 = mouse1_antennas[idx_m1_t1 + 1]
        t2_m1 = mouse1_times[idx_m1_t1 + 1]
        nexta_m2 = mouse2_antennas[idx_m2 + 1]
        t2_m2 = mouse2_times[idx_m2 + 1]
        if next_am1 == a_m1 and t2_m1 < t2_m2:  # mouse 1 backs out first
            return 'mouse 1'
        elif next_am2 == a_m2 and t2_m2 < t1_m1: #mouse 2 backs out first
            return 'mouse 2'

        
