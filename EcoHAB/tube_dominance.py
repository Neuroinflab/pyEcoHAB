from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utility_functions as utils
import numpy as np
import matplotlib.pyplot as plt
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot, single_heat_map
from numba import jit
from collections import OrderedDict
nbins = 10
homepath = os.path.expanduser("~/")

pipe_opposite_antenna = { 1:2,
                          2:1,
                          3:4,
                          4:3,
                          5:6,
                          6:5,
                          7:8,
                          8:7}
cage_opposite_antenna = {1:8,
                         2:3,
                         3:2,
                         4:5,
                         5:4,
                         6:7,
                         7:6,
                         8:1}


def get_idx_pre(t0, times):
    idxs = np.where(np.array(times) < t0)[0]
    if len(idxs):
        return idxs[0]
    return None

def get_idx_post(t0, t1, times):
    return  np.where((np.array(times) >= t0) & (np.array(times) < t1))[0]

#1 check the mouse 2 before the readout of mouse 1
def check_mouse2_antenna_pre_time1(antenna_m1,
                                   t_m1,
                                   mouse2_antennas,
                                   mouse2_times):
 
    idx = get_idx_pre(t_m1, mouse2_times)
    if idx is None:
        return False
    antenna_m2 = mouse2_antennas[idx]
    print('mouse 1:', antenna_m1, t_m1)
    print('mouse 2:', antenna_m2)
    if antenna_m2 == pipe_opposite_antenna[antenna_m1]:
        return True
    return False

def check_m1_entering_and_pushing_m2_out(idx_m1_t1,
                                         mouse1_antennas,
                                         mouse1_times,
                                         mouse2_antennas,
                                         mouse2_times):
    a_m1 = mouse1_antennas[idx_m1_t1]
    t1_m1 = mouse1_times[idx_m1_t1]
    if check_mouse2_antenna_pre_time1(a_m1, t1_m1,
                                      mouse2_antennas, mouse2_times):
        idx_m2 = get_idx_pre(t1_m1, mouse2_times)
        a_m2 = mouse2_antennas[idx_m2]
        t2_m1 = mouse1_times[idx_m1_t1 + 1]
        nexta_m2 = mouse2_antennas[idx_m2 + 1]
        t2_m2 = mouse2_times[idx_m2 + 1]
        print('mouse 1:', t2_m1)
        print('mouse 1:', nexta_m2, t2_m2)
        if nexta_m2 == a_m2 and t2_m2 < t2_m1: #mouse 2 backs out first
            print('True')
            return True
    return False

def check_mouse2_after_mouse1_reading(idx_m1_t1,
                                      mouse1_antennas,
                                      mouse1_times,
                                      mouse2_antennas,
                                      mouse2_times):
    antenna_m1 = mouse1_antennas[idx_m1_t1]
    t1_m1 = mouse1_times[idx_m1_t1]
    t2_m1 = mouse1_times[idx_m1_t1 + 1]
    m2_idxs = get_idx_post(t1_m1, t2_m1, mouse2_times)
    print('mouse 1:', antenna_m1, t1_m1, t2_m1)
    print('mouse 2:', m2_idxs)
    opposite_antenna = pipe_opposite_antenna[antenna_m1]
    opposite_cage = cage_opposite_antenna[opposite_antenna]
    for idx in m2_idxs:
        print(mouse2_antennas[idx], mouse2_times[idx])
        if mouse2_antennas[idx] == opposite_antenna:
            print('Opposite antenna')
            if idx+1 < len(mouse2_antennas):
                if mouse2_antennas[idx+1] == opposite_antenna or mouse2_antennas[idx+1] == opposite_cage:
                    print(opposite_antenna, opposite_cage)
                    return True
    

def tube_dominance_2_mice_single_phase(ehd, mouse1, mouse2, t_start, t_end):
    """We're checking here, how many times mouse1 dominates over mouse2
    between t_start and t_end.

    """
    domination_counter = 0
    ehd.mask_data(t_start, t_end)
    antennas1 = ehd.getantennas(mouse1)
    times1 = ehd.gettimes(mouse1)
    antennas2 = ehd.getantennas(mouse2)
    times2 = ehd.gettimes(mouse2)
    for idx in range(len(antennas1[:-1])):
        print('Check pre')
        if check_m1_entering_and_pushing_m2_out(idx,
                                                antennas1,
                                                times1,
                                                antennas2,
                                                times2):
            domination_counter += 1
        print('Check post')
        if check_mouse2_after_mouse1_reading(idx,
                                             antennas1,
                                             times1,
                                             antennas2,
                                             times2):
            domination_counter += 1
    print(domination_counter)
    return domination_counter

def tube_domination_single_phase(ehd, cf, phase, print_out=True):
    mice = ehd.mice
    st, en = cf.gettime(phase)
    domination =  np.zeros((len(mice), len(mice)))
    if print_out:
        print(phase)
    for i, mouse1 in enumerate(mice):
        for j, mouse2 in enumerate(mice):
            if i != j:
                domination[i, j] = tube_dominance_2_mice_single_phase(ehd,
                                                                      mouse1,
                                                                      mouse2,
                                                                      st,
                                                                      en)
    return domination


def tube_domination_whole_experiment(ehd, cf, main_directory, prefix, remove_mouse=None, print_out=True):
    phases = cf.sections()
    phases = utils.filter_dark(phases)
    mice = ehd.mice
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    domination = np.zeros((len(phases), len(mice), len(mice)))
    fname_ = 'tube_domination_%s%s.csv' % (prefix, add_info_mice)
    for i, phase in enumerate(phases):
        domination[i] = tube_domination_single_phase(ehd, cf, phase, print_out=print_out)
        save_single_histograms(domination[i],
                               'tube_domination_alternative',
                               mice,
                               phase,
                               main_directory,
                               'tube_domination_alternative/histograms',
                               prefix,
                               additional_info=add_info_mice)
        single_heat_map(domination[i],
                        'tube_domination_alternative',
                        main_directory,
                        mice,
                        prefix,
                        xlabels='domineering mouse',
                        ylabels='pushed out mouse',
                        subdirectory='tube_domination_alternative/figs',
                        vmax=None,
                        vmin=None,
                        xticks=mice,
                        yticks=mice)
    write_csv_rasters(mice,
                      phases,
                      domination,
                      main_directory,
                      'tube_domination_alternative/raster_plots',
                      fname_)
    make_RasterPlot(main_directory,
                    'tube_domination_alternative/raster_plots',
                    domination,
                    phases,
                    fname_,
                    mice,
                    title='# dominations')
    
if __name__ == '__main__':
    for new_path in datasets:
        path = os.path.join(homepath, new_path)
        prefix = utils.make_prefix(path)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
        else:
            remove_mouse = None
        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 500
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        prefix = utils.make_prefix(path)
        res_dir = utils.results_path(path)
        cf = ExperimentConfigFile(path)
        tube_domination_whole_experiment(ehd, cf, res_dir, prefix, remove_mouse=None, print_out=True)
