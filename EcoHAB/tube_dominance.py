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

mouse_attention_span = 10  # sec
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


def get_states_and_readouts(antennas, times, t1, t2):
    before = utils.get_idx_pre(t1, times)
    between = utils.get_idx_between(t1, t2, times)
    after = utils.get_idx_post(t2, times)
    states = []
    readouts = []
    if before is not None:
        states.append(antennas[before])
        readouts.append(times[before])
    for idx in between:
        states.append(antennas[idx])
        readouts.append(times[idx])
    assert(len(states) == len(readouts))
    return states, readouts


def mice_in_different_spots(states1, states2):
    for s1 in states1:
        if s1 in states2:
            return False
    return True


def get_times_antennas(ehd, mouse, t_1, t_2):
    ehd.mask_data(t_1, t_2)
    antennas, times = ehd.getantennas(mouse), ehd.gettimes(mouse)
    ehd.unmask_data()
    return times, antennas

def does_mouse1_push_out(m1_states, m1_times, antennas2, times2):
    first_antenna = m1_states[0]
    opposite_antenna = pipe_opposite_antenna[first_antenna]
    #m1 backs out
    
    if utils.skipped_antennas(m1_states):
        #print("Skipped antennas")
        return False
    if len(m1_states) == 1:
        #print("one antenna")
        return False
    #m1 is crossing the chamber
    if len(set(m1_states)) == 3:
        end_state = m1_states[-2]
        end_time = m1_times[-2]
    else:
        end_state = m1_states[-1]
        end_time = m1_times[-1]

    if utils.in_chamber(m1_states[0], end_state):
        #print("in chamber")
        return False
    # mouse1 is backing off not moving forward
    if utils.mouse_backing_off(m1_states):
        #print("backing off")
        return False

    
    m2_states, m2_readouts = get_states_and_readouts(antennas2, times2,
                                                     m1_times[0], end_time)
    between = utils.get_idx_between(m1_times[0], end_time, times2)

    # mouse2 does not move, when mouse 1 moves, skip
    
    if len(between) == 0 :
        #print("no activity of mouse 2")
        return False
    # there are errors in antenna readings, skip
    if utils.skipped_antennas(m2_states):
        #print("Skipped antennas")
        return False
    # mouse1 and mouse2 are in different parts of EcoHAB, skip
    if mice_in_different_spots(m1_states, m2_states): # mice in different parts of the system
        #print("mice in different spots")
        return False
        # start with mouse 2 standing at the opposite antenna
    opposite_idxs = np.where(np.array(m2_states) == opposite_antenna)[0]
    if not len(opposite_idxs):
        #print("mouse 2 never on the opposite side of the pipe")
        # mouse2 is never opposite mouse 1, skip
        return False
    # mouse2 starts at the same antenna
    if m1_states[0] in m2_states[:opposite_idxs[0]]:
        #print("mice start at the same antenna")
        # mice start at the same antenna (following not tube dominance)
        return False


    m2_m1_in_pipe = m2_states[opposite_idxs[0]:]
    m2_times_m1_in_pipe = m2_readouts[opposite_idxs[0]:]
    idx_after = utils.get_idx_post(end_time, times2)
    if idx_after is not None:
        m2_after = antennas2[idx_after]
    else:
        m2_after = m1_states[0]
    #print('mouse 1', m1_states, m1_times)
    #print('mouse 2', m2_states, m2_readouts, m2_after)
    if m2_states[0] == opposite_antenna:
        if m2_readouts[1] - m2_readouts[0] < 5:
            idx = utils.get_idx_pre(m2_readouts[0], times2)
            if idx is not None:
                if antennas2[idx] == m1_states[0] and m1_times[0] - times2[idx] < 5:
                    #print(antennas2[idx], times2[idx])
                    return False
    if np.all(np.array(m2_m1_in_pipe) == opposite_antenna) and m2_after != m1_states[0]: #
        return True

    return False


def get_more_states(antennas, times, midx):
    #save first antenna
    states = [antennas[midx]]
    readouts = [times[midx]]
    midx += 1
    idx = 1
    while True:
        if midx >= len(antennas):
            break
        #read in next antenna
        new_antenna = antennas[midx]
        new_readout = times[midx]
        #if pause too long break
        if new_readout > readouts[idx - 1] + mouse_attention_span:
            break

        states.append(new_antenna)
        readouts.append(new_readout)

        idx += 1
        #if more than 2 antennas, break
        if len(set(states)) == 3:
            # go back to the last readout of the opposite antenna not to loose it
            break
        midx += 1
        
    return states, readouts, midx


def check_mouse1_pushing_out_mouse2(antennas1, times1, antennas2, times2):
    idx = 0
    domination_counter = 0
    while True:
        m1_states, m1_readouts, idx = get_more_states(antennas1, times1, idx)
        if does_mouse1_push_out(m1_states, m1_readouts, antennas2, times2):
            domination_counter += 1
        if idx >= len(antennas1):
            break
        
    return domination_counter
    

def tube_dominance_2_mice_single_phase(ehd, mouse1, mouse2, t_start, t_end):
    """We're checking here, how many times mouse1 dominates over mouse2
    between t_start and t_end.

    """
      
    m1_times, m1_antennas = get_times_antennas(ehd, mouse1, t_start, t_end)
    m2_times, m2_antennas = get_times_antennas(ehd, mouse2, t_start, t_end)
    domination_counter = check_mouse1_pushing_out_mouse2(m1_antennas, m1_times, m2_antennas, m2_times)
        
    return domination_counter

def check_mice(t1_m1, t2_m1, a1_m1, a2_m1, mouse2_antennas, mouse2_times):
    #print('mouse 1', t1_m1, a1_m1, t2_m1, a2_m1)
   
    
    m2_before = get_idx_pre(t1_m1, mouse2_times)
    m2_idxs = get_idx_between(t1_m1, t2_m1, mouse2_times)
    m2_after = get_idx_post(t2_m1, mouse2_times)

    m2_states = []
    m2_times = []
    if m2_before:
        m2_states.append(mouse2_antennas[m2_before])
        m2_times.append(mouse2_times[m2_before])
    for idx in m2_idxs:
        m2_states.append(mouse2_antennas[idx])
        m2_times.append(mouse2_times[idx])
    if m2_after:
        m2_states.append(mouse2_antennas[m2_after])
        m2_times.append(mouse2_times[m2_after])
    for i, t in enumerate(m2_times):
        print('mouse2', t, m2_states[i],)


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
                               'tube_dominance_alternative',
                               mice,
                               phase,
                               main_directory,
                               'tube_dominance_alternative/histograms',
                               prefix,
                               additional_info=add_info_mice)
        single_heat_map(domination[i],
                        'tube_dominance_alternative',
                        main_directory,
                        mice,
                        prefix,
                        phase,
                        xlabel='domineering mouse',
                        ylabel='pushed out mouse',
                        subdirectory='tube_dominance_alternative/histograms',
                        vmax=None,
                        vmin=None,
                        xticks=mice,
                        yticks=mice)
    write_csv_rasters(mice,
                      phases,
                      domination,
                      main_directory,
                      'tube_dominance_alternative/raster_plots',
                      fname_)
    make_RasterPlot(main_directory,
                    'tube_dominance_alternative/raster_plots',
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
