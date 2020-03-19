# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import
import numpy as np
from . import utility_functions as utils
from .write_to_file import save_single_histograms, write_csv_rasters
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot, single_heat_map
from . import exec_functions as dispatch
from . import dominance_in_2_cages as dom2
how_many_antennas = 3

mas = 10  # sec
nbins = 10

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



def mice_in_different_spots(states1, states2):
    for s1 in states1:
        if s1 in states2:
            return False
    return True


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

    
    m2_states, m2_readouts = utils.get_states_and_readouts(antennas2, times2,
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
    idx_after = utils.get_idx_post(end_time, times2)
    if idx_after is not None:
        m2_after = antennas2[idx_after]
    else:
        m2_after = m1_states[0]
    if np.all(np.array(m2_m1_in_pipe) == opposite_antenna):
        if m2_after != m1_states[0]:
            return True
    if m2_readouts[opposite_idxs[0]] > m1_times[0]:
        if m1_states[0] not in m2_states[opposite_idxs[0]:]:
            return True
    return False


def check_mouse1_pushing_out_mouse2(antennas1, times1, antennas2, times2, normalization):
    idx = 0
    dominance_counter = 0
    while True:
        if idx == len(antennas1) or idx > len(antennas1):
            break
        m1_states, m1_readouts, idx = utils.get_more_states(antennas1,
                                                            times1,
                                                            idx,
                                                            mas,
                                                            how_many_antennas)
        if does_mouse1_push_out(m1_states, m1_readouts, antennas2, times2):
            dominance_counter += 1
    if normalization is None:
        return dominance_counter
    if normalization == "m1_activity":
        return dominance_counter/len(antennas1)
    if normalization == "m2_activity":
        return dominance_counter/len(antennas2)
    if normalization == "m1_m2_activity":
        return dominance_counter/len(antennas2)/len(antennas1)
    

def tube_dominance_single_phase(ehd, cf, phase, normalization):
    mice = ehd.mice
    t_start, t_end = cf.gettime(phase)
    dominance =  np.zeros((len(mice), len(mice)))
    for i, mouse1 in enumerate(mice):
        m1_times, m1_antennas = utils.get_times_antennas(ehd, mouse1,
                                                         t_start, t_end)
        for j, mouse2 in enumerate(mice):
            if i != j:
                m2_times, m2_antennas = utils.get_times_antennas(ehd, mouse2,
                                                                 t_start, t_end)

                dominance[i, j] = check_mouse1_pushing_out_mouse2(m1_antennas, m1_times,
                                                        m2_antennas, m2_times,
                                                        normalization)
    return dominance

def get_tube_dominance(ehd, cf, prefix="", res_dir="", normalization=None,
                       delimiter=";"):
    if normalization is None:
        fname = 'tube_dominance_no_normalization'
    else:
        fname = 'tube_dominance_%s' % normalization
    if prefix is "":
        prefix = ehd.prefix
    if res_dir is "":
        res_dir = ehd.res_dir

    if ehd.how_many_antennas == 2:
        dom2.get_tube_dominance_2_cages(ehd, cf, res_dir, prefix)
        dom2.get_subversion_evaluation(ehd, cf, res_dir, prefix)
        dom2.get_visits_to_stimulus_cage(ehd, cf, res_dir, prefix)
    dispatch.evaluate_whole_experiment(ehd, cf, res_dir, prefix,
                                       tube_dominance_single_phase,
                                       fname, 'dominating mouse',
                                       'pushed out mouse',
                                       '# dominances',
                                       args=[normalization], vmin=0, vmax=25,
                                       delimiter=delimiter)
