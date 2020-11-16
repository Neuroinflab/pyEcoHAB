# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import
import numpy as np
from . import utility_functions as utils
from .write_to_file import save_single_histograms, write_csv_rasters
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot
from .plotting_functions import single_heat_map
from . import exec_functions as dispatch
from . import dominance_in_2_cages as dom2


def mice_in_different_spots(states1, states2):
    for s1 in states1:
        if s1 in states2:
            return False
    return True


def does_mouse1_push_out(m1_states, m1_times, antennas2, times2, config):
    m2_states, m2_readouts = utils.get_states_and_readouts(antennas2, times2,
                                                           m1_times[0],
                                                           m1_times[-1])
    between = utils.get_idx_between(m1_times[0], m1_times[-1], times2)

    if len(between) == 0:
        return False

    if mice_in_different_spots(m1_states, m2_states):
        return False
    first_antenna = m1_states[0]
    other_antennas = config.other_tunnel_antenna(first_antenna)
    if not len(other_antennas):
        return False
    if len(other_antennas) == 1:
        opposite_antenna = other_antennas[0]
    else:
        for ant in other_antennas:
            if ant not in config.internal_antennas:
                other_antenna = ant
            else:
                # if there is only one entrance antenna to the tunnel,
                # there is no tube domination
                return False

    if opposite_antenna not in m2_states:
        return False
    opposite_idx = m2_states.index(opposite_antenna)
    if opposite_idx > 0:
        if m2_states[opposite_idx - 1] == first_antenna:
            return False
        elif m2_states[opposite_idx - 1] == config.address[first_antenna]:
            return False
    m2_m1_in_pipe = m2_states[opposite_idx:]
    idx_after = utils.get_idx_post(m1_times[-1], times2)
    if idx_after is not None:
        m2_after = antennas2[idx_after]
    else:
        m2_after = m1_states[0]
    if np.all(np.array(m2_m1_in_pipe) == opposite_antenna):
        if m2_after != m1_states[0]:
            return True
    if m2_readouts[opposite_idx] > m1_times[0]:
        if m1_states[0] not in m2_states[opposite_idx:]:
            return True
    return False


def check_mouse1_pushing(antennas1, times1, antennas2, times2,
                         config, normalization=None):
    if len(antennas1) < 2:
        return False
    idx = 1
    dominance_counter = 0
    while idx < len(antennas1):
        a1, a2 = antennas1[idx-1:idx+1]
        t1, t2 = times1[idx-1:idx+1]
        if a1 != a2:
            if config.same_tunnel[a1] == config.same_tunnel[a2]:
                temp_ants = [a1, a2]
                temp_times = [t1, t2]
                idx = idx + 1
                while idx < len(antennas1):
                    a3 = antennas1[idx]
                    t3 = times1[idx]
                    if a3 != a1 and a3 != a2:
                        break
                    temp_ants.append(a3)
                    temp_times.append(t3)
                    idx = idx + 1
                if idx == len(antennas1) or config.address[a1] != config.address[a3]:
                    dominance_counter += does_mouse1_push_out(temp_ants,
                                                              temp_times,
                                                              antennas2,
                                                              times2,
                                                              config)
        idx = idx + 1

    if normalization is None:
        return dominance_counter
    if normalization == "m1_activity":
        return dominance_counter/len(antennas1)
    if normalization == "m2_activity":
        return dominance_counter/len(antennas2)
    if normalization == "m1_m2_activity":
        return dominance_counter/len(antennas2)/len(antennas1)


def tube_dominance_single_phase(ehd, timeline, phase, normalization):
    mice = ehd.mice
    t_start, t_end = timeline.get_time_from_epoch(phase)
    dominance = np.zeros((len(mice), len(mice)))
    for i, mouse1 in enumerate(mice):
        m1_times, m1_antennas = utils.get_times_antennas(ehd, mouse1,
                                                         t_start, t_end)
        for j, mouse2 in enumerate(mice):
            if i != j:
                m2_times, m2_antennas = utils.get_times_antennas(ehd, mouse2,
                                                                 t_start,
                                                                 t_end)

                dominance[i, j] = check_mouse1_pushing(m1_antennas,
                                                       m1_times,
                                                       m2_antennas,
                                                       m2_times,
                                                       ehd.setup_config,
                                                       normalization)
    return dominance


def get_tube_dominance(ehd, timeline, prefix="", res_dir="", normalization=None,
                       delimiter=";"):
    if normalization is None:
        fname = 'tube_dominance_no_normalization'
    else:
        fname = 'tube_dominance_%s' % normalization
    if prefix == "":
        prefix = ehd.prefix
    if res_dir == "":
        res_dir = ehd.res_dir

    if len(ehd.setup_config.tunnels) == 1:
        dom2.get_tube_dominance_2_cages(ehd, timeline, res_dir, prefix)
        dom2.get_subversion_evaluation(ehd, timeline, res_dir, prefix)
        dom2.get_visits_to_stimulus_cage(ehd, timeline, res_dir, prefix)
    dispatch.evaluate_whole_experiment(ehd, timeline, res_dir, prefix,
                                       tube_dominance_single_phase,
                                       fname, 'dominating mouse',
                                       'pushed out mouse',
                                       '# dominances',
                                       args=[normalization], vmin=0, vmax=25,
                                       delimiter=delimiter)
