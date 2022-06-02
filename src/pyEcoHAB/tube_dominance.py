# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from .utils import general as utils
from .following import exec_fun
from . import dominance_in_2_cages as dom2


def tube_dominance_single_direction(intervals_m1, intervals_m2):
    t_starts_m1, t_ends_m1 = intervals_m1
    t_starts_m2, t_ends_m2 = intervals_m2
    counter = 0
    time_together = 0
    intervals = []
    #mouse2 enters the tunnel when mouse1 is already there
    for i in range(len(t_starts_m1)):
        indxs = utils.get_idx_between(t_starts_m1[i],
                                      t_ends_m1[i],
                                      t_starts_m2)
        for idx in indxs:
            if t_ends_m2[idx] < t_ends_m1[i]:
                counter += 1
                time_together += t_ends_m2[idx] - t_starts_m1[i]
                intervals.append(t_ends_m2[idx] - t_starts_m1[i])
                break
    #mouse1 enters the tunnel when mouse2 is already there
    for i in range(len(t_starts_m2)):
        indxs = utils.get_idx_between(t_starts_m2[i],
                                      t_ends_m2[i],
                                      t_starts_m1)
        for idx in indxs:
            if t_ends_m2[i] < t_ends_m1[idx]:
                counter += 1
                time_together += t_ends_m2[i] - t_starts_m1[idx]
                intervals.append(t_ends_m2[i] - t_starts_m1[idx])
                break

    return counter, time_together, intervals


def tube_dominance_single_pair(direction_m1, backing_m2):
    pushing = 0
    intervals = []
    time_together = 0

    for key_push in direction_m1.keys():
        key = key_push.split(" ")
        key_back = "%s %s" % (key[1], key[1])
        out = tube_dominance_single_direction(direction_m1[key_push],
                                              backing_m2[key_back])
        dom_single_dir, time_single_dir, ints_single_dir = out
        pushing += dom_single_dir
        time_together += time_single_dir
        intervals += ints_single_dir
    return pushing, time_together, intervals  


def tube_dominance_matrices(data, mice, t_start, t_stop):
    assert t_stop - t_start > 0
    durations = t_stop - t_start
    tube_dom = utils.make_results_dict(mice)
    time_dom = utils.make_results_dict(mice)
    labels = utils.all_mouse_pairs(mice)
    interval_details = {label: [] for label in labels}
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            print(data)
            out = tube_dominance_single_pair(data["directions"][mouse1],
                                             data["backing out"][mouse2])
            tube_dom[mouse1][mouse2], time_in_pipe, mouse_intervals = out
            time_dom[mouse1][mouse2] = time_in_pipe/durations
            key = "%s|%s" % (mouse1, mouse2)
            interval_details[key] += mouse_intervals
    return tube_dom, time_dom, interval_details
 

def get_tube_dominance(ecohab_data, timeline, N, binsize="whole_phase",
                       res_dir="", prefix="", remove_mouse=None,
                       save_distributions=True, save_figures=False,
                       return_median=False, delimiter=";",
                       save_times=False, seed=None,
                       full_dir_tree=True):
    data_prep = utils.prepare_for_tube_dominance
    exec_fun(ecohab_data, timeline, N, var_name="tube_dominance", action1_name="pushing",
             action2_name="pushed", data_prep=data_prep, function=tube_dominance_matrices,
             binsize=binsize, res_dir=res_dir, prefix=prefix, remove_mouse=remove_mouse,
             save_distributions=save_distributions, save_figures=save_figures,
             return_median=return_median, delimiter=delimiter, save_times=save_times,
             seed=seed, full_dir_tree=full_dir_tree)
    
