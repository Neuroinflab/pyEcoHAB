from __future__ import division, print_function
import EcoHab
import numpy as np
import os
import utility_functions as utils
from write_to_file import save_data_cvs
from ExperimentConfigFile import ExperimentConfigFile
from collections import OrderedDict

def visits_and_durations(intervals, t_start, t_stop):
    visits, durations = 0, 0
    interval_array = np.array(intervals)
    if not len(interval_array):
        return visits, durations
    idx_pre = utils.get_idx_pre(t_start, interval_array[:, 0])

    #check visit before t_start
    if idx_pre is not None:
        start_pre, stop_pre = intervals[int(idx_pre)]
        if  stop_pre > t_start:
            visits += 1
            if stop_pre > t_stop:
                durations += t_stop - t_start
                return visits, durations
            durations +=  stop_pre - t_start

    idx_between = utils.get_idx_between(t_start, t_stop, interval_array[:, 0])

    for idx in idx_between:
        i_start, i_stop = intervals[idx]
        if i_start >= t_stop:
            break
        visits += 1
        if i_stop < t_stop:
            durations += i_stop - i_start
        else:
            durations += t_stop - i_start
    return visits, durations


def visits_and_durations_bins(intervals, time_start,
                              time_end, binsize):
    length = utils.get_length(time_start, time_end, binsize)
    visits = np.zeros(length)
    durations = np.zeros(length)
    time = time_start
    for i in range(length):
        t_start = time_start + i*binsize
        t_end = time_start + (i+1)*binsize
        visits[i], durations[i] = visits_and_durations(intervals,
                                                       t_start, t_end)
    return visits, durations


# def get_cage_visits(ehs, cf,
#                     res_dir=None,
#                     prefix=None,
#                     which_phases=None,
#                     remove_mouse=None):
#     if prefix is None:
#         prefix = ehd.prefix
#     if res_dir is None:
#         res_dir = ehd.res_dir
#     mice = utils.get_mice(ehs.mice, remove_mouse)
#     add_info_mice = utils.add_info_mice_filename(remove_mouse)
#     full_results = np.zeros((len(phases), len(mice), len(mice)))
#     full_results_exp = np.zeros((len(phases), len(mice), len(mice)))
#Get intervals for the whole experiment and then loop?
#     for idx_phase, phase in enumerate(phases):
#             t_start, t_end = cf.gettime(phase)
#             data = utils.get_ehs_data(ehs, mouse, t_start, t_end)
#             for address in range(1, 5):
#                 intervals = utils.get_intervals(data, address)
