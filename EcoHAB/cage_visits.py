from __future__ import division, print_function
import EcoHab
import numpy as np
import os
import utility_functions as utils
from write_to_file import save_data_cvs
from ExperimentConfigFile import ExperimentConfigFile


def chamber_visits_and_durations(intervals, t_start, t_end):
    interval_array = np.array(intervals)

    pre = utils.get_idx_between(t_start, t_end,
                                interval_array[:, 0])
    post = utils.get_idx_between(t_start, t_end,
                                 interval_array[:, 1])
    if len(pre) == 0 or len(post) == 0:
        return 0, 0
    #check if a visit started before
    visits, durations = 0, 0
    if post[0] < pre[0]:
        durations += interval_array[0, 1] - t_start

    #check if a visit finished after
    if pre[-1] > post[-1]:
        durations += t_end - interval_array[-1, 0]
        visits += 1
    intervals_in_bin = interval_array[pre[0]:post[-1]+1, :]
    durations += sum(intervals_in_bin[:, 1] - intervals_in_bin[:, 0])
    visits += len(intervals_in_bin)
    return visits, durations


def get_one_phase(ehd, cf, phase, address, binsize):
    t_start, t_end = cf.gettime(phase)
    mice = ehd.mice
    length = int(np.ceil((t_end - t_start)/binsize/3600))
    visits = np.zeros(len(mice), length)
    durations = np.zeros(len(mice), length)
    data = utils.get_ehs_data(ehs, mouse, t_start, t_end)
    for address in range(1, 5):
        intervals = utils.get_intervals(data, address)
    


def get_cage_visits(ehs, cf,
                    res_dir=None,
                    prefix=None,
                    which_phases=None,
                    remove_mouse=None):
    if prefix is None:
        prefix = ehd.prefix
    if res_dir is None:
        res_dir = ehd.res_dir
    mice = utils.get_mice(ehs.mice, remove_mouse)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    full_results = np.zeros((len(phases), len(mice), len(mice)))
    full_results_exp = np.zeros((len(phases), len(mice), len(mice)))
    for idx_phase, phase in enumerate(phases):
        if get_data:
            data = utils.prepare_data(ehs, mice, cf.gettime(phase))
