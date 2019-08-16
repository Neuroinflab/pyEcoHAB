from __future__ import division, print_function, absolute_import
import numpy as np
import os
from collections import OrderedDict
from . import utility_functions as utils
from .write_to_file import save_data_cvs

def visits_and_durations(intervals, t_start, t_stop):
    visits, durations = 0, 0
    interval_array = np.array(intervals)
    if not len(interval_array):
        return visits, durations
    idx_between = utils.get_idx_between(t_start, t_stop, interval_array[:, 0])
    for idx in idx_between:
        i_start, i_stop = intervals[idx]
        if i_start >= t_stop:
            break
        visits += 1
        durations += i_stop - i_start
    return visits, durations


def visits_and_durations_bins(intervals, time_start,
                              time_end, binsize):
    length = utils.get_length(time_start, time_end, binsize)
    visits = np.zeros(length, dtype=int)
    durations = np.zeros(length)
    time = time_start
    for i in range(length):
        t_start = time_start + i*binsize
        t_end = time_start + (i+1)*binsize
        visits[i], durations[i] = visits_and_durations(intervals,
                                                       t_start, t_end)
    return visits, durations

def calculate_visits_and_durations(data, mice, address, t_start, t_end, binsize):
    visits = OrderedDict()
    durations = OrderedDict()
    for mouse in mice:
        intervals = utils.get_intervals(data[mouse], address)
        visits[mouse], durations[mouse] = visits_and_durations_bins(intervals,
                                                                    t_start,
                                                                    t_end,
                                                                    binsize)
    return visits, durations


def get_visits(ehs, cf, binsize, cages=None,
               res_dir=None, prefix=None,
               remove_mouse=None, headers=None):
    if prefix is None:
        prefix = ehs.prefix
    if res_dir is None:
        res_dir = ehs.res_dir
    if headers is None:
        basic = ['Number of visits to box %d\n',
                 'Total time in box %d, seconds\n']
    
    phases = utils.filter_dark_light(cf.sections())
    fname = '%scollective_results_all_chambers_binsize_%f_h.csv'%(prefix,
                                                                  binsize//3600)
    mice = utils.get_mice(ehs.mice, remove_mouse)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    if cages is None:
        cages = {'1': 1, '3': 3, '2': 2, '4': 4}
    headers = {str(i):basic for i in cages.keys()}
    data = {c:{0:{},1:{}} for c in cages.keys()}
    data['mice'] = mice
    data['phases'] = phases
    data['time'] = {}
    ehs_data = utils.prepare_data(ehs, mice)
    for idx_phase, phase in enumerate(phases):
        t_start, t_end = cf.gettime(phase)
        for address in cages.keys():
            visits, durations = calculate_visits_and_durations(ehs_data,
                                                               mice,
                                                               cages[address],
                                                               t_start,
                                                               t_end,
                                                               binsize)
            data[address][0][phase] = visits
            data[address][1][phase] = durations
            data['time'][phase] = utils.get_times(binsize)
    save_data_cvs(data, fname, res_dir,
                  cages,
                  headers)
