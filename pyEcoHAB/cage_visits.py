from __future__ import division, print_function, absolute_import
import numpy as np
import os
from collections import OrderedDict
from . import utility_functions as utils
from .write_to_file import save_data_cvs
from .plotting_functions import make_visit_interval_histogram


def get_visits(intervals, t_start, t_stop):
    interval_array = np.array(intervals)
    visit_list = []
    if not len(interval_array):
        return visit_list
    idx_between = utils.get_idx_between(t_start, t_stop, interval_array[:, 0])
    for idx in idx_between:
        i_start, i_stop = intervals[idx]
        if i_start >= t_stop:
            break
        visit_list.append(i_stop - i_start)
    return visit_list


def get_visits_in_bins(intervals, time_start,
                       time_end, binsize):
    length = utils.get_length(time_start, time_end, binsize)
    visits = []
    for i in range(length):
        t_start = time_start + i*binsize
        t_end = time_start + (i+1)*binsize
        visits.append(get_visits(intervals, t_start, t_end))
    return visits

def calculate_visits_and_durations(data, mice, address, t_start, t_end, binsize):
    visits = OrderedDict()
    durations = OrderedDict()
    all_visits = OrderedDict()
    for mouse in mice:
        intervals = utils.get_intervals(data[mouse], address)
        out = get_visits_in_bins(intervals, t_start,
                                        t_end, binsize)
        visits[mouse] = [float(len(o)) for o in out]
        durations[mouse] = [sum(o) for o in out]
        all_visits[mouse] = out
    return visits, durations, all_visits


def get_all_visits(ehs, cf, binsize, cages=None,
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
    histogram_fname = 'histograms_binsize_%f_h' % (binsize//3600)
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
        visits_in_cages = {}
        for address in cages.keys():
            visit_data = calculate_visits_and_durations(ehs_data,
                                                        mice,
                                                        cages[address],
                                                        t_start,
                                                        t_end,
                                                        binsize)
            data[address][0][phase] = visit_data[0]
            data[address][1][phase] = visit_data[1]
            data['time'][phase] = utils.get_times(binsize)
            visits_in_cages[address] = visit_data[2]
        if "dark" in phase or "DARK" in phase:
            make_visit_interval_histogram(visits_in_cages, phase, mice,
                                          histogram_fname, res_dir,
                                          "histograms_of_visits_binsize_%f" % (binsize//3600),
                                          prefix, add_info_mice)
    save_data_cvs(data, fname, res_dir,
                  cages,
                  headers)
