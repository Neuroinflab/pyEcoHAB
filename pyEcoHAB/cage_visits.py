from __future__ import division, print_function, absolute_import
import numpy as np
import os
from collections import OrderedDict
from . import utility_functions as utils
from .write_to_file import save_data_cvs, save_visit_duration
from .plotting_functions import make_visit_duration_histogram


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


def get_activity(ehs, cf, binsize, cages=None,
                   res_dir=None, prefix=None,
                   remove_mouse=None, headers=None):
    if prefix is None:
        prefix = ehs.prefix
    if res_dir is None:
        res_dir = ehs.res_dir
    if headers is None:
        basic = ['Number of visits to box %s\n',
                 'Total time in box %s, seconds\n']
    
    phases = utils.filter_dark_light(cf.sections())
    fname = '%sactivity_bin_%3.1f_h.csv'%(prefix,
                                                                  binsize//3600)
    histogram_fname = 'activity_histograms_bin_%3.1f_h' % (binsize//3600)
    mice = utils.get_mice(ehs.mice, remove_mouse)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    if cages is None:
        cages = OrderedDict([(4, "A"), (1, "B"), (2, "C"), (3, "D")])

    headers = {i:basic for i in cages.keys()}
    data = {c:{0:{},1:{}} for c in cages.keys()}
    
    ehs_data = utils.prepare_data(ehs, mice)
    bin_labels = {}
    for idx_phase, phase in enumerate(phases):
        t_start, t_end = cf.gettime(phase)
        visits_in_cages = {}
        for address in cages.keys():
            visit_data = calculate_visits_and_durations(ehs_data,
                                                        mice,
                                                        address,
                                                        t_start,
                                                        t_end,
                                                        binsize)
            data[address][0][phase] = visit_data[0]
            data[address][1][phase] = visit_data[1]
            bin_labels[phase] = utils.get_times(binsize)
            visits_in_cages[address] = visit_data[2]
        if "dark" in phase or "DARK" in phase:
             make_visit_duration_histogram(visits_in_cages,
                                           bin_labels[phase],
                                           phase, mice, cages,
                                           histogram_fname, res_dir,
                                           "other_variables/visit_histograms_binsize_%3.1f"
                                           % (binsize//3600),
                                           prefix, add_info_mice)
             save_visit_duration(visits_in_cages,
                                 bin_labels[phase],
                                 phase, mice, cages,
                                 histogram_fname, res_dir,
                                 "other_variables/visit_histograms_binsize_%3.1f"
                                 % (binsize//3600),
                                 prefix, add_info_mice)
    save_data_cvs(data, phases, mice, bin_labels, fname, res_dir,
                  cages, headers)
    save_data_cvs(data, phases, mice, bin_labels, fname, res_dir,
                  cages, headers, target_dir="social_approach")
