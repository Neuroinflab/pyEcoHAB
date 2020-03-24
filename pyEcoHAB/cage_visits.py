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
    added_interval = False

    if not len(interval_array):
        return visit_list, added_interval

    first = utils.get_idx_pre(t_start, interval_array[:, 0])
    if first is not None:
        if interval_array[first, 1] > t_start:
            added_interval = True
            if t_stop < interval_array[first, 1]:
                int_stop = t_stop
            else:
                int_stop = interval_array[first, 1]
            visit_list.append(int_stop - t_start)

    idx_between = utils.get_idx_between(t_start, t_stop, interval_array[:, 0])
    if len(idx_between) == 0:
        return visit_list, added_interval

    for idx in idx_between:
        i_start, i_stop = intervals[idx]
        if i_start >= t_stop:
            break
        if t_stop < i_stop:
            i_stop = t_stop
        visit_list.append(i_stop - i_start)

    return visit_list, added_interval


def get_visits_in_bins(intervals, time_start,
                       time_end, binsize):
    length = utils.get_length(time_start, time_end, binsize)
    visits = []
    added_visit = []
    for i in range(length):
        t_start = time_start + i*binsize
        t_end = time_start + (i+1)*binsize
        last_visits, outs = get_visits(intervals, t_start, t_end)
        visits.append(last_visits)
        added_visit.append(outs)
    return visits, added_visit

def calc_visit_per_mouse(intervals, t_start, t_end, binsize):
    visits_in_bins, added_visit = get_visits_in_bins(intervals,
                                                     t_start,
                                                     t_end,
                                                     binsize)
    visits = []
    for i, o in enumerate(visits_in_bins):
        visits.append(len(o) - added_visit[i])
    durations = [sum(o) for o in visits_in_bins]
    return visits, durations, visits_in_bins


def calculate_visits_and_durations(data, mice, address, t_start, t_end, binsize):
    visits = OrderedDict()
    durations = OrderedDict()
    all_visits = OrderedDict()
    for m in mice:
        ints = utils.get_intervals(data[m], address)
        visits[m], durations[m], all_visits[m] = calc_visit_per_mouse(ints,
                                                                      t_start,
                                                                      t_end,
                                                                      binsize)
    return visits, durations, all_visits


def get_activity(ehs, cf, binsize, res_dir="", prefix="", remove_mouse="",
                 save_histogram=False, delimiter=";",
                 headers=['Number of visits to box',
                          'Total time (sec) in box']):
    """
    Calculate activity of each mouse in time bins across the phases
    of the experiment.


    This function both counts visits of every mouse to each Eco-HAB compartment
    and calculates time spent in each compartment. It is based on visits
    calculated while reading in the experiment data based on antenna readings.
    Activity values are saved as comma separated values in a file with 
    activity_bin_{binsize}_h.csv in res_dir/activity directory.  This function
    also automatically saves activity to res_dir/approach_to_social to aid
    further analysis.

    Args:
        ehs : Loader or Loader_like
           Eco-HAB dataset.
        cf : ExperimentConfigFile
           timeline of the experiment.
        binsize : number 
           time bins for calculating activity.
           A number value specifies number of seconds in each bin, e.g. binsize
           equal 3600 results in 1 h bins.
        res_dir : string
           destination directory
           default value is the destination directory established for ehs.
        prefix : string
           string added to every name of the results file
        remove_mouse : string or list
           name of mouse or mice to be removed from the results file
        save_histogram :
           if True save visit durations to every bin and generate histograms
           of vist durations
        delimiter : str, optional
           String or character separating columns.
        headers : list of strings
           strings that will be written above activity parameters for each
           compartment.
           
    """
    if prefix == "":
        prefix = ehs.prefix
    if res_dir == "":
        res_dir = ehs.res_dir
    
    phases = utils.filter_dark_light(cf.sections())
    fname = '%sactivity_bin_%3.1f_h.csv'%(prefix,
                                          binsize//3600)
    histogram_fname = 'activity_histograms_bin_%3.1f_h' % (binsize//3600)
    mice = utils.get_mice(ehs.mice, remove_mouse)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    

    data = {c:{0:{},1:{}} for c in ehs.cages}
    ehs_data = utils.prepare_data(ehs, mice)
    bin_labels = {}
    for idx_phase, phase in enumerate(phases):
        t_start, t_end = cf.gettime(phase)
        visits_in_cages = {}
        for address in ehs.cages:
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
            
        if save_histogram:
            make_visit_duration_histogram(visits_in_cages,
                                          bin_labels[phase],
                                          phase, mice,
                                          histogram_fname, res_dir,
                                          "other_variables/visit_histograms_binsize_%3.1f"
                                          % (binsize//3600),
                                          prefix, add_info_mice)
            save_visit_duration(visits_in_cages,
                                bin_labels[phase],
                                phase, mice,
                                histogram_fname, res_dir,
                                "other_variables/visit_histograms_binsize_%3.1f"
                                % (binsize//3600),
                                prefix, add_info_mice)
    save_data_cvs(data, phases, mice, bin_labels, fname, res_dir, ehs.cages,
                  headers)
    save_data_cvs(data, phases, mice, bin_labels, fname, res_dir, ehs.cages,
                  headers, target_dir="social_approach")
