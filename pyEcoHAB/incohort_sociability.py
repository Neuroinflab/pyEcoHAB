# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import
import os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from . import utility_functions as utils
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot
from .write_to_file import write_binned_data, write_csv_rasters, write_csv_alone

def prepare_mice_intervals(data_mice, address):
    ints = {}
    for mouse in data_mice.keys():
        ints[mouse] = utils.intervals2lists(data_mice[mouse], address)
    return ints


def check_interval(intervals_mouse1, intervals_mouse2, idx, new_idx):
    original_s, original_e = intervals_mouse1[0][idx], intervals_mouse1[1][idx]
    other_s, other_e = intervals_mouse2[0][new_idx], intervals_mouse2[1][new_idx]
    if  other_s > original_e:
        return False
    if original_s >  other_s and original_e > other_e:
        intervals_mouse1[0][idx] = other_e
        return False
    if original_s < other_s and original_e < other_e:
        intervals_mouse1[1][idx] = other_s
        return False
    if original_s >= other_s and original_e <= other_e:
        # delete mouse1 interval
        intervals_mouse1[0].remove(original_s)
        intervals_mouse1[1].remove(original_e)
        return True
    #  cut the original interval in half
    intervals_mouse1[1][idx] = other_s
    if other_e < original_e:
        intervals_mouse1[0].insert(idx + 1, other_e)
        intervals_mouse1[1].insert(idx + 1, original_e)
    if other_s == original_s:
        intervals_mouse1[0].remove(original_s)
        intervals_mouse1[1].remove(other_s)
        return True
    return False

def remove_overlapping_intervals(intervals_mouse1, intervals_mouse2):
    """"
    Eliminate all the intervals, when mouse 1 is with mouse 2
    """
    if len(intervals_mouse1[0]) == 0:
        return
    if len(intervals_mouse2[0]) == 0:
        return
    i = 0
    while i < len(intervals_mouse1[0]):
        new_idx = utils.get_idx_pre(intervals_mouse1[0][i], intervals_mouse2[0])
        if new_idx is None:
            new_idx = utils.get_idx_pre(intervals_mouse1[1][i], intervals_mouse2[0])
            if new_idx is None:
                i = i + 1
                continue
        if intervals_mouse2[1][new_idx] < intervals_mouse1[0][i]:
            new_idx += 1
        if new_idx  >= len(intervals_mouse2[0]):
            break
        removed = check_interval(intervals_mouse1, intervals_mouse2, i, new_idx)
        if not removed:
            i = i + 1


def mouse_alone(data_mice, address):
    ints = prepare_mice_intervals(data_mice, address)
    new_intervals = {}
    result = {}
    for mouse in ints.keys():
        mouse_ints = []
        mouse_ints.append(ints[mouse][0][:])
        mouse_ints.append(ints[mouse][1][:])
        for other_mouse in ints.keys():
            if mouse != other_mouse:
                remove_overlapping_intervals(mouse_ints, ints[other_mouse])
        new_intervals[mouse] = mouse_ints
    for mouse in new_intervals:
        result[mouse] = utils.get_duration(new_intervals[mouse][0], new_intervals[mouse][1])
    return result

def make_solitude_output(addresses, mice):
    output = OrderedDict()
    for address in addresses:
        output[address] = OrderedDict()
        for mouse in mice:
            output[address][mouse] = OrderedDict()
    return output

def get_solitude(ehs, cf, res_dir="", prefix=""):
    if prefix is "":
        prefix = ehs.prefix
    if res_dir is "":
        res_dir = ehs.res_dir
    phases = utils.filter_dark(cf.sections())
    output = make_solitude_output(ehs.cages, ehs.mice)

    for phase in phases:
        times = cf.gettime(phase)
        data = utils.prepare_data(ehs, ehs.mice, times)
        for address in ehs.cages:
            alone = mouse_alone(data, address)
            for mouse in ehs.mice:
                output[address][mouse][phase] = alone[mouse]
    write_csv_alone(output, phases, res_dir, prefix, delimiter=delimiter)


def mice_overlap(ints1, ints2):
    """Return time overlap of mice m1 and m2 in cage <address>."""
    total_overlap = 0
    for int1 in ints1:
        for int2 in ints2:
            total_overlap += utils.interval_overlap(int1, int2)
    return total_overlap


def time_fraction_together_one_cage(ints1, ints2, total_time):
    assert total_time > 0
    return mice_overlap(ints1, ints2)/total_time


def expected_time_fraction_together_one_cage(ints1, ints2, total_time):
    durations_m1 = utils.calculate_total_duration(ints1)
    durations_m2 = utils.calculate_total_duration(ints2)
    return durations_m1/total_time*durations_m2/total_time


def mice_together(data_mice, m1, m2, addresses, total_time):
    """Return the time spent together by two mice and expected time 
    assuming independence."""
    time_together = 0
    exp_time_together = 0
    for address in addresses:
        ints1 = utils.get_intervals(data_mice[m1], address)
        ints2 = utils.get_intervals(data_mice[m2], address)
        time_together += time_fraction_together_one_cage(ints1,
                                                         ints2,
                                                         total_time)
        exp_time_together += expected_time_fraction_together_one_cage(ints1,
                                                                      ints2,
                                                                      total_time)
    return time_together, exp_time_together

def make_results_dict(mice):
    result = OrderedDict()
    for mouse1 in mice:
        result[mouse1] = OrderedDict()
        for mouse2 in mice:
            result[mouse1][mouse2] = 0

    return result

def single_phase_results(data, mice, addresses, total_time):
    res = make_results_dict(mice)
    res_exp = make_results_dict(mice)
    for ii, m1 in enumerate(mice):
        for jj in range(ii + 1, len(mice)):
            m2 = mice[jj]
            res[m1][m2], res_exp[m1][m2] = mice_together(data,
                                                         m1,
                                                         m2,
                                                         addresses,
                                                         total_time)
    return res, res_exp


def get_dark_light_data(phase, cf, ehs, mice):

    if phase == "dark" or phase == "DARK" or phase == "Dark":
        phases = utils.filter_dark(cf.sections())
    elif phase == "light" or phase == "LIGHT" or phase == "Light":
        phases = utils.filter_light(cf.sections())
    out_phases = [phase]
    data = {mouse:[] for mouse in mice}
    total_time = 0
    for i, ph in enumerate(phases):
        time = cf.gettime(ph)
        out = utils.prepare_data(ehs, mice, time)
        for mouse in mice:
            data[mouse].extend(out[mouse])
        total_time += (time[1] - time[0])
    out_data = {phase: {0: data}}
    return out_phases, {phase: {0: total_time}}, {phase: {0: data}}

def prepare_fnames_and_totals(ehs, cf, prefix, bins, mice, filter_dark=True):
    if bins in ["ALL", "all", "All"]:
        phases = ["ALL"]
        time = cf.gettime("ALL")
        total_time = {"ALL": {0: (time[1] - time[0])}}
        data = {"ALL": {0: utils.prepare_data(ehs, mice, time)}}
        keys = [["ALL"], [0]]
    elif bins in ['dark', "DARK", "Dark", "light", "LIGHT", "Light"]:
        phases, total_time, data = get_dark_light_data(bins, cf, ehs, mice)
        keys = [list(data.keys()), [0]]
    elif isinstance(bins, int) or isintance(bins, float):
        phases = []
        data = OrderedDict()
        total_time = OrderedDict()
        if filter_dark:
            all_phases = utils.filter_dark(cf.sections())
        else:
            all_phases = utils.filter_dark_light(cf.sections())
        bin_labels = utils.get_times(bins)
        for phase in all_phases:
            t_start, t_end = cf.gettime(phase)
            phases.append("%s_%5.2fh" % (phase, bins))
            data[phase] = OrderedDict()
            total_time[phase] = OrderedDict()
            j = 0
            while t_start < t_end:
                t_e = t_start + bins
                if t_e > t_end:
                    t_e = t_end
                time = [t_start, t_e]
                data[phase][bin_labels[j]] = utils.prepare_data(ehs, mice, time)
                total_time[phase][bin_labels[j]] = time[1] - time[0]
                t_start += bins
                j += 1
        keys = [all_phases, bin_labels]

   
    return phases, total_time, data, keys


def make_all_results_dict(phases, bins):
    result = OrderedDict()
    for phase in phases:
        result[phase] = OrderedDict()
        for bin1 in bins:
            result[phase][bin1] = 0

    return result


def get_incohort_sociability(ehs, cf, binsize, res_dir="",
                             prefix="", remove_mouse="", delimiter=";"):

    """
    Calculate in-cohort sociability for each pair of mice in time bins across
    the phases of the experiment.

    In-cohort sociability is a measure of sociability unique to the Eco-HAB
    system. It evaluates time spent together by each pair of mice in each
    of Eco-HAB compartments, taking into account expected time spent
    together by that pair of mice based on mice preference of Eco-HAB
    compartments.

    in-cohort sociability results for each dark phase of the experiment
    are saved in two formats: as csv files with in-cohort sociability values
    for each mouse pair, where the mouse pair is specified by row and column
    in results/incohort_sociability/bins_{bin_length}_h/histograms/data,
    and csv files where the pair of mice is specified by row and column
    and raster files with in-cohort sociability calculated for all mouse pairs
    in results/incohort_


    Args:
        ehs : Loader or Loader_like
           Eco-HAB dataset.
        cf : ExperimentConfigFile
           timeline of the experiment.
        binsize : string or number 
           time bins for calculating activity. Possible string values are:
           "ALL" -- calculate activity for the whole experiment,
           "dark" -- calculate activity for all dark phases,
           "light" -- calculate activity for all light phases.
           A number value specifies number of seconds in each bin, e.g. binsize
           equal 3600 results in 1 h bins.
        res_dir : string
           destination directory
           default value is the destination directory established for ehs.
        prefix : string
           string added to the name of every generated results file
           default value is the prefix established for ehs
        remove_mouse : string or list
           name of mouse or mice to be removed from the results file
           As a default activity will be established for every mouse registered
           in ehs.
        delimiter : str, optional
           String or character separating columns.
    """
    if prefix == "":
        prefix = ehs.prefix
    if res_dir == "":
        res_dir = ehs.res_dir
    mice = utils.get_mice(ehs.mice, remove_mouse)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)

   
    fname_measured_prefix = "incohort_sociability_measured_time_%s_%s" % (prefix,
                                                                     add_info_mice)
    fname_expected_prefix = "incohort_sociability_expected_time_%s_%s" % (prefix,
                                                                          add_info_mice)
    fname_excess_prefix = "incohort_sociability_excess_time_%s_%s" % (prefix,
                                                                      add_info_mice)
    phases, time, data, keys = prepare_fnames_and_totals(ehs,
                                                         cf,
                                                         prefix,
                                                         binsize,
                                                         mice,
                                                         filter_dark)
    if isinstance(binsize, int) or isinstance(binsize, float):
        binsize_name = "%3.2f_h" % (binsize/3600)
        if binsize == 43200:
            csv_results_incohort = np.zeros((len(phases), len(mice),
                                             len(mice)))
            csv_results_incohort_exp = np.zeros((len(phases), len(mice),
                                                 len(mice)))
    else:
        binsize_name = binsize
    if time == 0:
        return
    full_results = make_all_results_dict(*keys)
    full_results_exp = make_all_results_dict(*keys)
    out_dict_hist = os.path.join("incohort_sociability", "histograms",
                                 "bins_%s" % binsize_name)
    out_dict_rasters = os.path.join("incohort_sociability", "raster_plots",
                                    "bins_%s" % binsize_name)
    all_phases, bin_labels = keys[0], keys[1]
    for idx_phase, ph in enumerate(all_phases):
        new_phase = phases[idx_phase].replace(' ', '_')
        for lab in bin_labels:
            full_results[ph][lab],\
                full_results_exp[ph][lab] = single_phase_results(data[ph][lab],
                                                                 mice,
                                                                 ehs.cages,
                                                                 time[ph][lab])

        write_binned_data(full_results[ph],
                          'incohort_sociability_measured_time',
                          mice, bin_labels, new_phase, res_dir, 
                          out_dict_hist,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter)
        write_binned_data(full_results_exp[ph],
                          'incohort_sociability_expected_time',
                          mice, bin_labels, new_phase, res_dir, 
                          out_dict_hist,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter)
        excess_time = utils.calc_excess(full_results[ph],
                                        full_results_exp[ph])

        write_binned_data(excess_time,
                          'incohort_sociability_excess_time',
                          mice, bin_labels, new_phase, res_dir, 
                          out_dict_hist, prefix, additional_info=add_info_mice,
                          delimiter=delimiter)
        if isinstance(binsize, int) or isinstance(binsize, float):
            if int(binsize) == 12*3600:
                fname = "incohort_sociability_"
                res = utils.dict_to_array_2D(full_results[ph][0],
                                             mice, mice)
                exp_res = utils.dict_to_array_2D(full_results_exp[ph][0],
                                                 mice, mice)
                single_in_cohort_soc_plot(res,
                                          exp_res,
                                          mice,
                                          new_phase,
                                          fname,
                                          res_dir,
                                          out_dict_hist,
                                          prefix+add_info_mice)
                csv_results_incohort[idx_phase] = res
                csv_results_incohort_exp[idx_phase] = exp_res

        fname_measured = "%s_%s.csv" % (fname_measured_prefix, new_phase)
        fname_excess = "%s_%s.csv" % (fname_excess_prefix, new_phase)
        fname_expected = "%s_%s.csv" % (fname_expected_prefix, new_phase)
        raster_labels = [bin_label/3600 for bin_label in bin_labels]
        phase_full_results = utils.dict_to_array_3D(full_results[ph],
                                                    bin_labels,
                                                    mice, mice)
        phase_exp_full_results = utils.dict_to_array_3D(full_results_exp[ph],
                                                        bin_labels,
                                                        mice, mice)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results,
                          res_dir,
                          out_dict_rasters,
                          fname_measured,
                          delimiter=delimiter)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_exp_full_results,
                          res_dir,
                          out_dict_rasters,
                          fname_expected,
                          delimiter=delimiter)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results - phase_exp_full_results,
                          res_dir,
                          out_dict_rasters,
                          fname_excess, delimiter=delimiter)
      
    if isinstance(binsize, int) or isinstance(binsize, float):
        if binsize == 43200:
            write_csv_rasters(mice,
                              all_phases,
                              csv_results_incohort,
                              res_dir,
                              out_dict_rasters,
                              "incohort_sociability_measured_time_ALL_phases_binned.csv",
                              delimiter=delimiter)
            write_csv_rasters(mice,
                              all_phases,
                              csv_results_incohort_exp,
                              res_dir,
                              out_dict_rasters,
                              "incohort_sociability_expected_time_ALL_phases_binned.csv",
                              delimiter=delimiter)
            write_csv_rasters(mice,
                              all_phases,
                              csv_results_incohort - csv_results_incohort_exp,
                              res_dir,
                              out_dict_rasters,
                               "incohort_sociability_excess_time_ALL_phases_binned.csv",
                              delimiter=delimiter)
            make_RasterPlot(res_dir,
                            out_dict_rasters,
                            csv_results_incohort,
                            all_phases,
                            "incohort_sociability_measured_time_ALL_phases_binned",
                            mice,
                            prefix=prefix,
                            to_file=True,
                            vmin=-1,
                            vmax=1,
                            title="Measured in-cohort sociability",
                            symmetric=True)
            make_RasterPlot(res_dir,
                            out_dict_rasters,
                            csv_results_incohort_exp,
                            all_phases,
                            "incohort_sociability_expected_time_ALL_phases_binned",
                            mice,
                            prefix=prefix,
                            to_file=True,
                            vmin=-1,
                            vmax=1,
                            title="Expected in-cohort sociability",
                            symmetric=True)
            make_RasterPlot(res_dir,
                            out_dict_rasters,
                            csv_results_incohort-csv_results_incohort_exp,
                            all_phases,
                            "incohort_sociability_excess_time_ALL_phases_binned",
                            mice,
                            prefix=prefix,
                            to_file=True,
                            vmin=-1,
                            vmax=1,
                            title="Excess in-cohort sociability",
                            symmetric=True)
