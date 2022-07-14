# SPDX-License-Identifier: LGPL-2.1-or-later
# !/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
from collections import OrderedDict

from . import get_activity
from .utils import general as utils
from .utils import random_data_generation as rdg
from .write_to_file import write_csv_rasters
from .write_to_file import write_interpair_intervals
from .write_to_file import write_bootstrap_results
from .write_to_file import write_sum_data
from .write_to_file import write_binned_data
from .plotting_functions import single_in_cohort_soc_plot
from .plotting_functions import make_RasterPlot
from .plotting_functions import pooled_hists
from .plotting_functions import make_histograms_for_every_mouse
from .plotting_functions import pooled_hists_for_every_mouse
from .plotting_functions import single_histogram_figures


def bootstrap_single_phase(direction_dict, mice_list,
                           t_start, t_stop, function):

    action_count = utils.make_results_dict(mice_list, tolist=True)
    times_together = utils.make_results_dict(mice_list, tolist=True)
    for i, new_directions in enumerate(direction_dict):
        out = function(new_directions, mice_list,
                       t_start, t_stop)
        for m1 in mice_list:
            for m2 in mice_list:
                if m1 != m2:
                    action_count[m1][m2].append(out[0][m1][m2])
                    times_together[m1][m2].append(out[1][m1][m2])
    return action_count, times_together


def resample_single_phase(directions_dict, mice, t_start, t_stop, phase,
                          var_name, function, return_median=False,
                          save_figures=False,
                          save_distributions=True, res_dir=None, prefix=None,
                          stf=False, full_dir_tree=True):
    """If return_median is False, function returns mean value
    of the resampled following distribution

    stf: save times following"""
    N = len(directions_dict)
    if res_dir is None:
        res_dir = ecohab_data.res_dir
    if prefix is None:
        prefix = ecohab_data.prefix
    suff = (var_name.split("_")[0][0]+var_name.split("_")[1][0]).upper()
    name = var_name.replace("_", " ")
    actions, action_durations = bootstrap_single_phase(directions_dict,
                                                       mice,
                                                       t_start, t_stop,
                                                       function)
    binsize = (t_stop - t_start)/3600
    hist_dir = os.path.join("other_variables",
                            "%s_hists" % var_name,
                            "bin_%4.2f" % binsize)
    hist_time_dir = os.path.join("other_variables",
                                 "durations_%s_hists" % var_name,
                                 "bin_%4.2f" % binsize)
    fname_following = "%s_%s_count_dist_%d_%4.2f" % (prefix, suff, N, binsize)
    fname_times = "%s_%s_durations_dis_%d_%4.2f" % (prefix, suff, N, binsize)
    if save_figures:
        for mouse1 in mice:
            for mouse2 in mice:
                if mouse1 == mouse2:
                    continue
                key = "%s_%s" % (mouse1, mouse2)
                fname1 = "%s%s_hist_%s_%s_N_%d_%4.2f" % (prefix, suff,
                                                         phase.replace(' ',
                                                                       '_'),
                                                         key, N, binsize)
                fname2 = "%s%s_hist_%s_%s_N_%d_%4.2f" % (prefix,
                                                         "durations_%s" % suff,
                                                         phase.replace(' ',
                                                                       '_'),
                                                         key, N, binsize)
                single_histogram_figures(actions[mouse1][mouse2],
                                         fname1, res_dir,
                                         hist_dir,
                                         "%s count" % name,
                                         xlabel=name,
                                         ylabel="count",
                                         median_mean=True,
                                         full_dir_tree=full_dir_tree)
                if stf:
                    single_histogram_figures(action_durations[mouse1][mouse2],
                                             fname2,
                                             res_dir,
                                             hist_time_dir,
                                             "%s durations" % name,
                                             xlabel="duration",
                                             ylabel="count", nbins=10,
                                             median_mean=True,
                                             full_dir_tree=full_dir_tree)
    dist_dir_fol = os.path.join("other_variables",
                                "%s_hists" % var_name,
                                "bin_%4.2f" % binsize)
    dist_dir_time = os.path.join("other_variables",
                                 "durations_%s_hists" % var_name,
                                 "bin_%4.2f" % binsize)
    if save_distributions:
        write_bootstrap_results(actions, phase, mice,
                                fname_following, res_dir,
                                dist_dir_fol, prefix,
                                full_dir_tree=full_dir_tree)
        if stf:
            write_bootstrap_results(action_durations, phase, mice,
                                    fname_times, res_dir,
                                    dist_dir_time,
                                    prefix,
                                    full_dir_tree=full_dir_tree)
    out_actions = utils.make_results_dict(mice)
    out_times = utils.make_results_dict(mice)
    for m1 in mice:
        for m2 in mice:
            if m1 == m2:
                continue
            if return_median:
                out_actions[m1][m2] = np.median(actions[m1][m2])
                out_times[m1][m2] = np.median(action_durations[m1][m2])
            else:
                out_actions[m1][m2] = np.mean(actions[m1][m2])
                out_times[m1][m2] = np.mean(action_durations[m1][m2])

    return out_actions, out_times


def following_single_pair(direction_m1, direction_m2):
    followings = 0
    intervals = []
    time_together = 0
    assert direction_m1.keys() == direction_m2.keys()
    for key in direction_m1.keys():
        out = following_single_direction(direction_m1[key],
                                         direction_m2[key])
        f_single_dir, time_single_dir, ints_single_dir = out
        followings += f_single_dir
        time_together += time_single_dir
        intervals += ints_single_dir
    return followings, time_together, intervals


def following_single_direction(intervals_m1, intervals_m2):
    t_starts_m1, t_ends_m1 = intervals_m1
    t_starts_m2, t_ends_m2 = intervals_m2
    counter = 0
    time_together = 0
    intervals = []
    for i in range(len(t_starts_m1)):
        indxs = utils.get_idx_between(t_starts_m1[i],
                                      t_ends_m1[i],
                                      t_starts_m2)
        for idx in indxs:
            if t_ends_m2[idx] > t_ends_m1[i]:
                counter += 1
                time_together += t_ends_m1[i] - t_starts_m2[idx]
                intervals.append(t_ends_m2[idx] - t_starts_m1[i])
                break
    return counter, time_together, intervals


def add_intervals(all_intervals, phase_intervals):
    for mouse in phase_intervals.keys():
        all_intervals[mouse].extend(phase_intervals[mouse])


def following_matrices(directions_dict, mice, t_start, t_stop):
    assert t_stop - t_start > 0
    durations = t_stop - t_start
    followings = utils.make_results_dict(mice)
    time_together = utils.make_results_dict(mice)
    labels = utils.all_mouse_pairs(mice)
    interval_details = {label: [] for label in labels}
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            out = following_single_pair(directions_dict[mouse1],
                                        directions_dict[mouse2])
            followings[mouse1][mouse2], time_in_pipe, mouse_intervals = out
            time_together[mouse1][mouse2] = time_in_pipe/durations
            key = "%s|%s" % (mouse1, mouse2)
            interval_details[key] += mouse_intervals
    return followings, time_together, interval_details


def get_dynamic_interactions(ecohab_data, timeline, N, binsize="whole_phase",
                             res_dir="", prefix="", remove_mouse=None,
                             save_distributions=False, save_figures=False,
                             return_median=False, delimiter=";",
                             save_times=False, seed=None,
                             full_dir_tree=True):

    return exec_fun(ecohab_data, timeline, N, var_name="dynamic_interactions",
                    action1_name="leading", action2_name="following",
                    data_prep=utils.prepare_tube_data,
                    function=following_matrices, binsize=binsize,
                    res_dir=res_dir, prefix=prefix, remove_mouse=remove_mouse,
                    save_distributions=save_distributions,
                    save_figures=save_figures, return_median=return_median,
                    delimiter=delimiter, save_times=save_times, seed=seed,
                    full_dir_tree=full_dir_tree)


def exec_fun(ecohab_data, timeline, N, var_name, action1_name,
             action2_name, data_prep, function, binsize,
             res_dir, prefix, remove_mouse,
             save_distributions, save_figures,
             return_median, delimiter,
             save_times, seed,
             full_dir_tree):

    if isinstance(seed, int):
        np.random.seed(seed)
    if res_dir == "":
        res_dir = ecohab_data.res_dir
    if prefix == "":
        prefix = ecohab_data.prefix
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    name = var_name.replace("_", " ")
    mice = utils.get_mice(ecohab_data.mice, remove_mouse)

    phases, times, data, data_keys = utils.get_registrations_bins(ecohab_data,
                                                                  timeline,
                                                                  binsize,
                                                                  mice,
                                                                  data_prep)

    surrogate_data = rdg.generate_surrogate_data(ecohab_data, timeline,
                                                 binsize, mice, N,
                                                 data_prep)

    sur_data_list = rdg.reshape_surrogate_data(surrogate_data)
    all_phases, bin_labels = data_keys
    action_m = utils.make_all_results_dict(*data_keys)
    action_m_exp = utils.make_all_results_dict(*data_keys)
    time_together = utils.make_all_results_dict(*data_keys)
    time_together_exp = utils.make_all_results_dict(*data_keys)

    if isinstance(binsize, int) or isinstance(binsize, float):
        binsize_name = "%3.2f_h" % (binsize/3600)
        if int(binsize) == 24*3600:
            csv_results_action_m = np.zeros((len(phases), len(mice),
                                            len(mice)))
            csv_results_action_m_exp = np.zeros((len(phases), len(mice),
                                                len(mice)))
            if save_times:
                csv_results_time = np.zeros((len(phases), len(mice),
                                             len(mice)))
                csv_results_time_exp = np.zeros((len(phases), len(mice),
                                                 len(mice)))
    elif isinstance(binsize, str) and binsize.lower() in ["whole_phase",
                                                          "whole phase",
                                                          "whole_phases",
                                                          "whole phases"]:
        binsize_name = binsize.replace(" ", "_")
        csv_results_action_m = np.zeros((len(phases), len(mice),
                                        len(mice)))
        csv_results_action_m_exp = np.zeros((len(phases), len(mice),
                                            len(mice)))
        if save_times:
            csv_results_time = np.zeros((len(phases), len(mice),
                                         len(mice)))
            csv_results_time_exp = np.zeros((len(phases), len(mice),
                                             len(mice)))
    else:
        binsize_name = binsize.replace(" ", "_")
    if return_median:
        method = "median_N_%d" % N
    else:
        method = "mean_N_%d" % N
    fname = '%s_%s_%s' % (var_name, method, add_info_mice)
    fname_csv = '%s_%s_%s.csv' % (action2_name, prefix, add_info_mice)

    fname_beg = '%s_excess' % action2_name
    fname_rev = '%s_excess' % action1_name
    if full_dir_tree:
        fname_exc = '%s_%s_%s%s.csv' % (fname_beg,
                                        method,
                                        prefix,
                                        add_info_mice)
        fname_exc_rev = '%s_%s_%s%s.csv' % (fname_rev,
                                            method,
                                            prefix,
                                            add_info_mice)
        fname_rev_ = '%s_%s%s.csv' % (action1_name, prefix, add_info_mice)
    else:
        fname_exc = 'rasters_%s_%s_%s%s.csv' % (fname_beg,
                                                method,
                                                prefix,
                                                add_info_mice)
        fname_exc_rev = 'rasters_%s_%s_%s%s.csv' % (fname_rev,
                                                    method,
                                                    prefix,
                                                    add_info_mice)
        fname_rev_ = 'rasters_%s_%s_%s.csv' % (action1_name, prefix,
                                               add_info_mice)

    keys = utils.all_mouse_pairs(mice)
    interval_details = {key: [] for key in keys}
    if ecohab_data.how_many_antennas() > 2:
        vmax = 20
        vmin1 = -20
        vmax1 = 20
        vmaxt = 0.005
        vmin1t = -0.005
        vmax1t = 0.005
    else:
        vmax = 150
        vmin1 = -150
        vmax1 = 150
        vmaxt = 0.01
        vmin1t = -0.01
        vmax1t = 0.01
    if full_dir_tree:
        raster_dir = os.path.join(var_name,
                                  "raster_plots",
                                  "bins_%s" % binsize_name)
        raster_dir_add = os.path.join(var_name,
                                      'additionals',
                                      'raster_plots', "bins_%s" % binsize_name)
        hist_dir = os.path.join(var_name,
                                'histograms',
                                "bins_%s" % binsize_name)
        hist_dir_add = os.path.join(var_name,
                                    "additionals",
                                    "histograms", "bins_%s" % binsize_name)
        other_dir = os.path.join('other_variables',
                                 'durations_%s' % var_name,
                                 'histograms',
                                 "bins_%s" % binsize_name)
        other_hist = os.path.join("other_variables",
                                  "histograms_of_%s_intervals" % var_name,
                                  "bins_%s" % binsize_name)
        other_exc_hist = os.path.join('other_variables',
                                      '%s_excess_histograms' % var_name,
                                      "bins_%s" % binsize_name)
        other_raster_dir = os.path.join("other_variables",
                                        "durations_%s" % var_name,
                                        "rasters",
                                        "bin_%s" % binsize)
    else:
        raster_dir = var_name
        raster_dir_add = var_name
        hist_dir = var_name
        hist_dir_add = var_name
        suff = (var_name.split("_")[0][0]+var_name.split("_")[1][0]).upper()
        other_dir = os.path.join('other_variables',
                                 'durations_%s' % suff)
        other_hist = os.path.join("other_variables",
                                  "histograms_of_%s_intervals" % suff)
        other_exc_hist = os.path.join('other_variables',
                                      '%s_excess_histograms' % suff)
        other_raster_dir = os.path.join("other_variables",
                                        "durations_%s" % suff)
    meas_prefix = "measured_%s_%s_%s" % (var_name, prefix, add_info_mice)
    exp_prefix = "expected_%s_%s_%s" % (var_name, prefix, add_info_mice)
    exc_prefix = "excess_%s_%s_%s" % (var_name, prefix, add_info_mice)
    meas_prefix_dur = "durations_%s_%s_%s" % (var_name, prefix, add_info_mice)
    exp_prefix_dur = "exp_durations_%s_%s_%s" % (var_name, prefix,
                                                 add_info_mice)
    exc_prefix_dur = "excess_durations_%s_%s_%s" % (var_name, prefix,
                                                    add_info_mice)

    mouse_action2_sum = OrderedDict()
    mouse_action1_sum = OrderedDict()
    mouse_action2_sum_exc = OrderedDict()
    mouse_action1_sum_exc = OrderedDict()
    mouse_activity = OrderedDict()
    mouse_action2_sum_div_activ = OrderedDict()
    mouse_action1_sum_div_activ = OrderedDict()
    mouse_action2_sum_div_activ_exc = OrderedDict()
    mouse_action1_sum_div_activ_exc = OrderedDict()
    visits = get_activity(ecohab_data, timeline, binsize)
    mouse_activity = utils.sum_activity(visits, all_phases, mice, bin_labels)

    for idx_phase, ph in enumerate(all_phases):
        new_phase = phases[idx_phase]
        for i, lab in enumerate(bin_labels[ph]):
            t_start, t_stop = times[ph][lab]
            directions_dict = data[ph][lab]
            out = function(directions_dict, mice, t_start, t_stop)
            action_m[ph][lab], time_together[ph][lab], phase_intervals1 = out
            duration = t_stop - t_start
            out_expected = resample_single_phase(sur_data_list[ph][lab],
                                                 mice,
                                                 t_start,
                                                 t_stop,
                                                 new_phase,
                                                 var_name,
                                                 function,
                                                 res_dir=res_dir,
                                                 prefix=prefix,
                                                 save_distributions=save_distributions,
                                                 stf=save_times,
                                                 save_figures=save_figures,
                                                 full_dir_tree=full_dir_tree)
            action_m_exp[ph][lab], time_together_exp[ph][lab] = out_expected
            add_intervals(interval_details, phase_intervals1)
        mouse_action2_sum[ph] = utils.sum_per_mouse(action_m[ph], mice,
                                                    bin_labels[ph],
                                                    action2_name)
        mouse_action1_sum[ph] = utils.sum_per_mouse(action_m[ph], mice,
                                                    bin_labels[ph],
                                                    action1_name)
        mouse_action2_sum_div_activ[ph] = utils.divide_sum_activity(mouse_action2_sum[ph],
                                                                    mouse_activity[ph])
        mouse_action1_sum_div_activ[ph] = utils.divide_sum_activity(mouse_action1_sum[ph],
                                                                    mouse_activity[ph])
        if full_dir_tree:
            meas_fname = meas_prefix
            exp_fname = exp_prefix
            exc_fname = exc_prefix
        else:
            meas_fname = "histograms_%s" % meas_prefix
            exp_fname = "histograms_%s" % exp_prefix
            exc_fname = "histograms_%s" % exc_prefix
        write_binned_data(action_m[ph],
                          meas_fname,
                          mice, bin_labels[ph], new_phase, res_dir,
                          hist_dir_add,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter, full_dir_tree=full_dir_tree)
        write_binned_data(action_m_exp[ph],
                          '%s_%s' % (exp_fname, method),
                          mice, bin_labels[ph], new_phase, res_dir,
                          hist_dir_add,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter,
                          full_dir_tree=full_dir_tree)
        exc_action_m = utils.calc_excess(action_m[ph], action_m_exp[ph])

        mouse_action2_sum_exc[ph] = utils.sum_per_mouse(action_m_exp[ph],
                                                        mice,
                                                        bin_labels[ph],
                                                        action2_name)
        mouse_action1_sum_exc[ph] = utils.sum_per_mouse(action_m_exp[ph],
                                                        mice,
                                                        bin_labels[ph],
                                                        action1_name)
        mouse_action2_sum_div_activ_exc[ph] = utils.divide_sum_activity(mouse_action2_sum_exc[ph],
                                                                        mouse_activity[ph])
        mouse_action1_sum_div_activ_exc[ph] = utils.divide_sum_activity(mouse_action1_sum_exc[ph],
                                                                        mouse_activity[ph])

        write_binned_data(exc_action_m,
                          '%s_%s' % (exc_fname, method),
                          mice, bin_labels[ph], new_phase, res_dir,
                          hist_dir,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter, full_dir_tree=full_dir_tree)

        if isinstance(binsize, int) or isinstance(binsize, float):
            if int(binsize) == 24*3600:
                if full_dir_tree:
                    fname = "%s_N_%d_%s" % (var_name, N, method)
                else:
                    fname = "histograms_%s_N_%d_%s" % (var_name, N, method)
                res = utils.dict_to_array_2D(action_m[ph][0],
                                             mice, mice)
                exp_res = utils.dict_to_array_2D(action_m_exp[ph][0],
                                                 mice, mice)
                single_in_cohort_soc_plot(res,
                                          exp_res,
                                          mice,
                                          new_phase,
                                          fname,
                                          res_dir,
                                          hist_dir,
                                          prefix+add_info_mice,
                                          hist=False,
                                          vmin=0,
                                          vmax=vmax,
                                          vmin1=vmin1,
                                          vmax1=vmax1,
                                          titles=['# %s' % name,
                                                  '# expected %s' % name,
                                                  '# excess %s' % name,
                                                  'histogram of # excess %s'
                                                  % name],
                                          labels=['%s mouse' % action2_name,
                                                  '%s mouse' % action1_name],
                                          full_dir_tree=full_dir_tree)
                csv_results_action_m[idx_phase] = res
                csv_results_action_m_exp[idx_phase] = exp_res
        elif isinstance(binsize, str):
            if binsize.lower() in ["whole phase", "whole_phase",
                                   "whole phases", "whole_phases"]:
                fname = "%s_N_%d_%s" % (var_name, N, method)
                res = utils.dict_to_array_2D(action_m[ph][0],
                                             mice, mice)
                exp_res = utils.dict_to_array_2D(action_m_exp[ph][0],
                                                 mice, mice)
                single_in_cohort_soc_plot(res,
                                          exp_res,
                                          mice,
                                          new_phase,
                                          fname,
                                          res_dir,
                                          hist_dir,
                                          prefix+add_info_mice,
                                          hist=False,
                                          vmin=0,
                                          vmax=vmax,
                                          vmin1=vmin1,
                                          vmax1=vmax1,
                                          titles=['# %s' % name,
                                                  '# expected %s' % name,
                                                  '# excess %s' % name,
                                                  'histogram of # excess %s'
                                                  % name],
                                          labels=['%s mouse' % action2_name,
                                                  '%s mouse' % action1_name],
                                          full_dir_tree=full_dir_tree)
                csv_results_action_m[idx_phase] = res
                csv_results_action_m_exp[idx_phase] = exp_res
        if full_dir_tree:
            fname_measured = "%s_%s.csv" % (meas_prefix, new_phase)
            fname_exc = "%s_%s.csv" % (exc_prefix, new_phase)
            fname_expected = "%s_%s.csv" % (exp_prefix, new_phase)
            action2_fname = "excess_%s_%s" % (action2_name, new_phase)
            action1_fname = "excess_%s_%s" % (action1_name, new_phase)
        else:
            fname_measured = "raster_%s_%s.csv" % (meas_prefix, new_phase)
            fname_exc = "raster_%s_%s.csv" % (exc_prefix, new_phase)
            fname_expected = "raster_%s_%s.csv" % (exp_prefix, new_phase)
            action2_fname = "raster_excess_%s_%s" % (action2_name, new_phase)
            action1_fname = "raster_excess_%s_%s" % (action1_name, new_phase)

        raster_labels = [bin_label/3600 for bin_label in bin_labels[ph]]
        phase_full_results = utils.dict_to_array_3D(action_m[ph],
                                                    bin_labels[ph],
                                                    mice, mice)
        phase_exp_full_results = utils.dict_to_array_3D(action_m_exp[ph],
                                                        bin_labels[ph],
                                                        mice, mice)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results,
                          res_dir,
                          raster_dir_add,
                          fname_measured,
                          delimiter=delimiter,
                          symmetrical=False, prefix=prefix,
                          full_dir_tree=full_dir_tree)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_exp_full_results,
                          res_dir,
                          raster_dir_add,
                          fname_expected,
                          delimiter=delimiter,
                          symmetrical=False, prefix=prefix,
                          full_dir_tree=full_dir_tree)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results - phase_exp_full_results,
                          res_dir,
                          raster_dir,
                          action1_fname,
                          delimiter=delimiter,
                          symmetrical=False,
                          reverse=True, prefix=prefix,
                          full_dir_tree=full_dir_tree)

        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results - phase_exp_full_results,
                          res_dir,
                          raster_dir,
                          action2_fname,
                          delimiter=delimiter,
                          symmetrical=False, prefix=prefix,
                          full_dir_tree=full_dir_tree)

        if save_times:
            if full_dir_tree:
                fname_times = 'durations_%s' % var_name
            else:
                fname_times = 'histogram_durations_%s' % var_name
            write_binned_data(time_together[ph],
                              fname_times,
                              mice, bin_labels[ph], new_phase, res_dir,
                              other_dir,
                              prefix, additional_info=add_info_mice,
                              delimiter=delimiter,
                              full_dir_tree=full_dir_tree)
            write_binned_data(time_together_exp[ph],
                              '%s_expected_%s'
                              % (fname_times, method),
                              mice, bin_labels[ph], new_phase, res_dir,
                              other_dir,
                              prefix, additional_info=add_info_mice,
                              delimiter=delimiter,
                              full_dir_tree=full_dir_tree)
            exc_time = utils.calc_excess(time_together[ph],
                                         time_together[ph])
            write_binned_data(exc_time,
                              '%s_expected_%s'
                              % (fname_times, method),
                              mice, bin_labels[ph], new_phase, res_dir,
                              other_dir,
                              prefix, additional_info=add_info_mice,
                              delimiter=delimiter,
                              full_dir_tree=full_dir_tree)
            if isinstance(binsize, int) or isinstance(binsize, float):
                if int(binsize) == 24*3600:
                    fname = "%s_N_%d_%s" % (fname_times, N, method)
                    res = utils.dict_to_array_2D(time_together[ph][0],
                                                 mice, mice)
                    exp_res = utils.dict_to_array_2D(time_together_exp[ph][0],
                                                     mice, mice)
                    single_in_cohort_soc_plot(res,
                                              exp_res,
                                              mice,
                                              new_phase,
                                              fname,
                                              res_dir,
                                              other_dir,
                                              prefix+add_info_mice,
                                              hist=False,
                                              vmin=0,
                                              vmax=vmaxt,
                                              vmin1=vmin1t,
                                              vmax1=vmax1t,
                                              titles=['Fraction of duration %s'
                                                      % name,
                                                      '# expected duration',
                                                      '# excess duration',
                                                      'histogram of # excess duration %s'
                                                      % name],
                                              labels=['%s mouse' % action2_name,
                                                      '%s mouse'
                                                      % action1_name],
                                              full_dir_tree=full_dir_tree)
                    csv_results_time[idx_phase] = res
                    csv_results_time_exp[idx_phase] = exp_res
            elif isinstance(binsize, str):
                if binsize.lower() in ["whole phase", "whole_phase",
                                       "whole phases", "whole_phases"]:
                    fname = "%s_N_%d_%s" % (fname_times, N, method)
                    res = utils.dict_to_array_2D(time_together[ph][0],
                                                 mice, mice)
                    exp_res = utils.dict_to_array_2D(time_together_exp[ph][0],
                                                     mice, mice)
                    single_in_cohort_soc_plot(res,
                                              exp_res,
                                              mice,
                                              new_phase,
                                              fname,
                                              res_dir,
                                              other_dir,
                                              prefix+add_info_mice,
                                              hist=False,
                                              vmin=0,
                                              vmax=vmaxt,
                                              vmin1=vmin1t,
                                              vmax1=vmax1t,
                                              titles=['Fraction of duration %s'
                                                      % name,
                                                      '# expected duration',
                                                      '# excess duration',
                                                      'histogram of # excess duration %s'
                                                      % name],
                                              labels=['%s mouse' %
                                                      action2_name,
                                                      '%s mouse' %
                                                      action1_name],
                                              full_dir_tree=full_dir_tree)
                    csv_results_time[idx_phase] = res
                    csv_results_time_exp[idx_phase] = exp_res
            if full_dir_tree:
                r_fname_measured = "%s_%s.csv" % (meas_prefix_dur, new_phase)
                r_fname_expected = "%s_%s.csv" % (exp_prefix_dur, new_phase)
                r_fname_exc = "%s_%s.csv" % (exc_prefix_dur, new_phase)
            else:
                r_fname_measured = "raster_%s_%s.csv" % (meas_prefix_dur,
                                                         new_phase)
                r_fname_expected = "raster_%s_%s.csv" % (exp_prefix_dur,
                                                         new_phase)
                r_fname_exc = "raster_%s_%s.csv" % (exc_prefix_dur,
                                                    new_phase)

            raster_labels = [bin_label/3600 for bin_label in bin_labels[ph]]
            phase_full_results = utils.dict_to_array_3D(action_m[ph],
                                                        bin_labels[ph],
                                                        mice, mice)
            phase_exp_full_results = utils.dict_to_array_3D(action_m_exp[ph],
                                                            bin_labels[ph],
                                                            mice, mice)
            write_csv_rasters(mice,
                              raster_labels,
                              phase_full_results,
                              res_dir,
                              other_dir,
                              r_fname_measured,
                              delimiter=delimiter,
                              symmetrical=False, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            write_csv_rasters(mice,
                              raster_labels,
                              phase_exp_full_results,
                              res_dir,
                              other_dir,
                              r_fname_expected,
                              delimiter=delimiter,
                              symmetrical=False, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            write_csv_rasters(mice,
                              raster_labels,
                              phase_full_results - phase_exp_full_results,
                              res_dir,
                              other_dir,
                              r_fname_exc, delimiter=delimiter,
                              symmetrical=False, prefix=prefix,
                              full_dir_tree=full_dir_tree)

    write_sum_data(mouse_action1_sum,
                   "mouse_%s_sum_%s" % (action1_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action2_sum,
                   "mouse_%s_sum_%s" % (action1_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action1_sum_exc,
                   "mouse_%s_sum_excess_%s" % (action1_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action1_sum_exc,
                   "mouse_%s_sum_excess_%s" % (action2_name,
                                               binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action1_sum_div_activ,
                   "mouse_%s_activity_%s" % (action1_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action2_sum_div_activ,
                   "mouse_%s_activity_%s" % (action2_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action1_sum_div_activ_exc,
                   "mouse_%s_activity_excess_%s"
                   % (action1_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)
    write_sum_data(mouse_action2_sum_div_activ_exc,
                   "mouse_%s_activity_excess_%s"
                   % (action1_name, binsize_name),
                   mice, bin_labels, all_phases,
                   res_dir,
                   raster_dir_add,
                   prefix,
                   additional_info="ALL",
                   delimiter=";",
                   bool_bins=True,
                   full_dir_tree=full_dir_tree)

    if isinstance(binsize, str):
        if binsize.lower() in ["whole phase", "whole_phase",
                               "whole_phases", "whole phases"]:
            write_csv_rasters(mice,
                              phases,
                              csv_results_action_m -
                              csv_results_action_m_exp,
                              res_dir,
                              raster_dir,
                              fname_exc,
                              symmetrical=False,
                              delimiter=delimiter,
                              reverse=True, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            write_csv_rasters(mice,
                              phases,
                              csv_results_action_m -
                              csv_results_action_m_exp,
                              res_dir,
                              raster_dir,
                              fname_exc_rev,
                              symmetrical=False,
                              delimiter=delimiter, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            write_csv_rasters(mice,
                              phases,
                              csv_results_action_m,
                              res_dir,
                              raster_dir_add,
                              fname_csv,
                              symmetrical=False,
                              reverse=True,
                              delimiter=delimiter, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            write_csv_rasters(mice,
                              phases,
                              csv_results_action_m,
                              res_dir,
                              raster_dir_add,
                              fname_rev_,
                              symmetrical=False,
                              delimiter=delimiter, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            make_RasterPlot(res_dir,
                            raster_dir,
                            (csv_results_action_m -
                             csv_results_action_m_exp),
                            phases,
                            fname_exc,
                            mice,
                            title='excess %s' % action2_name,
                            symmetrical=False, prefix=prefix,
                            full_dir_tree=full_dir_tree)
            pooled_hists(action_m,
                         action_m_exp,
                         all_phases,
                         '%s_histogram' % var_name,
                         res_dir,
                         other_exc_hist,
                         prefix,
                         additional_info=add_info_mice,
                         full_dir_tree=full_dir_tree)

    if save_times:
        make_histograms_for_every_mouse(interval_details,
                                        "%s_intervals_hist" % var_name,
                                        mice,
                                        res_dir,
                                        other_hist,
                                        prefix,
                                        additional_info=add_info_mice,
                                        full_dir_tree=full_dir_tree)
        pooled_hists_for_every_mouse(interval_details,
                                     "%s_intervals_hist" % var_name,
                                     mice,
                                     res_dir,
                                     other_hist,
                                     prefix,
                                     additional_info=add_info_mice,
                                     full_dir_tree=full_dir_tree)
        write_interpair_intervals(interval_details,
                                  other_hist,
                                  res_dir, "%s_intervals" % var_name,
                                  prefix, additional_info=add_info_mice,
                                  delimiter=delimiter,
                                  full_dir_tree=full_dir_tree)
    if isinstance(binsize, str) and binsize.lower() in ["whole phase",
                                                        "whole_phase",
                                                        "whole_phases",
                                                        "whole phases"]:
        if save_times:
            write_csv_rasters(mice,
                              phases,
                              csv_results_time - csv_results_time_exp,
                              res_dir,
                              other_raster_dir,
                              "raster_excess_duration_%s" % action2_name,
                              symmetrical=False,
                              reverse=True,
                              delimiter=delimiter, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            write_csv_rasters(mice,
                              phases,
                              csv_results_time - csv_results_time_exp,
                              res_dir,
                              other_raster_dir,
                              "raster_excess_duration_%s" % action1_name,
                              symmetrical=False,
                              delimiter=delimiter, prefix=prefix,
                              full_dir_tree=full_dir_tree)
            make_RasterPlot(res_dir,
                            other_raster_dir,
                            (csv_results_time - csv_results_time_exp),
                            phases,
                            "raster_excess_durations_%s" % var_name,
                            mice,
                            title='excess duration %s' % name,
                            symmetrical=False, prefix=prefix,
                            full_dir_tree=full_dir_tree)
    return action_m, action_m_exp, phases, mice
