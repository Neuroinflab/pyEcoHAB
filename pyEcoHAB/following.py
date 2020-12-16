# -*- coding: utf-8 -*-
from __future__ import division, absolute_import
import random
import os
import numpy as np

from . import utility_functions as utils
from .write_to_file import save_single_histograms
from .write_to_file import write_csv_rasters
from .write_to_file import write_binned_data
from .write_to_file import write_interpair_intervals
from .write_to_file import write_bootstrap_results
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot
from .plotting_functions import make_pooled_histograms
from .plotting_functions import make_histograms_for_every_mouse
from .plotting_functions import make_pooled_histograms_for_every_mouse
from .plotting_functions import single_histogram_figures


def insert_interval(candidate_t_start, interval,
                    t_starts, t_ends, duration):

    if candidate_t_start in t_starts:
        return 0

    candidate_t_end = candidate_t_start + interval

    if candidate_t_end in t_ends:
        return 0

    if candidate_t_end > duration:
        return 0

    beg = utils.get_idx_pre(candidate_t_end, t_starts)
    if beg is None:
        t_starts.insert(0, candidate_t_start)
        t_ends.insert(0, candidate_t_end)
        return 1

    end = utils.get_idx_pre(candidate_t_start, t_ends)
    if end is None:
        t_starts.append(candidate_t_start)
        t_ends.append(candidate_t_end)
        return 1

    if beg != end:
        return 0

    t_starts.insert(beg + 1, candidate_t_start)
    t_ends.insert(beg + 1, candidate_t_end)
    return 1


def generate_intervals(t_starts, t_stops, duration):
    intervals = sorted(utils.get_interval_durations_2_lists(t_starts,
                                                            t_stops),
                       reverse=True)
    new_t_starts, new_t_stops = [], []
    ints_len = len(intervals)
    i = 0
    iterations = 0
    while i < ints_len:
        interval = intervals[i]
        can_t_start = random.randrange(0, duration)
        out = insert_interval(can_t_start, interval,
                              new_t_starts, new_t_stops,
                              duration)
        iterations += 1
        i += out
        if iterations > 2*ints_len:
            # start over
            i = 0
            iterations = 0
            new_t_starts, new_t_stops = [], []
    return new_t_starts, new_t_stops


def gen_directions_dict(directions_dict, duration, keys):
    new_dict = {}
    for key in keys:
        old_intervals = directions_dict[key]
        new_dict[key] = generate_intervals(old_intervals[0],
                                           old_intervals[1],
                                           duration)
    return new_dict


def bootstrap_single_phase(directions_dict, mice_list,
                           t_start, t_stop, keys, N=1000):
    followings = utils.make_results_dict(mice_list, tolist=True)
    times_together = utils.make_results_dict(mice_list, tolist=True)
    new_directions = {}
    for i in range(N):
        for mouse in mice_list:
            new_directions[mouse] = gen_directions_dict(directions_dict[mouse],
                                                        t_stop - t_start, keys)
        out = following_matrices(new_directions, mice_list,
                                 t_start, t_stop, keys)
        for m1 in mice_list:
            for m2 in mice_list:
                if m1 != m2:
                    followings[m1][m2].append(out[0][m1][m2])
                    times_together[m1][m2].append(out[1][m1][m2])
    return followings, times_together


def resample_single_phase(directions_dict, mice, t_start, t_stop, N, phase,
                          keys, return_median=False, save_figures=False,
                          save_distributions=True, res_dir=None, prefix=None,
                          stf=False):
    """If return_median is False, function returns mean value
    of the resampled following distribution

    stf: save times following"""

    if res_dir is None:
        res_dir = ecohab_data.res_dir
    if prefix is None:
        prefix = ecohab_data.prefix

    followings, times_following = bootstrap_single_phase(directions_dict,
                                                         mice,
                                                         t_start, t_stop, keys,
                                                         N=N)
    binsize = (t_stop - t_start)/3600
    hist_dir = os.path.join("other_variables",
                            "dynamic_interactions_hists",
                            "bin_%4.2f" % binsize)
    hist_time_dir = os.path.join("other_variables",
                                 "durations_dynamic_interactions_hists",
                                 "bin_%4.2f" % binsize)
    fname_following = "%s_DI_count_dist_%d_%4.2f" % (prefix, N, binsize)
    fname_times = "%s_DI_durations_dis_%d_%4.2f" % (prefix, N, binsize)
    if save_figures:
        for mouse1 in mice:
            for mouse2 in mice:
                if mouse1 == mouse2:
                    continue
                key = "%s_%s" % (mouse1, mouse2)
                fname1 = "%s%s_hist_%s_%s_N_%d_%4.2f" % (prefix, "DI",
                                                         phase.replace(' ',
                                                                       '_'),
                                                         key, N, binsize)
                fname2 = "%s%s_hist_%s_%s_N_%d_%4.2f" % (prefix,
                                                         "durations_DI",
                                                         phase.replace(' ',
                                                                       '_'),
                                                         key, N, binsize)
                single_histogram_figures(followings[mouse1][mouse2],
                                         fname1, res_dir,
                                         hist_dir,
                                         "Dynamic interactions count",
                                         xlabel="dynamic interactions",
                                         ylabel="count",
                                         median_mean=True)
                if stf:
                    single_histogram_figures(times_following[mouse1][mouse2],
                                             fname2,
                                             res_dir,
                                             hist_time_dir,
                                             "Dynamic interaction durations",
                                             xlabel="duration",
                                             ylabel="count", nbins=10,
                                             median_mean=True)
    dist_dir_fol = os.path.join("other_variables",
                                "dynamic_interactions_hists",
                                "bin_%4.2f" % binsize)
    dist_dir_time = os.path.join("other_variables",
                                 "durations_dynamic_interactions_hists",
                                 "bin_%4.2f" % binsize)
    if save_distributions:
        write_bootstrap_results(followings, phase, mice,
                                fname_following, res_dir,
                                dist_dir_fol, prefix)
        if stf:
            write_bootstrap_results(times_following, phase, mice,
                                    fname_times, res_dir,
                                    dist_dir_time,
                                    prefix)
    out_followings = utils.make_results_dict(mice)
    out_times = utils.make_results_dict(mice)
    for m1 in mice:
        for m2 in mice:
            if m1 == m2:
                continue
            if return_median:
                out_followings[m1][m2] = np.median(followings[m1][m2])
                out_times[m1][m2] = np.median(times_following[m1][m2])
            else:
                out_followings[m1][m2] = np.mean(followings[m1][m2])
                out_times[m1][m2] = np.mean(times_following[m1][m2])

    return out_followings, out_times


def following_single_pair(directions_m1, directions_m2, keys):
    followings = 0
    intervals = []
    time_together = 0

    for key in keys:
        out = following_single_direction(directions_m1[key],
                                         directions_m2[key])
        f_single_dir, time_single_dir, ints_single_dir = out
        followings += f_single_dir
        time_together += time_single_dir
        intervals += ints_single_dir
    return followings, time_together, intervals


def following_matrices(directions_dict, mice, t_start, t_stop, keys):
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
                                        directions_dict[mouse2], keys)
            followings[mouse1][mouse2], time_in_pipe, mouse_intervals = out
            time_together[mouse1][mouse2] = time_in_pipe/durations
            key = "%s|%s" % (mouse1, mouse2)
            interval_details[key] += mouse_intervals
    return followings, time_together, interval_details


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


def get_dynamic_interactions(ecohab_data, timeline, N, binsize=12*3600,
                             res_dir="", prefix="", remove_mouse=None,
                             save_distributions=True, save_figures=False,
                             return_median=False, delimiter=";",
                             save_times_following=False, seed=None):
    if res_dir == "":
        res_dir = ecohab_data.res_dir
    if prefix == "":
        prefix = ecohab_data.prefix
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    mice = utils.get_mice(ecohab_data.mice, remove_mouse)
    phases, times, data, data_keys = utils.prepare_binned_registrations(ecohab_data,
                                                                        timeline,
                                                                        binsize,
                                                                        mice)
    if isinstance(seed, int):
        random.seed(seed)
    all_phases, bin_labels = data_keys
    following = utils.make_all_results_dict(*data_keys)
    following_exp = utils.make_all_results_dict(*data_keys)
    time_together = utils.make_all_results_dict(*data_keys)
    time_together_exp = utils.make_all_results_dict(*data_keys)

    if isinstance(binsize, int) or isinstance(binsize, float):
        binsize_name = "%3.2f_h" % (binsize/3600)
        if int(binsize) == 43200 or int(binsize) == 24*3600:
            csv_results_following = np.zeros((len(phases), len(mice),
                                              len(mice)))
            csv_results_following_exp = np.zeros((len(phases), len(mice),
                                                  len(mice)))
            if save_times_following:
                csv_results_time = np.zeros((len(phases), len(mice),
                                             len(mice)))
                csv_results_time_exp = np.zeros((len(phases), len(mice),
                                                 len(mice)))
    else:
        binsize_name = binsize
    if return_median:
        method = "median_N_%d" % N
    else:
        method = "mean_N_%d" % N
    fname = 'dynamics_interactions_%s_%s' % (method, add_info_mice)
    fname_ = 'following_%s%s.csv' % (prefix, add_info_mice)
    fname_rev_ = 'leading_%s%s.csv' % (prefix, add_info_mice)
    fname_beg = 'following_excess'
    fname_rev = 'leading_excess'
    fname_excess = '%s_%s_%s%s.csv' % (fname_beg,
                                       method,
                                       prefix,
                                       add_info_mice)
    fname_excess_rev = '%s_%s_%s%s.csv' % (fname_rev,
                                           method,
                                           prefix,
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
    raster_dir = os.path.join("dynamic_interactions", "raster_plots",
                              "bins_%s" % binsize_name)
    raster_dir_add = os.path.join('dynamic_interactions', 'additionals',
                                  'raster_plots', "bins_%s" % binsize_name)
    hist_dir = os.path.join('dynamic_interactions', 'histograms',
                            "bins_%s" % binsize_name)
    hist_dir_add = os.path.join("dynamic_interactions", "additionals",
                                "histograms", "bins_%s" % binsize_name)
    other_dir = os.path.join('other_variables',
                             'durations_dynamic_interactions', 'histograms',
                             "bins_%s" % binsize_name)
    other_hist = os.path.join("other_variables",
                              "histograms_of_dynamic_interactions_intervals",
                              "bins_%s" % binsize_name)
    other_excess_hist = os.path.join('other_variables',
                                     'dynamic_interactions_excess_histograms',
                                     "bins_%s" % binsize_name)
    meas_prefix = "measured_dynamic_interactions_%s_%s" % (prefix,
                                                           add_info_mice)
    exp_prefix = "expected_dynamic_interactions_%s_%s" % (prefix,
                                                          add_info_mice)
    excess_prefix = "excess_dynamic_interactions_%s_%s" % (prefix,
                                                           add_info_mice)

    meas_prefix_dur = "durations_dynamic_interactions_%s_%s" % (prefix,
                                                                add_info_mice)
    exp_prefix_dur = "exp_durations_dynamic_interactions_%s_%s" % (prefix,
                                                                   add_info_mice)
    excess_prefix_dur = "excess_durations_dynamic_interactions_%s_%s" % (prefix,
                                                                         add_info_mice)

    other_raster_dir = os.path.join("other_variables",
                                    "durations_dynamic_interactions",
                                    "rasters",
                                    "bin_%s" % binsize)
    for idx_phase, ph in enumerate(all_phases):
        new_phase = phases[idx_phase]
        for i, lab in enumerate(bin_labels):
            t_start, t_stop = times[ph][lab]
            directions_dict = data[ph][lab]
            out = following_matrices(directions_dict, mice, t_start, t_stop,
                                     ecohab_data.directions)
            following[ph][lab], time_together[ph][lab], phase_intervals1 = out
            duration = t_stop - t_start
            out_expected = resample_single_phase(directions_dict,
                                                 mice,
                                                 t_start,
                                                 t_stop,
                                                 N,
                                                 new_phase,
                                                 ecohab_data.directions,
                                                 res_dir=res_dir,
                                                 prefix=prefix,
                                                 stf=save_times_following,
                                                 save_figures=save_figures)
            following_exp[ph][lab], time_together_exp[ph][lab] = out_expected
            add_intervals(interval_details, phase_intervals1)

        write_binned_data(following[ph],
                          'dynamic_interactions',
                          mice, bin_labels, new_phase, res_dir,
                          hist_dir_add,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter)
        write_binned_data(following_exp[ph],
                          'dynamic_interactions_expected_%s' % method,
                          mice, bin_labels, new_phase, res_dir,
                          hist_dir_add,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter)
        excess_following = utils.calc_excess(following[ph],
                                             following_exp[ph])
        write_binned_data(excess_following,
                          'dynamic_interactions_excess_%s' % method,
                          mice, bin_labels, new_phase, res_dir,
                          hist_dir,
                          prefix, additional_info=add_info_mice,
                          delimiter=delimiter)

        if isinstance(binsize, int) or isinstance(binsize, float):
            if int(binsize) == 12*3600 or int(binsize) == 24*3600:
                fname = "dynamic_interactions_N_%d_%s" % (N, method)
                res = utils.dict_to_array_2D(following[ph][0],
                                             mice, mice)
                exp_res = utils.dict_to_array_2D(following_exp[ph][0],
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
                                          titles=['# dynamic interactions',
                                                  '# expected dynamic interactions',
                                                  '# excess dynamic interactions',
                                                  'histogram of # excess dynamic interactions', ],
                                          labels=['following mouse',
                                                  'followed mouse'])
                csv_results_following[idx_phase] = res
                csv_results_following_exp[idx_phase] = exp_res
        fname_measured = "%s_%s.csv" % (meas_prefix, new_phase)
        fname_excess = "%s_%s.csv" % (excess_prefix, new_phase)
        fname_expected = "%s_%s.csv" % (exp_prefix, new_phase)
        raster_labels = [bin_label/3600 for bin_label in bin_labels]
        phase_full_results = utils.dict_to_array_3D(following[ph],
                                                    bin_labels,
                                                    mice, mice)
        phase_exp_full_results = utils.dict_to_array_3D(following_exp[ph],
                                                        bin_labels,
                                                        mice, mice)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results,
                          res_dir,
                          raster_dir_add,
                          fname_measured,
                          delimiter=delimiter,
                          symmetric=False)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_exp_full_results,
                          res_dir,
                          raster_dir_add,
                          fname_expected,
                          delimiter=delimiter,
                          symmetric=False)
        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results - phase_exp_full_results,
                          res_dir,
                          raster_dir,
                          "excess_following_%s.csv" % new_phase,
                          delimiter=delimiter,
                          symmetric=False,
                          reverse=True)

        write_csv_rasters(mice,
                          raster_labels,
                          phase_full_results - phase_exp_full_results,
                          res_dir,
                          raster_dir,
                          "excess_leading_%s.csv" % new_phase,
                          delimiter=delimiter,
                          symmetric=False)

        if save_times_following:
            write_binned_data(time_together[ph],
                              'durations_dynamic_interactions',
                              mice, bin_labels, new_phase, res_dir,
                              other_dir,
                              prefix, additional_info=add_info_mice,
                              delimiter=delimiter)
            write_binned_data(time_together_exp[ph],
                              'durations_dynamic_interactions_expected_%s' % method,
                              mice, bin_labels, new_phase, res_dir,
                              other_dir,
                              prefix, additional_info=add_info_mice,
                              delimiter=delimiter)
            excess_time = utils.calc_excess(time_together[ph],
                                            time_together[ph])
            write_binned_data(excess_time,
                              'durations_dynamic_interactions_expected_%s' % method,
                              mice, bin_labels, new_phase, res_dir,
                              other_dir,
                              prefix, additional_info=add_info_mice,
                              delimiter=delimiter)
            if isinstance(binsize, int) or isinstance(binsize, float):
                if int(binsize) == 12*3600 or int(binsize) == 24*3600:
                    fname = "durations_dynamic_interactions_N_%d_%s" % (N,
                                                                        method)
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
                                              titles=['Fraction of duration dynamics interation',
                                                      '# expected duration',
                                                      '# excess duration',
                                                      'histogram of # excess duration dynamic interactions', ],
                                              labels=['following mouse',
                                                      'followed mouse'])
                    csv_results_time[idx_phase] = res
                    csv_results_time_exp[idx_phase] = exp_res
            fname_measured = "%s_%s.csv" % (meas_prefix_dur, new_phase)
            fname_expected = "%s_%s.csv" % (exp_prefix_dur, new_phase)
            raster_labels = [bin_label/3600 for bin_label in bin_labels]
            phase_full_results = utils.dict_to_array_3D(following[ph],
                                                        bin_labels,
                                                        mice, mice)
            phase_exp_full_results = utils.dict_to_array_3D(following_exp[ph],
                                                            bin_labels,
                                                            mice, mice)
            write_csv_rasters(mice,
                              raster_labels,
                              phase_full_results,
                              res_dir,
                              other_dir,
                              fname_measured,
                              delimiter=delimiter,
                              symmetric=False)
            write_csv_rasters(mice,
                              raster_labels,
                              phase_exp_full_results,
                              res_dir,
                              other_dir,
                              fname_expected,
                              delimiter=delimiter,
                              symmetric=False)
            write_csv_rasters(mice,
                              raster_labels,
                              phase_full_results - phase_exp_full_results,
                              res_dir,
                              other_dir,
                              fname_excess, delimiter=delimiter,
                              symmetric=False)

    if isinstance(binsize, int) or isinstance(binsize, float):
        if binsize == 43200:
            write_csv_rasters(mice,
                              phases,
                              csv_results_following
                              - csv_results_following_exp,
                              res_dir,
                              raster_dir,
                              fname_excess,
                              symmetric=False,
                              delimiter=delimiter,
                              reverse=True)
            write_csv_rasters(mice,
                              phases,
                              csv_results_following
                              - csv_results_following_exp,
                              res_dir,
                              raster_dir,
                              fname_excess_rev,
                              symmetric=False,
                              delimiter=delimiter)
            write_csv_rasters(mice,
                              phases,
                              csv_results_following,
                              res_dir,
                              raster_dir_add,
                              fname_,
                              symmetric=False,
                              reverse=True,
                              delimiter=delimiter)
            write_csv_rasters(mice,
                              phases,
                              csv_results_following,
                              res_dir,
                              raster_dir_add,
                              fname_rev_,
                              symmetric=False,
                              delimiter=delimiter)
            make_RasterPlot(res_dir,
                            raster_dir,
                            (csv_results_following - csv_results_following_exp),
                            phases,
                            fname_excess,
                            mice,
                            title='% excess following',
                            symmetric=False)
            make_pooled_histograms(following,
                                   following_exp,
                                   all_phases,
                                   'Dynamic_interactions_histogram',
                                   res_dir,
                                   other_excess_hist,
                                   prefix,
                                   additional_info=add_info_mice)

    if save_times_following:
        make_histograms_for_every_mouse(interval_details,
                                        "dynamic_interactions_intervals_hist",
                                        mice,
                                        res_dir,
                                        other_hist,
                                        prefix,
                                        additional_info=add_info_mice)
        make_pooled_histograms_for_every_mouse(interval_details,
                                               "dynamic_interactions_intervals_hist",
                                               mice,
                                               res_dir,
                                               other_hist,
                                               prefix,
                                               additional_info=add_info_mice)
        write_interpair_intervals(interval_details,
                                  other_hist,
                                  res_dir, "dynamic_interactions_intervals",
                                  prefix, additional_info=add_info_mice,
                                  delimiter=delimiter)
        if binsize == 43200:
            write_csv_rasters(mice,
                              phases,
                              csv_results_time - csv_results_time_exp,
                              res_dir,
                              other_raster_dir,
                              "excess_duration_following",
                              symmetric=False,
                              reverse=True,
                              delimiter=delimiter)
            write_csv_rasters(mice,
                              phases,
                              csv_results_time - csv_results_time_exp,
                              res_dir,
                              other_raster_dir,
                              "excess_duration_leading",
                              symmetric=False,
                              delimiter=delimiter)
            make_RasterPlot(res_dir,
                            other_raster_dir,
                            (csv_results_time - csv_results_time_exp),
                            phases,
                            "excess_durations_dynamic_interactions",
                            mice,
                            title='% excess duration dynamic interactions',
                            symmetric=False)
    return following, following_exp, phases, mice
