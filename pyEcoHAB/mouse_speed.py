from __future__ import print_function, division, absolute_import
import numpy as np

from . import utility_functions as utils
from .write_to_file import save_single_histograms
from .write_to_file import write_csv_rasters
from .write_to_file import write_csv_tables
from .write_to_file import write_csv_alone
from .write_to_file import write_interpair_intervals
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot
from .plotting_functions import make_pooled_histograms
from .plotting_functions import make_histograms_for_every_mouse
from .plotting_functions import make_pooled_histograms_for_every_mouse

phase_duration = 12*3600

KEY_DICT = {
    "12": 0,
    "21": 0,
    "34": 0,
    "43": 0,
    "56": 0,
    "65": 0,
    "78": 0,
    "87": 0,
}


def calculate_expected(p1, p2):
    expected = 0
    for key in p1:
        expected += p1[key]*p2[key]
    return expected


def calculate_expected_matrices(time_dict_m1,
                                time_dict_m2,
                                tot_time_dict_m1,
                                count_dict_m2,
                                duration,
                                mice_list):
    assert len(mice_list) > 1
    assert duration > 0
    followings = np.zeros((len(mice_list), len(mice_list)))
    time = np.zeros((len(mice_list), len(mice_list)))
    for j, mouse1 in enumerate(mice_list):
        for k, mouse2 in enumerate(mice_list):
            if mouse1 != mouse2:
                assert mouse1 in time_dict_m1
                assert mouse2 in time_dict_m2
                assert mouse1 in tot_time_dict_m1
                assert mouse2 in count_dict_m2
                time_m1 = time_dict_m1[mouse1]
                time_m2 = time_dict_m2[mouse2]
                tot_time_m1 = tot_time_dict_m1[mouse1]
                count_m2 = count_dict_m2[mouse2]

                are_mouse_dicts_correct(time_m2, count_m2)
                time[j, k] = calculate_expected(time_m1, time_m2)
                time[j, k] /=  (duration**2)
                followings[j, k] = calculate_expected(tot_time_m1, count_m2)
                followings[j, k] /= duration
    return followings, time


def are_mouse_dicts_correct(dict1, dict2):
    keys = []
    for key in dict1.keys():
        if key not in dict2:
            raise KeyError
        if dict1[key] == 0:
            if dict2[key] != 0:
                raise ValueError
        keys.append(key)
    assert sorted(keys) == sorted(dict2.keys())


def check_2nd_mouse(antenna1, next_antenna1, t1,
                    threshold, antennas2, times2):
    idxs = utils.get_idx_between(t1, t1 + threshold, times2)
    for ci in idxs:
        if ci + 1 >= len(antennas2):
            return 0
        if antennas2[ci] == antenna1 and antennas2[ci+1] == next_antenna1:
            if times2[ci+1] > t1 + threshold:
                return 1
    return 0


def get_intervals(antenna1, next_antenna1,
                  t1, threshold,
                  antennas2, times2):
    idxs = utils.get_idx_between(t1, t1 + threshold, times2)
    for ci in idxs:
        if ci + 1 >= len(antennas2):
            return 0
        if antennas2[ci] == antenna1 and antennas2[ci+1] == next_antenna1:
            if times2[ci+1] > t1 + threshold:
                #time spent together in tunnel, whole time, time spent in tunnel by mouse 2
                return  times2[ci + 1] -  t1, t1 + threshold - times2[ci], times2[ci + 1] - times2[ci]
    return 0


def following_2_mice_in_pipe(antennas1, times1,
                             antennas2, times2):
    change_indices = utils.change_state(antennas1)
    followings = 0
    intervals = []
    deltas_t1 = KEY_DICT.copy() #time spent in tunnels by mouse1, while being followed
    deltas_t2 = KEY_DICT.copy() #time spent in tunnels by mouse2, while following mouse1
    followings_mouse2 = KEY_DICT.copy() #count of instances of mouse2 following
    tot_time_tunnel = KEY_DICT.copy()
    time_together = 0

    for idx in change_indices:
        antenna1, next_antenna1 = antennas1[idx:idx+2]
        delta_t1 = times1[idx+1] - times1[idx]
        key = utils.get_key_for_frequencies(antenna1, next_antenna1)
        if key:
            tot_time_tunnel[key] += delta_t1
        if utils.in_tube(antenna1, next_antenna1) and delta_t1 < phase_duration:
            out = check_2nd_mouse(antenna1, next_antenna1,
                                  times1[idx], delta_t1,
                                  antennas2, times2)
            if out:
                followings += out
                int_combined, int_overlap, delta_t2 = get_intervals(antenna1,
                                                                      next_antenna1,
                                                                      times1[idx],
                                                                      delta_t1,
                                                                      antennas2,
                                                                      times2)
                intervals += [int_combined]
                time_together += int_overlap
                key = "%d%d" % (antenna1, next_antenna1)
                deltas_t1[key] += delta_t1
                deltas_t2[key] += delta_t2
                followings_mouse2[key] += 1
    assert followings == sum(followings_mouse2.values())
    return followings, time_together, intervals, deltas_t1, deltas_t2, tot_time_tunnel, followings_mouse2


def add_values_to_dict(dictionary, values):
    for key in values.keys():
        dictionary[key] += values[key]


def initialize_dict(mice):
    return {mouse:KEY_DICT.copy() for mouse in mice}


def following_matrices(times_antennas, mice, durations):
    followings = np.zeros((len(mice), len(mice)))
    time_together = np.zeros((len(mice), len(mice)))
    labels = utils.all_pairs(mice)
    interval_details = {label:[] for label in labels}
    intervals_mouse1 = initialize_dict(mice)
    intervals_mouse2 = initialize_dict(mice)
    time_tunnel_mouse1 = initialize_dict(mice)
    followings_mouse2 = initialize_dict(mice)

    for i, mouse1 in enumerate(mice):
        for j, mouse2 in enumerate(mice):
            if mouse1 == mouse2:
                continue
            times1, antennas1 = times_antennas[mouse1]
            times2, antennas2 = times_antennas[mouse2]
            out = following_2_mice_in_pipe(antennas1,
                                           times1,
                                           antennas2,
                                           times2)
            followings[i, j], time_in_pipe, mouse_intervals, int_m1, int_m2, time_m1, f_m2 = out
            time_together[i, j] = time_in_pipe/durations
            key =  "%s|%s" % (mouse1, mouse2)
            interval_details[key] += mouse_intervals
            add_values_to_dict(intervals_mouse1[mouse1], int_m1)
            add_values_to_dict(intervals_mouse2[mouse2], int_m2)
            add_values_to_dict(time_tunnel_mouse1[mouse1], time_m1)
            add_values_to_dict(followings_mouse2[mouse2], f_m2)
    return followings, time_together, interval_details, intervals_mouse1, intervals_mouse2, time_tunnel_mouse1, followings_mouse2


def prepare_data(ehd, st, en):
    times_antennas = {}
    for j, mouse1 in enumerate(ehd.mice):
        times_antennas[mouse1] = utils.get_times_antennas(ehd,
                                                          mouse1,
                                                          st,
                                                          en)
    return times_antennas


def following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase):
    st, en = cf.gettime(phase) 
    duration = en - st
    assert duration > 0
    times_readings = prepare_data(ehd, st, en)
    return following_matrices(times_readings, ehd.mice, duration)


def add_intervals(all_intervals, phase_intervals):
    for mouse in phase_intervals.keys():
        all_intervals[mouse].extend(phase_intervals[mouse])


def get_following(ehd, cf, res_dir=None, prefix=None,
              remove_mouse=None):
    if res_dir is None:
        res_dir = ehd.res_dir
    if prefix is None:
        prefix = ehd.prefix
    phases = utils.filter_dark(cf.sections())
    mice = [mouse[-4:] for mouse in ehd.mice]
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    following = np.zeros((len(phases), len(mice), len(mice)))
    following_exp = np.zeros((len(phases), len(mice), len(mice)))
    time_together = np.zeros((len(phases), len(mice), len(mice)))
    time_together_exp = np.zeros((len(phases), len(mice),
                                       len(mice)))
    fname = 'following_in_pipe_%s' % (add_info_mice)
    fname_ = 'following_in_pipe_%s%s' % (prefix,
                                         add_info_mice)
    fname_beg = 'relative_following_in_pipe_excess'
    fname_exp = '%s_%s%s.csv' % (fname_beg, prefix,
                                 add_info_mice)

    keys = utils.all_pairs(ehd.mice)
    interval_details = {key:[] for key in keys}

    for i, phase in enumerate(phases):
        out = following_2nd_mouse_in_pipe_single_phase(ehd,
                                                       cf,
                                                       phase)
        following[i], time_together[i], phase_intervals, ints_m1, ints_m2, time_m1, foll_m2  = out
        start, end = cf.gettime(phase)
        duration = end - start
        assert duration > 0
        out_expected = calculate_expected_matrices(ints_m1,
                                                   ints_m2,
                                                   time_m1,
                                                   foll_m2,
                                                   duration,
                                                   ehd.mice)
        following_exp[i], time_together_exp[i] = out_expected
        add_intervals(interval_details, phase_intervals)
        save_single_histograms(following[i],
                               'following_in_pipe',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(following_exp[i],
                               'following_in_pipe_expected_time',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_in_pipe_relative_excess_following',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        vmin1 = (following[i] - following_exp[i]).min()
        vmax1 = (following[i] - following_exp[i]).max()
        single_in_cohort_soc_plot(following[i],
                                  following_exp[i],
                                  mice,
                                  phase,
                                  fname,
                                  res_dir,
                                  'following_in_pipe/histograms',
                                  prefix+add_info_mice,
                                  hist=False,
                                  vmin=0,
                                  vmax=max(following[i].max(), following_exp[i].max()),
                                  vmin1=vmin1,
                                  vmax1=vmax1,
                                  titles=['# followings',
                                          '# expected followings',
                                          '# excess followings',
                                          'histogram of # excess followings',],
                                  labels=['following mouse', 'followed mouse'])


        save_single_histograms(time_together[i],
                               'time_together_in_pipe',
                               ehd.mice,
                               phase,
                               res_dir,
                               'time_together_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(time_together_exp[i],
                               'time_together_in_pipe_expected_time',
                               ehd.mice,
                               phase,
                               res_dir,
                               'time_together_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((time_together[i]-time_together_exp[i]),
                               'time_together_in_pipe_relative_excess_time_together',
                               ehd.mice,
                               phase,
                               res_dir,
                               'time_together_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        vmin1 = (time_together[i] - time_together_exp[i]).min()
        vmax1 = (time_together[i] - time_together_exp[i]).max()
        single_in_cohort_soc_plot(time_together[i],
                                  time_together_exp[i],
                                  mice,
                                  phase,
                                  "fraction_time_following",
                                  res_dir,
                                  'time_together_in_pipe/histograms',
                                  prefix+add_info_mice,
                                  hist=False,
                                  vmin=0,
                                  vmax=max(time_together[i].max(), time_together_exp[i].max()),
                                  vmin1=vmin1,
                                  vmax1=vmax1,
                                  titles=['Fraction of time following',
                                          '# expected time',
                                          '# excess time',
                                          'histogram of # excess time following',],
                                  labels=['following mouse', 'followed mouse'])

    write_csv_rasters(ehd.mice,
                      phases,
                      following,
                      res_dir,
                      'following_in_pipe/raster_plots',
                      fname_)
    write_csv_rasters(ehd.mice,
                      phases,
                      (following-following_exp),
                      res_dir,
                      'following_in_pipe/raster_plots',
                      fname_exp)
    

    make_RasterPlot(res_dir,
                    'following_in_pipe/raster_plots',
                    following,
                    phases,
                    fname_,
                    mice,
                    title='# followings')
    make_RasterPlot(res_dir,
                    'following_in_pipe/raster_plots',
                    (following-following_exp),
                    phases,
                    fname_exp,
                    mice,
                    title='% excess following')

    make_pooled_histograms(following,
                           following_exp,
                           phases,
                           'Following_histogram',
                           res_dir,
                           'following_in_pipe/raster_plots',
                           prefix,
                           additional_info=add_info_mice)

    make_histograms_for_every_mouse(interval_details,
                                    "followings_intervals_histogram",
                                    ehd.mice,
                                    res_dir,
                                    "following_in_pipe/histograms_of_following_intervals",
                                    prefix,
                                    additional_info=add_info_mice)
    make_pooled_histograms_for_every_mouse(interval_details,
                                           "followings_intervals_histogram",
                                           ehd.mice,
                                           res_dir,
                                           "following_in_pipe/histograms_of_following_intervals",
                                           prefix,
                                           additional_info=add_info_mice)
    write_interpair_intervals(interval_details, "following_in_pipe/histograms_of_following_intervals",
                              res_dir, "following_intervals", prefix,
                              additional_info=add_info_mice)
    return following, following_exp, phases, mice
