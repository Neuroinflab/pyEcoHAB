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
    3: 0,
    7: 0,
    11: 0,
    15: 0,
}


def calculate_time_in_tunnel(times, antennas):
    tot_time_tunnel = KEY_DICT.copy()
    change_indices = utils.change_state(antennas)
    for idx in change_indices:
        antenna1, next_antenna1 = antennas[idx:idx+2]
        key = utils.get_key_for_frequencies(antenna1, next_antenna1)
        if key:
            tot_time_tunnel[key] += times[idx+1] - times[idx]
    return tot_time_tunnel


def calculate_expected_per_tunnel(times, antennas, key, bins, t_stop):
    out_follows = np.zeros(len(bins))
    out_times = np.zeros(len(bins))

    if len(bins) > 1:
        binsize = bins[1] - bins[0]
    else:
        binsize = 0

    for i, t_beg in enumerate(bins):
        if i == len(bins) - 1:
            t_end = t_stop
        else:
            t_end = t_beg + binsize
        idxs = utils.get_idx_between(t_beg, t_end, times)
        for idx in idxs:
            try:
                antenna2, next_antenna2 = antennas[idx], antennas[idx+1]
            except IndexError:
                break
            delta_t2 = times[idx+1] - times[idx]
            if utils.in_tube(antenna2,
                             next_antenna2) and delta_t2 < phase_duration:
                if utils.get_key_for_frequencies(antenna2,
                                                 next_antenna2) == key:
                    out_follows[i] += 1
                    out_times[i] += delta_t2
    return out_follows.tolist(), out_times.tolist()


def expected_for_one_pair(times_antennas, tot_time_tunnels,
                          t_start, t_stop):
    times2, antennas2 = times_antennas
    mean_follows = np.zeros((len(tot_time_tunnels.keys(),)))
    std_follows = np.zeros((len(tot_time_tunnels.keys(),)))
    mean_times = np.zeros((len(tot_time_tunnels.keys(),)))
    std_times = np.zeros((len(tot_time_tunnels.keys(),)))
    for i, key in enumerate(tot_time_tunnels.keys()):
        binsize = tot_time_tunnels[key]
        try:
            t_begs = utils.get_times(binsize, t_start, t_stop)
        except ZeroDivisionError:
            continue
        follows, times = calculate_expected_per_tunnel(times2,
                                                       antennas2,
                                                       key,
                                                       t_begs,
                                                       t_stop)
        mean_follows[i] = np.mean(follows)
        mean_times[i] =  np.mean(times)
        std_follows[i] = np.std(follows)
        std_times[i] =  np.mean(times)
    return mean_follows, std_follows, mean_times, std_times

def expected_matrices(times_antennas, mice_list,
                      t_start, t_stop):
    assert len(mice_list) > 1
    assert t_stop - t_start > 0
    exp_followings = np.zeros((len(mice_list), len(mice_list)))
    exp_time = np.zeros((len(mice_list), len(mice_list)))
    for j, mouse1 in enumerate(mice_list):
        tot_time_tunnels = calculate_time_in_tunnel(*times_antennas[mouse1])
        for k, mouse2 in enumerate(mice_list):
            if mouse1 != mouse2:
                times2, antennas2 = times_antennas[mouse2]
                out = expected_for_one_pair(times_antennas[mouse2],
                                            tot_time_tunnels,
                                            t_start, t_stop)
                mean_follows, std_follows, mean_times, std_times = out
                exp_followings[j, k] = (mean_follows/std_follows).sum()
                exp_time[j, k] = (mean_times/std_times).sum()
                exp_time[j, k] /=  (t_stop - t_start)
    return exp_followings, exp_time


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
                #time spent together in tunnel, whole time
                return  times2[ci + 1] -  t1, t1 + threshold - times2[ci]
    return 0


def following_2_mice_in_pipe(antennas1, times1,
                             antennas2, times2):
    change_indices = utils.change_state(antennas1)
    followings = 0
    intervals = []
    time_together = 0

    for idx in change_indices:
        antenna1, next_antenna1 = antennas1[idx:idx+2]
        delta_t1 = times1[idx+1] - times1[idx]
        if utils.in_tube(antenna1, next_antenna1) and delta_t1 < phase_duration:
            out = check_2nd_mouse(antenna1, next_antenna1,
                                  times1[idx], delta_t1,
                                  antennas2, times2)
            if out:
                followings += out
                int_combined, int_overlap = get_intervals(antenna1,
                                                          next_antenna1,
                                                          times1[idx],
                                                          delta_t1,
                                                          antennas2,
                                                          times2)
                intervals += [int_combined]
                time_together += int_overlap
    return followings, time_together, intervals


def add_values_to_dict(dictionary, values):
    for key in values.keys():
        dictionary[key] += values[key]


def initialize_dict(mice):
    return {mouse:KEY_DICT.copy() for mouse in mice}


def following_matrices(times_antennas, mice, t_start, t_stop):
    assert t_stop - t_start > 0
    durations = t_stop - t_start
    followings = np.zeros((len(mice), len(mice)))
    time_together = np.zeros((len(mice), len(mice)))
    labels = utils.all_pairs(mice)
    interval_details = {label:[] for label in labels}

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
            followings[i, j], time_in_pipe, mouse_intervals = out
            time_together[i, j] = time_in_pipe/durations
            key =  "%s|%s" % (mouse1, mouse2)
            interval_details[key] += mouse_intervals
    return followings, time_together, interval_details


def prepare_data(ehd, st, en):
    times_antennas = {}
    for j, mouse1 in enumerate(ehd.mice):
        times_antennas[mouse1] = utils.get_times_antennas(ehd,
                                                          mouse1,
                                                          st,
                                                          en)
    return times_antennas


def get_matrices_single_phase(ehd, cf, phase, function):
    t_start, t_stop = cf.gettime(phase)
    assert  t_stop - t_start > 0
    times_readings = prepare_data(ehd, t_start, t_stop)
    return function(times_readings, ehd.mice,
                    t_start, t_stop)


def add_intervals(all_intervals, phase_intervals):
    for mouse in phase_intervals.keys():
        all_intervals[mouse].extend(phase_intervals[mouse])


def get_following(ehd, cf, res_dir=None, prefix=None,
              remove_mouse=None):
    if res_dir is None:
        res_dir = ehd.res_dir
    if prefix is None:
        prefix = ehd.prefix
    phases = utils.filter_dark_light(cf.sections())
    mice = [mouse[-4:] for mouse in ehd.mice]
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    following = np.zeros((len(phases), len(mice), len(mice)))
    following_exp = np.zeros((len(phases), len(mice), len(mice)))
    time_together = np.zeros((len(phases), len(mice), len(mice)))
    time_together_exp = np.zeros((len(phases), len(mice),
                                       len(mice)))
    fname = 'following_weighted_sum_in_pipe_%s' % (add_info_mice)
    fname_ = 'following_weighted_sum_in_pipe_%s%s' % (prefix,
                                         add_info_mice)
    fname_beg = 'relative_following_weighted_sum_in_pipe_excess'
    fname_exp = '%s_%s%s.csv' % (fname_beg, prefix,
                                 add_info_mice)

    keys = utils.all_pairs(ehd.mice)
    interval_details = {key:[] for key in keys}
    if ehd.how_many_antennas() > 2:
        vmax = 40
        vmin1 = -40
        vmax1 = 40
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

    for i, phase in enumerate(phases):
        out = get_matrices_single_phase(ehd,
                                        cf,
                                        phase,
                                        following_matrices)
        following[i], time_together[i], phase_intervals  = out
        start, end = cf.gettime(phase)
        duration = end - start
        assert duration > 0
        out_expected = get_matrices_single_phase(ehd,
                                                 cf,
                                                 phase,
                                                 expected_matrices)
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
                               'following_in_pipe_expected_time_weighted_sum',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_in_pipe_relative_excess_following_weighted_sum',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
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
                                  vmax=vmax,
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
                               'time_together_in_pipe_expected_time_weighted_sum',
                               ehd.mice,
                               phase,
                               res_dir,
                               'time_together_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((time_together[i]-time_together_exp[i]),
                               'time_together_in_pipe_relative_excess_time_together_weighted_sum',
                               ehd.mice,
                               phase,
                               res_dir,
                               'time_together_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)

        single_in_cohort_soc_plot(time_together[i],
                                  time_together_exp[i],
                                  mice,
                                  phase,
                                  "fraction_time_following_weighted_sum",
                                  res_dir,
                                  'time_together_in_pipe/histograms',
                                  prefix+add_info_mice,
                                  hist=False,
                                  vmin=0,
                                  vmax=vmaxt,
                                  vmin1=vmin1t,
                                  vmax1=vmax1t,
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
