from __future__ import print_function, division, absolute_import
import random
import numpy as np

from . import utility_functions as utils
from .write_to_file import save_single_histograms
from .write_to_file import write_csv_rasters
from .write_to_file import write_csv_tables
from .write_to_file import write_csv_alone
from .write_to_file import write_interpair_intervals
from .write_to_file import write_bootstrap_results
from .plotting_functions import single_in_cohort_soc_plot, make_RasterPlot
from .plotting_functions import make_pooled_histograms
from .plotting_functions import make_histograms_for_every_mouse
from .plotting_functions import make_pooled_histograms_for_every_mouse
from .plotting_functions import single_histogram_figures

phase_duration = 12*3600
keys = ['12', '21', '34', '43', '56', '65', '78', '87']
KEY_DICT = {
    '12': 0,
    '21': 0,
    '34': 0,
    '43': 0,
    '56': 0,
    '65': 0,
    '78': 0,
    '87': 0,
}

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
            #start over
            i = 0
            iterations = 0
            new_t_starts, new_t_stops = [], []
    return new_t_starts, new_t_stops


def generate_directions_dict(directions_dict, duration):
    new_dict = {}
    for key in keys:
        old_intervals = directions_dict[key]
        new_dict[key] = generate_intervals(old_intervals[0],
                                           old_intervals[1],
                                           duration)
    return new_dict

def bootstrap_single_phase(directions_dict, mice_list,
                           t_start, t_stop, N=1000):
    followings = np.zeros((len(mice_list), len(mice_list), N), dtype=int)
    times_together = np.zeros((len(mice_list), len(mice_list), N))
    new_directions = {}
    for i in range(N):
        for mouse in mice_list:
            new_directions[mouse] = generate_directions_dict(directions_dict[mouse],
                                                             t_stop - t_start)
        out = following_matrices(new_directions, mice_list,
                                 t_start, t_stop)
        followings[:, :, i] = out[0]
        times_together[:, :, i] = out[1]
    return followings, times_together

def resample_single_phase(ehd, cf, phase, N,
                          return_median=False,
                          save_figures=False,
                          save_distributions=True,
                          res_dir=None, prefix=None,
                          stf=False):
    """If return_median is False, function returns mean value
    of the resampled following distribution

    stf: save times following"""
    t_start, t_stop = cf.gettime(phase)
    mice = ehd.mice
    assert  t_stop - t_start > 0
    if res_dir is None:
        res_dir = ehd.res_dir
    if prefix is None:
        prefix = ehd.prefix
    directions_dict = prepare_data(ehd, t_start, t_stop)
    followings, times_following = bootstrap_single_phase(directions_dict,
                                                         mice,
                                                         t_start, t_stop,
                                                         N=N)
    if save_figures:
        fname_following = "following_count_distribution_%d" % N
        fname_times = "following_times_distribution_%d" % N
        for i, mouse1 in enumerate(mice):
            for j, mouse2 in enumerate(mice):
                if mouse1 != mouse2:
                    key = "%s|%s" % (mouse1, mouse2)
                    fname1 = "%s_histogram_%s_%s_N_%d" % ("following",
                                                          phase.replace(' ',
                                                                        '_'),
                                                          key, N)
                    fname2 = "%s_histogram_%s_%s_N_%d" % ("time_together",
                                                          phase.replace(' ',
                                                                        '_'),
                                                          key, N)
                    single_histogram_figures(followings[i, j],
                                             fname1, res_dir,
                                             "other_variables/following_hists",
                                             "Following count distribution",
                                             xlabel="followings",
                                             ylabel="count",
                                             median_mean=True)
                    if stf:
                        single_histogram_figures(times_following[i, j], fname2,
                                                 "Following times distribution",
                                                 "other_variables/time_following_hists",
                                                 add_text_2,
                                                 xlabel="time_together",
                                                 ylabel="count", nbins=10,
                                                 median_mean=True)
    if save_distributions:
        fname_following = "following_count_distribution_%d" % N
        fname_times = "following_times_distribution_%d" % N
        write_bootstrap_results(followings, phase, mice,
                                fname_following, res_dir,
                                "other_variables/following_hists", prefix)
        if stf:
            write_bootstrap_results(times_following, phase, mice,
                                    fname_times, res_dir,
                                    "other_variables/time_following_hists",
                                    prefix)

    if return_median:
        return np.median(followings, axis=2), np.median(times_following, axis=2)
    return followings.mean(axis=2), times_following.mean(axis=2)


def following_single_pair(directions_m1, directions_m2):
    
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


def following_matrices(directions_dict, mice, t_start, t_stop):
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
            out = following_single_pair(directions_dict[mouse1],
                                        directions_dict[mouse2])
            followings[i, j], time_in_pipe, mouse_intervals = out
            time_together[i, j] = time_in_pipe/durations
            key =  "%s|%s" % (mouse1, mouse2)
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


def extract_directions(times, antennas, last_antenna):
    direction_dict = {key:[[], []] for key in keys}
    change_indices = utils.change_state(antennas)
    for c_idx in change_indices:
        if c_idx + 1 >= len(antennas):
            break
        ant, next_ant = antennas[c_idx], antennas[c_idx + 1]
        key = utils.get_key_for_frequencies(ant, next_ant)
        if key is not None:
            try:
                third_antenna = antennas[c_idx + 2]
            except IndexError:
               third_antenna = last_antenna
            if third_antenna == ant:
                continue
            direction_dict[key][0].append(times[c_idx])
            direction_dict[key][1].append(times[c_idx + 1])
    return direction_dict


def prepare_data(ehd, st, en):
    directions = {}
    for j, mouse1 in enumerate(ehd.mice):
        times_antennas = utils.get_times_antennas(ehd,
                                                  mouse1,
                                                  st,
                                                  en)
        last_times, last_antennas = utils.get_times_antennas(ehd,
                                                             mouse1,
                                                             en,
                                                             en+(en-st))
        try:
            last_antenna = last_antennas[0]
        except IndexError:
            last_antenna = None
        directions[mouse1] = extract_directions(times_antennas[0],
                                                times_antennas[1],
                                                last_antenna)
    return directions

def get_matrices_single_phase(ehd, cf, phase, function):
    t_start, t_stop = cf.gettime(phase)
    assert  t_stop - t_start > 0
    directions_dict = prepare_data(ehd, t_start, t_stop)
    return function(directions_dict, ehd.mice,
                    t_start, t_stop)


def add_intervals(all_intervals, phase_intervals):
    for mouse in phase_intervals.keys():
        all_intervals[mouse].extend(phase_intervals[mouse])


def get_following(ehd, cf, N, res_dir="", prefix="",
                  remove_mouse=None, save_distributions=True,
                  save_figures=False, return_median=False, delimiter=";",
                 save_times_following=False):
    if res_dir == "":
        res_dir = ehd.res_dir
    if prefix == "":
        prefix = ehd.prefix
    phases = utils.filter_dark_light(cf.sections())
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    following = np.zeros((len(phases), len(ehd.mice), len(ehd.mice)))
    following_exp = np.zeros((len(phases), len(ehd.mice), len(ehd.mice)))
    time_together = np.zeros((len(phases), len(ehd.mice), len(ehd.mice)))
    time_together_exp = np.zeros((len(phases), len(ehd.mice),
                                  len(ehd.mice)))
    if return_median:
        method = "median_N_%d" % N
    else:
        method = "mean_N_%d" % N
    fname = 'following_%s_%s' % (method, add_info_mice)
    fname_ = 'following_%s%s.csv' % (prefix, add_info_mice)
    fname_rev_ = 'leading_%s%s.csv' % (prefix, add_info_mice)
    fname_beg = 'following_excess'
    fname_rev = 'leading_excess'
    fname_exp = '%s_%s_%s%s.csv' % (fname_beg,
                                    method,
                                    prefix,
                                    add_info_mice)
    fname_exp_rev = '%s_%s_%s%s.csv' % (fname_rev,
                                    method,
                                    prefix,
                                    add_info_mice)
    keys = utils.all_pairs(ehd.mice)
    interval_details = {key:[] for key in keys}
    if ehd.how_many_antennas() > 2:
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

    for i, phase in enumerate(phases):
        out = get_matrices_single_phase(ehd,
                                        cf,
                                        phase,
                                        following_matrices)
        following[i], time_together[i], phase_intervals  = out
        start, end = cf.gettime(phase)
        duration = end - start
        assert duration > 0
        out_expected = resample_single_phase(ehd,
                                             cf,
                                             phase,
                                             N,
                                             res_dir=res_dir,
                                             prefix=prefix,
                                             stf=save_times_following,
                                             save_figures=save_figures)
        following_exp[i], time_together_exp[i] = out_expected
        add_intervals(interval_details, phase_intervals)
        save_single_histograms(following[i],
                               'following',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following/histograms',
                               prefix,
                               additional_info=add_info_mice,
                               delimiter=delimiter)
        save_single_histograms(following_exp[i],
                               'following_expected_time_%s' % method,
                               ehd.mice,
                               phase,
                               res_dir,
                               'following/histograms',
                               prefix,
                               additional_info=add_info_mice,
                               delimiter=delimiter)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_excess_%s' %method,
                               ehd.mice,
                               phase,
                               res_dir,
                               'following/histograms',
                               prefix,
                               additional_info=add_info_mice,
                               delimiter=delimiter)
        single_in_cohort_soc_plot(following[i],
                                  following_exp[i],
                                  ehd.mice,
                                  phase,
                                  fname,
                                  res_dir,
                                  'following/histograms',
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
        if save_times_following:
            save_single_histograms(time_together[i],
                                   'time_following',
                                   ehd.mice,
                                   phase,
                                   res_dir,
                                   'other_variables/time_following/histograms',
                                   prefix,
                                   additional_info=add_info_mice,
                                   delimiter=delimiter)
            save_single_histograms(time_together_exp[i],
                                   'time_following_expected_%s' % method,
                                   ehd.mice,
                                   phase,
                                   res_dir,
                                   'other_variables/time_following/histograms',
                                   prefix,
                                   additional_info=add_info_mice,
                                   delimiter=delimiter)
            save_single_histograms((time_together[i]-time_together_exp[i]),
                                   'time_following_excess_%s' % method,
                                   ehd.mice,
                                   phase,
                                   res_dir,
                                   'other_variables/time_following/histograms',
                                   prefix,
                                   additional_info=add_info_mice,
                                   delimiter=delimiter)

            single_in_cohort_soc_plot(time_together[i],
                                      time_together_exp[i],
                                      ehd.mice,
                                      phase,
                                      "time_following_%s" % method,
                                      res_dir,
                                      'other_variables/time_following/histograms',
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
                                      labels=['following mouse',
                                              'followed mouse'])

    write_csv_rasters(ehd.mice,
                      phases,
                      following,
                      res_dir,
                      'following/raster_plots',
                      fname_,
                      symmetric=False,
                      delimiter=delimiter)
    write_csv_rasters(ehd.mice,
                      phases,
                      following,
                      res_dir,
                      'following/raster_plots',
                      fname_rev_,
                      symmetric=False,
                      reverse_order=True)
    write_csv_rasters(ehd.mice,
                      phases,
                      (following-following_exp),
                      res_dir,
                      'following/raster_plots',
                      fname_exp,
                      symmetric=False)
    write_csv_rasters(ehd.mice,
                      phases,
                      (following-following_exp),
                      res_dir,
                      'following/raster_plots',
                      fname_exp_rev,
                      symmetric=False,
                      reverse_order=True,
                      delimiter=delimiter)


    make_RasterPlot(res_dir,
                    'following/raster_plots',
                    (following-following_exp),
                    phases,
                    fname_exp,
                    ehd.mice,
                    title='% excess following',
                    symmetric=False)

    make_pooled_histograms(following,
                           following_exp,
                           phases,
                           'Following_histogram',
                           res_dir,
                           'other_variables/following_excess_histograms',
                           prefix,
                           additional_info=add_info_mice)

    if save_times_following:
        make_histograms_for_every_mouse(interval_details,
                                        "followings_intervals_histogram",
                                        ehd.mice,
                                        res_dir,
                                        "other_variables/histograms_of_following_intervals",
                                        prefix,
                                        additional_info=add_info_mice)
        make_pooled_histograms_for_every_mouse(interval_details,
                                               "followings_intervals_histogram",
                                               ehd.mice,
                                               res_dir,
                                               "other_variables/histograms_of_following_intervals",
                                               prefix,
                                               additional_info=add_info_mice)
        write_interpair_intervals(interval_details,
                                  "other_variables/histograms_of_following_intervals",
                                  res_dir, "following_intervals", prefix,
                                  additional_info=add_info_mice,
                                  delimiter=delimiter)
    return following, following_exp, phases, ehd.mice
