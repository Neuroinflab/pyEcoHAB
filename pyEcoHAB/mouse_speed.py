from __future__ import print_function, division, absolute_import
import timeit
import random
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
    intervals = utils.shuffle_intervals(t_starts, t_stops)
    new_t_starts, new_t_stops = [], []
    ints_len = len(intervals)
    i = 0
    while i < ints_len:
        interval = intervals[i]
        can_t_start = random.randrange(0, duration)
        out = insert_interval(can_t_start, interval,
                              new_t_starts, new_t_stops,
                              duration)
        i += out
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
    tstart = timeit.default_timer()
    followings = np.zeros((len(mice_list), len(mice_list), N))
    times_together = np.zeros((len(mice_list), len(mice_list), N))
    new_directions = {}
    assert t_stop - t_start > 0
    for i in range(N):
        for mouse in mice_list:
            new_directions[mouse] = generate_directions_dict(directions_dict[mouse],
                                                             t_stop - t_start)
        out = following_matrices(new_directions, mice_list,
                                 t_start, t_stop)
        followings[:, :, i] = out[0]
        times_together[:, :, i] = out[1]
    tstop = timeit.default_timer()
    print("Took", tstop - tstart)
    return followings, times_together


def compare_single_phase(ehd, cf, phase, N=100):
    t_start, t_stop = cf.gettime(phase)
    mice = ehd.mice
    assert  t_stop - t_start > 0
    res_dir = ehd.res_dir
    prefix = ehd.prefix
    fname_following = "following_count_distribution_%d" % N
    fname_times = "following_times_distribution_%d" % N
    out = get_matrices_single_phase(ehd,
                                    cf,
                                    phase,
                                    following_matrices)
    following, time_together, intervals = out
    directions_dict = prepare_data(ehd, t_start, t_stop)
    out = bootstrap_single_phase(directions_dict,
                                 mice,
                                 t_start, t_stop,
                                 N=N)
    for i, mouse1 in enumerate(mice):
        for j, mouse2 in enumerate(mice):
            if mouse1 != mouse2:
                key = "%s|%s" % (mouse1, mouse2)
                print(mouse1, mouse2,
                      "following", following[i, j],
                      "expected mean", out[0][i, j].mean(),
                      "expected median", np.median(out[0][i, j]))
                print(mouse1, mouse2,
                      "time together", time_together[i, j],
                      "expected mean", out[1][i, j].mean(),
                      "expected median", np.median(out[1][i, j]))
                fname1 = "%s_histogram_%s_%s_N_%d" % ("following",
                                                      phase.replace(' ', '_'),
                                                      key, N)
                fname2 = "%s_histogram_%s_%s_N_%d" % ("time_together",
                                                      phase.replace(' ', '_'),
                                                      key, N)
                add_text_1 = "measured following count %d" %  following[i, j]
                add_text_2 = "measured time_together %f" %  time_together[i, j]
                single_histogram_figures(out[0][i, j], fname1, res_dir,
                                         "following_hists", add_text_1,
                                         xlabel="followings", ylabel="count",
                                         median_mean=True)
                single_histogram_figures(out[1][i, j], fname2, res_dir,
                                         "time_following_hists", add_text_2,
                                         xlabel="time_together", ylabel="count", nbins=10,
                                         median_mean=True)

    write_bootstrap_results(out[0], phase, mice,
                            fname_following, res_dir,
                            "following_hists", prefix)
    write_bootstrap_results(out[1], phase, mice,
                            fname_times, res_dir,
                            "time_following_hists", prefix)

def bootstrap_2_mice(directions_dict1, directions_dict2, duration=12*3600.,
                     N=1000):
    pass


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


def extract_directions(times, antennas):
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
                if third_antenna == ant:
                    continue
            except IndexError:
                pass
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
        directions[mouse1] = extract_directions(*times_antennas)
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
    fname = 'following_resampling_in_pipe_%s' % (add_info_mice)
    fname_ = 'following_resampling_in_pipe_%s%s' % (prefix,
                                         add_info_mice)
    fname_beg = 'relative_following_resampling_in_pipe_excess'
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
        start = timeit.default_timer()
        out = get_matrices_single_phase(ehd,
                                        cf,
                                        phase,
                                        following_matrices)
        stop = timeit.default_timer()
        print("Elapsed time", stop - start)
        following[i], time_together[i], phase_intervals  = out
        start, end = cf.gettime(phase)
        duration = end - start
        assert duration > 0
        # out_expected = get_matrices_single_phase(ehd,
        #                                          cf,
        #                                          phase,
        #                                          expected_matrices)
        # following_exp[i], time_together_exp[i] = out_expected
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
                               'following_in_pipe_expected_time_resampling',
                               ehd.mice,
                               phase,
                               res_dir,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_in_pipe_relative_excess_following_resampling',
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
                               'time_together_in_pipe_expected_time_resampling',
                               ehd.mice,
                               phase,
                               res_dir,
                               'time_together_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((time_together[i]-time_together_exp[i]),
                               'time_together_in_pipe_relative_excess_time_together_resampling',
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
                                  "fraction_time_following_resampling",
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
