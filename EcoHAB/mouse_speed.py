from __future__ import print_function, division, absolute_import
import numpy as np

from . import utility_functions as utils
from .write_to_file import save_single_histograms
from .write_to_file import write_csv_rasters
from .write_to_file import write_csv_tables
from .write_to_file import write_csv_alone
from .write_to_file import write_interpair_intervals
from .plotfunctions import single_in_cohort_soc_plot, make_RasterPlot
from .plotfunctions import make_pooled_histograms
from .plotfunctions import make_histograms_for_every_mouse
from .plotfunctions import make_pooled_histograms_for_every_mouse

titles = {
    3: '12',
    7: '34',
    11: '56',
    15: '78',
}
threshold = 12*3600


def frequency_mouse_in_tube(times, antennas, period):
    change_indices = utils.change_state(antennas)
    frequency = {
        3: 0,
        7: 0,
        11: 0,
        15: 0,
    }
    window = {
        3: 0,
        7: 0,
        11: 0,
        15: 0,
    }
    for idx in change_indices:
        antenna, next_antenna = antennas[idx:idx + 2]
        delta_t = times[idx+1] - times[idx]
        key = utils.get_key_for_frequencies(antenna, next_antenna)
        if key and delta_t < threshold:
            frequency[key] += 1
            window[key] += delta_t
    for key in frequency:
        frequency[key] = frequency[key]/period
    return frequency, window


def frequencies_for_all(ehd, cf, phase):
    t_st, t_en = cf.gettime(phase)
    period = t_en - t_st
    frequency = {}
    window = {}
    for mouse1 in ehd.mice:
        times, antennas = utils.get_times_antennas(ehd, mouse1, t_st, t_en)
        frequency[mouse1], window[mouse1] = frequency_mouse_in_tube(times,
                                                                    antennas,
                                                                    period)
    return frequency, window


def calculate_expected_followings(window_mouse1, frequency_mouse2):
    expected_followings = 0
    for key in window_mouse1:
        expected_followings += window_mouse1[key]*frequency_mouse2[key]
    return expected_followings


def expected_following_in_pipe_single_phase(ehd, cf, phase):
    out = np.zeros((len(ehd.mice), len(ehd.mice)))
    frequency, window = frequencies_for_all(ehd, cf, phase)
    for j, mouse1 in enumerate(ehd.mice):
        for k, mouse2 in enumerate(ehd.mice):
            if mouse1 != mouse2:
                out[j, k] = calculate_expected_followings(window[mouse1],
                                                          frequency[mouse2])
    return out


def check_2nd_mouse(antenna1, next_antenna1, t1, threshold, antennas2, times2):
    idxs = utils.get_idx_between(t1, t1 + threshold, times2)
    for ci in idxs:
        if ci + 1 >= len(antennas2):
            return 0
        if antennas2[ci] == antenna1 and antennas2[ci+1] == next_antenna1:
            if times2[ci+1] > t1 + threshold:
                return 1
    return 0

def following_interval(antenna1, next_antenna1, t1, threshold, antennas2, times2):
    idxs = utils.get_idx_between(t1, t1 + threshold, times2)
    for ci in idxs:
        if ci + 1 >= len(antennas2):
            return 0
        if antennas2[ci] == antenna1 and antennas2[ci+1] == next_antenna1:
            if times2[ci+1] > t1 + threshold:
                return  times2[ci+1] -  times2[ci]
    return 0

def following_2_mice_in_pipe(antennas1, times1,
                             antennas2, times2):
    change_indices = utils.change_state(antennas1)
    followings = 0
    intervals = []
    for idx in change_indices:
        antenna1, next_antenna1 = antennas1[idx:idx+2]
        delta_t1 = times1[idx+1] - times1[idx]
        if utils.in_tube(antenna1, next_antenna1) and delta_t1 < threshold:
            out = check_2nd_mouse(antenna1, next_antenna1,
                                  times1[idx], delta_t1,
                                  antennas2, times2)
            if out:
                followings += out
                intervals += [following_interval(antenna1,
                                                 next_antenna1,
                                                 times1[idx],
                                                 delta_t1,
                                                 antennas2,
                                                 times2)]
    return followings, intervals


def following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase):
    st, en = cf.gettime(phase) 
    followings = np.zeros((len(ehd.mice), len(ehd.mice)))
    interval_details = {}
    for mouse1 in ehd.mice:
        for mouse2 in ehd.mice:
            if mouse1 != mouse2:
                key = "%s_%s" % (mouse1, mouse2)
                interval_details[key] = []
    for j, mouse1 in enumerate(ehd.mice):
        times1, antennas1 = utils.get_times_antennas(ehd, mouse1,
                                                     st, en)
        for k, mouse2 in enumerate(ehd.mice):
            if mouse2 != mouse1:
                times2, antennas2 = utils.get_times_antennas(ehd, mouse2,
                                                             st, en)
                mouse_followings, mouse_intervals = following_2_mice_in_pipe(antennas1,
                                                                             times1,
                                                                             antennas2,
                                                                             times2)
                followings[j, k] = mouse_followings
                key =  "%s_%s" % (mouse1, mouse2)
                interval_details[key] += mouse_intervals
    return followings, interval_details

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
    fname = 'following_in_pipe_%s' % (add_info_mice)
    fname_ = 'following_in_pipe_%s%s' % (prefix,
                                                           add_info_mice)
    fname_beg = 'relative_following_in_pipe_excess'
    fname_exp = '%s_%s%s.csv' % (fname_beg, prefix, add_info_mice)
    interval_details = {}
    for mouse1 in ehd.mice:
        for mouse2 in ehd.mice:
            if mouse1 != mouse2:
                key = "%s_%s" % (mouse1, mouse2)
                interval_details[key] = []

    for i, phase in enumerate(phases):
        following[i], phase_intervals = following_2nd_mouse_in_pipe_single_phase(ehd,
                                                                                 cf,
                                                                                 phase)
        following_exp[i] = expected_following_in_pipe_single_phase(ehd,
                                                                   cf,
                                                                   phase)
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
                                    "following_in_pipe",
                                    prefix,
                                    additional_info=add_info_mice)
    make_pooled_histograms_for_every_mouse(interval_details,
                                           "followings_intervals_histogram",
                                           ehd.mice,
                                           res_dir,
                                           "following_in_pipe",
                                           prefix,
                                           additional_info=add_info_mice)
    write_interpair_intervals(interval_details, "following_in_pipe",
                              res_dir, "following_intervals", prefix,
                              additional_info=add_info_mice)
    return following, following_exp, phases, mice
