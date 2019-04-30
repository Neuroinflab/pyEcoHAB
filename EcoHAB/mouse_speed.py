from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
import utility_functions as utils
import numpy as np
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot, make_pooled_histograms
from numba import jit
titles = {
    3: '12',
    7: '34',
    11: '56',
    15: '78',
}
threshold = 2


@jit
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


@jit
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

@jit
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


def following_2_mice_in_pipe(antennas1, times1,
                             antennas2, times2):
    change_indices = utils.change_state(antennas1)
    followings = 0
    for idx in change_indices:
        antenna1, next_antenna1 = antennas1[idx:idx+2]
        delta_t1 = times1[idx+1] - times1[idx]
        if utils.in_tube(antenna1, next_antenna1) and delta_t1 < threshold:
            followings += check_2nd_mouse(antenna1, next_antenna1,
                                          times1[idx], delta_t1,
                                          antennas2, times2)
    return followings


def following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase):
    st, en = cf.gettime(phase) 
    followings = np.zeros((len(ehd.mice), len(ehd.mice)))
    for j, mouse1 in enumerate(ehd.mice):
        times1, antennas1 = utils.get_times_antennas(ehd, mouse1,
                                                     st, en)
        for k, mouse2 in enumerate(ehd.mice):
            if mouse2 != mouse1:
                times2, antennas2 = utils.get_times_antennas(ehd, mouse2,
                                                             st, en)
                followings[j, k] = following_2_mice_in_pipe(antennas1,
                                                            times1,
                                                            antennas2,
                                                            times2)
    return followings


def following_for_all_2nd_mouse_in_pipe(ehd, cf, main_directory,
                                        prefix, remove_mouse=None):
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
    for i, phase in enumerate(phases):
        following[i] = following_2nd_mouse_in_pipe_single_phase(ehd,
                                                                cf,
                                                                phase)
        following_exp[i] = expected_following_in_pipe_single_phase(ehd,
                                                                   cf,
                                                                   phase)
        save_single_histograms(following[i],
                               'following_in_pipe',
                               ehd.mice,
                               phase,
                               main_directory,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(following_exp[i],
                               'following_in_pipe_expected_time',
                               ehd.mice,
                               phase,
                               main_directory,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_in_pipe_relative_excess_following',
                               ehd.mice,
                               phase,
                               main_directory,
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
                                  main_directory,
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
                      main_directory,
                      'following_in_pipe/raster_plots',
                      fname_)
    write_csv_rasters(ehd.mice,
                      phases,
                      (following-following_exp),
                      main_directory,
                      'following_in_pipe/raster_plots',
                      fname_exp)
    

    make_RasterPlot(main_directory,
                    'following_in_pipe/raster_plots',
                    following,
                    phases,
                    fname_,
                    mice,
                    title='# followings')
    make_RasterPlot(main_directory,
                    'following_in_pipe/raster_plots',
                    (following-following_exp),
                    phases,
                    fname_exp,
                    mice,
                    title='% excess following')

    make_pooled_histograms(following,
                           following_exp,
                           phases,
                           'All_phases_histogram',
                           main_directory,
                           '',
                           prefix,
                           additional_info=add_info_mice)

    return following, following_exp, phases, mice


if __name__ == '__main__':
    nbins = 10
    import matplotlib.pyplot as plt
    from data_info import *
    import os
    homepath = os.path.expanduser("~/")
    followings_list = []
    followings_exp_list = []
    phases_list = []
    prefixes = []
    max_len = 0
    datasets = ['EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/',
                'EcoHAB_data_November/C57 13-24.04 long/',]
    for new_path in datasets:
       
        path = os.path.join(homepath, new_path)
        prefix = utils.make_prefix(path)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
        else:
            remove_mouse = None
        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 500
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        prefix = utils.make_prefix(path)
        res_dir = utils.results_path(path)

        cf = ExperimentConfigFile(path)
        following, following_exp, phases, mice = following_for_all_2nd_mouse_in_pipe(ehd, cf, res_dir, prefix)
        followings_list.append(following)
        followings_exp_list.append(following_exp)
        phases_list.append(phases)
        prefixes.append(prefix)
        if len(phases) > max_len:
            max_len = len(phases)
    if len(datasets) > 1:
        shape = (len(datasets), max_len)
        fig, axes = plt.subplots(*shape, figsize=(shape[1]*10,shape[0]*10))
        max_bins = 0
        min_bins = 200
        min_count = 200
        max_count = 0
        for i in range(shape[0]):
            plt.text(0.2, (0.8+i)/shape[0], datasets[i], fontsize=34, transform=fig.transFigure)
            for j in range(shape[1]):
                if j < followings_list[i].shape[0]:
                    results = followings_list[i][j]
                    results_exp = followings_exp_list[i][j]
                    what = results[results > 0] - results_exp[results > 0]
                    n, bins, patches = axes[i][j].hist(what, bins=nbins)

                    axes[i][j].set_title(phases_list[i][j], fontsize=34)
                    if bins.min() < min_bins:
                        min_bins = bins.min()
                    if bins.max() > max_bins:
                        max_bins = bins.max()
                    if max(n) > max_count:
                        max_count = max(n)
                    if min(n) < min_count:
                        min_count = min(n)
        print(min_bins, max_bins, min_count, max_count)
        for i in range(shape[0]):
            for j in range(shape[1]):
                axes[i][j].set_xlim([min_bins, max_bins])
                axes[i][j].set_ylim([min_count, max_count+4])
                axes[i][j].plot([0, 0],
                                [min_count, max_count+4],
                                color='r',
                                linewidth=2)
                if j:
                    axes[i][j].set_yticklabels([])
                else:
                    for tick in axes[i][j].yaxis.get_major_ticks():
                        tick.label.set_fontsize(34) 
                if i < shape[0]-1:
                    axes[i][j].set_xticklabels([])
                else:
                    for tick in axes[i][j].xaxis.get_major_ticks():
                        tick.label.set_fontsize(34) 

        fig.subplots_adjust(wspace=0.2)
        fig.subplots_adjust(hspace=0.2)
        fig.savefig('Long_data_hist.png', dpi=300)
        #plt.show()
