from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utility_functions as utils
import numpy as np
import matplotlib.pyplot as plt
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot
from numba import jit
from collections import OrderedDict
nbins = 10
homepath = os.path.expanduser("~/")
titles = {
    3: '12',
    7: '34',
    11: '56',
    15: '78',
}
threshold = 2

#@jit
def check_2nd_mouse(antenna1, next_antenna1, t1, threshold, antennas2, times2, out):
    idxs1 = np.where(np.array(times2) >= t1)[0]
    idxs2 = np.where(np.array(times2) <= t1 + threshold)[0]
    common_idxs = list(set(idxs1)&set(idxs2))
    for ci in common_idxs:
        antenna2 = antennas2[ci]
        if ci + 1 < len(antennas2) and antenna2 == antenna1:
            next_antenna2 = antennas2[ci+1]
            delta_t2 = times2[ci+1] - times2[ci]
            if utils.in_tube(antenna2, next_antenna2):
                delta_t = times2[ci] - t1
                out.append(delta_t)
                return 1
    return 0
                           

#@jit            
def following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en, print_out=False):
    intervals = []
    ehd.mask_data(st, en)
    antennas1 = ehd.getantennas(mouse1)
    times1 = ehd.gettimes(mouse1)
    antennas2 = ehd.getantennas(mouse2)
    times2 = ehd.gettimes(mouse2)
    change_indices = np.where((np.array(antennas1[1:]) - np.array(antennas1[:-1])) != 0)[0]
    followings = 0
    for idx in change_indices:
        antenna1 = antennas1[idx]
        next_antenna1 = antennas1[idx+1]
        delta_t1 = times1[idx+1] - times1[idx]
        if utils.in_tube(antenna1, next_antenna1) and delta_t1 < threshold:
            followings += check_2nd_mouse(antenna1, next_antenna1, times1[idx], delta_t1,
                                          antennas2, times2, out=intervals)
    if print_out:
        print('%s following %s %d times'%(mouse2, mouse1, followings))
    return followings, intervals
    

def following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase, print_out=False):
    mice = ehd.mice
    st, en = cf.gettime(phase) 
    followings = np.zeros((len(mice), len(mice)))
    all_mice_dict = {}
    if print_out:
        print(phase)
    for j, mouse1 in enumerate(mice):
        for k, mouse2 in enumerate(mice):
            if mouse2 != mouse1:
                followings[j, k], intervals = following_2_mice_in_pipe_condition(ehd,
                                                                                 mouse1, mouse2, st, en,
                                                                                 print_out=print_out)
                key = '%s_%s' % (mouse1, mouse2)
                all_mice_dict[key] = intervals
    return followings, all_mice_dict
    
def following_for_all_2nd_mouse_in_pipe(ehd, cf, main_directory, prefix, remove_mouse=None, print_out=False):

    phases = cf.sections()
    phases = utils.filter_dark(phases)
    mice = ehd.mice
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    following = np.zeros((len(phases), len(mice), len(mice)))
    following_exp = np.zeros((len(phases), len(mice), len(mice)))
    fname = 'following_in_pipe_2nd_mouse_in_pipe_%s' % (add_info_mice)
    fname_ = 'following_in_pipe_2nd_mouse_in_pipe_%s%s' % (prefix, add_info_mice)
    fname_exp = 'relative_following_in_pipe_excess_2nd_mouse_in_pipe_%s%s.csv' % (prefix, add_info_mice)
    interval_dict = OrderedDict()
    for i, phase in enumerate(phases):
        following[i], intervals = following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase, print_out=print_out)
        following_exp[i] = expected_following_in_pipe_single_phase(ehd, cf, phase)
        interval_dict[phase] = intervals
        save_single_histograms(following[i],
                               'following_in_pipe_2nd_mouse_in_pipe',
                               mice,
                               phase,
                               main_directory,
                               'following_in_pipe_2nd_mouse_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(following_exp[i],
                               'following_in_pipe_2nd_mouse_in_pipe_expected_time',
                               mice,
                               phase,
                               main_directory,
                               'following_in_pipe_2nd_mouse_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_in_pipe_2nd_mouse_in_pipe_relative_excess_following',
                               mice,
                               phase,
                               main_directory,
                               'following_in_pipe_2nd_mouse_in_pipe/histograms',
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
                                  'following_in_pipe_2nd_mouse_in_pipe/histograms',
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
    write_csv_rasters(mice,
                      phases,
                      following,
                      main_directory,
                      'following_in_pipe_2nd_mouse_in_pipe/raster_plots',
                      fname_)
    write_csv_rasters(mice,
                      phases,
                      (following-following_exp),
                      main_directory,
                      'following_in_pipe_2nd_mouse_in_pipe/raster_plots',
                      fname_exp)
    

    make_RasterPlot(main_directory,
                    'following_in_pipe_2nd_mouse_in_pipe/raster_plots',
                    following,
                    phases,
                    fname_,
                    mice,
                    title='# followings')
    make_RasterPlot(main_directory,
                    'following_in_pipe_2nd_mouse_in_pipe/raster_plots',
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

    return following, following_exp, phases, mice, interval_dict


@jit
def frequency_mouse_in_tube(ehd, mouse, st, en):
    ehd.mask_data(st, en)
    antennas = ehd.getantennas(mouse)
    times = ehd.gettimes(mouse)
    change_indices = np.where((np.array(antennas[1:]) - np.array(antennas[:-1])) != 0)[0]
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
        antenna = antennas[idx]
        next_antenna = antennas[idx+1]
        delta_t = times[idx+1] - times[idx]
        key = False
        if antenna % 2 and next_antenna == antenna + 1:
            key = next_antenna + antenna
        elif next_antenna % 2 and antenna == next_antenna + 1:
            key = next_antenna + antenna
            
        if key and delta_t < threshold:
            frequency[key] += 1
            window[key] += delta_t
    for key in frequency:
        frequency[key] /= (en-st)
    ehd.unmask_data()
    return frequency, window
   

@jit
def frequencies_for_all(phase, cf, ehd):
    mice = ehd.mice
    st, en = cf.gettime(phase)
    frequency = {}
    window = {}
    for mouse1 in mice:
        frequency[mouse1], window[mouse1] = frequency_mouse_in_tube(ehd, mouse1, st, en)
    return frequency, window
   
@jit
def expected_following_in_pipe_single_phase(ehd, cf, phase):
    mice = ehd.mice
    expected_followings = np.zeros((len(mice), len(mice)))
    frequency, window = frequencies_for_all(phase, cf, ehd)
    for j, mouse1 in enumerate(mice):
        for k, mouse2 in enumerate(mice):
            if mouse1 != mouse2:
                for key in window[mouse1]:
                    expected_followings[j, k] += window[mouse1][key]*frequency[mouse2][key]
                expected_followings[j, k] = np.round(expected_followings[j, k])
    return expected_followings
            

@jit
def following_for_all_in_pipe_condition(ehd, cf):
    mice = ehd.mice
    phases = cf.sections()
    phases = utils.filter_dark(phases)
    followings = np.zeros((len(phases), len(mice), len(mice)))
    for i, phase in enumerate(phases):
        st, en = cf.gettime(phase)
        for j, mouse1 in enumerate(mice):
            for k, mouse2 in enumerate(mice):
                if mouse2 != mouse1:
                    followings[i, j, k] = following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en, collect_intervals)
    return followings, phases, mice

def save_all_followings_data_to_file(directory, prefix, followings, ending=''):
    dir_ = utils.check_directory(directory, 'following_in_pipe')
    fname_base = "%s_followings_wyniki_%s_%s.csv"
    header = 'mouse pair; measurement\n'
    measurement_keys = ['antennas', 'mouse1_time', 'mouse1_delta_t', 'mouse2_time', 'mouse2_delta_t']
    phases = followings.keys()
    for phase in phases:
        fname = os.path.join(dir_, fname_base % (ending, prefix, phase))
        f = open(fname, 'w')
        f.write(header)
        for key in followings[phase]:
            for mk in measurement_keys:
                f.write(key +';'+mk)
                for meas in followings[phase][key][mk]:
                    f.write(';'+str(meas))
                f.write('\n')
        f.close()

def save_all_intervals_to_file(directory, prefix, followings, ending=''):
    dir_ = utils.check_directory(directory, 'intervals_following_in_pipe')
    fname_base = "%s_followings_intervals_wyniki_%s_%s.csv"
    header = 'mouse pair; following interval\n'
    phases = followings.keys()
    for phase in phases:
        fname = os.path.join(dir_, fname_base % (ending, prefix, phase))
        f = open(fname, 'w')
        f.write(header)
        for key in followings[phase]:
            f.write(key)
            for meas in followings[phase][key]:
                f.write(';'+str(meas))
            f.write('\n')
        f.close()

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
