from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utils
import numpy as np
import matplotlib.pyplot as plt
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot

bins = 2000
homepath = os.path.expanduser("~/")
threshold = 2
from numba import jit
titles = {
    3: '12',
    7: '34',
    11: '56',
    15: '78',
}
def in_tube(antenna, next_antenna):
    if antenna % 2:
        if next_antenna  == antenna  + 1:
            return True
    else:
        if next_antenna == antenna - 1:
            return True
    return False

def make_out():
    return {'antennas':[],
           'mouse1_time':[],
           'mouse2_time':[],
           'mouse1_delta_t':[],
           'mouse2_delta_t':[]}

def add_to_out(out, key, t1, t2, delta_t1, delta_t2):
     out['antennas'].append(key)
     out['mouse1_time'].append(t1)
     out['mouse2_time'].append(t2)
     out['mouse1_delta_t'].append(delta_t1)
     out['mouse2_delta_t'].append(delta_t2)

@jit
def check_2nd_mouse(antenna1, next_antenna1, t1, threshold, antennas2, times2):
    idxs1 = np.where(np.array(times2) >= t1)[0]
    idxs2 = np.where(np.array(times2) <= t1 + threshold)[0]
    common_idxs = list(set(idxs1)&set(idxs2))
    followings = 0
    for ci in common_idxs:
        antenna2 = antennas2[ci]
        if ci + 1 < len(antennas2) and antenna2 == antenna1:
            next_antenna2 = antennas2[ci+1]
            delta_t2 = times2[ci+1] - times2[ci]
            if in_tube(antenna2, next_antenna2):
                followings = 1
    return followings
                           
@jit
def following_2_mice(ehd, mouse1, mouse2, st, en):
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
        if in_tube(antenna1, next_antenna1):
            followings += check_2nd_mouse(antenna1, next_antenna1, times1[idx], threshold, antennas2, times2)
    return followings

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
            
        if key:
            frequency[key] += 1
            window[key] += delta_t
    for key in frequency:
        frequency[key] /= (en-st)
    ehd.unmask_data()
    return frequency, window


def histo():
    figure, ax = plt.subplots(2, 2)
    n, bins1, patches = ax[0][0].hist(out[3], bins=bins)
    ax[0][0].set_title('12')
    ax[0][0].set_xscale("log", nonposx='clip')
    n, bins2, patches = ax[0][1].hist(out[7], bins=bins)
    ax[0][1].set_title('34')
    ax[0][1].set_xscale("log", nonposx='clip')
    n, bins3, patches = ax[1][0].hist(out[11], bins=bins)
    ax[1][0].set_title('56')
    ax[1][0].set_xscale("log", nonposx='clip')
    n, bins4, patches = ax[1][1].hist(out[15], bins=bins)
    ax[1][1].set_title('78')
    ax[1][1].set_xscale("log", nonposx='clip')

    max_lim = max(max(bins1[-1], bins2[-1]), max(bins3[-1], bins4[-1]))
    
    ax[0][0].set_xlim([0, max_lim])
    ax[0][1].set_xlim([0, max_lim])
    ax[1][0].set_xlim([0, max_lim])
    ax[1][1].set_xlim([0, max_lim])

    
@jit            
def following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en):
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
        if in_tube(antenna1, next_antenna1):
            followings += check_2nd_mouse(antenna1, next_antenna1, times1[idx], delta_t1, antennas2, times2)
    return  followings



def following_single_phase(ehd, cf, phase):
    mice = ehd.mice
    st, en = cf.gettime(phase) 
    followings = np.zeros((len(mice), len(mice)))
    for j, mouse1 in enumerate(mice):
        for k, mouse2 in enumerate(mice):
            if mouse2 != mouse1:
                followings[j, k] = following_2_mice(ehd, mouse1, mouse2, st, en)
    return followings
    
@jit
def following_for_all(ehd, cf, main_directory, prefix, remove_mouse=None, threshold=threshold):
    phases = cf.sections()
    phases = utils.filter_dark_light(phases)
    mice = ehd.mice
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    following = np.zeros((len(phases), len(mice), len(mice)))
    following_exp = np.zeros((len(phases), len(mice), len(mice)))
    fname = 'following_in_pipe_threshold_%f_%s' % (threshold, add_info_mice)
    name_ = 'following_in_pipe_threshold_%f_%s%s' % (threshold, prefix, add_info_mice)
    name_exp_ = 'relative_following_in_pipe_excess_threshold_%f_%s%s.csv' % (threshold, prefix, add_info_mice)
    for i, phase in enumerate(phases):
        following[i] = following_single_phase(ehd, cf, phase)
        following_exp[i] = expected_following_in_pipe_single_phase_threshold(ehd, cf, phase, threshold)
        save_single_histograms(following[i],
                               'following_in_pipe_threshold_%f'%threshold,
                               mice,
                               phase,
                               main_directory,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(following_exp[i],
                               'following_in_pipe_threshold_%f_expected_time'%threshold,
                               mice,
                               phase,
                               main_directory,
                               'following_in_pipe/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms((following[i]-following_exp[i]),
                               'following_in_pipe_threshold_%f_relative_excess_following'%threshold,
                               mice,
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
                                          'histogram of # excess followings',])
    write_csv_rasters(mice,
                      phases,
                      following,
                      main_directory,
                      'following_in_pipe/raster_plots',
                      name_)
    write_csv_rasters(mice,
                      phases,
                      (following-following_exp),
                      main_directory,
                      'following_in_pipe/raster_plots',
                      name_exp_)
    

    make_RasterPlot(main_directory,
                    'following_in_pipe/raster_plots',
                    following,
                    phases,
                    name_,
                    mice,
                    title='# followings')
    make_RasterPlot(main_directory,
                    'following_in_pipe/raster_plots',
                    (following-following_exp),
                    phases,
                    name_exp_,
                    mice,
                    title='relative excess following')                     

def following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase):
    mice = ehd.mice
    st, en = cf.gettime(phase) 
    followings = np.zeros((len(mice), len(mice)))
    for j, mouse1 in enumerate(mice):
        for k, mouse2 in enumerate(mice):
            if mouse2 != mouse1:
                followings[j, k] = following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en)
    return followings
    
def following_for_all_2nd_mouse_in_pipe(ehd, cf, main_directory, prefix, remove_mouse=None):
    phases = cf.sections()
    phases = utils.filter_dark_light(phases)
    mice = ehd.mice
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    following = np.zeros((len(phases), len(mice), len(mice)))
    following_exp = np.zeros((len(phases), len(mice), len(mice)))
    fname = 'following_in_pipe_2nd_mouse_in_pipe_%s' % (add_info_mice)
    fname_ = 'following_in_pipe_2nd_mouse_in_pipe_%s%s' % (prefix, add_info_mice)
    fname_exp = 'relative_following_in_pipe_excess_2nd_mouse_in_pipe_%s%s.csv' % (prefix, add_info_mice)
    for i, phase in enumerate(phases):
        following[i] = following_2nd_mouse_in_pipe_single_phase(ehd, cf, phase)
        following_exp[i] = expected_following_in_pipe_single_phase(ehd, cf, phase)
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
                                          'histogram of # excess followings',])
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
def expected_following_in_pipe_single_phase_threshold(ehd, cf, phase, threshold=threshold):
    mice = ehd.mice
    st, en = cf.gettime(phase) 
    expected = np.zeros((len(mice), len(mice)))
    frequency, window = frequencies_for_all(phase, cf, ehd)
    for j, mouse1 in enumerate(mice):
        for k, mouse2 in enumerate(mice):
            if mouse2 != mouse1:
                for key in window[mouse1]:
                    expected[j, k] += frequency[mouse1][key]*threshold*frequency[mouse2][key]*(en-st)
    return expected
    
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
    phases = utils.filter_dark_light(phases)
    followings = np.zeros((len(phases), len(mice), len(mice)))
    for i, phase in enumerate(phases):
        st, en = cf.gettime(phase)
        for j, mouse1 in enumerate(mice):
            for k, mouse2 in enumerate(mice):
                if mouse2 != mouse1:
                    followings[i, j, k] = following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en)
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


if __name__ == '__main__':

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
        following_for_all_2nd_mouse_in_pipe(ehd, cf, res_dir, prefix, remove_mouse=None)
        following_for_all(ehd, cf, res_dir, prefix, remove_mouse=None)
        
