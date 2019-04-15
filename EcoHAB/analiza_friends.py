# -*- coding: utf-8 -*-
from __future__ import print_function, division
import EcoHab
import numpy as np
import matplotlib.pyplot as plt
import os
import utility_functions as utils
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
from numba import jit
from ExperimentConfigFile import ExperimentConfigFile


def prepare_mice_intervals(data_mice, address):
    ints = {}
    for mouse in data_mice.keys():
        ints[mouse] = utils.get_intervals_2_lists(data_mice[mouse], address)
    return ints


def check_interval(intervals_mouse1, intervals_mouse2, idx, new_idx):
    original_s, original_e = intervals_mouse1[0][idx], intervals_mouse1[1][idx]
    other_s, other_e = intervals_mouse2[0][new_idx], intervals_mouse2[1][new_idx]
    if  other_s > original_e:
        return False
    if original_s >  other_s and original_e > other_e:
        intervals_mouse1[0][idx] = other_e
        return False
    if original_s < other_s and original_e < other_e:
        intervals_mouse1[1][idx] = other_s
        return False
    if original_s >= other_s and original_e <= other_e:
        # delete mouse1 interval
        intervals_mouse1[0].remove(original_s)
        intervals_mouse1[1].remove(original_e)
        return True
    #  cut the original interval in half
    intervals_mouse1[1][idx] = other_s
    if other_e < original_e:
        intervals_mouse1[0].insert(idx + 1, other_e)
        intervals_mouse1[1].insert(idx + 1, original_e)
    if other_s == original_s:
        intervals_mouse1[0].remove(original_s)
        intervals_mouse1[1].remove(other_s)
        return True
    return False

def remove_overlapping_intervals(intervals_mouse1, intervals_mouse2):
    """"
    Eliminate all the intervals, when mouse 1 is with mouse 2
    """
    if len(intervals_mouse1[0]) == 0:
        return
    if len(intervals_mouse2[0]) == 0:
        return
    i = 0
    while i < len(intervals_mouse1[0]):
        new_idx = utils.get_idx_pre(intervals_mouse1[0][i], intervals_mouse2[0])
        if new_idx is None:
            new_idx = utils.get_idx_pre(intervals_mouse1[1][i], intervals_mouse2[0])
            if new_idx is None:
                i = i + 1
                continue
        if intervals_mouse2[1][new_idx] < intervals_mouse1[0][i]:
            new_idx += 1
        if new_idx  >= len(intervals_mouse2[0]):
            break
        removed = check_interval(intervals_mouse1, intervals_mouse2, i, new_idx)
        if not removed:
            i = i + 1


def mouse_alone(data_mice, address):
    ints = prepare_mice_intervals(data_mice, address)
    new_intervals = {}
    result = {}
    for mouse in ints.keys():
        mouse_ints = []
        mouse_ints.append(ints[mouse][0][:])
        mouse_ints.append(ints[mouse][1][:])
        for other_mouse in ints.keys():
            if mouse != other_mouse:
                remove_overlapping_intervals(mouse_ints, ints[other_mouse])
        new_intervals[mouse] = mouse_ints
    for mouse in new_intervals:
        result[mouse] = utils.get_duration(new_intervals[mouse][0], new_intervals[mouse][1])
    return result

@jit
def mice_overlap(ints1, ints2):
    """Return time overlap of mice m1 and m2 in cage <address>."""
    total_overlap = 0
    for int1 in ints1:
        for int2 in ints2:
            total_overlap += utils.interval_overlap(int1, int2)
    return total_overlap


def time_fraction_together_one_cage(ints1, ints2, total_time):
    assert total_time > 0
    return mice_overlap(ints1, ints2)/total_time

def total_time_results(mice_data, mice):
    result = np.zeros((4, len(mice)))
    for address in [1, 2, 3, 4]:
        for i,mouse in enumerate(mice):
            ints = utils.get_intervals(mice_data[mouse], address)
            result[address-1,i] = calculate_total_time(ints)
    return result


@jit
def mouse_alone_ehs(ehs, cf, main_directory, prefix):
    phases = utils.filter_dark(cf.sections())
    mice = ehs.mice
    output = np.zeros((4, len(mice), len(phases)+1))
    for phase, sec in enumerate(phases):
        times = cf.gettime(sec)
        data = utils.prepare_data(ehs, mice, times)
        for i in range(1,5):
            alone = mouse_alone(data, i)
            for j, mouse in enumerate(mice):
                output[i-1, j, phase] = alone[mouse]
    phases.append('ALL DARK')
    output[:,:,-1] = output[:,:,:-1].sum(axis=2)  # last column -- sum of activity in all dark phases
    write_csv_alone(output, phases, mice, main_directory, prefix)

def get_dark_light_data(phase, cf, ehs):
    print(phase)
    if phase == "dark" or phase == "DARK":
        phases = utils.filter_dark(cf.sections())
    elif phase == "light" or phase == "LIGHT":
        phases = utils.filter_light(cf.sections())
    
    
    mice = ehs.mice
    datas = []
    total_time = 0
    for i, ph in enumerate(phases):
        time = cf.gettime(ph)
        datas.append(utils.prepare_data(ehs, mice, time))
        total_time += time[1] - time[0]
    
    data = datas[0].copy()
    for dat in datas[1:]:
        for mm in dat.keys():
            data[mm].extend(dat[mm])
            
    return phases, total_time, mice, data

@jit
def in_cohort_sociability_all_dark_light(ehs, cf, main_directory, prefix, remove_mouse=None, phase="dark"):
    
    phases, total_time, mice, data = get_dark_light_data(phase, cf, ehs)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    fname = 'incohort_sociability_%s_ALL_%s' % (add_info_mice, phase)
    name_ = 'incohort_sociability_measured_time_%s_%s_ALL_%s.csv' % (prefix, add_info_mice, phase)
    name_exp_ = 'incohort_sociability_excess_time_%s_%s_ALL_%s.csv' % (prefix, add_info_mice, phase)

    full_results = np.zeros((1, len(mice), len(mice)))
    full_results_exp = np.zeros((1, len(mice), len(mice)))
    
    phase = "ALL_"+phase
    full_results[0], full_results_exp[0] = single_phase_results(data,
                                                                mice,
                                                                total_time=total_time)
    save_single_histograms(full_results[0],
                           'incohort_sociability_measured_time',
                           mice,
                           phase,
                           main_directory,
                           'in_cohort_sociability/histograms',
                           prefix,
                           additional_info=add_info_mice)
    save_single_histograms(full_results_exp[0],
                           'incohort_sociability_expected_time',
                           mice,
                           phase,
                           main_directory,
                           'in_cohort_sociability/histograms',
                           prefix,
                           additional_info=add_info_mice)
    save_single_histograms(full_results[0]-full_results_exp[0],
                           'incohort_sociability_excess_time',
                           mice,
                           phase,
                           main_directory,
                           'in_cohort_sociability/histograms',
                           prefix,
                           additional_info=add_info_mice)
    
    single_in_cohort_soc_plot(full_results[0],
                              full_results_exp[0],
                              mice,
                              phase,
                              fname,
                              main_directory,
                              'in_cohort_sociability/histograms',
                              prefix+add_info_mice)
    write_csv_rasters(mice,
                      [phase],
                      full_results,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_)
    write_csv_rasters(mice,
                      [phase],
                      full_results-full_results_exp,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_exp_)
    

    make_RasterPlot(main_directory,
                    'in_cohort_sociability/raster_plots',
                    full_results,
                    [phase],
                    name_,
                    mice,
                    vmin=0,
                    vmax=0.5,
                    title='% time together')
    make_RasterPlot(main_directory,
                    'in_cohort_sociability/raster_plots',
                    full_results-full_results_exp,
                    [phase],
                    name_exp_,
                    mice,
                    title='% time together',
                    vmin=-.25,
                    vmax=.25)
    
@jit
def in_cohort_sociability_all_phases(ehs, cf, main_directory, prefix, remove_mouse=None):
    phase="ALL"
    phases = [phase]
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    fname = 'incohort_sociability_%s_ALL' % add_info_mice
    name_ = 'incohort_sociability_measured_time_%s_%s_ALL.csv' % (prefix, add_info_mice)
    name_exp_ = 'incohort_sociability_excess_time_%s_%s_ALL.csv' % (prefix, add_info_mice)
    
    mice = ehs.mice
    time = cf.gettime("ALL")
    data = utils.prepare_data(ehs, mice, time)
    full_results = np.zeros((1, len(mice), len(mice)))
    full_results_exp = np.zeros((1, len(mice), len(mice)))
    total_time = time[1] - time[0]
    full_results[0], full_results_exp[0] = single_phase_results(data,
                                                                mice,
                                                                total_time=total_time)

    save_single_histograms(full_results[0],
                           'incohort_sociability_measured_time',
                           mice,
                           phase,
                           main_directory,
                           'in_cohort_sociability/histograms',
                           prefix,
                           additional_info=add_info_mice)
    save_single_histograms(full_results_exp[0],
                           'incohort_sociability_expected_time',
                           mice,
                           phase,
                           main_directory,
                           'in_cohort_sociability/histograms',
                           prefix,
                           additional_info=add_info_mice)
    save_single_histograms(full_results[0]-full_results_exp[0],
                           'incohort_sociability_excess_time',
                           mice,
                           phase,
                           main_directory,
                           'in_cohort_sociability/histograms',
                           prefix,
                           additional_info=add_info_mice)
    
    single_in_cohort_soc_plot(full_results[0],
                              full_results_exp[0],
                              mice,
                              phase,
                              fname,
                              main_directory,
                              'in_cohort_sociability/histograms',
                              prefix+add_info_mice)
    write_csv_rasters(mice,
                      phases,
                      full_results,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_)
    write_csv_rasters(mice,
                      phases,
                      full_results-full_results_exp,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_exp_)
    

    make_RasterPlot(main_directory,
                    'in_cohort_sociability/raster_plots',
                    full_results,
                    phases,
                    name_,
                    mice,
                    vmin=0,
                    vmax=0.5,
                    title='% time together')
    make_RasterPlot(main_directory,
                    'in_cohort_sociability/raster_plots',
                    full_results-full_results_exp,
                    phases,
                    name_exp_,
                    mice,
                    title='% time together',
                    vmin=-.25,
                    vmax=.25)
    
def single_phase_results(data, mice, total_time=43200.):
    results = np.zeros((len(mice), len(mice)))
    results_exp = np.zeros((len(mice), len(mice)))

    for ii, mouse1 in enumerate(mice):
        for jj, mouse2 in enumerate(mice):
            if ii < jj:
                res = mice_together(data, mouse1, mouse2, total_time)
                results[ii, jj] = res[0]
                results_exp[ii, jj] = res[1]
    return results, results_exp

    
def in_cohort_sociability(ehs, cf, main_directory, prefix, remove_mouse=None):
   
    mice = ehs.mice
    phases = utils.filter_dark(cf.sections())
    full_results = np.zeros((len(phases), len(mice), len(mice)))
    full_results_exp = np.zeros((len(phases), len(mice), len(mice)))
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    fname = 'incohort_sociability_%s' % add_info_mice
    name_ = 'incohort_sociability_measured_time_%s%s.csv' % (prefix, add_info_mice)
    name_exp_ = 'incohort_sociability_excess_time_%s%s.csv' % (prefix, add_info_mice)
    
    for idx_phase, phase in enumerate(phases):
        data = utils.prepare_data(ehs, mice, cf.gettime(phase))
        full_results[idx_phase], full_results_exp[idx_phase] = single_phase_results(data, mice)
        
        save_single_histograms(full_results[idx_phase],
                               'incohort_sociability_measured_time',
                               mice,
                               phase,
                               main_directory,
                               'in_cohort_sociability/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(full_results_exp[idx_phase],
                               'incohort_sociability_expected_time',
                               mice,
                               phase,
                               main_directory,
                               'in_cohort_sociability/histograms',
                               prefix,
                               additional_info=add_info_mice)
        save_single_histograms(full_results[idx_phase]-full_results_exp[idx_phase],
                               'incohort_sociability_excess_time',
                               mice,
                               phase,
                               main_directory,
                               'in_cohort_sociability/histograms',
                               prefix,
                               additional_info=add_info_mice)
        
        single_in_cohort_soc_plot(full_results[idx_phase],
                              full_results_exp[idx_phase],
                              mice,
                              phase,
                              fname,
                              main_directory,
                                  'in_cohort_sociability/histograms',
                                  prefix+add_info_mice)

    write_csv_rasters(mice,
                      phases,
                      full_results,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_)
    write_csv_rasters(mice,
                      phases,
                      full_results-full_results_exp,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_exp_)
    

    make_RasterPlot(main_directory,
                    'in_cohort_sociability/raster_plots',
                    full_results,
                    phases,
                    name_,
                    mice,
                    vmin=0,
                    vmax=0.5,
                    title='% time together')
    make_RasterPlot(main_directory,
                    'in_cohort_sociability/raster_plots',
                    full_results-full_results_exp,
                    phases,
                    name_exp_,
                    mice,
                    title='% time together',
                    vmin=-.25,
                    vmax=.25)


if __name__ == '__main__':
    homepath = os.path.expanduser("~/")
    threshold = 3
    new_path = "EcoHAB_data_November/Maciek_01_30_2018"
    path = os.path.join(homepath, new_path)
    prefix = utils.make_prefix(path)
    ehd = EcoHab.EcoHabData(path=path,
                            _ant_pos=None,
                            how_many_appearances=10)
    ehs = EcoHab.EcoHabSessions(ehd)
    cf = ExperimentConfigFile(path)
    tstart, tend = cf.gettime('ALL')
    directory = utils.results_path(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    mouse_alone_ehs(ehs, cf, directory, prefix)
    in_cohort_sociability(ehs, cf, directory, prefix)
    in_cohort_sociability_all_phases(ehs, cf, directory, prefix)
    in_cohort_sociability_all_dark_light(ehs, cf, directory, prefix, phase="light")
    in_cohort_sociability_all_dark_light(ehs, cf, directory, prefix, phase="dark")
