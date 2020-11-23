from __future__ import division, print_function, absolute_import
import os

import numpy as np
from pyEcoHAB import utility_functions as uf
from pyEcoHAB.plotting_functions import single_histogram_figures
from pyEcoHAB.plotting_functions import histograms_antenna_transitions
from pyEcoHAB.plotting_functions import histograms_transitions_cages_tunnels
from pyEcoHAB.utils.for_loading import save_mismatches

directory = "antenna_transitions"

def save_antenna_transitions(transition_times, fname, res_dir, directory):
    dir_correct = os.path.join(res_dir, directory)
    out_dir = uf.check_directory(dir_correct, "data")
    fname = os.path.join(out_dir, fname)
    f = open(fname, "w")
    for key in transition_times.keys():
        f.write("%s;" % key)
        for duration in transition_times[key]:
            f.write("%f;" % duration)
        f.write("\n")
    f.close()

def single_mouse_antenna_transitions(antennas1, times1):
    out = {}
    for i, a1 in enumerate(antennas1[:-1]):
        a2 = antennas1[i+1]
        key = "%s %s" % (a1, a2)
        if key not in out:
            out[key] = [] 
        out[key].append(times1[i+1]-times1[i])
    return out


def antenna_transtions_in_phases(data, phase_bounds, phases,
                                 data_keys, setup_config,
                                 res_dir):
    transition_times = {}
    all_phases, bin_labels = data_keys
    for idx_phase, ph in enumerate(all_phases):
        new_phase = phases[idx_phase]
        transition_times[ph] = {}
        for i, lab in enumerate(bin_labels):
            t_start, t_stop = phase_bounds[ph][lab]
            tunnels_antennas_dict = data[ph][lab]
            transition_times[ph][lab] = {}
            for key in setup_config.all_pairs:
                transition_times[ph][lab][key] = []
            for mouse in tunnels_antennas_dict.keys():
                antennas = tunnels_antennas_dict[mouse]["antennas"]
                times = tunnels_antennas_dict[mouse]["times"]
                out = single_mouse_antenna_transitions(antennas, times)
                for key in out:
                    transition_times[ph][lab][key].extend(out[key])
            save_antenna_transitions(transition_times[ph][lab],
                                     "transition_durations_%s%s.csv" % (new_phase, lab),
                                     res_dir, "antenna_transitions")

    histograms_antenna_transitions(transition_times, setup_config,
                                   res_dir, "antenna_transitions")
    return transition_times


def get_antenna_transitions(ecohab_data, timeline, bins=12*3600):
    """Save and plot histograms of consecutive tag registrations
    by pairs of antennas
    All - all phases
    filter_dark, filter_light
    """

    phases, tot_times, data, data_keys = uf.prepare_binned_registrations(ecohab_data,
                                                                         timeline,
                                                                         bins,
                                                                         ecohab_data.mice,
                                                                         uf.get_times_antennas_list_of_mice)
    transitions = antenna_transtions_in_phases(data, tot_times, phases,
                                               data_keys, ecohab_data.setup_config,
                                               ecohab_data.res_dir)
    return transitions


def get_registration_trains(ecohab_data):
    title = "Series of registrations by "
    fname_duration = "total_duration_of_registration_trains"
    fname_count = "total_count_of_registration_trains"
    directory = "trains_of_registrations"
    registration_trains = {}
    counts_in_trains = {}
    for antenna in ecohab_data.all_antennas:
        registration_trains[antenna] = []
        counts_in_trains[antenna] = []
    for mouse in ecohab_data.mice:
        times = ecohab_data.get_times(mouse)
        antennas = ecohab_data.get_antennas(mouse)
        previous_antenna = antennas[0]
        previous_t_start = times[0]
        count = 1
        i = 1
        for i, a in enumerate(antennas[1:]):
            if a == previous_antenna:
                count += 1
            else:
                if count > 2:
                    duration = times[i] - previous_t_start
                    registration_trains[previous_antenna].append(duration)
                    counts_in_trains[previous_antenna].append(count)
                count = 1
                previous_antenna = a
                previous_t_start = times[i+1]
           
    histograms_registration_trains(registration_trains, ecohab_data.setup_config,
                                   fname_duration, ecohab_data.res_dir, directory,
                                   title=title,
                                   xlabel="Duration (s)")
    histograms_registration_trains(counts_in_trains, ecohab_data.setup_config,
                                   fname_count, ecohab_data.res_dir, directory,
                                   title=title,
                                   xlabel="#registrations")
    save_antenna_transitions(registration_trains, "train_durations.csv",
                             ecohab_data.res_dir, directory)
    save_antenna_transitions(counts_in_trains, "counts_in_trains.csv",
                             ecohab_data.res_dir, directory)
    return registration_trains, counts_in_trains


def histograms_registration_trains(data_dict, config, fname, res_dir, directory,
                                   title, xlabel=""):
    
    titles = {}
    fnames = {}
    xmin = 1000
    xmax = 0
    max_count = 0
    nbins = 30
    xlogscale = True
    for key in data_dict.keys():
        if not len(data_dict[key]):
            continue
        titles[key] = "%s %s" % (title, key)
        fnames[key] = "%s_%s" % (fname, key)
        hist, bins = np.histogram(data_dict[key], nbins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),
                                  len(bins))
        hist, bins = np.histogram(data_dict[key], bins=logbins)
        if max(hist) > max_count:
            max_count = max(hist) + 1
        if xmin > min(data_dict[key]):
            xmin =  min(data_dict[key]) - 0.5
        if xmax < max(data_dict[key]):
            xmax = max(data_dict[key]) + 0.5
        len(data_dict[key])
    for key in data_dict.keys():
        if not len(data_dict[key]):
            continue
        single_histogram_figures(data_dict[key], fnames[key],
                                 res_dir, directory, titles[key],
                                 nbins=nbins, xlogscale=xlogscale,
                                 xlabel=xlabel,
                                 ylabel="count", xmin=xmin, xmax=xmax,
                                 ymin=0, ymax=max_count,
                                 fontsize=14, median_mean=True)

