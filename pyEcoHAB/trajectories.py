from __future__ import division, print_function, absolute_import
import os

import numpy as np
from pyEcoHAB.utility_functions import check_directory
from pyEcoHAB.plotting_functions import single_histogram_figures
from pyEcoHAB.plotting_functions import histograms_antenna_transitions

directory = "antenna_transitions"

def save_antenna_transitions(transition_times, config, res_dir,
                             directory):
    dir_correct = os.path.join(res_dir, directory)
    out_dir = check_directory(dir_correct, "data")
    fname = os.path.join(out_dir, "transition_durations.csv")
    f = open(fname, "w")
    for key in transition_times.keys():
        f.write("%s;" % key)
        for duration in transition_times[key]:
            f.write("%f;" % duration)
        f.write("\n")
    f.close()


def get_antenna_transitions(ehd):
    transition_times = {}
    for mouse in ehd.mice:
        antennas = ehd.get_antennas(mouse)
        times = ehd.get_times(mouse)
        for i, a1 in enumerate(antennas[:-1]):
            a2 = antennas[i+1]
            key = "%s %s" % (a1, a2)
            if key not in transition_times:
                transition_times[key] = []
            transition_times[key].append(times[i+1]-times[i])
    histograms_antenna_transitions(transition_times, ehd.setup_config,
                                   ehd.res_dir, directory)
    save_antenna_transitions(transition_times, ehd.setup_config, ehd.res_dir,
                             directory)
def get_registration_trains(data):
    title = "Series of registrations by "
    fname_duration = "total_duration_of_registration_trains"
    fname_count = "total_count_of_registration_trains"
    directory = "trains_of_registrations"
    registration_trains = {}
    counts_in_trains = {}
    for antenna in data.all_antennas:
        registration_trains[antenna] = []
        counts_in_trains[antenna] = []
    for mouse in data.mice:
        times = data.get_times(mouse)
        antennas = data.get_antennas(mouse)
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
           
    histograms_registration_trains(registration_trains, data.setup_config,
                                   fname_duration, data.res_dir, directory,
                                   title=title,
                                   xlabel="Duration (s)")
    histograms_registration_trains(counts_in_trains, data.setup_config,
                                   fname_count, data.res_dir, directory,
                                   title=title,
                                   xlabel="#registrations")
    #save_registration_trains(registration_trains, data.setup_config, data.res_dir)
    #save_registration_trains(counts_in_trains, data.setup_config, data.res_dir)


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

def mouse_registered_in_tunnel_by_2_antennas(ehd):
    pass
