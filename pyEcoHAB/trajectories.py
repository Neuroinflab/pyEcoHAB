from __future__ import division, print_function, absolute_import
import os

import numpy as np
from pyEcoHAB.utility_functions import check_directory
from pyEcoHAB.plotting_functions import single_histogram_figures
from pyEcoHAB.plotting_functions import histograms_antenna_transitions
from pyEcoHAB.utils.for_loading import save_mismatches

directory = "antenna_transitions"

def save_antenna_transitions(transition_times, fname, res_dir, directory):
    dir_correct = os.path.join(res_dir, directory)
    out_dir = check_directory(dir_correct, "data")
    fname = os.path.join(out_dir, fname)
    f = open(fname, "w")
    for key in transition_times.keys():
        f.write("%s;" % key)
        for duration in transition_times[key]:
            f.write("%f;" % duration)
        f.write("\n")
    f.close()

def single_mouse_antenna_transitions(ants, ts):
    out = {}
    for i, a1 in enumerate(ants[:-1]):
        a2 = ants[i+1]
        key = "%s %s" % (a1, a2)
        if key not in out:
            out[key] = [] 
        out[key].append(ts[i+1]-ts[i])
    return out


def get_antenna_transitions(ehd):
    transition_times = {}
    for antenna1 in ehd.setup_config.all_antennas:
        for antenna2 in ehd.setup_config.all_antennas:
            key = "%s %s" % (antenna1, antenna2)
            transition_times[key] = []
        
    for mouse in ehd.mice:
        antennas = ehd.get_antennas(mouse)
        times = ehd.get_times(mouse)
        out = single_mouse_antenna_transitions(antennas, times)
        for key in out:
            transition_times[key].extend(out[key])
    histograms_antenna_transitions(transition_times, ehd.setup_config,
                                   ehd.res_dir, directory)
    save_antenna_transitions(transition_times, "transition_durations.csv",
                             ehd.res_dir, directory)
    return transition_times

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
    save_antenna_transitions(registration_trains, "train_durations.csv",
                             data.res_dir, directory)
    save_antenna_transitions(counts_in_trains, "counts_in_trains.csv",
                             data.res_dir, directory)
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

def incorrect_tunnel_single_mouse(keys, antennas, times, durations):
    count = {}
    total_count = {}
    for key in keys:
        count[key] = 0
        total_count[key] = 0
    for i, a1 in enumerate(antennas[:-1]):
        a2 = antennas[i+1]
        key = "%s %s" % (min(a1, a2), max(a1, a2))
        if key in keys:
            t1, t2 = times[i], times[i+1]
            d1 = durations[i]/1000
            total_count[key] += 1
            if t2 <= t1 + d1:
                count[key] += 1
    return count, total_count


def get_incorrect_tunnel_registrations(ehd):
    count = {}
    directions = ehd.setup_config.directions
    total_count = {}
    for direction in directions:
        a1, a2 = direction.split(" ")
        key = "%s %s" % (min(a1, a2), max(a1, a2))
        count[key] = 0
        total_count[key] = 0

    for mouse in ehd.mice:
        antennas1 = ehd.get_antennas(mouse)
        times1 = ehd.get_times(mouse)
        durations1 = ehd.get_durations(mouse)
        out_count, out_tot_count = incorrect_tunnel_single_mouse(count.keys(),
                                                                 antennas1,
                                                                 times1,
                                                                 durations1)
        for key in count:
            count[key] += out_count[key]
            total_count[key] += out_tot_count[key]

    fname = "incorrect_tunnel_registrations.csv"
    header = "tunnel, count, percentage of all passings through the tunnel\n"
    save_mismatches(count, total_count, ehd.res_dir, fname=fname, header=header)
    return count, total_count
                    
        
