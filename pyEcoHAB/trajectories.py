from __future__ import division, print_function, absolute_import
import os

import numpy as np
from pyEcoHAB.utility_functions import check_directory
from pyEcoHAB.plotting_functions import single_histogram_figures
directory = "antenna_transitions"


def histograms_antenna_transitions(transition_times, config, res_dir):    
    dir_correct = os.path.join(directory, "correct_antenna_transition")
    dir_incorrect = os.path.join(directory,
                                 "incorrect_antenna_transition")
    title_double_a = "consecutive crossing of antenna %s entrance to %s"
    max_count = 0
    nbins = {}
    title = {}
    fname = {}
    dir_name = {}
    xlogscale = {}
    xmin = 1000
    xmax = 0
    for key in transition_times.keys():
        first, last = key.split(" ")
        gen_key = "%s %s" % (min(first, last), max(first, last))
        if gen_key in config.mismatched_pairs:
            dir_name[key] = dir_incorrect
        else:
            dir_name[key] = dir_correct
        if first == last:
            title[key] = title_double_a % (first,
                                          config.address[first])
        else:
            if key in config.directions:
                title[key] = "%s (tunnel)" % key
            elif gen_key in config.mismatched_pairs:
                title[key] = key
            else:
                title[key] = "%s cage %s" %(key,
                                           config.address[first])
        fname[key] = "transition_times_antennas_%s" % key
        if len(transition_times[key]) > 1000:
            nbins[key] = 40
            if max(transition_times[key]) > 1000*min(transition_times[key]):
                xlogscale[key] = True
            else:
                xlogscale[key] = False
        else:
            nbins[key] = 10
            xlogscale[key] = False
        while True:
            if 0 in transition_times[key]:
                transition_times[key].remove(0)
            else:
                break
        hist, bins = np.histogram(transition_times[key], nbins[key])
        if xlogscale[key]:
            logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),
                                  len(bins))
            hist, bins = np.histogram(transition_times[key], bins=logbins)

 
        if max(hist) > max_count:
            max_count = max(hist) + 5
        if xmin > min(transition_times[key]):
            xmin =  min(transition_times[key]) - 0.5
        if xmax < max(transition_times[key]):
            xmax = max(transition_times[key]) + 0.5
    for key in transition_times.keys():
        if dir_name[key] == dir_correct:
            single_histogram_figures(transition_times[key], fname[key],
                                     res_dir, dir_name[key], title[key],
                                     nbins=nbins[key],
                                     xlogscale=xlogscale[key],
                                     xlabel="Transition times (s)",
                                     ylabel="count", xmin=xmin, xmax=xmax,
                                     ymin=0, ymax=max_count,
                                     fontsize=14, median_mean=True)
        else:
            single_histogram_figures(transition_times[key], fname[key],
                                     res_dir, dir_name[key], title[key],
                                     nbins=nbins[key],
                                     xlogscale=xlogscale[key],
                                     xlabel="Transition times (s)",
                                     ylabel="count",
                                     fontsize=14, median_mean=True)


def save_antenna_transitions(transition_times, config, res_dir):    
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
                                   ehd.res_dir)
    save_antenna_transitions(transition_times, ehd.setup_config, ehd.res_dir)
