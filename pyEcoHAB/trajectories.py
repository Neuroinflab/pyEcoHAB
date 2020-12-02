from __future__ import division, print_function, absolute_import
import os
from collections import OrderedDict
import numpy as np
from pyEcoHAB import utility_functions as utils
from pyEcoHAB.plotting_functions import single_histogram_figures
from pyEcoHAB.plotting_functions import histograms_antenna_transitions
from pyEcoHAB.utils.for_loading import save_mismatches
from pyEcoHAB.write_to_file import save_antenna_transitions

directory = "antenna_transitions"


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
                                 res_dir, prefix, delimiter):
    transition_times = {}
    all_phases, bin_labels = data_keys
    for idx_phase, ph in enumerate(all_phases):
        transition_times[ph] = {}
        for i, lab in enumerate(bin_labels):
            t_start, t_stop = phase_bounds[ph][lab]
            tunnels_antennas_dict = data[ph][lab]
            transition_times[ph][lab] = {}
            for key in setup_config.all_pairs:
                transition_times[ph][lab][key] = []
            for mouse in sorted(tunnels_antennas_dict.keys()):
                antennas = tunnels_antennas_dict[mouse]["antennas"]
                times = tunnels_antennas_dict[mouse]["times"]
                out = single_mouse_antenna_transitions(antennas, times)
                for key in out:
                    transition_times[ph][lab][key].extend(out[key])
    save_antenna_transitions(transition_times,
                             "transition_durations",
                             res_dir, prefix, directory, delimiter=delimiter)
    histograms_antenna_transitions(transition_times, setup_config,
                                   res_dir, directory,
                                   "transition_times_antennas", prefix)
    cages_tunnels = get_cage_tunnel_transitions(transition_times,
                                                setup_config)
    save_antenna_transitions(cages_tunnels,
                             "transition_durations",
                             res_dir, prefix, directory, delimiter=delimiter)
    
    histograms_antenna_transitions(cages_tunnels, setup_config,
                                   res_dir, directory,
                                   "transitions_all", prefix)
    light_dark = get_light_dark_transitions(transition_times)
    save_antenna_transitions(light_dark, 
                             "transition_durations",
                             res_dir, prefix, directory, delimiter=delimiter)
    histograms_antenna_transitions(light_dark, setup_config,
                                   res_dir, directory,
                                   "transition_times_antennas", prefix)
    return transition_times




def get_light_dark_transitions(transitions):
    if transitions.keys() == ["ALL"]:
        return {}
    out = {"dark": {0: {}}, "light": {0: {}}}
    for phase in transitions.keys():
        for label in transitions[phase].keys():
            for key in transitions[phase][label].keys():
                if "dark" in phase or "Dark" in phase or "DARK" in phase:
                    if key in out["dark"][0]:
                        out["dark"][0][key] += transitions[phase][label][key]
                    else:
                        out["dark"][0][key] = transitions[phase][label][key]
                elif  "light" in phase or "Light" in phase or "LIGHT" in phase:
                    if key in out["light"][0]:
                        out["light"][0][key] += transitions[phase][label][key]
                    else:
                        out["light"][0][key] = transitions[phase][label][key]
    return out


def get_cage_tunnel_transitions(transitions, setup_config):
    out = {}
    tunnel_pairs = setup_config.tunnel_pairs()
    cage_pairs = setup_config.cage_pairs()
    for phase in transitions.keys():
        out[phase] = {}
        for label in transitions[phase]:
            out[phase][label] = {"cages": [], "tunnels": []}
            for key in transitions[phase][label].keys():
                if key in tunnel_pairs:
                    out[phase][label]["tunnels"] += transitions[phase][label][key]
                elif key in cage_pairs:
                    out[phase][label]["cages"] += transitions[phase][label][key]
    return out


def get_antenna_transition_durations(ecohab_data, timeline, bins=12*3600,
                                      res_dir="", prefix="", remove_mouse="",
                                      delimiter=";"):
    """Save and plot histograms of durations between consecutive tag
    registrations by antenna pairs.

    Args:
        ecohab_data : Loader or Loader_like
           Eco-HAB dataset.
        timeline : Timeline
           timeline of the experiment.
        binsize : string or number
           time bins for calculating activity. Possible string values are:
           "ALL" -- calculate activity for the whole experiment,
           A number value specifies number of seconds in each bin, e.g. binsize
           equal 3600 results in 1 h bins.
        res_dir : string
           destination directory
           default value is the destination directory established
           for ecohab_data.
        prefix : string
           string added to the name of every generated results file
           default value is the prefix established for ecohab_data
        remove_mouse : string or list
           name of mouse or mice to be removed from the results file
           As a default activity will be established for every mouse registered
           in ecohab_data.
        delimiter : str, optional
           String or character separating columns

    """
    if prefix == "":
        prefix = ecohab_data.prefix
    if res_dir == "":
        res_dir = ecohab_data.res_dir
    mice = utils.get_mice(ecohab_data.mice, remove_mouse)
    function = utils.get_times_antennas_list_of_mice
    phases, times, data, keys = utils.prepare_binned_registrations(ecohab_data,
                                                                   timeline,
                                                                   bins,
                                                                   mice,
                                                                   function)
    transitions = antenna_transtions_in_phases(data, times, phases,
                                               keys, ecohab_data.setup_config,
                                               res_dir, prefix, delimiter=";")
    return transitions


def get_registration_trains(ecohab_data):
    title = "Series of registrations by "
    fname_duration = "total_duration_of_registration_trains"
    fname_count = "total_count_of_registration_trains"
    directory = "trains_of_registrations"
    registration_trains = {"ALL":{0:{}}}
    counts_in_trains = {"ALL":{0:{}}}
    for antenna in ecohab_data.all_antennas:
        registration_trains["ALL"][0][antenna] = []
        counts_in_trains["ALL"][0][antenna] = []
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
                    registration_trains["ALL"][0][previous_antenna].append(duration)
                    counts_in_trains["ALL"][0][previous_antenna].append(count)
                count = 1
                previous_antenna = a
                previous_t_start = times[i+1]
           
    histograms_registration_trains(registration_trains["ALL"][0], ecohab_data.setup_config,
                                   fname_duration, ecohab_data.res_dir, directory,
                                   title=title,
                                   xlabel="Duration (s)")
    histograms_registration_trains(counts_in_trains["ALL"][0], ecohab_data.setup_config,
                                   fname_count, ecohab_data.res_dir, directory,
                                   title=title,
                                   xlabel="#registrations")
    save_antenna_transitions(registration_trains, ["ALL"],
                             "train_durations.csv",
                             ecohab_data.res_dir, "", directory)
    save_antenna_transitions(counts_in_trains, ["ALL"], "counts_in_trains.csv",
                             ecohab_data.res_dir, "", directory)
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

