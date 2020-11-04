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


