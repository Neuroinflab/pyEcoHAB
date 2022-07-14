# SPDX-License-Identifier: LGPL-2.1-or-later
# !/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
from collections import OrderedDict

from .utils import general as utils
from .write_to_file import write_registrations_stats
from .plotting_functions import single_timeline_heat_map


def get_single_antenna_stats(ecohab_data, timeline, binsize, antennas="ALL",
                             res_dir="", prefix="", remove_mouse="",
                             delimiter=";"):
    """
    Count number and combined durations of registrations of each mouse tag
    by specified antennas in bins of size binsize for tags
    registered in every phase.

    Args:
        ecohab_data : Loader or Loader_like
           Eco-HAB dataset.
        timeline : Timeline
           timeline of the experiment.
        binsize : number
           time bins for calculating activity
           A number value specifies number of seconds in each bin, e.g. binsize
           equal 3600 results in 1 h bins.
        antennas: string, int or list of ints
           Ids of registering antennas.
           Default value is all antennas
        res_dir : string
           destination directory
           default value is the destination directory established for
           ecohab_data.
        prefix : string
           string added to the name of every generated results file
           default value is the prefix established for ecohab_data
        remove_mouse : string or list
           name of mouse or mice to be removed from the results file
           As a default activity will be established for every mouse registered
           in ecohab_data.
        delimiter : str, optional
           String or character separating columns.
    """
    if prefix == "":
        prefix = ecohab_data.prefix
    if res_dir == "":
        res_dir = ecohab_data.res_dir
    mice = utils.get_mice(ecohab_data.mice, remove_mouse)
    add_info_mice = utils.add_info_mice_filename(remove_mouse)
    out_dir = os.path.join("other_variables", "registration_stats")
    bin_ = binsize/3600
    fname_durations = "registration_duration_%4.2fh" % bin_
    fname_count = "registration_count_%4.2fh" % bin_
    shortest_phase = utils.get_shortest_phase_duration(timeline)
    if binsize <= shortest_phase:
        phases = timeline.sections()
        times = [timeline.get_time_from_epoch(phase) for phase in phases]
    else:
        t_start, t_end = timeline.get_time_from_epoch("ALL")
        phases = []
        times = []
        i = 0
        while t_start < t_end:
            phases.append("%dxbin_%5.2fh" % (i, bin_))
            times.append((t_start, t_start+binsize))
            t_start += binsize
            i += 1
    if antennas == "ALL":
        antennas = sorted(set(ecohab_data.get_antennas(ecohab_data.mice)))
    if antennas in ecohab_data.all_antennas:
        antennas = [antennas]
    if not isinstance(antennas, list):
        raise Exception("""Incorrect antenna format.
        You should either provide a list of ints, an int or 'ALL'""")

    for i, phase in enumerate(phases):
        t_start, t_end = times[i]
        count = OrderedDict()
        durations = OrderedDict()
        for antenna in antennas:
            count[antenna] = OrderedDict()
            durations[antenna] = OrderedDict()
            for mouse in mice:
                results = ecohab_data.get_registration_stats(mouse,
                                                             t_start,
                                                             t_end,
                                                             antenna,
                                                             binsize)
                count[antenna][mouse], durations[antenna][mouse] = results

            single_timeline_heat_map(durations[antenna],
                                     res_dir,
                                     mice,
                                     prefix,
                                     phase,
                                     binsize,
                                     antenna,
                                     out_dir)
        new_fname_count = "%s_%s" % (fname_count, antenna)
        new_fname_durations = "%s_%s" % (fname_durations, antenna)
        write_registrations_stats(count, phase, mice, binsize,
                                  new_fname_count, res_dir,
                                  out_dir, prefix,
                                  add_info=add_info_mice, delimiter=";")
        write_registrations_stats(durations, phase, mice, binsize,
                                  new_fname_durations, res_dir,
                                  out_dir, prefix,
                                  add_info=add_info_mice, delimiter=";")
