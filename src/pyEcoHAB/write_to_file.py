# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
from .utils import general as utils


def make_header_for_activity(phases, delimiter):
    header = 'mouse%s\"time [h]\"' % delimiter
    for phase in phases:
        header += delimiter + '\"' + phase + "\""
    return header


def write_single_chamber(f, header, phases, mice, time, data_stim,
                         delimiter, floats=False):

    f.write(header+'\n')
    for i, mouse in enumerate(mice):
        longest = 0
        for phase in phases:
            if len(data_stim[phase][mouse]) > longest:
                longest = len(data_stim[phase][mouse])
        lines = [mouse for i in range(longest)]
        for phase in phases:
            for k, t in enumerate(time[phase]):
                if phase == phases[0]:
                    lines[k] += '%s%3.2f' % (delimiter, t/3600)
                try:
                    if floats:
                        lines[k] += "%s%7.3f" % (delimiter,
                                                 data_stim[phase][mouse][k])
                    else:
                        lines[k] += delimiter + str(data_stim[phase][mouse][k])
                except IndexError:
                    lines[k] += delimiter

        for line in lines:
            f.write(line + '\n')


def save_data_cvs(data, phases, mice, bin_labels, fname,
                  path, which, headers, target_dir="activity",
                  delimiter=";"):
    new_path = os.path.join(path, target_dir)
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    fname = os.path.join(new_path, fname)
    print(fname)
    f = open(fname, "w")
    head = make_header_for_activity(phases, delimiter)
    for stim in which:
        for j, h in enumerate(headers):
            f.write("%s %s\n" % (h, stim))
            write_single_chamber(f, head, phases, mice, bin_labels,
                                 data[stim][j], delimiter, floats=j)
    f.close()


def write_binned_data(data_stim, fname, mice, bin_labels, phase,
                      path, target_dir, prefix, additional_info="",
                      delimiter=";", full_dir_tree=True):
    if full_dir_tree:
        new_path = os.path.join(path, target_dir, "data")
    else:
        new_path = os.path.join(path, target_dir)
    fname = os.path.join(new_path, '%s_%s_%s_%s.csv' % (fname,
                                                        prefix,
                                                        phase,
                                                        additional_info))
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    print(fname)
    f = open(fname, "w")
    header = make_header_for_activity(mice, delimiter)

    f.write(header+'\n')
    for mouse1 in mice:
        lines = [mouse1 for l in range(len(bin_labels))]
        for j, mouse2 in enumerate(mice):
            for k, t in enumerate(bin_labels):
                if not j:
                    lines[k] += '%s%3.2f' % (delimiter, t/3600)
                try:
                    lines[k] += delimiter + str(data_stim[t][mouse1][mouse2])
                except TypeError:
                    lines[k] += delimiter
        for line in lines:
            f.write(line + '\n')
    f.close()


def save_single_histograms(result, fname, mice, phase, main_directory,
                           directory, prefix, additional_info="",
                           delimiter=";", full_dir_tree=True):
    if full_dir_tree:
        new_name = os.path.join(directory, 'data')
    else:
        new_name = directory
    directory = utils.check_directory(main_directory, new_name)
    fname = os.path.join(directory, '%s_%s_%s_%s.csv' % (fname,
                                                         prefix,
                                                         phase,
                                                         additional_info))
    print(fname)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
        return None
    for i, mouse in enumerate(mice):
        f.write(delimiter + mouse)
    f.write('\n')
    for i, mouse in enumerate(mice):
        f.write(mouse)
        for j, mouse in enumerate(mice):
            f.write(delimiter + str(result[i, j]))
        f.write('\n')


def write_csv_rasters(mice, phases, output, directory,
                      dirname, fname, symmetrical=True,
                      reverse=False, delimiter=";", prefix="",
                      full_dir_tree=True):
    if full_dir_tree:
        new_name = os.path.join(dirname, 'data')
    else:
        new_name = dirname
    directory = utils.check_directory(directory, new_name)
    if prefix:
        fname = "%s_%s" % (prefix, fname)
    fname = os.path.join(directory, fname)
    
    print(fname)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
    header = 'mouse pair'
    for phase in phases:
        header += delimiter + str(phase) + " h"
    header += '\n'
    f.write(header)
    if symmetrical:
        new_output, pairs = utils.make_table_of_pairs(output, phases, mice)
    else:
        new_output, pairs = utils.make_table_of_all_mouse_pairs(output, phases,
                                                                mice, reverse)
    for i, pair in enumerate(pairs):
        f.write(pair)
        for j in range(len(phases)):
            f.write(delimiter + str(new_output[i, j]))
        f.write('\n')
    f.close()


def write_csv_tables(results, phases, mice, main_directory,
                     dirname, fname, prefix,
                     delimiter=";"):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(main_directory, new_name)
    fname = os.path.join(directory, '%s_%s.csv' % (fname, prefix))
    print(fname)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return

    phase_pairs = [phases[2*i] + (len(mice) + 2)*delimiter +
                   phases[2*i + 1]+'\n' for i in range(len(phases)//2)]
    new_results = [(results[2*i],
                    results[2*i + 1])for i in range(len(phases)//2)]
    mice_header = ""
    mice_len = 2
    if len(phases) == 1:
        mice_len = 1
    for i in range(mice_len):
        for mouse in mice:
            mice_header += delimiter + mouse
        if not i:
            mice_header += 2*delimiter
    mice_header += '\n'

    if len(phases) % 2:
        phase_pairs.append((phases[-1]+'\n'))
        new_results.append((results[-1],))
    for i, new_pair in enumerate(phase_pairs):
        f.write(new_pair)
        if len(new_results[i]) > 1:
            f.write(mice_header)
        else:
            single = 1+len(mice)*(len(mice[0])+1)
            f.write(mice_header[:single]+'\n')
        for j, mouse in enumerate(mice):
            for l, nr in enumerate(new_results[i]):
                f.write(mouse)
                for k, mouse in enumerate(mice):
                    f.write(delimiter + str(nr[j, k]))
                if not l:
                    f.write(delimiter*2)
            f.write('\n')
        f.write('\n')
    f.close()


def write_csv_alone(alone, phases, main_directory, prefix,
                    header='Mice alone in %s\n',
                    fname='mouse_alone_%s.csv',
                    directory="solitude",
                    delimiter=";"):
    directory = utils.check_directory(main_directory, directory)
    fname = os.path.join(directory, fname % prefix)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return
    phases_header = 'Mouse'
    for phase in phases:
        phases_header += delimiter + phase
    phases_header += '\n'
    for address in alone.keys():
        f.write(header % address)
        f.write(phases_header)
        for mouse in alone[address].keys():
            f.write(mouse)
            for phase in phases:
                f.write(delimiter + str(alone[address][mouse][phase]))
            f.write('\n')
    f.close()


def write_interpair_intervals(results, main_directory,
                              directory, fname, prefix,
                              additional_info="",
                              delimiter=";",
                              full_dir_tree=True):
    if full_dir_tree:
        new_name = os.path.join(main_directory, 'data')
    else:
        new_name = main_directory
    directory = utils.check_directory(directory, new_name)
    fname = os.path.join(directory, '%s_%s_%s.csv' % (fname,
                                                      prefix,
                                                      additional_info))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
        return None
    print(fname)
    f.write("followed mouse %s following mouse %s intervals\n" % (delimiter,
                                                                  delimiter))
    keys = sorted(results.keys())
    for key in keys:
        mouse1, mouse2 = key.split('|')
        f.write("%s%s%s%s" % (mouse1, delimiter, mouse2, delimiter))
        for interval in results[key]:
            f.write("%f%s" % (interval, delimiter))
        f.write("\n")
    f.close()


def save_visit_duration(results, time, phase, mice,
                        fname, main_directory, directory,
                        prefix, add_info="", delimiter=";"):
    new_dir = os.path.join(main_directory, directory)
    new_dir = utils.check_directory(new_dir, "data")
    for mouse in mice:
        new_name = os.path.join(new_dir, '%s_%s_%s_%s_%s.csv' % (fname,
                                                                 mouse,
                                                                 phase,
                                                                 prefix,
                                                                 add_info))
        print(new_name)
        f = open(new_name, "w")
        for address in results.keys():
            f.write("Visit durations to %s " % address)
            f.write("time%s durations\n" % delimiter)
            for j, out in enumerate(results[address][mouse]):
                f.write("%2.2f" % time[j])
                for single in out:
                    f.write("%s%2.2f" % (delimiter, single))
                f.write("\n")
        f.close()


def write_bootstrap_results(results, phase, mice_list,
                            fname, main_directory,
                            directory, prefix,
                            add_info="", delimiter=";",
                            full_dir_tree=True):
    new_dir = os.path.join(main_directory, directory)
    if full_dir_tree:
        new_dir = utils.check_directory(new_dir, "data")
    else:
        new_dir = utils.check_directory(new_dir)
    new_name = os.path.join(new_dir, '%s_%s_%s_%s.csv' % (fname,
                                                          phase.replace(' ',
                                                                        '_'),
                                                          prefix, add_info))

    f = open(new_name, "w")
    for mouse1 in mice_list:
        for mouse2 in mice_list:
            if mouse1 != mouse2:
                key = "%s|%s" % (mouse1, mouse2)
                f.write(key)
                for value in results[mouse1][mouse2]:
                    f.write(delimiter + str(value))
                f.write("\n")
    f.close()


def write_registrations_stats(crossings, phase, mice_list,
                              binsize, fname, main_directory,
                              directory, prefix,
                              add_info="", delimiter=";"):
    header = ""
    new_dir = os.path.join(main_directory, directory)
    new_dir = utils.check_directory(new_dir, "data")
    new_name = os.path.join(new_dir, '%s_%s_%s_%s.csv' % (fname,
                                                          phase.replace(' ',
                                                                        '_'),
                                                          prefix, add_info))

    antennas = sorted(crossings.keys())
    n_rows = len(crossings[antennas[0]][mice_list[0]])

    f = open(new_name, "w")
    for i in range(n_rows):
        header += "%s%4.2f" % (delimiter, i*binsize/3600)
    f.write(delimiter)
    for antenna in antennas:
        f.write("Antenna %s" % antenna + (n_rows)*delimiter)
    f.write('\n')
    for antenna in antennas:
        f.write(header)
    f.write("\n")
    for mouse in mice_list:
        f.write(mouse)
        for antenna in antennas:
            for i, row in enumerate(range(n_rows)):
                f.write(delimiter+str(crossings[antenna][mouse][i]))
        f.write("\n")
    f.close()


def save_antenna_transitions(transition_times,
                             fname, res_dir, prefix, directory,
                             delimiter=";"):
    dir_correct = os.path.join(res_dir, directory)
    out_dir = utils.check_directory(dir_correct, "data")
    for idx_phase, phase in enumerate(transition_times.keys()):
        new_phase = phase.replace(" ", "_")
        for label in transition_times[phase].keys():
            new_fname = "%s_%s_%s.csv" % (fname, new_phase, label)
            new_path = os.path.join(out_dir, new_fname)
            print(new_path)
            f = open(new_path, "w")
            for key in transition_times[phase][label].keys():
                f.write("%s%s" % (key, delimiter))
                for duration in transition_times[phase][label][key]:
                    f.write("%f%s" % (duration, delimiter))
                f.write("\n")
            f.close()


def write_sum_data(data, fname, mice, bin_labels, phases,
                   path, target_dir, prefix, additional_info="",
                   delimiter=";", bool_bins=bool,
                   full_dir_tree=True):
    if full_dir_tree:
        new_path = os.path.join(path, target_dir, "data")
    else:
        new_path = os.path.join(path, target_dir)
    fname = os.path.join(new_path, '%s_%s_%s.csv' % (fname, prefix,
                                                     additional_info))
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    print(fname)
    f = open(fname, "w")
    header = 'mouse'

    for phase in phases:
        if bool_bins == True:
            for bin in bin_labels[phase]:
                header += delimiter + str(bin/3600) + "h "+ str(phase)
        else:
            header += delimiter + str(phase)
    header += '\n'
    f.write(header)
    phase = phases[0]

    for mouse_label in mice:
        f.write(mouse_label + delimiter)
        for phase in phases:
            for bi in bin_labels[phase]:
                for mouse in mice:
                    if mouse == mouse_label:
                        try:
                            f.write(str(data[phase][bi][mouse]) + delimiter)
                        except KeyError:
                            f.write(delimiter)
                    else:
                        continue
        f.write("\n")
    f.close()

def write_two_values(data1, data2, list_of_param, fname, mice, bin_labels,
                     phases, path, target_dir, prefix, additional_info="",
                     delimiter=";", full_dir_tree=True):
    if full_dir_tree:
        new_path = os.path.join(path, target_dir, "data")
    else:
        new_path = os.path.join(path, target_dir)
    fname = os.path.join(new_path, '%s_%s_%s.csv' % (fname, prefix,
                                                     additional_info))
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    print(fname)
    f = open(fname, "w")
    header = 'mouse'

    for phase in phases:
        for bi in bin_labels[phase]:
            for param in list_of_param:
                header += delimiter + str(param)+ " "
                header += str(bi / 3600) + "h " + str(phase)
    header += '\n'
    f.write(header)

    phase = phases[0]

    for mouse_label in mice:
        f.write(mouse_label + delimiter)
        for phase in phases:
            for bi in bin_labels[phase]:
                for param in list_of_param:
                    for mouse in mice:
                        if mouse == mouse_label:
                            if param == list_of_param[0]:
                                f.write(str(data1[phase][bi][mouse])
                                        + delimiter)
                            elif param == list_of_param[1]:
                                f.write(str(data2[phase][bi][mouse])
                                        + delimiter)
                        else:
                            continue
        f.write("\n")
    f.close()


