from __future__ import division, print_function, absolute_import
import os
import numpy as np
from . import utility_functions as utils


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
            if len(data_stim[phase][mouse])> longest:
                longest = len(data_stim[phase][mouse])
        lines = [mouse for i in range(longest)]
        for phase in phases:
            for k, t in enumerate(time[phase]):
                if phase == phases[0]:
                    lines[k] += '%s%3.2f'%(delimiter, t/3600)
                try:
                    if floats:
                        lines[k] += "%s%7.3f" % (delimiter,
                                                 data_stim[phase][mouse][k])
                    else:
                        lines[k] += delimiter + str(data_stim[phase][mouse][k])
                except IndexError:
                    print("Phase too short", phase)
                    pass

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
            f.write("%s %s\n" %(h, stim))
            write_single_chamber(f, head, phases, mice, bin_labels,
                                 data[stim][j], delimiter, floats=j)
        


def write_binned_data(data_stim, fname, mice, bin_labels, phase,
                      path, target_dir, prefix, additional_info="",
                      delimiter=";"):
    new_path = os.path.join(path, target_dir, "data")
    fname =  os.path.join(new_path, '%s_%s_%s_%s.csv'% (fname,
                                                        prefix,
                                                        phase,
                                                        additional_info))
    if not os.path.exists(new_path):
        print(new_path)
        os.makedirs(new_path)
    f = open(fname, "w")
    header = make_header_for_activity(mice, delimiter)

    f.write(header+'\n')

    assert len(bin_labels) == len(data_stim.keys())
    for mouse1 in mice:
        lines = [mouse1 for l in range(len(bin_labels))]
        for j, mouse2 in enumerate(mice):
            for k, t in enumerate(bin_labels):
                if not j:
                    lines[k] += '%s%3.2f'%(delimiter, t/3600)
                lines[k] += delimiter + str(data_stim[t][mouse1][mouse2])
        for line in lines:
            f.write(line + '\n')
    f.close()


def save_single_histograms(result, fname, mice, phase, main_directory,
                           directory, prefix, additional_info="",
                           delimiter=";"):
    new_name = os.path.join(directory, 'data')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s_%s_%s.csv'% (fname,
                                                         prefix,
                                                         phase,
                                                         additional_info))
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
                      dirname, fname, symmetric=True,
                      reverse_order=False, delimiter=";"):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(directory, new_name)
    fname = os.path.join(directory, fname)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
   
    header = 'mouse pair'
    for phase in phases:
        header += delimiter + str(phase)+ " h"
        
    header += '\n'
    f.write(header)
    if symmetric:
        new_output, pairs = utils.make_table_of_pairs(output, phases, mice)
    else:
        new_output, pairs = utils.make_table_of_all_pairs(output, phases,
                                                          mice, reverse_order)
    for i, pair in enumerate(pairs):
        f.write(pair)
        for j in range(len(phases)):
            f.write(delimiter + str(new_output[i,j]))
        f.write('\n')
    f.close()

def write_csv_tables(results, phases, mice, main_directory,
                     dirname, fname, prefix,
                     delimiter=";"):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s.csv' % (fname, prefix))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return
    
    phase_pairs = [phases[2*i]+(len(mice)+2)*delimiter\
                   + phases[2*i+1]+'\n' for i in range(len(phases)//2)]
    new_results = [(results[2*i],
                    results[2*i+1])for i in range(len(phases)//2)]
    
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
                    header='Mice alone in chamber %s\n',
                    fname='mouse_alone_%s.csv',
                    directory="solitude",
                    delimiter=";"):
    directory = utils.check_directory(main_directory, directory)
    fname =  os.path.join(directory, fname % prefix)
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
                              delimiter=";"):
    new_name = os.path.join(main_directory, 'data')
    directory = utils.check_directory(directory, new_name)
    fname =  os.path.join(directory, '%s_%s_%s.csv'% (fname,
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

def save_visit_duration(results, time, phase, mice,
                        fname, main_directory, directory,
                        prefix, add_info="", delimiter=";"):
    new_dir = os.path.join(main_directory, directory)
    new_dir = utils.check_directory(new_dir, "data")
    for mouse in mice:
        new_name =  os.path.join(new_dir, '%s_%s_%s_%s_%s.csv'%(fname,
                                                                mouse,
                                                                phase,
                                                                prefix,
                                                                add_info))
        print(new_name)
        f = open(new_name, "w")
        for address in results.keys():
            f.write("Visit durations to chamber %s " % address)
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
                            add_info="", delimiter=";"):
    new_dir = os.path.join(main_directory, directory)
    new_dir = utils.check_directory(new_dir, "data")
    new_name =  os.path.join(new_dir, '%s_%s_%s_%s.csv'%(fname,
                                                         phase.replace(' ',
                                                                       '_'),
                                                         prefix,
                                                         add_info))

    f = open(new_name, "w")
    for i, mouse1 in enumerate(mice_list):
        for j, mouse2 in enumerate(mice_list):
            if mouse1 != mouse2:
                key = "%s|%s" % (mouse1, mouse2)
                f.write(key)
                for value in results[i, j]:
                    f.write(delimiter + str(value))
                f.write("\n")
    f.close()
