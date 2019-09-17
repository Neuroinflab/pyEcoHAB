from __future__ import division, print_function, absolute_import
import os
import numpy as np
from . import utility_functions as utils


def write_single_chamber(f, header, heads, address, mice, phases, time, data_stim ):
    longest = 0
    for j, h in enumerate(heads):
        f.write(h%address)
        f.write(header+'\n')
        for i, mouse in enumerate(mice):
            longest = 0
            for phase in phases:
                if len(data_stim[j][phase][mouse])> longest:
                    longest = len(data_stim[j][phase][mouse])
            lines = [mouse for i in range(longest)]

            for phase in phases:
                for k, t in enumerate(time[phase]):
                    if phase == phases[0]:
                        lines[k] += ';%3.2f'%(t/3600)
                    try:
                        lines[k] += ';'+str(data_stim[j][phase][mouse][k])
                    except IndexError:
                        print("Phase to short", phase)
                        return

            for line in lines:
                f.write(line+'\n')

def save_data_cvs(data, fname, path, which, headers,
                  target_dir="time_in_chambers"):
    new_path = os.path.join(path, target_dir)
    if not os.path.exists(new_path):
        print(new_path)
        os.makedirs(new_path)

    fname = os.path.join(new_path, fname)
    f = open(fname,'w')
    phases = data['phases']
    mice = data['mice']
    header = 'mouse;\"time [h]\"'
    for phase in phases:
        header += ';\"'+phase+"\""
  
    for stim in which.keys():
        heads = headers[stim]
        write_single_chamber(f,header, heads, which[stim], mice, phases, data['time'], data[stim])


def save_single_histograms(result, fname, mice, phase, main_directory, directory, prefix, additional_info=None):
    new_name = os.path.join(directory, 'data')
    directory = utils.check_directory(main_directory, new_name)
    if isinstance(additional_info, str):
        fname =  os.path.join(directory, '%s_%s_%s_%s.csv'% (fname, prefix, phase, additional_info))
        
    else:
        fname =  os.path.join(directory, '%s_%s_%s.csv'% (fname, prefix, phase))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
        return None
    for i, mouse in enumerate(mice):
        f.write(';'+mouse)
    f.write(';\n')
    for i, mouse in enumerate(mice):
        f.write(mouse + ';')
        for j, mouse in enumerate(mice):
            f.write(str(result[i, j])+';')
        f.write('\n')


def write_csv_rasters(mice, phases, output, directory, dirname, fname):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(directory, new_name)
    fname = os.path.join(directory, fname)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
   
    header = 'mouse pair'
    for phase in phases:
        header += ';' + phase
        
    header += '\n'
    f.write(header)
    new_output, pairs = utils.make_table_of_pairs(output, phases, mice)
    for i, pair in enumerate(pairs):
        f.write(pair)
        for j in range(len(phases)):
            f.write(';')
            f.write(str(new_output[i,j]))
        f.write('\n')
    f.close()

def write_csv_tables(results, phases, mice, main_directory, dirname, fname, prefix):
    new_name = os.path.join(dirname, 'data')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s.csv' % (fname, prefix))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return
    
    phase_pairs = [phases[2*i]+(len(mice)+2)*';'+ phases[2*i+1]+'\n' for i in range(len(phases)//2)]
    new_results = [(results[2*i], results[2*i+1])for i in range(len(phases)//2)]
    
    mice_header = ';'
    mice_len = 2
    if len(phases) == 1:
        mice_len = 1
    for i in range(mice_len):
        for mouse in mice:
            mice_header += mouse + ';'
        if not i:
            mice_header += ';;'
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
                f.write(mouse+';')
                for k, mouse in enumerate(mice):
                    f.write(str(nr[j, k])+';')
                if not l:
                    f.write(';')
            f.write('\n')
        f.write('\n')
    f.close()

def write_csv_alone(alone, phases, mice, main_directory, prefix, labels=["1", "2", "3", "4"],  header='Mice alone in chamber %s\n', fname='mouse_alone_%s.csv', directory="mouse_alone"):
    directory = utils.check_directory(main_directory, directory)
    fname =  os.path.join(directory, fname % prefix)
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to ', fname)
        return
    phases_header = 'Mouse;'
    for phase in phases:
        phases_header += phase + ';'
    phases_header += '\n'
    for i in range(1, alone.shape[0]+1):
        f.write(header % labels[i-1])
        f.write(phases_header)
        for j, mouse in enumerate(mice):
            f.write(mouse+';')
            for k, phase in enumerate(phases):
                f.write(str(alone[i-1, j, k])+';')
            f.write('\n')
    f.close()

def write_interpair_intervals(results, main_directory,
                              directory, fname, prefix,
                              additional_info=None):
    new_name = os.path.join(main_directory, 'data')
    directory = utils.check_directory(directory, new_name)
    if isinstance(additional_info, str):
        fname =  os.path.join(directory, '%s_%s_%s.csv'% (fname,
                                                          prefix,
                                                          additional_info))
    else:
        fname =  os.path.join(directory, '%s_%s.csv'% (fname,
                                                       prefix))
    try:
        f = open(fname, 'w')
    except IOError:
        print('Could not write to file', fname)
        return None
    print(fname)
    f.write("followed mouse, following mouse, intervals\n")
    keys = sorted(results.keys())
    for key in keys:
        mouse1, mouse2 = key.split('|')
        f.write("%s,%s," % (mouse1, mouse2))
        for interval in results[key]:
            f.write("%f," % interval)
        f.write("\n")

def save_visit_duration(results, time, phase, mice,
                        fname, main_directory,
                        directory, prefix,
                        add_info=None):
    new_dir = os.path.join(main_directory, directory)
    new_dir = utils.check_directory(new_dir, "data")
    for mouse in mice:
        if isinstance(add_info, str):
            new_name =  os.path.join(new_dir, '%s_%s_%s_%s_%s.csv'%(fname,
                                                                      mouse,
                                                                      phase,
                                                                      prefix,
                                                                      add_info))
        else:
            new_name =  os.path.join(new_dir, '%s_%s_%s_%s.csv'%(fname,
                                                                   mouse,
                                                                   phase,
                                                                   prefix))
        print(new_name)
        f = open(new_name, "w")
        for address in results.keys():
            f.write("Visit durations to chamber %s" % address)
            f.write("time, durations\n")
            for j, out in enumerate(results[address][mouse]):
                f.write("%2.2f" % time[j])
                for single in out:
                    f.write(",%2.2f" % single)
                f.write("\n")
        f.close()
