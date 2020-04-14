from __future__ import division, print_function, absolute_import
import os
import time
import sys
from collections import OrderedDict
import numpy as np

def results_path(path):
    return os.path.join(path, 'Results')


def make_prefix(path):
    """
    Read-in info.txt and make a prefix for all files for results.
    Parameters
    ----------
    path: str
    """
    key_list = [
        'genotype',
        "Strain",
        'sex',
        'gender', 
        'Experimentator',
        'type of experiment',
        "Type of Experiment",
        'date of experiment',
        "Start date and hour",
        'social odor',
        'no social odor',
    ]

    fname = os.path.join(path, 'info.txt')
    try:
        f = open(fname)
    except IOError:
        return ''

    info_dict = {}
    for line in f:
        try:
            key, info = line.split(':')
        except ValueError:
            continue
        info = info.strip()
        info = info.replace(' ', '_')
        info_dict[key] = info
    prefix = ''
    for key in key_list:
        if key not in info_dict:
            continue
        if key == 'social odor' or key == 'non-social odor':
            if info_dict[key] == 'none' or info_dict[key] == 'None':
                key = key.replace(' ', '_')
                prefix += '_no_' + key.replace(' ', '_') + '_'
            else:
                prefix += key.replace(' ', '_') + '_' + info_dict[key] + '_'
        else:
            prefix += info_dict[key] + '_'
    return prefix


def parse_fname(fname):
    """"Extracts time and date from data filename"""
    try:
        date, hour = fname.split("_")
    except ValueError:
        parts = fname.split("_")
        if len(parts) == 3:
            date, hour = parts[:2]
        else:
            print("Unnkown filename format %s.")
            raise
    hour = hour.split(".")[0]
    date_in_sec = time.mktime(time.strptime(date, '%Y%m%d'))
    datenext = time.strftime('%Y%m%d', time.localtime(
        date_in_sec + 24*3600.))

    return hour, date, datenext

def print_human_time(tt):
    """convert seconds to date and time since epoch """
    st = time.localtime(tt)
    return time.asctime(st)


def time_to_sec(tt):
    """Convert date and time to seconds since epoch"""
    try:
        more_than_sec, less_than_sec = tt.split('.')
    except ValueError:
        return time.mktime(time.strptime(tt,
                                             '%Y%m%d %H:%M:%S'))

    seconds = time.mktime(time.strptime(more_than_sec,
                                        '%Y%m%d %H:%M:%S'))
    return seconds + float(less_than_sec)/1000.


def reformat_date_time(date, time):
    return "%s %s" %(date.replace('.',''), time)


def process_line_more_elements(elements):
    """remove dot from 2nd column of new data files"""
    date_time = reformat_date_time(elements[1], elements[2])
    return [elements[0], date_time] + elements[3:]


def process_line_5_elements(elements, date):
    """Add date to data (old data files)"""
    elements[1] = ' '.join([date, elements[1]])
    return elements


def get_filenames(path):
    try:
        f_list = os.listdir(path)
    except FileNotFoundError:
        return []
    out = []
    for f_name in f_list:
        if f_name.endswith("0000.txt"):
            out.append(f_name)
        else:
            split = f_name.split("_")
            if len(split) < 3:
                continue
            if split[-1].endswith(".txt") and split[1].endswith("0000"):
                out.append(f_name)
    return out


def _read_single_file(dir_path, fname):
    """Reads in a single data file"""
    hour, date, datenext = parse_fname(fname)
    raw_data = []
    f = open(os.path.join(dir_path, fname),'r')
    for line in f:
        elements = line.split()
        if len(elements) == 5:
            if hour == '23' and elements[1][:2] == '00':
                line = process_line_5_elements(elements, datenext)
            else:
                line = process_line_5_elements(elements, date)
        elif len(elements) > 5:
            line = process_line_more_elements(elements)
        else:
            raise(IOError('Unknown data format in file %s' %f))
        raw_data += [line]
    return raw_data


def remove_antennas(data, antennas):
    if isinstance(antennas, int):
        antennas = [antennas]
    if isinstance(antennas, list):
        keys = data.keys()
        new_data = {key:[] for key in keys}
        for i, antenna in enumerate(data['Antenna']):
            if antenna not in antennas:
                for key in keys:
                    new_data[key].append(data[key][i])
            else:
                print("Removing record",
                      [data[key][i] for key in keys])
        return new_data
    return data


def remove_ghost_tags(raw_data, how_many_appearances,
                       how_many_days, tags=[]):
    """
    Remove animal tag registrations that are untrustworthy.

    This method removes all animal tag registration, when the Eco-HAB
    system registered the animal tag, if:
    1. less times than how_many days,
    2. during less than how_many_days of the experiment,
    3. the tag was provided in tags.

    Args:
    raw_data: a list of list or an 2D array
        raw data read by Loader._read_in_raw_data
    how_many_appearances: int
        minimum number of tag registration
    how_many_days: float
        minimum number of days, on which the animal tag was registred
    tags: list
        animal tags to be removed from raw_data
    """
    new_data = []
    ghost_mice = []
    counters = {}
    dates = {}
    if len(tags):
        for tag in tags:
            ghost_mice.append(tag)
    for d in raw_data:
        mouse = d[4]
        if mouse not in counters:
            counters[mouse] = 0
        if mouse not in dates:
            dates[mouse] = set()
        counters[mouse] += 1
        dates[mouse].add(d[1].split()[0])

    for mouse in counters:
        if counters[mouse] < how_many_appearances or len(dates[mouse]) <= how_many_days:
            if mouse not in ghost_mice:
                ghost_mice.append(mouse)
    for d in raw_data:
        mouse = d[4]
        if mouse not in ghost_mice:
            new_data.append(d)

    return new_data[:]



def check_antenna_presence(raw_data, max_break):
    t_start = raw_data['Time'][0]
    all_times = np.array(raw_data['Time'])
    breaks = {}

    for antenna in range(1, 9):
        antenna_idx = np.where(np.array(raw_data['Antenna']) == antenna)[0]
        times = all_times[antenna_idx]
        breaks[antenna] = []
        if len(times):
            if times[0] - t_start > max_break:
                breaks[antenna].append([0, np.round(times[0])])
            intervals = times[1:] - times[0:-1]
            where_breaks = np.where(intervals > max_break)[0]
            if len(where_breaks):
                for i in where_breaks:
                    breaks[antenna].append([np.round(times[i]), np.round(times[i+1])])
        else:
            breaks[antenna].append([np.round(t_start),
                                    raw_data['Time'][-1]])
    return breaks

def antenna_mismatch(raw_data):
    t_start = raw_data['Time'][0]
    all_times = np.array(raw_data['Time'])
    weird_transit = [[], []]
    mice = set(raw_data['Tag'])
    for mouse in mice:
        mouse_idx = np.where(np.array(raw_data['Tag']) == mouse)[0]
        times = all_times[mouse_idx]
        antennas = np.array(raw_data['Antenna'])[mouse_idx]
        for i, a in enumerate(antennas[:-1]):
            if abs(a - antennas[i+1]) not in [0,1,7]:
                weird_transit[0].append(times[i])
                if a < antennas[i+1]:
                    weird_transit[1].append("\t    %d,\t\t\t %d" % (a,
                                                                    antennas[i+1]))
                else:
                    weird_transit[1].append("\t    %d,\t\t\t %d" % (antennas[i+1],
                                                                    a))
    pairs = list(set(weird_transit[1]))

    mismatches = OrderedDict()
    print("Mismatched antenna readings\n")
    print("First reading, consecutive reading,  count, percentage\n")
    for pair in pairs:
        mismatches[pair] = weird_transit[1].count(pair)
        print("%s,\t%d, %3.2f per 100"% (pair, mismatches[pair],
                                         np.round(100*mismatches[pair]/len(raw_data['Antenna']))))
    return weird_transit


def run_diagnostics(raw_data, max_break):
    antenna_breaks = check_antenna_presence(raw_data, max_break)
    if antenna_breaks:
        print('No registrations on antennas:')
        for antenna in antenna_breaks:
            print(antenna, ':')
            for breaks in antenna_breaks[antenna]:
                print(print_human_time(breaks[0]),
                      print_human_time(breaks[1]))
                print((breaks[1] - breaks[0])/3600, 'h')
    antenna_mismatch(raw_data)


def _from_raw_data(raw_data, antenna_positions, remove_antennas=[]):
    data = {}
    data['Id'] = [d[0] for d in raw_data]
    data['Time'] = [ufl.time_to_sec(d[1]) for d in raw_data]
    data['Antenna'] = [antenna_positions[d[2]] for d in raw_data]
    data['Duration'] = [d[3] for d in raw_data]
    data['Tag'] = [d[4] for d in raw_data]

    return remove_antennas(data, remove_antennas)
