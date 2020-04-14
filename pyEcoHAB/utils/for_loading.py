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
