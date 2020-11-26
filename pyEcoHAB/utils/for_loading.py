from __future__ import division, print_function, absolute_import
import os
import time
import calendar
import sys
from collections import OrderedDict, Counter
import numpy as np
from pyEcoHAB.utility_functions import check_directory


try:
    basestring
except NameError:
    basestring = str

PAIRS = ["1 3", "1 4", "1 5", "1 6", "1 7", "2 4", "2 5", "2 6", "2 7", "2 8",
         "3 5", "3 6", "3 7", "3 8", "4 6", "4 7", "4 8", "5 7", "5 8", "6 8"]


def results_path(path, res_dir):
    return os.path.join(path, res_dir)


def make_prefix(path):
    """
    Read-in info.txt and make a prefix for all files for results.
    Parameters
    ----------
    path: str
    """
    key_list = [
        'genotype',
        "strain",
        'sex',
        'gender',
        'experimentator',
        'type of experiment',
        'date of experiment',
        "start date and hour",
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
        new_key = key.strip().lower()
        info = info.strip()
        info = info.replace(' ', '_')
        info_dict[new_key] = info
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
    date_in_sec = calendar.timegm(time.strptime(date, '%Y%m%d'))
    datenext = time.strftime('%Y%m%d', time.localtime(
        date_in_sec + 24*3600.))

    return hour, date, datenext


def print_human_time(tt):
    """convert seconds to date and time since epoch """
    st = time.gmtime(tt)
    return time.asctime(st)


def time_to_sec(tt):
    """Convert date and time to seconds since epoch"""
    try:
        more_than_sec, less_than_sec = tt.split('.')
    except ValueError:
        return calendar.timegm(time.strptime(tt + " UTC",
                                             '%Y%m%d %H:%M:%S %Z'))

    seconds = calendar.timegm(time.strptime(more_than_sec + " UTC",
                                            '%Y%m%d %H:%M:%S %Z'))
    return seconds + float(less_than_sec)/1000


def reformat_date_time(date, time):
    return "%s %s" % (date.replace('.', ''), time)


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


def read_single_file(dir_path, fname):
    """Reads in a single data file"""
    hour, date, datenext = parse_fname(fname)
    raw_data = []
    f = open(os.path.join(dir_path, fname), 'r')
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
            raise(IOError('Unknown data format in file %s' % f))
        raw_data += [line]
    f.close()
    return raw_data


def remove_one_antenna(data, antenna):
    """
    Remove animal tags registered by a specified antenna from 2D data array
    """
    where = []
    for i, a in enumerate(data["Antenna"]):
        if a == antenna:
            where.append(i)
    if len(where) == 0:
        return data
    return np.delete(data, where, axis=0)


def remove_antennas(data, antennas):
    """
    Remove all animal tag registrations by specified antennas from data.

    Args:
    data: a list of lista or an 2D array
        data dictionary (raw_data transformed by from_raw_data)
    antennas: int or list
        either a single antenna or a list of antennas to remove


    Returns:
       data dictionary (same as data)
    """
    if not isinstance(antennas, list):
        antennas = [antennas]
    new_data = data.copy()
    if isinstance(antennas, list):
        for antenna in antennas:
            new_data = remove_one_antenna(new_data, antenna)
    return new_data


def remove_ghost_tags(raw_data, legal_tags="ALL"):
    """
    Leave animal tag registrations that are trustworthy.

    If a list of correct tags is provided this method removes all registrations
    of incorrect  tags.
  
    Args:
    raw_data: a list of lists or an 2D array
        raw_data read by Loader._read_in_raw_data
    legal_tags: list
        animal tags to be kept in raw_data
        Default "ALL". Keep all tags.

    Returns:
       a list of lists or an 2D array (the same type as raw_data)
    """
    if legal_tags == "ALL":
        return raw_data

    new_data = []
    if isinstance(legal_tags, basestring):
        legal_tags = [legal_tags]
    
    for d in raw_data:
        mouse = d[4]
        if mouse in legal_tags:
            new_data.append(d)
    return new_data[:]


def check_antenna_presence(raw_data, setup_config, max_break):
    if not len(raw_data):
        raise Exception("Empty dataset")
    all_times = raw_data['Time']
    t_start = raw_data['Time'][0]
    breaks = {}
    t_end = raw_data['Time'][-1]
    for antenna in setup_config.all_antennas:
        antenna_idx = []
        for i, a in enumerate(raw_data['Antenna']):
            if antenna == a:
                antenna_idx.append(i)
        times = all_times[antenna_idx]
        breaks[antenna] = []
        if len(times):
            if times[0] - t_start > max_break:
                breaks[antenna].append([t_start, np.round(times[0])])
            intervals = times[1:] - times[0:-1]
            where_breaks = np.where(intervals > max_break)[0]
            if len(where_breaks):
                for i in where_breaks:
                    breaks[antenna].append([np.round(times[i]),
                                            np.round(times[i+1])])
            if t_end - times[-1] > max_break:
                breaks[antenna].append([np.round(times[-1]), t_end])
        else:
            breaks[antenna].append([np.round(t_start),
                                    raw_data['Time'][-1]])
    return breaks


def antenna_mismatch(raw_data, setup_config):
    if not len(raw_data):
        raise Exception("Empty dataset")
    pairs = setup_config.mismatched_pairs
    all_times = raw_data['Time']
    mice = set(raw_data['Tag'])
    mismatches = OrderedDict()
    for pair in pairs:
        mismatches[pair] = 0
    for mouse in mice:
        mouse_idx = np.where(np.array(raw_data['Tag']) == mouse)[0]
        times = all_times[mouse_idx]
        ant = raw_data['Antenna'][mouse_idx]
        for i, a in enumerate(ant[:-1]):
            key = "%s %s" % (min(a, ant[i+1]), max(a, ant[i+1]))
            if key in pairs:
                mismatches[key] += 1
    return mismatches


def total_mismatches(mismatches):
    all_antennas = set()
    for pair in mismatches.keys():
        a1, a2 = pair.split()
        all_antennas.add(a1)
        all_antennas.add(a2)
    all_antennas = sorted(all_antennas)
    out = {}
    for antenna in all_antennas:
        out[antenna] = 0
        for pair in mismatches:
            if antenna in pair:
                out[antenna] += mismatches[pair]
    return out


def skipped_registrations(raw_data, setup_config):
    if not len(raw_data):
        raise Exception("Empty dataset")
    one_skipped = setup_config.skipped_one()
    skipped_two = setup_config.skipped_two()
    skipped_more = setup_config.skipped_more()

    mice = set(raw_data['Tag'])
    mismatches = OrderedDict([("skipped one", 0), ("skipped two", 0), ("skipped more", 0)])
    for mouse in mice:
        mouse_idx = np.where(np.array(raw_data['Tag']) == mouse)[0]
        ant = raw_data['Antenna'][mouse_idx]
        for i, a in enumerate(ant[:-1]):
            key = "%s %s" % (a, ant[i+1])
            if key in one_skipped:
                mismatches["skipped one"] += 1
            elif key in skipped_more:
                mismatches["skipped more"] += 1
            elif key in skipped_two:
                mismatches["skipped two"] += 1
    return mismatches

def save_skipped_registrations(skipped, tot_registrations, res_dir,
                               fname="skipped_registrations.csv",
                               header=u"type, count, percentage\n"):
    out_f1 = header
    for key in skipped.keys():
        try:
            exact_mis = np.round(100*skipped[key]/tot_registrations)
        except ZeroDivisionError:
            exact_mis = 1
        out_f1 += u"%s, %d, %3.2f per 100\n" % (key, skipped[key],
                                                 exact_mis)
    new_path = check_directory(res_dir, "diagnostics")
    fpath1 = os.path.join(new_path, fname)
    f1 = open(fpath1, "w")
    f1.write(out_f1)
    f1.close()
    return out_f1


def save_mismatches(mismatches, tot_registrations, res_dir,
                    fname="antenna_mismatches.csv",
                    header=u"antenna pair,  count, percentage\n"):
    out_f1 = header
    for pair in mismatches.keys():
        a1, a2 = pair.split(" ")
        if isinstance(tot_registrations, dict):
            try:
                exact_mis = np.round(100*mismatches[pair]/tot_registrations[pair])
            except ZeroDivisionError:
                exact_mis = 0
        else:
            exact_mis = np.round(100*mismatches[pair]/tot_registrations)
        out_f1 += u"%s, %d, %3.2f per 100\n" % (pair, mismatches[pair],
                                                exact_mis)
    new_path = check_directory(res_dir, "diagnostics")
    fpath1 = os.path.join(new_path, fname)
    f1 = open(fpath1, "w")
    f1.write(out_f1)
    f1.close()
    return out_f1


def save_total_mismatches(tot_mismatches, counters, res_dir):
    out_f1 = u"antenna, incorrect transitions count, percentage of antenna recordings\n"
    for a1 in tot_mismatches.keys():
        if counters[a1]:
            exact_mis = 100*tot_mismatches[a1]/counters[a1]
        else:
            exact_mis = 0
        out_f1 += u"%s, %d, %3.2f per 100\n" % (a1, tot_mismatches[a1],
                                                 exact_mis)
    new_path = check_directory(res_dir, "diagnostics")
    fpath1 = os.path.join(new_path, "incorrect_antenna_transitions.csv")
    f1 = open(fpath1, "w")
    f1.write(out_f1)
    f1.close()
    return out_f1


def save_antenna_breaks(antenna_breaks, res_dir):
    out_f2 = u'Breaks in registrations on antennas:\n'
    for antenna in antenna_breaks:
        out_f2 += u"%s:\n" % antenna
        for breaks in antenna_breaks[antenna]:
            out_f2 += u"%s %s, %4.2f h\n" % (print_human_time(breaks[0]),
                                             print_human_time(breaks[1]),
                                             (breaks[1] - breaks[0])/3600)

    new_path = check_directory(res_dir, "diagnostics")
    fpath2 = os.path.join(new_path, "breaks_in_registrations.csv")
    f2 = open(fpath2, "w")
    f2.write(out_f2)
    f2.close()
    return out_f2

def run_diagnostics(raw_data, max_break, res_dir, setup_config):
    """
    Calculate parameters showing fidelity of obtained antenna
    registrations and overall experimental errors during the experiment.
    These parameters will be saved in "diagnostics" directory.

    Args:
    raw_data:  dictionary of sequences
       data read in from data files
    max_break: float
       maximum brea k in single antenna registrations (in sec)
    res_dir: string
       path to results directory
    setup_config: SetupConfig or ExperimentalSetupConfig
      object describing geometry of the (modular) experimental setup
    
    returns:
      string_1: text
        text showing count and percentage of pairs of mismatched antenna
        registrations
      string_2: text
        text describing periods when any antenna was silent for more
        than max_break
      string_3: text
        text showing count and percentage of mismatched antenna registrations
        per single antenna
      string_4: text
        text showing count and percentage of cases, when for 2 consecutive
        antenna registration one antenna, two antennas and more than two
        antennas were skipped
      string_5: text
        text showing count and percentage of cases, when two entrance antennas
        to the same tunnel registered an animal simultaneously
    """
    mismatches = antenna_mismatch(raw_data, setup_config)
    string_1 = save_mismatches(mismatches, len(raw_data["Antenna"]),
                               res_dir)
    antenna_breaks = check_antenna_presence(raw_data, setup_config, max_break)
    string_2 = save_antenna_breaks(antenna_breaks, res_dir)
    tot_mismatches = total_mismatches(mismatches)
    counters = Counter(raw_data["Antenna"])
    string_3 = save_total_mismatches(tot_mismatches, counters, res_dir)
    skip = skipped_registrations(raw_data, setup_config)
    string_4 = save_skipped_registrations(skip, len(raw_data["Tag"]), res_dir)
    count, total_count = incorrect_tunnel_registrations(raw_data, setup_config)
    header = u"tunnel, count, percentage of all passings through the tunnel\n"
    string_5 = save_mismatches(count, total_count, res_dir,
                               fname="incorrect_tunnel_registrations.csv",
                               header=header)
    return string_1, string_2, string_3, string_4, string_5


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


def incorrect_tunnel_registrations(raw_data, setup_config):
    count = OrderedDict()
    directions = setup_config.directions
    total_count = {}
    for direction in sorted(directions):
        a1, a2 = direction.split(" ")
        key = "%s %s" % (min(a1, a2), max(a1, a2))
        count[key] = 0
        total_count[key] = 0
    mice = set(raw_data['Tag'])

    for mouse in mice:
        mouse_idx = np.where(np.array(raw_data['Tag']) == mouse)[0]
        times1 = raw_data["Time"][mouse_idx]
        antennas1 = raw_data['Antenna'][mouse_idx]
        durations1 = raw_data["Duration"][mouse_idx]
        out_count, out_tot_count = incorrect_tunnel_single_mouse(count.keys(),
                                                                 antennas1,
                                                                 times1,
                                                                 durations1)
        for key in count:
            count[key] += out_count[key]
            total_count[key] += out_tot_count[key]
    return count, total_count

def transform_raw(row):
    return (int(row[0]), time_to_sec(row[1]),
            row[2], int(row[3]), row[4])

def from_raw_data(raw_data):
    """
    Transform raw data, which is in the form of a double list
    to a 2D structured array with registrations of animal tags.
    """
    new_data = []
    for row in raw_data:
        new_data.append(transform_raw(row))
    data_type = [("Id", int),
                 ("Time", float),
                 ("Antenna", "U15"),
                 ("Duration", int),
                 ("Tag", "U15")]
    return np.array(new_data, dtype=data_type)


def transform_visits(data):
    data_type = [("Address", "U30"),
                 ("Tag", "U15"),
                 ("AbsStartTimecode", float),
                 ("AbsEndTimecode", float),
                 ("VisitDuration", float),
                 ("ValidVisitSolution", bool)]
    return np.array(data, dtype=data_type)


def rename_antennas(name, dataset):
    new_data = dataset.copy()
    for row in new_data:
        row["Antenna"] = "%s_%s" % (row["Antenna"], name)
    return new_data


def append_data_sources(data_sets):
    new_data = np.concatenate(data_sets)
    new_data.sort(order="Time")
    return new_data


class NamedDict(dict):
    """Creates a python dict with a name and attribute access of keys.

    Usage: mydict = NamedDict(name,**kwargs)
    where **kwargs are used to create dictionary key/value pairs.
    e.g.: params = NamedDict('modelParams',x=15,y=0)

    dict keys can be accessed and written as keys or attributes:
        myNamedDict['k'] is equivalent to myNamedDict.k, and
        myNamedDict['k'] = newvalue is equivalent to myNamedDict.k=newvalue.

    New entries/attributes can be created:
        myNamedDict.newkey = newvalue OR myNamedDict['newkey']= newvalue.

    Note: Dict ignores attributes beginning with underscore, so
    myNamedDict.__name__ returns the NamedDict name, but there is no dict key
    == "__name__"

    Note: all dict keys must be valid attribute names: that is, strings with
    first character in a-z/A-Z. This could be changed to allow all valid python
    dict keys as keys, but these keys would not have attribute access.

    """

    def __init__(self, name, **kwargs):
        super(NamedDict, self).__init__(**kwargs)
        self.__dict__ = dict(**kwargs)
        self.__name__ = name

    def __repr__(self):
        items = ('{}={}'.format(k, v) for (k, v) in self.items())
        length = len(self.__name__) + 1
        sep = ',\n' + ' '*length
        return '{}({})'.format(self.__name__, sep.join(items))

    def __setitem__(self, k, v):
        super(NamedDict, self).__setitem__(k, v)
        setattr(self, k, v)

    def __getattribute__(self, k):
        # attributes have higher priority
        try:
            return super(NamedDict, self).__getattribute__(k)
        except AttributeError:
            return super(NamedDict, self).__getitem__(k)

    def __setattr__(self, k, v):
        super(NamedDict, self).__setattr__(k, v)
        if not k.startswith('_'):
            super(NamedDict, self).__setitem__(k, v)

    def __dir__(self):
        dirlist = super(NamedDict, self).__dir__()
        return dirlist
