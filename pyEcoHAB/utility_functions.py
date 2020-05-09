from __future__ import division, print_function, absolute_import
import os
import time
import sys
from collections import OrderedDict
import numpy as np


same_pipe = {1: [1, 2],
             2: [1, 2],
             3: [3, 4],
             4: [3, 4],
             5: [5, 6],
             6: [5, 6],
             7: [7, 8],
             8: [7, 8]}

opposite_pipe = {1: [5, 6],
                 2: [5, 6],
                 3: [7, 8],
                 4: [7, 8],
                 5: [1, 2],
                 6: [1, 2],
                 7: [3, 4],
                 8: [3, 4]}

address = {1: "A", #4
           2: "B", #1,
           3: "B", #1,
           4: "C", #2,
           5: "C", #2,
           6: "D", #3,
           7: "D", #3,
           8: "A", #4
}

address_not_adjacent = {1: "B", #1,
                        2: "A", #4,
                        3: "C", #2,
                        4: "B", #1,
                        5: "D", #3,
                        6: "C", #2,
                        7: "A", #4,
                        8: "D", #3
}
# Surrounding: difference between antennas only 2 or 6 -- skipped one antenna
surrounding = {(1, 3): "B", #1,
               (1, 7): "A", #4,
               (2, 4): "B", #1,
               (2, 8): "A", #4,
               (3, 5): "C", #2,
               (4, 6): "C", #2,
               (5, 7): "D", #3,
               (6, 8): "D", #3
}

def check_directory(directory, subdirectory=None):
    if subdirectory:
        new_path = os.path.join(directory, subdirectory)
    else:
        new_path = directory
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    return new_path


def results_path(path):
    return os.path.join(path, 'Results')


def make_figure(title):
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, aspect='equal')
    fig.suptitle('%s'%(title), fontsize=14, fontweight='bold')
    return fig, ax


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

def list_of_pairs(mice):
    pair_labels = []
    for j, mouse in enumerate(mice):
        for k in range(j+1, len(mice)):
            pair_labels.append(mice[j]+'|'+mice[k])
    return pair_labels

def all_pairs(mice, reverse_order=False):
    pair_labels = []
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            if reverse_order is True:
                pair_labels.append('%s|%s' % (mouse2, mouse1))
            else:
                pair_labels.append('%s|%s' % (mouse1, mouse2))
    return pair_labels


def make_table_of_pairs(FAM, phases, mice):
    new_shape = (len(mice)*(len(mice)-1)//2, len(phases))
    output = np.zeros(new_shape)
    pair_labels = list_of_pairs(mice)
    for i in range(len(phases)):
        l = 0
        for j in range(len(mice)):
            for k in range(j + 1, len(mice)):
                output[l, i] = FAM[i, j, k]
                l += 1

    return output, pair_labels

def make_table_of_all_pairs(FAM, phases, mice, reverse_order=False):
    new_shape = (len(mice)*(len(mice)-1), len(phases))
    output = np.zeros(new_shape)
    pair_labels = all_pairs(mice,
                            reverse_order=reverse_order)
    for i, phase in enumerate(phases):
        l = 0
        for j in range(len(mice)):
            for k in range(len(mice)):
                if j == k:
                    continue
                if reverse_order is True:
                    output[l, i] = FAM[i, k, j]
                else:
                    output[l, i] = FAM[i, j, k]
                l += 1
    return output, pair_labels


def filter_dark(phases):
    out = []
    for phase in phases:
        if 'dark' in phase:
            out.append(phase)
        elif 'DARK' in phase:
            out.append(phase)
        elif 'Dark' in phase:
            out.append(phase)
    return out


def filter_light(phases):
    out = []
    for phase in phases:
        if 'light' in phase:
            out.append(phase)
        elif 'LIGHT' in phase:
            out.append(phase)
        elif 'Light' in phase:
            out.append(phase)
    return out


def filter_dark_light(phases):
    out = []
    for phase in phases:
        if 'light' in phase:
            out.append(phase)
        elif 'LIGHT' in phase:
            out.append(phase)
        elif 'Light' in phase:
            out.append(phase)
        elif 'dark' in phase:
            out.append(phase)
        elif 'DARK' in phase:
            out.append(phase)
        elif 'Dark' in phase:
            out.append(phase)

    return out


def get_mice(mouse_list, remove_mouse):
    if remove_mouse is None:
        return mouse_list

    if isinstance(remove_mouse, str):
        remove_mouse = [remove_mouse]

    if isinstance(remove_mouse, list):
        for mouse in mouse_list:
            if mouse in remove_mouse:
                mouse_list.remove(mouse)
    return mouse_list


def add_info_mice_filename(remove_mouse):
    if remove_mouse is None:
        return ''
    if isinstance(remove_mouse, str):
        remove_mouse = [remove_mouse]
    add_info_mice = ''
    if isinstance(remove_mouse, list):
        add_info_mice = 'remove'
        for mouse in remove_mouse:
            add_info_mice += '_' + mouse 
    return add_info_mice


def get_idx_pre(t0, times):
    idxs = np.where(np.array(times) < t0)[0]
    if len(idxs):
        return idxs[-1]
    return None

def get_idx_between(t0, t1, times):
    return  np.where((np.array(times) >= t0) & (np.array(times) < t1))[0]


def get_idx_post(t1, times):
    idxs = np.where(np.array(times) > t1)[0]
    if len(idxs):
        return idxs[0]
    return None


def in_tube(antenna, next_antenna):
    if antenna % 2:
        if next_antenna  == antenna  + 1:
            return True
    else:
        if next_antenna == antenna - 1:
            return True
    return False


def in_chamber(antenna, next_antenna):
    antenna = antenna % 8
    next_antenna = next_antenna % 8
    if antenna % 2:
        if next_antenna  == antenna - 1:
            return True
    else:
        if next_antenna == antenna + 1:
            return True
    return False


def change_state(antennas):
    return np.where(abs(np.array(antennas[:-1]) - np.array(antennas[1:])) !=0)[0]

def mouse_going_forward(antennas):
    assert len(antennas) > 2
    first_antenna, last_antenna = antennas[0], antennas[-1]
    if first_antenna % 2 and last_antenna % 2:
        return True
    if not first_antenna % 2 and not last_antenna % 2:
        return True
    return False
    
def mouse_backing_off(antennas):
    first_antenna, last_antenna = antennas[0], antennas[-1]
    how_many = len(set(antennas))
    if how_many == 3:
        if first_antenna % 2 and not last_antenna % 2:
            return True
        if not first_antenna % 2 and last_antenna % 2:
            return True
        return False

    if first_antenna == last_antenna and len(antennas) > 2:
        return True
    return False

def skipped_antennas(antennas):
    change = abs(np.array(antennas[:-1]) - np.array(antennas[1:]))
    if len(np.intersect1d(np.where(change>=2)[0], np.where(change<=6)[0])):
        return True
    return False

    
def change_seven_to_one(change):
    if isinstance(change, list):
        change = np.array(change)
    seven = np.where(change == 7)[0]
    if len(seven):
        change[seven] = -1
    minus_seven = np.where(change == -7)[0]
    if len(minus_seven):
        change[minus_seven] = 1
    return change


def mouse_going_clockwise(antennas):
    change = np.array(antennas[:-1]) - np.array(antennas[1:])
    change = change_seven_to_one(change)
    if sum(change) < 0:
        return True
    return False


def mouse_going_counterclockwise(antennas):
    change = np.array(antennas[:-1]) - np.array(antennas[1:])
    change = change_seven_to_one(change)
    if sum(change) > 0:
        return True
    return False


def get_times_antennas(ehd, mouse, t_1, t_2):
    if t_1 == 0 and t_2 == -1:
        return ehd.gettimes(mouse), ehd.getantennas(mouse)
    ehd.mask_data(t_1, t_2)
    antennas, times = ehd.getantennas(mouse), ehd.gettimes(mouse)
    ehd.unmask_data()
    return times, antennas


def get_states_and_readouts(antennas, times, t1, t2):
    before = get_idx_pre(t1, times)
    between = get_idx_between(t1, t2, times)
    after = get_idx_post(t2, times)
    states = []
    readouts = []
    if before is not None:
        states.append(antennas[before])
        readouts.append(times[before])
    for idx in between:
        states.append(antennas[idx])
        readouts.append(times[idx])
    assert(len(states) == len(readouts))
    return states, readouts


def get_more_states(antennas, times, midx,
                    mouse_attention_span,
                    how_many_antennas):
    #save first antenna
    states = [antennas[midx]]
    readouts = [times[midx]]
    midx += 1
    idx = 1
    while True:
        if midx >= len(antennas):
            break
        #read in next antenna
        new_antenna = antennas[midx]
        new_readout = times[midx]
        #if pause too long break
        if new_readout > readouts[idx - 1] + mouse_attention_span:
            break

        states.append(new_antenna)
        readouts.append(new_readout)

        idx += 1
        #if more than 2 antennas, break
        if len(set(states)) == how_many_antennas:
            # go back to the last readout of the opposite antenna not to loose it
            break
        midx += 1
        
    return states, readouts, midx

def get_antennas(idxs, antennas):
    antenna_slice = []
    for new_idx in idxs:
        antenna_slice.append(antennas[new_idx])
    return antenna_slice


def get_timestamp(t_start, t_end, dt):
    return int(round((t_end - t_start)/dt))

def get_key_for_frequencies(antenna, next_antenna):
    if antenna % 2 and next_antenna == antenna + 1:
        return "%d%d" % (antenna, next_antenna)
    elif next_antenna % 2 and antenna == next_antenna + 1:
        return "%d%d" % (antenna, next_antenna)


def interval_overlap(int1, int2):
    """Return overlap between two intervals."""
    if int1[1] < int1[0]:
        int1 = int1[::-1]
    if int2[1] < int2[0]:
        int2 = int2[::-1]
    ints = sorted([int1, int2], key=lambda x: x[0])
    if ints[1][0] > ints[0][1]:
        return 0

    return min(ints[0][1], ints[1][1]) - ints[1][0]


def get_duration(starts, ends):
    return sum([abs(ends[i] - start) for i, start in enumerate(starts)])


def get_interval_durations(ints):
    return [x[1] - x[0] for x in ints]


def get_interval_durations_2_lists(starts, ends):
    return [abs(ends[i] - start) for i, start in enumerate(starts)]


def calculate_total_duration(intervals):
    return sum(get_interval_durations(intervals))


def get_intervals(data, address):
    return [[s, e] for a, s, e in data if a == address]


def intervals2lists(data, address):
    out = get_intervals(data, address)
    intervals_table = [[], []]
    for st, en in out:
        intervals_table[0] += [st]
        intervals_table[1] += [en]
    return intervals_table


def get_indices(t_start, t_end, starts, ends):
    idx_start = get_idx_between(t_start, t_end, starts).tolist()
    idx_end = get_idx_between(t_start, t_end, ends).tolist()
    return sorted(list(set(idx_start +  idx_end)))


def get_ehs_data_with_margin(ehs, mouse, t_start, t_end,
                             margin=12*3600):
    if t_start == 0 and t_end == -1:
        return ehs.getaddresses(mouse),\
            ehs.getstarttimes(mouse),\
            ehs.getendtimes(mouse)

    ehs.mask_data(t_start - margin, t_end +  margin)
    adresses = ehs.getaddresses(mouse)
    starts = ehs.getstarttimes(mouse)
    ends = ehs.getendtimes(mouse)
    ehs.unmask_data()
    return adresses, starts, ends


def prepare_data(ehs, mice, times=None):
    """Prepare masked data."""
    data = {}
    if not isinstance(mice, list):
        mice = [mice]
    if times is None:
        ehs.unmask_data()
        times = (ehs.getstarttimes(mice)[0],
                 ehs.getendtimes(mice)[-1])
    t_start, t_end = times
    for mouse in mice:
        data[mouse] = []
        ads, sts, ens = get_ehs_data_with_margin(ehs, mouse, t_start, t_end)
        idxs = get_indices(t_start, t_end, sts, ens)
        for i in idxs:
            data[mouse].append((ads[i],
                                max(sts[i], t_start),
                                min(ens[i], t_end)))
    return data


def get_animal_position(times, antennas, mouse, threshold):
    out = []
    for t_start, t_end, an_start, an_end in zip(times[:-1], times[1:],
                                                antennas[:-1],
                                                antennas[1:]):
        delta_t = t_end - t_start
        # Workflow as agreed on 14 May 2014
        if delta_t < threshold:
            continue
        delta_an = np.abs(an_end - an_start)
        if delta_an == 0:
            out.append((address[an_start], mouse,
                        t_start, t_end, t_end-t_start, True))
        elif delta_an in [1, 7]:
            if an_end in same_pipe[an_start]:
                continue
            else:
                out.append((address[an_start], mouse,
                            t_start, t_end, t_end-t_start, True))
        elif delta_an in [2, 6]:
            out.append((surrounding[(min(an_start, an_end),
                                     max(an_start, an_end))],
                        mouse, t_start, t_end, t_end-t_start, False))
        elif delta_an in [3, 4, 5]:
            if an_end in opposite_pipe[an_start]:
                continue
            else:
                out.append((address_not_adjacent[an_start],
                            mouse, t_start, t_end, t_end-t_start, False))
    return out


def get_length(time_start, time_end, binsize):
    return int(np.ceil((time_end - time_start)/binsize))

def get_times(binsize, time_start=None, time_end=None):
    if time_start is None:
        time_start = 0
    if time_end is None:
        time_end = 43200
    length = get_length(time_start, time_end, binsize)
    out = np.linspace(time_start, time_end - binsize, length)
    return out.tolist()

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
    f_list = os.listdir(path)
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


def dict_to_array_2D(dictionary, keys1, keys2):
    shape = (len(keys1), len(keys2))
    out = np.zeros(shape)
    for i, key1 in enumerate(keys1):
        for j, key2 in enumerate(keys2):
            out[i, j] = dictionary[key1][key2]
    return out

def dict_to_array_3D(dictionary, keys1, keys2, keys3):
    shape = (len(keys1), len(keys2), len(keys3))
    out = np.zeros(shape)
    for i, key1 in enumerate(keys1):
        for j, key2 in enumerate(keys2):
            for k, key3 in enumerate(keys3):
                out[i, j, k] = dictionary[key1][key2][key3]
    return out


def calc_excess(res, exp_res):
    excess = OrderedDict()
    for key1 in res.keys():
        excess[key1] = OrderedDict()
        for key2 in res[key1].keys():
            excess[key1][key2] = OrderedDict()
            for key3 in res[key1][key2].keys():
                excess[key1][key2][key3] = res[key1][key2][key3] - exp_res[key1][key2][key3]
    return excess
