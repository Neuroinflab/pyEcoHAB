# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
import time
import sys
from collections import OrderedDict
import numpy as np


# NamedDict class was originally written by Zbyszek Jędrzejewski-Szmek
# and Avrama Blackwell for moose_nerp https://github.com/neurord/moose_nerp


def check_directory(directory, subdirectory=None):
    if subdirectory:
        new_path = os.path.join(directory, subdirectory)
    else:
        new_path = directory
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    return new_path


def make_figure(title):
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, aspect='equal')
    fig.suptitle('%s' % (title), fontsize=14, fontweight='bold')
    return fig, ax


def list_of_pairs(mice):
    pair_labels = []
    for j, mouse in enumerate(mice):
        for k in range(j+1, len(mice)):
            pair_labels.append(mice[j] + '|' + mice[k])
    return pair_labels


def all_mouse_pairs(mice, reverse=False):
    pair_labels = []
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            if reverse is True:
                pair_labels.append('%s|%s' % (mouse2, mouse1))
            else:
                pair_labels.append('%s|%s' % (mouse1, mouse2))
    return pair_labels


def make_table_of_pairs(FAM, phases, mice):
    new_shape = (len(mice)*(len(mice)-1)//2, len(phases))
    output = np.zeros(new_shape)
    pair_labels = list_of_pairs(mice)
    for i in range(len(phases)):
        count = 0
        for j in range(len(mice)):
            for k in range(j + 1, len(mice)):
                output[count, i] = FAM[i, j, k]
                count += 1

    return output, pair_labels


def make_table_of_all_mouse_pairs(FAM, phases, mice, reverse=False):
    new_shape = (len(mice)*(len(mice)-1), len(phases))
    output = np.zeros(new_shape)
    pair_labels = all_mouse_pairs(mice,
                                  reverse=reverse)
    for i, phase in enumerate(phases):
        count = 0
        for j in range(len(mice)):
            for k in range(len(mice)):
                if j == k:
                    continue
                if reverse is True:
                    output[count, i] = FAM[i, k, j]
                else:
                    output[count, i] = FAM[i, j, k]
                count += 1
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
    if remove_mouse is None or remove_mouse == '' or len(remove_mouse) == 0:
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
    return np.where((np.array(times) >= t0) & (np.array(times) < t1))[0]


def get_idx_post(t1, times):
    idxs = np.where(np.array(times) > t1)[0]
    if len(idxs):
        return idxs[0]
    return None


def change_state(antennas):
    indx = []
    for i, a in enumerate(antennas[:-1]):
        if a != antennas[i+1]:
            indx.append(i)
    return indx


def get_times_antennas(e_data, mouse, t_1, t_2):
    if t_1 == 0 and t_2 == -1:
        e_data.unmask_data()
        return e_data.get_times(mouse), e_data.get_antennas(mouse)
    e_data.mask_data(t_1, t_2)
    antennas, times = e_data.get_antennas(mouse), e_data.get_times(mouse)
    e_data.unmask_data()
    return times, antennas


def get_times_antennas_list_of_mice(ecohab_data, mice, t_1, t_2):
    out = {}
    for mouse in mice:
        out[mouse] = {}
        times, antennas = get_times_antennas(ecohab_data, mouse, t_1, t_2)
        out[mouse]["times"] = times
        out[mouse]["antennas"] = antennas
    return out


def get_antennas(idxs, antennas):
    antenna_slice = []
    for new_idx in idxs:
        antenna_slice.append(antennas[new_idx])
    return antenna_slice


def get_timestamp(t_start, t_end, dt):
    return int(round((t_end - t_start)/dt))


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
    return sorted(list(set(idx_start + idx_end)))


def get_ecohab_data_with_margin(ecohab_data, mouse, t_start, t_end,
                                margin=12*3600):
    if t_start == 0 and t_end == -1:
        return ecohab_data.get_visit_addresses(mouse),\
            ecohab_data.get_starttimes(mouse),\
            ecohab_data.get_endtimes(mouse)
    ecohab_data.mask_data(t_start - margin, t_end + margin)
    adresses = ecohab_data.get_visit_addresses(mouse)
    starts = ecohab_data.get_starttimes(mouse)
    ends = ecohab_data.get_endtimes(mouse)
    ecohab_data.unmask_data()
    return adresses, starts, ends


def prepare_data(ecohab_data, mice, times=None):
    """Prepare masked data."""
    data = {}
    if not isinstance(mice, list):
        mice = [mice]
    if times is None:
        ecohab_data.unmask_data()
        times = (ecohab_data.get_starttimes(mice)[0],
                 ecohab_data.get_endtimes(mice)[-1])
    t_start, t_end = times
    for mouse in mice:
        data[mouse] = []
        ads, sts, ens = get_ecohab_data_with_margin(ecohab_data, mouse,
                                                    t_start, t_end)
        idxs = get_indices(t_start, t_end, sts, ens)
        for i in idxs:
            data[mouse].append((ads[i],
                                max(sts[i], t_start),
                                min(ens[i], t_end)))
    return data


def get_animal_position(times, antennas, mouse, threshold, same_pipe,
                        same_address, opposite_pipe, address, surrounding,
                        address_not_adjacent, internal_antennas):
    out = []
    i = 0
    if len(times) < 2:
        return []
    t_start, an_start = times[0], antennas[0]
    t_end, an_end = times[1], antennas[1]
    while i < len(times) - 1:
        delta_t = t_end - t_start
        if an_start in internal_antennas:
            while an_end == an_start:
                i = i + 1
                try:
                    t_end, an_end = times[i+1], antennas[i+1]
                except IndexError:
                    out.append((address[an_start], mouse,
                                t_start, t_end, t_end-t_start, True))
                    return out
            out.append((address[an_start], mouse,
                        t_start, t_end, t_end-t_start, True))
        elif an_end in internal_antennas:
            an_old_end = an_end
            while an_end == an_old_end:
                i = i + 1
                try:
                    an_end = antennas[i+1]
                    t_end = times[i+1]
                except IndexError:
                    out.append((address[an_old_end], mouse,
                                t_start, t_end, t_end-t_start, True))
                    return out
            out.append((address[an_old_end], mouse,
                        t_start, t_end, t_end-t_start, True))

        elif delta_t < threshold:
            pass
        elif an_end == an_start:
            out.append((address[an_start], mouse,
                        t_start, t_end, delta_t, True))
        elif an_start in same_pipe and an_end in same_pipe[an_start]:
            pass
        elif an_end in same_address[an_start]:
            out.append((address[an_start], mouse,
                        t_start, t_end, delta_t, True))
        elif (min(an_start, an_end), max(an_start, an_end)) in surrounding:
            out.append((surrounding[(min(an_start, an_end),
                                     max(an_start, an_end))],
                        mouse, t_start, t_end, delta_t, False))
        elif an_start in opposite_pipe and an_end in opposite_pipe[an_start]:
            pass
        else:
            out.append((address_not_adjacent[an_start],
                        mouse, t_start, t_end, delta_t, False))

        i = i + 1
        try:
            an_start, an_end = antennas[i], antennas[i+1]
            t_start, t_end = times[i], times[i+1]
        except IndexError:
            return out
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
                try:
                    out[i, j, k] = dictionary[key1][key2][key3]
                except TypeError:
                    out[i, j, k] = 0
    return out


def calc_excess(res, exp_res):
    excess = OrderedDict()
    for key1 in res.keys():
        excess[key1] = OrderedDict()
        for key2 in res[key1].keys():
            excess[key1][key2] = OrderedDict()
            for key3 in res[key1][key2].keys():
                excess[key1][key2][key3] = res[key1][key2][key3]\
                    - exp_res[key1][key2][key3]
    return excess


def get_dark_light_data(phase, timeline, ecohab_data, mice):
    if phase == "dark" or phase == "DARK" or phase == "Dark":
        phases = filter_dark(timeline.sections())
    elif phase == "light" or phase == "LIGHT" or phase == "Light":
        phases = filter_light(timeline.sections())
    out_phases = [phase]
    data = {mouse: [] for mouse in mice}
    total_time = 0
    for i, ph in enumerate(phases):
        time = timeline.get_time_from_epoch(ph)
        out = prepare_data(ecohab_data, mice, time)
        for mouse in mice:
            data[mouse].extend(out[mouse])
        total_time += (time[1] - time[0])
    out_data = {phase: {0: data}}
    return out_phases, {phase: {0: total_time}}, out_data

def get_custom_range_data(phase_def, timeline, ecohab_data, mice):
    out_phases = []
    out_data = {}   
    times = {}
    for i, phases in enumerate(phase_def):
        total_time = 0
        phase_name = ""
        data = {mouse: [] for mouse in mice}
        if isinstance(phases, str):
            phases = [phases]
        assert isinstance(phases, list)
        for ph in phases:
            phase_name += ph.replace(" ", "_") + "_"
            time = timeline.get_time_from_epoch(ph)
            out = prepare_data(ecohab_data, mice, time)
            for mouse in mice:
                data[mouse].extend(out[mouse])
            total_time += (time[1] - time[0])
        out_phases.append(phase_name)
        out_data[phase_name] = {0: data}
        times[phase_name] = {0: total_time}
    return out_phases, times, out_data

            
def prepare_binned_data(ecohab_data, timeline, bins, mice):
    total_time = OrderedDict()
    data = OrderedDict()
    if bins in ["ALL", "all", "All"]:
        phases = ["ALL"]
        time = timeline.get_time_from_epoch("ALL")
        total_time["ALL"] = {0: (time[1] - time[0])}
        data["ALL"] = {0: prepare_data(ecohab_data, mice, time)}
        keys = [["ALL"], {"ALL": [0]}]
    elif bins in ['dark', "DARK", "Dark", "light", "LIGHT", "Light"]:
        phases, total_time, data = get_dark_light_data(bins, timeline,
                                                       ecohab_data, mice)
        labels = {}
        for key in data.keys():
            labels[key] = [0]
        keys = [list(data.keys()), labels]
    elif isinstance(bins, list):
        phases, total_time, data = get_custom_range_data(bins, timeline,
                                                         ecohab_data, mice)
        labels = {}
        for key in data.keys():
            labels[key] = [0]
        keys = [list(data.keys()), labels]
    elif (isinstance(bins, str) and
          bins.lower() in ["whole_phase", "whole phase"]):
        phases = []
        all_phases = filter_dark_light(timeline.sections())
        phases = [phase.replace(" ", "_") for phase in all_phases]
        times = [timeline.get_time_from_epoch(phase)
                 for phase in all_phases]
        data = OrderedDict()
        total_time = OrderedDict()
        bin_labels = {}
        for phase in all_phases:
            bin_labels[phase] = [0]
            time = timeline.get_time_from_epoch(phase)
            total_time[phase] = {}
            total_time[phase][0] = time[-1] - time[0]
            data[phase] = {}
            data[phase][0] = prepare_data(ecohab_data, mice,
                                          time)
        keys = [all_phases, bin_labels]
    elif isinstance(bins, int) or isinstance(bins, float):
        phases = []
        all_phases = filter_dark_light(timeline.sections())
        shortest_phase = get_shortest_phase_duration(timeline)
        bin_labels = {}
        # you can not iterate by phases, if bins are longer than phases
        if bins > shortest_phase:
            t_start = timeline.get_time_from_epoch(all_phases[0])[0]
            t_end = timeline.get_time_from_epoch(all_phases[-1])[-1]
            all_phases = []
            times = []
            i = 1
            while t_start < t_end:
                times.append((t_start, t_start + bins))
                all_phases.append("%d_x" % i)
                bin_labels["%d_x" % i] = [0]
                i += 1
                t_start += bins
        else:
            times = [timeline.get_time_from_epoch(phase)
                     for phase in all_phases]
        for i, phase in enumerate(all_phases):
            t_start, t_end = times[i]
            phases.append("%s_%4.2fh" % (phase.replace(" ", "_"), bins/3600))
            bin_labels[phase] = get_times(bins, time_end=t_end-t_start)
            data[phase] = OrderedDict()
            total_time[phase] = OrderedDict()
            j = 0
            while t_start < t_end:
                t_e = t_start + bins
                if t_e > t_end:
                    t_e = t_end
                time = [t_start, t_e]
                data[phase][bin_labels[phase][j]] = prepare_data(ecohab_data,
                                                                 mice,
                                                                 time)
                total_time[phase][bin_labels[phase][j]] = time[1] - time[0]
                t_start += bins
                j += 1
        keys = [all_phases, bin_labels]
    return phases, total_time, data, keys


def extract_directions(times, antennas, last_antenna, keys):
    direction_dict = {key: [[], []] for key in keys}
    change_indices = change_state(antennas)
    for c_idx in change_indices:
        if c_idx + 1 >= len(antennas):
            break
        ant, next_ant = antennas[c_idx], antennas[c_idx + 1]
        key = "%s %s" % (ant, next_ant)
        if key in keys:
            try:
                third_antenna = antennas[c_idx + 2]
            except IndexError:
                third_antenna = last_antenna
            if third_antenna == ant:
                continue
            direction_dict[key][0].append(times[c_idx])
            direction_dict[key][1].append(times[c_idx + 1])
    return direction_dict


def extract_backing(times, antennas, last_antenna, setup):
    direction_dict = {key: [[], []] for key in setup.backing}
    internal = setup.internal_antennas
    i = 1
    while i < len(antennas) - 1:
        ant = antennas[i]
        prev_ant, next_ant = antennas[i-1], antennas[i+1]
        opposite_antenna = setup.other_tunnel_antenna(ant)[0]
        if ant not in internal:
            if prev_ant == opposite_antenna and next_ant == opposite_antenna:
                key = "%s %s" % (antennas[i-1],
                                 antennas[i+1])
                direction_dict[key][0].append(times[i-1])
                direction_dict[key][1].append(times[i+1])
                i = i+1
            elif ant == next_ant and prev_ant != opposite_antenna:
                key = "%s %s" % (antennas[i],
                                 antennas[i+1])
                direction_dict[key][0].append(times[i])
                direction_dict[key][1].append(times[i+1])
                i = i+1
        i = i+1

    return direction_dict


def prepare_for_tube_dominance(ecohab_data, mice, st, en):
    moves = {
        "directions": {},
        "backing out": {}
    }
    moves["directions"] = prepare_tube_data(ecohab_data, mice, st, en,
                                            "directions")
    moves["backing out"] = prepare_tube_data(ecohab_data, mice, st, en,
                                             "backing_out")
    return moves


def prepare_tube_data(ecohab_data, mice, st, en, what="directions"):
    directions = {}
    for j, mouse1 in enumerate(mice):
        times_antennas = get_times_antennas(ecohab_data,
                                            mouse1,
                                            st,
                                            en)
        last_times, last_antennas = get_times_antennas(ecohab_data,
                                                       mouse1,
                                                       en,
                                                       en+(en-st))

        try:
            last_antenna = last_antennas[0]
        except IndexError:
            last_antenna = None
        if what == "directions":
            directions[mouse1] = extract_directions(times_antennas[0],
                                                    times_antennas[1],
                                                    last_antenna,
                                                    ecohab_data.directions)
        elif what == "backing_out":
            directions[mouse1] = extract_backing(times_antennas[0],
                                                 times_antennas[1],
                                                 last_antenna,
                                                 ecohab_data.setup_config)
    return directions


def get_registrations_bins(ecohab_data, timeline, bins, mice,
                           function):
    total_time = OrderedDict()
    data = OrderedDict()
    if bins in ["ALL", "all", "All"]:
        phases = ["ALL"]
        time = timeline.get_time_from_epoch("ALL")
        data["ALL"] = {0: function(ecohab_data, mice, *time)}
        data_keys = [["ALL"], {"ALL": [0.0]}]
        total_time["ALL"] = {0: time}
    elif (isinstance(bins, str)
          and bins.lower() in ["whole_phase", "whole phase",
                               "whole phases",
                               "whole_phases"]):
        phases = []
        all_phases = filter_dark_light(timeline.sections())
        phases = [phase.replace(" ", "_") for phase in all_phases]
        times = [timeline.get_time_from_epoch(phase)
                 for phase in all_phases]
        data = OrderedDict()
        total_time = OrderedDict()
        bin_labels = {}
        for phase in all_phases:
            bin_labels[phase] = [0]
            time = timeline.get_time_from_epoch(phase)
            total_time[phase] = {}
            total_time[phase][0] = (time[0], time[-1])
            data[phase] = {}
            data[phase][0] = function(ecohab_data, mice,
                                      *time)
        data_keys = [all_phases, bin_labels]
    elif isinstance(bins, int) or isinstance(bins, float):
        phases = []
        all_phases = filter_dark_light(timeline.sections())
        min_phase = int(get_shortest_phase_duration(timeline))
        bin_labels = {}
        # you can not iterate by phases, if bins are longer than phases
        if bins > min_phase:
            t_start = timeline.get_time_from_epoch(all_phases[0])[0]
            t_end = timeline.get_time_from_epoch(all_phases[-1])[-1]
            all_phases = []
            times = []
            i = 1
            while t_start < t_end:
                times.append((t_start, t_start + bins))
                all_phases.append("%d_x" % i)
                bin_labels["%d_x" % i] = [0]
                i += 1
                t_start += bins
        else:
            times = [timeline.get_time_from_epoch(phase)
                     for phase in all_phases]
        for i, phase in enumerate(all_phases):
            t_start, t_end = times[i]
            phases.append("%s_%4.2fh" % (phase.replace(" ", "_"),
                                         bins/3600))
            data[phase] = OrderedDict()
            total_time[phase] = OrderedDict()
            bin_labels[phase] = get_times(bins, time_end=t_end-t_start)
            j = 0
            while t_start < t_end:
                t_e = t_start + bins
                if t_e > t_end:
                    t_e = t_end
                time = (t_start, t_e)
                data[phase][bin_labels[phase][j]] = function(ecohab_data,
                                                             mice,
                                                             *time)
                total_time[phase][bin_labels[phase][j]] = time
                t_start += bins
                j += 1
            data_keys = [all_phases, bin_labels]
    else:
        phases = []
        all_phases = filter_dark_light(timeline.sections())
        phases = [phase.replace(" ", "_") for phase in all_phases]
        times = [timeline.get_time_from_epoch(phase)
                 for phase in all_phases]
        data = OrderedDict()
        total_time = OrderedDict()
        bin_labels = {}
        for phase in all_phases:
            bin_labels[phase] = [0]
            time = timeline.get_time_from_epoch(phase)
            total_time[phase] = {}
            total_time[phase][0] = (time[0], time[-1])
            data[phase] = {}
            data[phase][0] = function(ecohab_data, mice,
                                      *time)
        data_keys = [all_phases, bin_labels]
    return phases, total_time, data, data_keys


def make_results_dict(mice, tolist=False):
    result = OrderedDict()
    for mouse1 in mice:
        result[mouse1] = OrderedDict()
        for mouse2 in mice:
            if not tolist:
                result[mouse1][mouse2] = 0
            else:
                result[mouse1][mouse2] = []

    return result


def make_all_results_dict(phases, bins):
    result = OrderedDict()
    for phase in phases:
        result[phase] = OrderedDict()
        for bin1 in bins[phase]:
            result[phase][bin1] = 0
    return result


def get_shortest_phase_duration(timeline):
    durs = []
    for phase in timeline.sections():
        time = timeline.get_time_from_epoch(phase)
        durs.append(time[1] - time[0])
    return min(durs)


def to_struck(string, fname=""):
    new_str = string + ' UTC'
    try:
        return time.strptime(new_str, '%d.%m.%Y%H:%M %Z')
    except ValueError:
        try:
            return time.strptime(new_str, '%d.%m.%Y%H:%M:%S %Z')
        except ValueError:
            raise Exception('Wrong date format in %s' % fname)


def diagonal_reflection_3D(matrix_data):
    result = OrderedDict()
    for key1 in matrix_data.keys():
        result[key1] = OrderedDict()
        for key2 in matrix_data[key1].keys():
            result[key1][key2] = OrderedDict()
            for key3 in matrix_data[key1][key2].keys():
                if key2 == key3:
                    result[key1][key2][key3] = matrix_data[key1][key2][key3]
                elif matrix_data[key1][key2][key3] == 0:
                    result[key1][key2][key3] = matrix_data[key1][key3][key2]
                else:
                    result[key1][key2][key3] = matrix_data[key1][key2][key3]
    return result


def sum_per_mouse(data, mice, binlabels, position="leading"):
    """
    position can either be leading or following.
    """
    sum_value = OrderedDict()
    for bi in binlabels:
        sum_value[bi] = OrderedDict()
        for mouse1 in mice:
            sum_value[bi][mouse1] = 0
            for mouse2 in mice:
                if mouse1 == mouse2:
                    continue
                if position == "leading":
                    sum_value[bi][mouse1] += data[bi][mouse1][mouse2]
                elif position == "following":
                    sum_value[bi][mouse1] += data[bi][mouse2][mouse1]
    return sum_value


def sum_activity(activity, phases, mice, bin_labels):
    """Activity data is in the dictionary [address][phase][mouse][bin_no]
    Following dicts are: [phase][bin_label][mouse_lead][mouse_follow]
    we want activity in the form [phase][bin_label][mouse],
    so dividing is easy
    """
    visits = OrderedDict()
    for phase in phases:
        visits[phase] = OrderedDict()
        for i, lab in enumerate(bin_labels[phase]):
            visits[phase][lab] = OrderedDict()
            for m in mice:
                visits[phase][lab][m] = 0
    for A in activity.keys():
        for phase in activity[A][0].keys():
            for i, lab in enumerate(bin_labels[phase]):
                for m in mice:
                    visits[phase][lab][m] += activity[A][0][phase][m][i]
    return visits


def divide_sum_activity(data_sum, data_activ):
    """
    Following/Leading (data_sum) dicts are: [phase][bin_label][mouse]
    activity (data) dict is [phase][bin_label][mouse]
    """
    result = OrderedDict()
    for label in data_sum.keys():
        result[label] = OrderedDict()
        for mouse in data_sum[label].keys():
            fol = data_sum[label][mouse]
            act = data_activ[label][mouse]
            try:
                result[label][mouse] = fol/act
            except ZeroDivisionError:
                result[label][mouse] = 0
    return result


def mean(numerator, N):
    result = OrderedDict()
    for key1 in numerator.keys():
        result[key1] = OrderedDict()
        for key2 in numerator[key1].keys():
            result[key1][key2] = numerator[key1][key2] / N
    return result


def standard_error(data, mean, N):
    result = OrderedDict()
    for key1 in data.keys():
        result[key1] = OrderedDict()
        for mouse1 in data[key1].keys():
            result[key1][mouse1] = 0
            for mouse2 in data[key1][mouse1]:
                if mouse1 != mouse2:
                    result[key1][mouse1] += (data[key1][mouse1][mouse2]
                                             - mean[key1][mouse1]) ** (2)
                else:
                    continue
            if N > 1:
                denom = N*(N-1)
                result[key1][mouse1] = (result[key1][mouse1] / denom) ** (0.5)
            else:
                result[bi][mouse1] = 0
    return result
