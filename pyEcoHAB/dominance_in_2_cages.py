#-*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import
import numpy as np
from . import utility_functions as utils
from . import exec_functions as dispatch
from .write_to_file import write_csv_alone

mouse_attention_span = 10  # sec
nbins = 10

opposite_antenna_dict = {1: 2,
                         2: 1,
                         3: 4,
                         4: 3,}
dt = 0.05
color_list = ['indianred', 'darkred', 'salmon',
              'darkolivegreen', 'forestgreen', 'darkslategrey',
              'darkslateblue', 'lightskyblue', 'navy',
              'mediumpurple', 'grey', 'gold',
              'orange', 'saddlebrown', 'magenta']


def get_states_and_home_cage_antenna(ehd, cf, dt):
    home_cage_antenna = ehd.get_home_cage_antenna()
    states = {}
    for mouse1 in ehd.mice:
        states[mouse1] = get_states(ehd, cf, mouse1, home_cage_antenna, dt)
    is_home_cage_antenna_ok = are_cages_correctly_assigned(states)
    if not is_home_cage_antenna_ok:
        home_cage_antenna = opposite_antenna_dict[home_cage_antenna]
        for mouse1 in ehd.mice:
            states[mouse1] = get_states(ehd, cf, mouse1, home_cage_antenna, dt)
    return states, home_cage_antenna

def get_time_spent(states, cage):
    if isinstance(states, list):
        return len(np.where(np.array(states) == cage)[0])
    if isinstance(states, np.ndarray):
        return len(np.where(states == cage)[0])

def are_cages_correctly_assigned(datas):
    out = []
    for mouse in datas:
        t_home = get_time_spent(datas[mouse], 0)
        t_chamber = get_time_spent(datas[mouse], 3)
        if t_home > t_chamber:
            out.append(0)
        else:
            out.append(3)

    if get_time_spent(out, 0) > get_time_spent(out, 3):
        return True
    return False

def get_states_mouse(antennas, times, t_start, t_end, home_antenna, dt):
    length = utils.get_timestamp(t_start, t_end, dt)
    states = np.zeros((length), dtype=int)
    if antennas[0] != home_antenna:
        timestamp = utils.get_timestamp(t_start, times[0], dt)
        states[:timestamp] = 3
        previous = 2
    else:
        previous = 1

    for i, a in enumerate(antennas[1:-1]):
        timestamp = utils.get_timestamp(t_start, times[i], dt)
        next_timestamp = utils.get_timestamp(t_start, times[i+1], dt)
        next_a = antennas[i + 1]

        if a != next_a:  # easy, the mouse is crossing the pipe
            states[timestamp:next_timestamp] = 1
            previous = 1
        elif a == next_a:
            if previous == 1:
                if a != home_antenna:
                    states[timestamp:next_timestamp] = 3

            else:
                if times[i+1] - times[i] > 2:
                    if a != home_antenna:
                        states[timestamp:next_timestamp] = 3
                else:
                    states[timestamp:next_timestamp] = 1

        previous = states[next_timestamp-1]
    if previous == 1:
        if next_a != home_antenna:
            states[next_timestamp:] = 3
    else:
        if  t_end - times[-1] < 2:
            states[next_timestamp:] = 1
        else:
            if next_a != home_antenna:
                states[next_timestamp:] = 3

    return states


def get_states(ehd, cf, mouse, home_antenna, dt=dt):
    """
    0 -- home cage, 1 -- pipe, 2 -- conditioning compartment
    """
    t_start, t_end = cf.gettime('ALL')
    times, antennas = utils.get_times_antennas(ehd, mouse,
                                               t_start, t_end)
    return get_states_mouse(antennas, times, t_start,
                            t_end, home_antenna, dt)

def find_stimulus_cage_mice(states, t_start, t_stop, beginning, dt):
    start = int(round((t_start - beginning)/dt))
    end = int(round((t_stop - beginning)/dt))
    mice = []
    for mouse in states:
        if np.any(states[mouse][start:end+1] == 3):
            mice.append(mouse)
    return mice


def get_dominating_mice(ehd, cf, phase, mouse, states, home_cage_antenna, dt):
    results = np.zeros((len(ehd.mice)))
    t_start, t_end = cf.gettime(phase)
    T_START, T_END = cf.gettime('ALL')
    time, antennas = utils.get_times_antennas(ehd, mouse, t_start, t_end)
    idx = 1
    mice = ehd.mice
    while True:
        if idx >= len(antennas):
            break
        if antennas[idx] == home_cage_antenna and antennas[idx-1] == home_cage_antenna:
            mice_list = find_stimulus_cage_mice(states, time[idx-1], time[idx], T_START, dt)
            for mouse in mice_list:
                results[mice.index(mouse)] += 1
            idx += 2
        idx += 1
    return results


def dominating_mice(ehd, cf, phase, states, home_cage_antenna, dt=dt):
    results = np.zeros((len(ehd.mice), len(ehd.mice)))
    for i, mouse in enumerate(ehd.mice):
        results[:, i] = get_dominating_mice(ehd, cf, phase, mouse, states,
                                            home_cage_antenna, dt=dt)
    return results


def tube_dominance_2_mice_single_phase(ehd, mouse1, mouse2, t_start, t_end, home_cage_antenna):
    """We're checking here, how many times mouse1 dominates over mouse2
    between t_start and t_end.

    """
      
    m1_times, m1_antennas = utils.get_times_antennas(ehd, mouse1,
                                                     t_start, t_end)
    m2_times, m2_antennas = utils.get_times_antennas(ehd, mouse2,
                                                     t_start, t_end)
    domination_counter = check_mouse1_defending(m1_antennas, m1_times, m2_antennas, m2_times, home_cage_antenna)
        
    return domination_counter


def tube_dominance_2_cages(ehd, cf, phase, home_cage_antenna):
    mice = ehd.mice
    st, en = cf.gettime(phase)
    dominance =  np.zeros((len(mice), len(mice)))
    for i, mouse1 in enumerate(mice):
        for j, mouse2 in enumerate(mice):
            if i != j:
                dominance[i, j] = tube_dominance_2_mice_single_phase(ehd,
                                                                     mouse1,
                                                                     mouse2,
                                                                     st,
                                                                     en,
                                                                     home_cage_antenna)
    return dominance


def count_attempts(tstamp1, tstamp2, times, antennas, home_cage_antenna):
    mouse_between = utils.get_idx_between(tstamp1, tstamp2, times)
    in_between_antennas = utils.get_antennas(mouse_between, antennas)
    opposite_antenna = opposite_antenna_dict[home_cage_antenna]
    mouse_after = mouse_between[-1] + 1
    if  len(mouse_between) == 1:
        if in_between_antennas[0] != opposite_antenna:
            if mouse_after < len(antennas):
                if antennas[mouse_after] != opposite_antenna:
                    return 1
        return 0

    counter = 0
    i = 1
    while True:
        if in_between_antennas[i] == home_cage_antenna and in_between_antennas[i-1] == home_cage_antenna:
            counter += 1
            i += 1
        i += 1
        if i >= len(in_between_antennas):
            break
    return counter

def check_mouse1_not_valid(mouse_previous_antenna,
                           mouse_antenna,
                           home_cage_antenna):

    if mouse_antenna != mouse_previous_antenna:
        return True
    if mouse_antenna == home_cage_antenna:
        return True # mouse1 is trying to enter the pipe

    return False


def check_mouse2_not_valid(mouse1_previous_timestamp, mouse1_timestamp,
                       antennas2, times2,
                       home_cage_antenna):

    mouse2_pre = utils.get_idx_pre(mouse1_previous_timestamp, times2)

    if mouse2_pre is None:
        return True
    mouse2_between = utils.get_idx_between(mouse1_previous_timestamp,
                                           mouse1_timestamp, times2)
    if len(mouse2_between) == 0:
        return True # mouse2 is not moving during antenna readouts

    mouse2_after = utils.get_idx_post(mouse1_timestamp, times2)

    if antennas2[mouse2_pre] != home_cage_antenna:
        return True #mouse 2 didn't start at the home cage
    return False


def check_mouse1_defending(antennas1, times1, antennas2, times2, home_cage_antenna):
    dominance_counter = 0
    for idx in range(1, len(antennas1)):
        mouse1_previous_antenna, mouse1_antenna = antennas1[idx-1:idx+1]
        mouse1_previous_timestamp, mouse1_timestamp = times1[idx-1:idx+1]
        if  check_mouse1_not_valid(mouse1_previous_antenna,
                                   mouse1_antenna,
                                   home_cage_antenna):
            continue
        if check_mouse2_not_valid(mouse1_previous_timestamp,
                                  mouse1_timestamp,
                                  antennas2, times2,
                                  home_cage_antenna):
            continue
        dominance_counter += count_attempts(mouse1_previous_timestamp,
                                            mouse1_timestamp, times2,
                                            antennas2,
                                            home_cage_antenna)
    return dominance_counter


def get_tube_dominance_2_cages(ehd, cf, res_dir=None, prefix=None, dt=dt,
                               delimiter=";"):
    states, home_cage_antenna = get_states_and_home_cage_antenna(ehd,
                                                                 cf,
                                                                 dt)
    if res_dir is None:
        res_dir = ehd.res_dir
    if prefix is None:
        prefix = ehd.prefix
    dispatch.evaluate_whole_experiment(ehd, cf, res_dir, prefix,
                                       tube_dominance_2_cages,
                                       'mouse_pushing_out_conditioning_compartment',
                                       'dominating mouse',
                                       'pushed out mouse',
                                       '# pushes',
                                       args=[home_cage_antenna],
                                       vmin=0, vmax=200,
                                       delimiter=delimiter)


def get_subversion_evaluation(ehd, cf, res_dir=None, prefix=None, dt=dt,
                              delimiter=";"):
    if res_dir is None:
        res_dir = ehd.res_dir
    if prefix is None:
        prefix = ehd.prefix
    states, home_cage_antenna = get_states_and_home_cage_antenna(ehd,
                                                                 cf,
                                                                 dt)
    dispatch.evaluate_whole_experiment(ehd, cf, res_dir, prefix,
                                       dominating_mice,
                                       'subversion_evaluation',
                                       'dominating mouse',
                                       'subversive mouse',
                                       '# times in small cage',
                                       args=[states, home_cage_antenna,
                                             dt], vmin=0, vmax=200,
                                       delimiter=delimiter)


def how_many_visits(states, t_start, t_end, T_0, dt):
    t1_stamp = utils.get_timestamp(T_0, t_start, dt)
    t2_stamp = utils.get_timestamp(T_0, t_end, dt)
    new_states = states[t1_stamp:t2_stamp]
    transitions = new_states[1:] - new_states[:-1]
    where = np.where(transitions == 2)[0]
    return len(where)


def get_visits_to_stimulus_cage(ehd, cf, res_dir="", prefix="", dt=dt,
                                delimiter=";"):
    if res_dir is "":
        res_dir = ehd.res_dir
    if prefix is "":
        prefix = ehd.prefix
    states, home_cage_antenna = get_states_and_home_cage_antenna(ehd,
                                                                 cf,
                                                                 dt)
    phases = utils.filter_dark_light(cf.sections())
    results = np.zeros((1,  len(ehd.mice), len(phases)))
    cumulative = np.zeros((1, len(ehd.mice)))
    T_0, T_1 = cf.gettime("ALL")
    for i, phase in enumerate(phases):
        t_start, t_end = cf.gettime(phase)
        for j, mouse in enumerate(ehd.mice):
            results[0, j, i] = how_many_visits(states[mouse], t_start, t_end, T_0, dt)
    for j, mouse in enumerate(ehd.mice):
        cumulative[0, j] = how_many_visits(states[mouse], T_0, T_1, T_0, dt)
        assert cumulative[0, j] == sum(results[0, j, :])
    write_csv_alone(results, phases, ehd.mice,
                    res_dir, prefix, labels=["conditioning compartment"],
                    header='Number of visits %s\n',
                    fname='visits_to_conditioning_compartment_%s.csv',
                    directory="mice_in_conditioning_compartment",
                    delimiter=delimiter)
