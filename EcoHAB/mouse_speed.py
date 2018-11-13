from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utils
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
bins = 2000
homepath = os.path.expanduser("~/")
threshold = 2
from numba import jit

def in_tube(antenna, next_antenna):
    if antenna % 2:
        if next_antenna  == antenna  + 1:
            return True
    else:
        if next_antenna == antenna - 1:
            return True
    return False

def make_out():
    return {'antennas':[],
           'mouse1_time':[],
           'mouse2_time':[],
           'mouse1_delta_t':[],
           'mouse2_delta_t':[]}

def add_to_out(out, key, t1, t2, delta_t1, delta_t2):
     out['antennas'].append(key)
     out['mouse1_time'].append(t1)
     out['mouse2_time'].append(t2)
     out['mouse1_delta_t'].append(delta_t1)
     out['mouse2_delta_t'].append(delta_t2)

@jit
def following_2_mice(ehd, mouse1, mouse2, st, en):
    ehd.mask_data(st, en)
    antennas1 = ehd.getantennas(mouse1)
    times1 = ehd.gettimes(mouse1)
    antennas2 = ehd.getantennas(mouse2)
    times2 = ehd.gettimes(mouse2)
    change_indices = np.where((np.array(antennas1[1:]) - np.array(antennas1[:-1])) != 0)[0]
    out = make_out()
    for idx in change_indices:
        antenna1 = antennas1[idx]
        next_antenna1 = antennas1[idx+1]
        delta_t1 = times1[idx+1] - times1[idx]
        key1 = str(antenna1) + str(next_antenna1)
        if in_tube(antenna1, next_antenna1):
            idxs1 = np.where(np.array(times2) >= times1[idx])[0]
            idxs2 = np.where(np.array(times2) <= times1[idx]+threshold)[0]
            common_idxs = list(set(idxs1)&set(idxs2))
            for ci in common_idxs:
                antenna2 = antennas2[ci]
                if ci + 1 < len(antennas2):
                    next_antenna2 = antennas2[ci+1]
                    delta_t2 = times2[ci+1] - times2[ci]
                    if antenna2 == antenna1:
                        if in_tube(antenna2, next_antenna2):
                           add_to_out(out, key1, times1[idx], times2[ci], delta_t1, delta_t2)
    return out


def following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en):
    ehd.mask_data(st, en)
    antennas1 = ehd.getantennas(mouse1)
    times1 = ehd.gettimes(mouse1)
    antennas2 = ehd.getantennas(mouse2)
    times2 = ehd.gettimes(mouse2)
    change_indices = np.where((np.array(antennas1[1:]) - np.array(antennas1[:-1])) != 0)[0]
    out = make_out()
     
    for idx in change_indices:
        antenna1 = antennas1[idx]
        next_antenna1 = antennas1[idx+1]
        delta_t1 = times1[idx+1] - times1[idx]
        key1 = str(antenna1) + str(next_antenna1)
        if in_tube(antenna1, next_antenna1):
            idxs1 = np.where(np.array(times2) >= times1[idx])[0]
            idxs2 = np.where(np.array(times2) <= times1[idx]+threshold)[0]
            common_idxs = list(set(idxs1)&set(idxs2))
            for ci in common_idxs:
                antenna2 = antennas2[ci]
                if ci + 1 < len(antennas2):
                    next_antenna2 = antennas2[ci+1]
                    delta_t2 = times2[ci+1] - times2[ci]
                    if antenna2 == antenna1 and times1[idx+1]>times2[ci]:
                        if in_tube(antenna2, next_antenna2):
                           add_to_out(out, key1, times1[idx], times2[ci], delta_t1, delta_t2)
    return out


@jit
def following_for_all(ehd, cf):
    mice = ehd.mice
    phases = cf.sections()
    phases.remove('ALL')
    followings = {}
    for phase in phases:
            print(phase)
            st, en = cf.gettime(phase)
            followings[phase] = {}
            for mouse1 in mice:
                antennas1 = ehd.getantennas(mouse1)
                times1 = ehd.gettimes(mouse1)
                for mouse2 in mice:
                    if mouse2 != mouse1:
                        key = mouse1 + '_' + mouse2
                        followings[phase][key] = following_2_mice(ehd, mouse1, mouse2, st, en)
    return followings

@jit
def following_for_all_in_pipe_condition(ehd, cf):
    mice = ehd.mice
    phases = cf.sections()
    phases.remove('ALL')
    followings = {}
    for phase in phases:
            print(phase)
            st, en = cf.gettime(phase)
            followings[phase] = {}
            for mouse1 in mice:
                antennas1 = ehd.getantennas(mouse1)
                times1 = ehd.gettimes(mouse1)
                for mouse2 in mice:
                    if mouse2 != mouse1:
                        key = mouse1 + '_' + mouse2
                        followings[phase][key] = following_2_mice_in_pipe_condition(ehd, mouse1, mouse2, st, en)
    return followings

def save_followings_to_file(directory, prefix, followings, ending=''):
    dir_ = utils.check_directory(directory, 'following_in_pipe')
    fname_base = "%s_followings_wyniki_%s_%s.csv"
    header = 'mouse pair; measurement\n'
    measurement_keys = ['antennas', 'mouse1_time', 'mouse1_delta_t', 'mouse2_time', 'mouse2_delta_t']
    phases = followings.keys()
    for phase in phases:
        fname = os.path.join(dir_, fname_base % (ending, prefix, phase))
        f = open(fname, 'w')
        f.write(header)
        for key in followings[phase]:
            for mk in measurement_keys:
                f.write(key +';'+mk)
                for meas in followings[phase][key][mk]:
                    f.write(';'+str(meas))
                f.write('\n')
        f.close()


if __name__ == '__main__':

    for new_path in datasets:
       
        path = os.path.join(homepath, new_path)
        prefix = utils.make_prefix(path)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
        else:
            remove_mouse = None
        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 500
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        cf = ExperimentConfigFile(path)
        followings = following_for_all(ehd, cf)
        prefix = utils.make_prefix(path)
        res_dir = utils.results_path(path)
        save_followings_to_file(res_dir, prefix, followings)
        followings = following_for_all_in_pipe_condition(ehd, cf)
        prefix = utils.make_prefix(path)
        res_dir = utils.results_path(path)
        save_followings_to_file(res_dir, prefix, followings, ending='first_mouse_in_pipe')
