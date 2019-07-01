# -*- coding: utf-8 -*-
from __future__ import division, print_function

import numpy as np
import os
import analiza_friends as af
import time_spent_in_cage as ts
import mouse_speed as ms
import utility_functions as utils
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
from write_to_file import save_data_cvs
import interactions
import tube_dominance as td
import dispatch

homepath = os.path.expanduser("~/")
threshold = 3
if __name__ == '__main__':
    

    for new_path in datasets:
       
        path = os.path.join(homepath, new_path)
        prefix = utils.make_prefix(path)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
        else:
            remove_mouse = None
        if new_path == "EcoHAB_data_November/C57_social_dominance_28-31.05/":
            remove_antennas = [3, 4, 5, 6, 7, 8]
        else:
            remove_antennas = []
        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 10
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path], remove_antennas=remove_antennas)
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path],
                                    remove_antennas=remove_antennas)

        ehs = EcoHab.EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path)
        tstart, tend = cf.gettime('ALL')
        directory = utils.results_path(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        af.mouse_alone_ehs(ehs, cf, directory, prefix)
        af.get_in_cohort_sociability(ehs, cf, directory, prefix,
                                 remove_mouse=remove_mouse)
        af.get_in_cohort_sociability(ehs, cf, directory, prefix,
                                 which_phases="ALL",
                                 remove_mouse=remove_mouse)
        af.get_in_cohort_sociability(ehs, cf, directory, prefix,
                                 which_phases="dark",
                                 remove_mouse=remove_mouse)
        af.get_in_cohort_sociability(ehs, cf, directory, prefix,
                                 which_phases="light",
                                 remove_mouse=remove_mouse)
        ms.get_following(ehd, cf, directory,
                         prefix,
                         remove_mouse=remove_mouse)
        for binsize in binsizes:
            print('Binsize ',binsize/3600)
            results_path = utils.results_path(path)
            fname_all_chambers = '%scollective_results_all_chambers_binsize_%f_h.csv'%(prefix, binsize//3600)
            try:
                cages = non_standard_cages[path]
            except KeyError:
                cages = standard_cages
            try:
                headers = non_standard_headers[path]
            except KeyError:
                headers = standard_headers
                
            data = ts.get_time_spent_in_each_chamber(ehs, cf, cages, binsize=binsize)
            data = ts.sum_data(data)
            save_data_cvs(data, fname_all_chambers, results_path, cages, headers)

       
        

            
        dispatch.evaluate_whole_experiment(ehd, cf, directory, prefix,
                                           td.tube_dominance_single_phase,
                                           'tube_dominance',
                                           'dominating mouse',
                                           'pushed out mouse',
                                           '# dominances', args=[None])
