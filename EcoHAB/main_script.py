from __future__ import division, print_function

import numpy as np
import os
import analiza_friends as af
import time_spent_in_cage as ts
import utils
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
from write_to_file import save_data_cvs

if __name__ == '__main__':
    

    for path in datasets:
        prefix = utils.make_prefix(path)
        if path in remove_tags:
            remove_mouse = remove_tags[path]
        else:
            remove_mouse = None
        if path not in antenna_positions:
            antenna_positions[path] = None
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[path])

        ehs = EcoHab.EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path)
        
        directory = utils.results_path(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        tstart, tend = cf.gettime('ALL')
        for binsize in binsizes:
            print('Binsize ',binsize/3600)
            results_path = utils.results_path(path)
            fname_all_chambers = 'collective_results_all_chambers_binsize_%f_h.csv'%(binsize//3600)
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

        af.mouse_alone_ehs(ehs, cf, directory, prefix)
        af.in_cohort_sociability(ehs, cf, directory, prefix, remove_mouse=remove_mouse)

        
