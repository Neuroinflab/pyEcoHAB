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
import interactions

homepath = os.path.expanduser("~/")
threshold = 3
if __name__ == '__main__':
    remove_tags = []
    #path1 = 'EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark'
    #path2 = 'EcoHAB_data_November/long_experiment_KO_from_phase_SNIFF_10_dark'
    path = 'EcoHAB_data_November/long_experiment_KO'
    antenna_pos = {
        "EcoHAB_data_November/long_experiment_KO":{'1':1,'2':5,'3':3,'4':6,'5':4,'6':2,'7':7,'8':8}
    }
    new_path = os.path.join(homepath, path)
    remove_mouse = []
    prefix = utils.make_prefix(path)

    ehd = EcoHab.EcoHabData(path=new_path,_ant_pos=antenna_pos[path],
                             which_phase="WRONG_ANTENNAS",
                             how_many_appearances=500)
    ehd2 =  EcoHab.EcoHabData(path=new_path,_ant_pos=None,
                             which_phase="CORRECT_ANTENNAS",
                              how_many_appearances=500)
    ehd.merge_experiment(ehd2)
    ehs = EcoHab.EcoHabSessions(ehd)
    cf = ExperimentConfigFile(new_path)
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

       
        
        directory = utils.results_path(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        af.mouse_alone_ehs(ehs, cf, directory, prefix)
        af.in_cohort_sociability(ehs, cf, directory, prefix, remove_mouse=remove_mouse)

            
