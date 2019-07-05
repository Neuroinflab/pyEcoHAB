# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import os
from EcoHab import * 
from .data_info import *

homepath = os.path.expanduser("~/")
threshold = 3
if __name__ == '__main__':
    

    for new_path in datasets:
       
        path = os.path.join(homepath, new_path)
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
            ehd = EcoHabData(path=path,
                             _ant_pos=antenna_positions[new_path],
                             remove_mice=remove_mouse,
                             how_many_appearances=how_many_appearances[new_path],
                             remove_antennas=remove_antennas)
        else:
            ehd = EcoHabData(path=path,
                             _ant_pos=antenna_positions[new_path],
                             how_many_appearances=how_many_appearances[new_path],
                             remove_antennas=remove_antennas)

        ehs = EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path)
        get_mouse_alone(ehs, cf)
        get_in_cohort_sociability(ehs, cf, remove_mouse=remove_mouse)
        get_in_cohort_sociability(ehs, cf,
                                  which_phases="ALL",
                                  remove_mouse=remove_mouse)
        get_in_cohort_sociability(ehs, cf,
                                  which_phases="dark",
                                  remove_mouse=remove_mouse)
        get_in_cohort_sociability(ehs, cf,
                                  which_phases="light",
                                  remove_mouse=remove_mouse)
        get_following(ehd, cf, remove_mouse=remove_mouse)
        for binsize in binsizes:
            get_visits(ehs, cf, binsize)
        get_tube_dominance(ehd, cf)
