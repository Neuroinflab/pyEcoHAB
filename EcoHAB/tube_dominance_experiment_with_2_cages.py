from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utility_functions as utils
import numpy as np
from write_to_file import save_single_histograms, write_csv_rasters
from plotfunctions import single_in_cohort_soc_plot, make_RasterPlot, single_heat_map
from numba import jit
import matplotlib.pyplot as plt
mouse_attention_span = 10  # sec
nbins = 10
homepath = os.path.expanduser("~/")
opposite_antenna_dict = {1: 2,
                    2: 1}
dt = 0.05
color_list = ['indianred', 'darkred', 'salmon',
              'darkolivegreen', 'forestgreen', 'darkslategrey',
              'darkslateblue', 'lightskyblue', 'navy',
              'mediumpurple', 'grey', 'gold',
              'orange', 'saddlebrown', 'magenta']
def get_states(ehd, cf, mouse, dt = 0.05):
    """
    0 -- home cage, 1 -- pipe, 2 -- cage with stimulus
    """
    t_start, t_end = cf.gettime('ALL')
    length = int((t_end - t_start)/dt)
    times, antennas = utils.get_times_antennas(ehd, mouse1,
                                               t_start, t_end)
    home_antenna = antennas[0]
    states = np.zeros((length, 1), dtype=int)
    for i, a in enumerate(antennas[:-1]):
        timestamp = int((times[i] - t_start)/dt)
        next_a = antennas[i + 1]
        next_timestamp = int((times[i+1] - t_start)/dt)
        if a != next_a:  # easy, the mouse is crossing the pipe
            states[timestamp:next_timestamp] = 1
        if a == next_a and a != home_antenna:
            states[timestamp:next_timestamp] = 2
    return states


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

def check_mouse1_defending(m1_states, m1_times, antennas2, times2, home_cage_antenna):
    return 0


if __name__ == '__main__':
    datasets = [#'EcoHAB_data_November/social_structure_swiss_webster_ctrl_05.02.18',]
                 'EcoHAB_data_November/Social structure males 02.03/',
                # 'EcoHAB_data_November/social_dominance_swiss_webster_dominant_remove_12.02.18',
                # 'EcoHAB_data_November/social_structure_16.01',
                # 'EcoHAB_data_November/social_structure_19.01.18_rep_II',
                 ]
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
            ehd1 = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd1 = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        prefix = utils.make_prefix(path)
        res_dir = utils.results_path(path)
        cf1 = ExperimentConfigFile(path)
        home_cage_antenna = ehd1.get_home_cage_antenna()
        states = {}
        for mouse1 in ehd1.mice:
            states[mouse1] = get_states(ehd1, cf1, mouse1, dt)
        T_START, T_END = cf1.gettime('ALL')
        for phase in cf1.sections():
            if phase == 'ALL':
                continue
            t_start, t_end = cf1.gettime(phase)
            print(phase, t_start, t_end)
                        
            idx_start = int((t_start-T_START)/dt)
            idx_end =  int((t_end-T_START)/dt)
            time = np.arange(t_start, t_end, dt)
            ehd1.mask_data(t_start, t_end)
            if phase != 'ALL':
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title(phase)
                for i, mouse in enumerate(ehd1.mice):
                    print(mouse)
                    print(ehd1.getantennas(mouse))

                    ax.plot(time/1000, 10*states[mouse][idx_start:idx_end]+0.1*i, color=color_list[i], label=mouse)
                ax.legend()
                ax.set_xlabel(time)
            ehd1.unmask_data()

        utils.evaluate_whole_experiment(ehd1, cf1, res_dir, prefix,
                                        tube_dominance_2_cages,
                                        'tube_dominance',
                                        'domineering mouse',
                                        'pushed out mouse',
                                        '# dominances',
                                        args=[home_cage_antenna])
