# -*- coding: utf-8 -*-
from __future__ import print_function,division
import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
import locale
from ExperimentConfigFile import ExperimentConfigFile
from analiza1 import smells, antenna_positions
import os
import utils
import plotfunctions
from write_to_file import save_single_histograms, write_csv_rasters, write_csv_tables, write_csv_alone
### How much time mice spend with each other

datarange = slice(0, 10, None)
datasets = [
    # "/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/Social structure males 02.03/",
    
    #   '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_dominance_swiss_webster_dominant_remove_12.02.18',
    #   '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_structure_16.01',
    #   '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_structure_19.01.18_rep_II',
    #   '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_structure_swiss_webster_ctrl_05.02.18',
    #"/home/jszmek/EcoHAB_data_November/mice K Wisniewska",
    #'/home/jszmek/EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/',
    #'/home/jszmek/EcoHAB_data_November/C57 13-24.04 long/',
    #"/home/jszmek/EcoHAB_data_November/C57 males long 11-22.05.18/",

    # "/home/jszmek/EcoHAB_data_November/C57 TIMP rep 2/",
    # "/home/jszmek/EcoHAB_data_November/C57 males rep 2/",
    # "/home/jszmek/EcoHAB_data_November/C57 males TIMP/",
    # "/home/jszmek/EcoHAB_data_November/BTBR males/",
    #"/home/jszmek/EcoHAB_data_November/long_experiment_WT",
    "/home/jszmek/EcoHAB_data_November/long_experiment_WT",
    "/home/jszmek/EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark",
    "/home/jszmek/EcoHAB_data_November/long_experiment_KO_from_phase_SNIFF_10_dark",
    ]

# # address = {1: 4, 2: 1, 3: 1, 4: 2, 5: 2, 6: 3, 7: 3, 8: 4}
# smells = {  'C57 CTRL EH (3) 25.07.14': {'soc': 1, 'nsoc': 3},
#             'C57 CTRL EH (4) 30.07.14': {'soc': 3, 'nsoc': 1},
#             'C57 CTRL EH (5) 20.11.14': {'soc': 3, 'nsoc': 1},
#             }
# antenna_positions = {'C57 CTRL EH (3) 25.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
#             'C57 CTRL EH (4) 30.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
#             'C57 CTRL EH (5) 20.11.14': None,
#             }

remove_tags = {
    "/home/jszmek/EcoHAB_data_November/C57 males rep 2/":["0065-0161984735"],
    "/home/jszmek/EcoHAB_data_November/BTBR males/":["0065-0136658439",
                                                             "0065-0141855614"]
}
how_many_appearances = {
    "/home/jszmek/EcoHAB_data_November/C57 males rep 2/":1000,
    "/home/jszmek/EcoHAB_data_November/BTBR males/":500,
    '/home/jszmek/EcoHAB_data_November/C57 30.04-11.05 LONG TIMP/':200
}
antenna_positions = {
    "/home/jszmek/EcoHAB_data_November/long_experiment_KO_mismatched_antennas_to_phase_SNIFF_10_dark":{'1': 1,
                                                                                                       '2': 5,
                                                                                                       '3': 3,
                                                                                                       '4': 6,
                                                                                                       '5': 4,
                                                                                                       '6': 2,
                                                                                                       '7': 7,
                                                                                                       '8': 8}}
def intervals(data_mice, mouse, address):
    return [[s, e] for a, s, e in data_mice[mouse] if a == address]

def mouse_alone(data_mice, address):
    ints = {}
    empty = 0
    for mouse in data_mice.keys():
        ints[mouse] = intervals(data_mice, mouse, address)
        if not len(ints[mouse]):
            del ints[mouse]
    result = {}
    for mouse in data_mice.keys():
        result[mouse] = 0
    if ints == {}:
        return result
    for mouse in ints.keys():
        other_mice = ints.keys()
        other_mice.remove(mouse)
        i = 0

        while True:
            for other_mouse in other_mice:
                j = 0
                remove = False
                if not len(ints[other_mouse]):
                    continue
                while True:
                    original_s, original_e = ints[mouse][i]
                    other_s, other_e = ints[other_mouse][j]
                    if other_e <= original_s:
                        j = j+1
                    elif other_s >= original_e:
                        break
                    elif original_s <= other_s:
                        ints[mouse][i][1] = other_s
                        if original_e >= other_e:
                            ints[other_mouse].remove([other_s, other_e])
                            ints[mouse].insert(i+1,[other_e, original_e])
                            break
                        else:
                            ints[other_mouse][j][0] = original_e
                            j = j+1
                    else:
                        ints[other_mouse][j][1] = original_s
                        if original_e <= other_e:
                            ints[mouse].remove([original_s, original_e])
                            remove = True
                            ints[other_mouse].insert(j+1, [original_e, other_e])
                            break
                        else:
                            ints[mouse][i][0] = other_e
                            j = j+1
                    if j >= len(ints[other_mouse]):
                        break
                if remove:
                    break
            if not remove:
                i = i+1
            if i >= len(ints[mouse]):
                break
        
   
    for mouse in ints.keys():
        result[mouse] = sum([e-s for s,e in ints[mouse]])
        
    return result

def interval_overlap(int1, int2):
    """Return overlap between two intervals."""
    ints = sorted([int1, int2], key=lambda x: x[0])
    if ints[1][0] > ints[0][1]:
        return 0
    else:
        return min(ints[0][1], ints[1][1]) - ints[1][0]

def mice_overlap(data_mice, m1, m2, address):
    """Return time overlap of mice m1 and m2 in cage <address>."""
    ints1 = intervals(data_mice, m1,address)
    ints2 = intervals(data_mice, m2,address)
    durs1 = [x[1] - x[0] for x in ints1]
    durs2 = [x[1] - x[0] for x in ints2]
    total = 0.
    for int1 in ints1:
        for int2 in ints2:
            total += interval_overlap(int1, int2)
    return total, sum(durs1), sum(durs2)
    
def mice_together(data_mice, m1, m2, total_time=43200.):
    """Return the time spent together by two mice and expected time 
    assuming independence."""
    result = np.zeros((4, 3))
    for address in [1, 2, 3, 4]:
        result[address-1] = mice_overlap(data_mice, m1, m2, address)
    fracs = result[:, 1:] / total_time
    time_together = result[:, 0].sum()
    exp_time_together = (fracs[:, 0] * fracs[:, 1]).sum()
    return time_together/total_time, exp_time_together

def prepare_data(ehs, mice, times, margin=12*3600.): 
    """Prepare masked data."""
    t1, t2 = times
    ehs.mask_data(t1 - margin, t2)
    data = {}
    for mm in mice:
        ads = ehs.getaddresses(mm)
        sts = ehs.getstarttimes(mm)
        ens = ehs.getendtimes(mm)
        data_mm = []
        for ad, st, en in zip(ads, sts, ens):
            if en > t1:
                data_mm.append((ad, max(st, t1), min(en, t2)))
        data[mm] = data_mm
    return data

def calculate_total_time(intervals):
    return sum([e-s for s, e in intervals])

def total_time_results(mice_data, mice):
    result = np.zeros((4, len(mice)))
    for address in [1, 2, 3, 4]:
        for i,mouse in enumerate(mice):
            ints = intervals(mice_data, mouse, address)
            result[address-1,i] = calculate_total_time(ints)
    return result

def test_results(data_mice,mice):
    #For surrogate data
    total_time_spent = total_time_results(data_mice,mice)
    alone = np.zeros((4,len(mice)))
    time_together = np.zeros((4,len(mice),len(mice)))
    for i,address in enumerate([1, 2, 3, 4]):
        mouse_alone_dictionary = mouse_alone(data_mice,address)
        for j, mouse in enumerate(mice):
            alone[i,j] = mouse_alone_dictionary[mouse]
            for k, other_mouse in enumerate(mice[j+1:]):
                time_together[i,j,j+1+k] = mice_overlap(data_mice, mouse, other_mouse, address)[0]
                time_together[i,j+1+k,j] = time_together[i,j,j+1+k]

def mouse_alone_ehs(ehs, cf, main_directory, prefix):
    phases = filter(lambda x: x.endswith('dark') or x.endswith('DARK'), cf.sections())
    mice = ehs.mice
    output = np.zeros((4, len(mice), len(phases)+1))
    for phase, sec in enumerate(phases):
        data = prepare_data(ehs, mice, cf.gettime(sec))
        for i in range(1,5):
            alone = mouse_alone(data, i)
            for j, mouse in enumerate(mice):
                output[i-1, j, phase] = alone[mouse]
    phases.append('ALL DARK')
    output[:,:,-1] = output[:,:,:-1].sum(axis=2)  # last column -- sum of activity in all dark phases
    write_csv_alone(output, phases, mice, main_directory, prefix)
       
def in_cohort_sociability(ehs, cf, main_directory, prefix, remove_mouse=None):

    mice = ehs.mice

    phases = filter(lambda x: x.endswith('dark') or x.endswith('DARK'), cf.sections())
        
    full_results = np.zeros((len(phases), len(mice), len(mice)))
    full_results_exp = np.zeros((len(phases), len(mice), len(mice)))
    if remove_mouse:
        fname = 'incohort_sociability_remove_%s' % remove_mouse
        name_ = 'incohort_sociability_measured_time_%sremove_%s.csv' % (prefix, remove_mouse)
        name_exp_ = 'incohort_sociability_excess_time_%sremove_%s.csv' % (prefix, remove_mouse)
    else:
        fname = 'incohort_sociability'
        name_ = 'incohort_sociability_measured_time_%s.csv' % prefix
        name_exp_ = 'incohort_sociability_excess_time_%s.csv' % prefix

    for idx_phase, phase in enumerate(phases):
        print(phase)
        data = prepare_data(ehs, mice, cf.gettime(phase))
        results = np.zeros((len(mice), len(mice)))
        results_exp = np.zeros((len(mice), len(mice)))

        for ii in range(len(mice)):
            for jj in range(len(mice)):
                if ii < jj:
                    res = mice_together(data, mice[ii], mice[jj])
                    results[ii, jj] = res[0]
                    results_exp[ii, jj] = res[1]
                    
        full_results[idx_phase] = results
        full_results_exp[idx_phase] = results_exp
        
        deltas = results[results > 0] - results_exp[results > 0]
        save_single_histograms(results,
                               'incohort_sociability_measured_time',
                               mice,
                               phase,
                               main_directory,
                               'in_cohort_sociability/histograms',
                               prefix,
                               additional_info=remove_mouse)
        save_single_histograms(results_exp,
                               'incohort_sociability_expected_time',
                               mice,
                               phase,
                               main_directory,
                               'in_cohort_sociability/histograms',
                               prefix,
                               additional_info=remove_mouse)
        save_single_histograms(results-results_exp,
                               'incohort_sociability_excess_time',
                               mice,
                               phase,
                               main_directory,
                               'in_cohort_sociability/histograms',
                               prefix,
                               additional_info=remove_mouse)
        
        single_in_cohort_soc_plot(results,
                              results_exp,
                              mice,
                              phase,
                              fname,
                              main_directory,
                                  'in_cohort_sociability/histograms',
                              prefix)

    write_csv_rasters(mice,
                      phases,
                      full_results,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_)
    write_csv_rasters(mice,
                      phases,
                      full_results-full_results_exp,
                      main_directory,
                      'in_cohort_sociability/raster_plots',
                      name_exp_)
    

    plotfunctions.make_RasterPlot(main_directory,
                                  'in_cohort_sociability/raster_plots',
                                  full_results,
                                  phases,
                                  name_,
                                  mice,
                                  vmin=0,
                                  vmax=0.5,
                                  title='% time together')
    plotfunctions.make_RasterPlot(main_directory,
                                  'in_cohort_sociability/raster_plots',
                                  full_results-full_results_exp,
                                  phases,
                                  name_exp_,
                                  mice,
                                  title='% time together',
                                  vmin=-.25,
                                  vmax=.25)

def single_in_cohort_soc_plot(results,
                              results_exp,
                              mice,
                              phase,
                              fname,
                              main_directory,
                              directory,
                              prefix):
    new_name = os.path.join(directory, 'figs')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s_%s'% (fname, prefix, phase))
    label_mice = [mouse[-4:] for mouse in mice]
    fig = plt.figure(figsize=(10, 6))
    ax = []
    for i in range(1,5):
        ax.append(fig.add_subplot(2,2,i))

    ax[0].imshow(results, vmin=0, vmax=0.5, interpolation='none')
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_title('% time together')
    ax[1].imshow(results_exp, vmin=0, vmax=0.5, interpolation='none')
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_title('Expected % time together')

    deltas = results[results > 0] - results_exp[results > 0]

    plt.subplot(223)
    try:
        im = ax[2].imshow(results - results_exp, vmin=-0.25, vmax=0.25,interpolation='none')
        ax[2].set_title('Excess % time together')
        cbar = fig.colorbar(im)
        ax[2].yaxis.tick_left()
        ax[2].get_yaxis().set_ticks([i for i,x in enumerate(mice)])
        ax[2].get_xaxis().set_ticks([i for i,x in enumerate(mice)])
        ax[2].set_xticklabels(label_mice)
        ax[2].set_yticklabels(label_mice)
        for label in ax[2].xaxis.get_ticklabels():
            label.set_rotation(90)

    except ValueError:
        pass

    try:
        ax[3].hist(deltas)
    except ValueError:
        pass
    ax[3].set_title('Histogram of excess % time together')
    ax[3].set_xlim([-0.1, 0.3])
    ax[3].get_xaxis().set_ticks([-0.1, 0., 0.1, 0.2, 0.3])
    ax[3].set_xticklabels([-0.1, 0., 0.1, 0.2, 0.3])
    fig.suptitle(phase)
    fig.subplots_adjust(left=0.3)
    fig.subplots_adjust(bottom=0.3)
    fig.subplots_adjust(wspace=0.25)
    fig.subplots_adjust(hspace=0.3)
    
    fig.savefig(fname+'.pdf')
    fig.savefig(fname+'.png', dpi=300)
            


if __name__ == '__main__':

    for path in datasets[datarange]:
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
        directory = utils.results_path(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        cf = ExperimentConfigFile(path)      
        mouse_alone_ehs(ehs, cf, directory, prefix)
        in_cohort_sociability(ehs, cf, directory, prefix, remove_mouse=remove_mouse)
  
   
