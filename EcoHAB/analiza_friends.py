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
### How much time mice spend with each other


datarange = slice(0, 10, None)
datasets = [
    "/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/Social structure males 02.03/",
    #' /home/jszmek/EcoHAB_data_November/Maciek_01_30_2018',
    # '/home/jszmek/EcoHAB_data_November/Maciek_social_structure_16.01',
    # '/home/jszmek/EcoHAB_data_November/Maciek_social_structure_19.01.18_rep_II',
     # '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_dominance_swiss_webster_dominant_remove_12.02.18',
     # '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_structure_16.01',
     # '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_structure_19.01.18_rep_II',
     # '/home/jszmek/Results_EcoHAB_data_November/do_analizy_in_z_cohort_z_sociability_z_numerami_transponderow/social_structure_swiss_webster_ctrl_05.02.18',

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
def intervals(data_mice, mouse,address):
    return [[s, e] for a, s, e in data_mice[mouse] if a == address]

def mouse_alone(data_mice,address):
    ints = {}
    empty = 0
    for mouse in data_mice.keys():
        ints[mouse] = intervals(data_mice,mouse,address)
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
                            ints[other_mouse].remove([other_s,other_e])
                            ints[mouse].insert(i+1,[other_e, original_e])
                            break
                        else:
                            ints[other_mouse][j][0] = original_e
                            j = j+1
                    else:
                        ints[other_mouse][j][1] = original_s
                        if original_e <= other_e:
                            ints[mouse].remove([original_s,original_e])
                            remove = True
                            ints[other_mouse].insert(j+1,[original_e, other_e])
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
    return sum([e-s for s,e in intervals])

def total_time_results(mice_data,mice):
    result = np.zeros((4,len(mice)))
    for address in [1,2,3,4]:
        for i,mouse in enumerate(mice):
            ints = intervals(mice_data,mouse,address)
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
            print(alone[i,j],total_time_spent[i,j])
def write_cvs_alone(alone,address,directory,phase):
    fname =  os.path.join(directory,'Mice_alone_%d_%s.csv'%(address,phase))
    f = open(fname,'w')
    for mouse in alone.keys():
        f.write(mouse+';'+str(alone[mouse])+'\n')
    f.close()


if __name__ == '__main__':


    alldeltas = {}
    allratios = {}
    allresults = {}

    # plt.set_cmap('option_d')
    for remove_mouse in [None,'0065-013665705']:
        for path in datasets[datarange]:
            alldeltas[path] = {}
            allratios[path] = {}
            allresults[path] = {}
            if path not in antenna_positions:
                antenna_positions[path] = None
            if remove_mouse:
                ehd = EcoHab.EcoHabData(path=path, _ant_pos=antenna_positions[path],remove_mice=[remove_mouse])
            else:
                ehd = EcoHab.EcoHabData(path=path, _ant_pos=antenna_positions[path])

            ehs = EcoHab.EcoHabSessions(ehd)

            cf = ExperimentConfigFile(path)
            tstart, tend = cf.gettime('ALL')
            mice = list(ehd.mice)
            mice = filter(lambda x: len(ehs.getstarttimes(x)) > 30, mice)
    

    
            phases = filter(lambda x: x.endswith('dark'), cf.sections())
            # phases = ['SNIFF 1 dark']
            # phases = filter(lambda x: x.endswith('dark') or x.endswith('light'), cf.sections())
            directory = utils.results_path(path)
            if not os.path.exists(directory):
                os.makedirs(directory)
            for sec in phases:
                data = prepare_data(ehs, mice, cf.gettime(sec))

                results = np.zeros((len(mice), len(mice)))
                results_exp = np.zeros((len(mice), len(mice)))
                for ii in range(len(mice)):
                    for jj in range(len(mice)):
                        if ii < jj:
                            res = mice_together(data, mice[ii], mice[jj])
                            results[ii, jj] = res[0]
                            results_exp[ii, jj] = res[1]
                fig = plt.figure(figsize=(10, 6))
                ax = []
                single_times = {}
                for i in range(1,5):
                    ax.append(fig.add_subplot(2,2,i))
                    alone = mouse_alone(data,i)
                    write_cvs_alone(alone,i,directory,sec)
     
                test_results(data,mice)
                ax[0].imshow(results, vmin=0, vmax=0.5,interpolation='none')
                ax[0].set_xticks([])
                ax[0].set_yticks([])
                ax[0].set_title('% time together')
                ax[1].imshow(results_exp, vmin=0, vmax=0.5,interpolation='none')
                ax[1].set_xticks([])
                ax[1].set_yticks([])
                ax[1].set_title('Expected % time together')
                deltas = results[results > 0] - results_exp[results > 0]
                alldeltas[path][sec] = deltas
                ratios = results[results > 0] / results_exp[results > 0]
                allratios[path][sec] = ratios
                allresults[path][sec] = results[results > 0]
                plt.subplot(223)
                try:
                    # plt.imshow(results - results_exp, vmin=-max(np.abs(deltas)),
                    #  vmax=max(np.abs(deltas)))
                    im = ax[2].imshow(results - results_exp, vmin=-0.25, vmax=0.25,interpolation='none')
                    ax[2].set_title('Excess % time together')
                    cbar = fig.colorbar(im)
                    ax[2].yaxis.tick_left()
                    ax[2].get_yaxis().set_ticks([i for i,x in enumerate(mice)])
                    ax[2].get_xaxis().set_ticks([i for i,x in enumerate(mice)])
                    ax[2].set_xticklabels(mice)
                    ax[2].set_yticklabels(mice)
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
                ax[3].get_xaxis().set_ticks([-0.1,0.,0.1,0.2,0.3])
                ax[3].set_xticklabels([-0.1,0.,0.1,0.2,0.3])
                
                
                fig.suptitle(sec)
                fig.subplots_adjust(left=0.3)
                fig.subplots_adjust(bottom=0.3)
                fig.subplots_adjust(wspace=0.25)
                fig.subplots_adjust(hspace=0.3)
                
                header = ''
                for mouse in mice:
                    header+= mouse+';'
                fig.savefig('%s/friends_%s_remove_%s.pdf' %(directory, sec,remove_mouse))
                fig.savefig('%s/friends_%s_remove_%s.png' %(directory, sec,remove_mouse), dpi=300)
                np.savetxt('%s/results_removed_mouse_%s_remove_%s.csv' %(directory, sec, remove_mouse), results, fmt='%.6f', delimiter=';',header=header,comments='')
                np.savetxt('%s/results_exp_removed_mouse_%s_remove_%s.csv' %(directory, sec, remove_mouse), results_exp,
                           fmt='%.6f', delimiter=';',header=header,comments='')
                np.savetxt('%s/results_final_removed_mouse_%s_remove_%s.csv' %(directory, sec, remove_mouse), results-results_exp, fmt='%.6f', delimiter=';',header=header,comments='')
