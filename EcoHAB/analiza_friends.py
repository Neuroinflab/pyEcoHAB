import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
import locale
from ExperimentConfigFile import ExperimentConfigFile
from analiza1 import smells, antenna_positions

### How much time mice spend with each other


datarange = slice(0, 10, None)
datasets = [
    '/home/jszmek/EcoHAB_data_November/Maciek_01_30_2018'
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

def interval_overlap(int1, int2):
    """Return overlap between two intervals."""
    ints = sorted([int1, int2], key=lambda x: x[0])
    if ints[1][0] > ints[0][1]:
        return 0
    else:
        return min(ints[0][1], ints[1][1]) - ints[1][0]

def mice_overlap(data_mice, m1, m2, address):
    """Return time overlap of mice m1 and m2 in cage <address>."""
    ints1 = [(s, e) for a, s, e in data_mice[m1] if a == address]
    ints2 = [(s, e) for a, s, e in data_mice[m2] if a == address]
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

if __name__ == '__main__':


    alldeltas = {}
    allratios = {}
    allresults = {}

    # plt.set_cmap('option_d')

    for path in datasets[datarange]:
        alldeltas[path] = {}
        allratios[path] = {}
        allresults[path] = {}
        ehd = EcoHab.EcoHabData(path, _ant_pos=antenna_positions[path],remove_mice=['0065-0136657055'])
        ehs = EcoHab.EcoHabSessions(ehd)

        cf = ExperimentConfigFile(path)
        tstart, tend = cf.gettime('ALL')
        mice = list(ehd.mice)
        mice = filter(lambda x: len(ehs.getstarttimes(x)) > 30, mice)
    

    
        phases = filter(lambda x: x.endswith('dark'), cf.sections())
        # phases = ['SNIFF 1 dark']
        # phases = filter(lambda x: x.endswith('dark') or x.endswith('light'), cf.sections())
        for sec in phases:
            data = prepare_data(ehs, mice, cf.gettime(sec))
        
            results = np.zeros((len(mice), len(mice)))
            results_exp = np.zeros((len(mice), len(mice)))
            for ii in range(len(mice)):
                for jj in range(len(mice)):
                    if ii < jj:
                        print ii, jj
                        res = mice_together(data, mice[ii], mice[jj])
                        results[ii, jj] = res[0]
                        results_exp[ii, jj] = res[1]

            plt.figure(figsize=(10, 6))
            plt.subplot(221)
            plt.imshow(results, vmin=0, vmax=0.5,interpolation='none')
            plt.xticks([])
            plt.yticks([])
            plt.title('% time together')
            plt.subplot(222)
            plt.imshow(results_exp, vmin=0, vmax=0.5,interpolation='none')
            plt.xticks([])
            plt.yticks([])
            plt.title('Expected % time together')
            deltas = results[results > 0] - results_exp[results > 0]
            alldeltas[path][sec] = deltas
            ratios = results[results > 0] / results_exp[results > 0]
            allratios[path][sec] = ratios
            allresults[path][sec] = results[results > 0]
            plt.subplot(223)
            try:
                # plt.imshow(results - results_exp, vmin=-max(np.abs(deltas)),
                #  vmax=max(np.abs(deltas)))
                plt.imshow(results - results_exp, vmin=-0.25, vmax=0.25,interpolation='none')
                plt.xticks([])
                plt.yticks([])
            except ValueError:
                pass
            plt.title('Excess % time together')
            plt.colorbar()
            plt.subplot(224)
            try:
                plt.hist(deltas)
            except ValueError:
                pass
            plt.title('Histogram of excess % time together')
            plt.xlim([-0.1, 0.3])
            plt.suptitle(path)
            
            plt.savefig('friends_%s_%s.pdf' %(path.translate(None, '/'), sec))
            plt.savefig('friends_%s_%s.png' %(path.translate(None, '/'), sec), dpi=300)
            np.savetxt('%s_results_removed_mouse_%s.csv' %(path, sec), results, fmt='%.6f', delimiter=';')
            np.savetxt('%s_results_exp_removed_mouse_%s.csv' %(path, sec), results_exp,
                    fmt='%.6f', delimiter=';')
            np.savetxt('%s_results_final_removed_mouse_%s.csv' %(path, sec), results-results_exp, 
                    fmt='%.6f', delimiter=';')
        # np.save('alldeltas_EH_new_%s.npy' %path, alldeltas[path])
        # np.save('allratios_EH_new_%s.npy' %path, allratios[path])
        # np.save('allesults_EH_new_%s.npy' %path, allresults[path])
            

    # np.save('alldeltas.npy', alldeltas)
    # np.save('allratios.npy', allratios)
    # np.save('allesults.npy', allresults)
