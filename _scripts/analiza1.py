import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
# import IntelliCage_tools as ict
import locale
from ExperimentConfigFile import ExperimentConfigFile

### How much time mice spend in the 'social' compartment

datarange = slice(10, 11, None)
datasets = ['BALB EH (1) 06.03.14', # 0
            'BALB EH (2) po raz 2 10.04.14',
            'C57 EH (1) 23.04.14',
            'FX KO F (3) EH (1) 15.05.14',
            'BALB VPA EH (1) 05.05.14', # ***
            'BALB VPA EH (2) 09.05.14', # 5
            'C57 CTRL EH (2) 28.04.14',
            'BTBR 1st time - stary EH',
            'EH BTBR po raz 2', #8
            'BALB CTRL EH (3) 12.06.14 - 15.06.14',
            'BALB VPA EH (4) 16.06.14 - 19.06.14', # 10 #
            'C57 VPA EH (1) 20.06.14',
            'C57 CTRL EH (3) 25.07.14', # 12
            'C57 CTRL EH (4) 30.07.14',
            'FX KO males (1) 08.08.14',
            'C57 VPA EH (1Biol) 28.10.2014', # 15, nowy Eco-HAB
            'BALB CTRL EH (4) 09.09.14',
            'BALB CTRL EH (5) 15.09.14',
            'FX WT males EH (1) 24.09.14',
            'FX KO females EH 2 29.09.14',
            'FX KO males EH (2) 19.09.14', # 20
            'C57 CTRL EH (5) 20.11.14',
            'FX KO males EH (3) 14.01.15', # 22
            'FX WT males EH (2) 07.01.15', # 23
            'FX WT females EH (1) 16.02.15', # 24
            'FX WT females EH (2) 20.02.15',
            'FX KO females EH (3) 11.02.15',
            'BTBR females (2) EH - 13.03.2015',
            'BALB VPA EH Biol (1) 23.03.15',
            'C57 VPA EH Biol (2) 27.03.15',
            'BTBR males EH (2) 31.03.15', # 30
            'C57 NaCl Biol (1) 26.06.15',
            'C57 NaCl Biol (2) 01.07.15',
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep1',
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep2',
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep1',
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep2',
            'BALB NaCl (1) Rep (1) 09.09.15',
            'BALB NaCl (2) Rep (1) 02.10.15',
            'BALB NaCl (2) Rep (2) 11.10.15', # 39
            'FX WT females EH Nen (1) 07.12.15', 
            'FAMB KO males (1) rep 2 26.02.16',
            'FAMB WT males (1) rep 1 22.02.16',
            'FAMB WT males (1) rep 2 02.03.16',
            'FX WT females EH (B) 05.05.16', # 44
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16',
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16',
            'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16',
            'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16',
            'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16',
            'PV Cre males HT EH (1) 09.06.16', # 50
            ]
# address = {1: 4, 2: 1, 3: 1, 4: 2, 5: 2, 6: 3, 7: 3, 8: 4}
smells = {'BALB EH (1) 06.03.14': {'soc': 3, 'nsoc': 1},
            'BALB EH (2) po raz 2 10.04.14': {'soc': 3, 'nsoc': 1},
            'C57 EH (1) 23.04.14': {'soc': 3, 'nsoc': 1},
            'FX KO F (3) EH (1) 15.05.14': {'soc': 3, 'nsoc': 1},
            'BALB VPA EH (1) 05.05.14': {'soc': 4 , 'nsoc': 2}, 
            'BALB VPA EH (2) 09.05.14': {'soc': 4, 'nsoc': 2},
            'C57 CTRL EH (2) 28.04.14': {'soc': 4, 'nsoc': 2},
            'BTBR 1st time - stary EH': {'soc': 3, 'nsoc': 1},
            'EH BTBR po raz 2': {'soc': 3, 'nsoc': 1},
            'BALB CTRL EH (3) 12.06.14 - 15.06.14': {'soc': 1, 'nsoc': 2}, # soc miedzy 2, 3; nsoc 4, 5
            'BALB VPA EH (4) 16.06.14 - 19.06.14': {'soc': 3, 'nsoc': 1}, # soc miedzy 6, 7; nsoc 2, 3
            'C57 VPA EH (1) 20.06.14': {'soc': 3, 'nsoc': 4},
            'C57 CTRL EH (3) 25.07.14': {'soc': 1, 'nsoc': 3},
            'C57 CTRL EH (4) 30.07.14': {'soc': 3, 'nsoc': 1},
            'FX KO males (1) 08.08.14': {'soc': 1, 'nsoc': 3},
            'C57 VPA EH (1Biol) 28.10.2014': {'soc': 1, 'nsoc': 3},
            'BALB CTRL EH (4) 09.09.14': {'soc': 1, 'nsoc': 3},
            'BALB CTRL EH (5) 15.09.14': {'soc': 3, 'nsoc': 1},
            'FX WT males EH (1) 24.09.14': {'soc': 1, 'nsoc': 3},
            'FX KO females EH 2 29.09.14': {'soc': 3, 'nsoc': 1},
            'FX KO males EH (2) 19.09.14': {'soc': 3, 'nsoc': 1},
            'C57 CTRL EH (5) 20.11.14': {'soc': 3, 'nsoc': 1},
            'FX KO males EH (3) 14.01.15': {'soc': 1, 'nsoc': 2},
            'FX WT males EH (2) 07.01.15': {'soc': 1, 'nsoc': 2}, # 23
            'FX WT females EH (1) 16.02.15': {'soc': 3, 'nsoc': 1}, # 24
            'FX WT females EH (2) 20.02.15': {'soc': 1, 'nsoc': 3},
            'FX KO females EH (3) 11.02.15': {'soc': 1, 'nsoc': 3},
            'BTBR females (2) EH - 13.03.2015': {'soc': 1, 'nsoc': 3},
            'BALB VPA EH Biol (1) 23.03.15': {'soc': 3, 'nsoc': 1},
            'C57 VPA EH Biol (2) 27.03.15': {'soc': 3, 'nsoc': 1},
            'BTBR males EH (2) 31.03.15': {'soc': 3, 'nsoc': 1},
            'C57 NaCl Biol (1) 26.06.15': {'soc': 1, 'nsoc': 3},
            'C57 NaCl Biol (2) 01.07.15': {'soc': 3, 'nsoc': 1},
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep1': {'soc': 3, 'nsoc': 1},
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep2': {'soc': 1, 'nsoc': 3},
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep1': {'soc': 1, 'nsoc': 3},
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep2': {'soc': 3, 'nsoc': 1},
            'BALB NaCl (1) Rep (1) 09.09.15': {'soc': 3, 'nsoc': 1},
            'BALB NaCl (2) Rep (1) 02.10.15': {'soc': 1, 'nsoc': 3},
            'BALB NaCl (2) Rep (2) 11.10.15': {'soc': 1, 'nsoc': 3},
            'FX WT females EH Nen (1) 07.12.15': {'soc': 1, 'nsoc': 3},
            'FAMB KO males (1) rep 2 26.02.16': {'soc': 1, 'nsoc': 3},
            'FAMB WT males (1) rep 1 22.02.16': {'soc': 1, 'nsoc': 3},
            'FAMB WT males (1) rep 2 02.03.16': {'soc': 3, 'nsoc': 1},
            'FX WT females EH (B) 05.05.16': {'soc': 3, 'nsoc': 1},
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16': {'soc': 3, 'nsoc': 1},
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16': {'soc': 1, 'nsoc': 3},
            'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16': {'soc': 3, 'nsoc': 1},
            'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16': {'soc': 3, 'nsoc': 1},
            'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16': {'soc': 1, 'nsoc': 3},
            'PV Cre males HT EH (1) 09.06.16': {'soc': 3, 'nsoc': 1},
            }
antenna_positions = {'BALB EH (1) 06.03.14': None,
            'BALB EH (2) po raz 2 10.04.14': None,
            'C57 EH (1) 23.04.14': {'6': 1, '5': 2, '7': 3, '8': 4, '2': 5, '1': 6, '4': 7, '3': 8},
            'FX KO F (3) EH (1) 15.05.14': None,
            'BALB VPA EH (1) 05.05.14': None,
            'BALB VPA EH (2) 09.05.14': None,
            'C57 CTRL EH (2) 28.04.14': None, 
            'BTBR 1st time - stary EH': None,
            'EH BTBR po raz 2': None,
            'BALB CTRL EH (3) 12.06.14 - 15.06.14': None,
            'BALB VPA EH (4) 16.06.14 - 19.06.14': None,
            'C57 VPA EH (1) 20.06.14': None,
            'C57 CTRL EH (3) 25.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
            'C57 CTRL EH (4) 30.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
            'FX KO males (1) 08.08.14': {'5': 1, '6': 2, '7': 3, '8': 4, '1': 5, '2': 6, '3': 7, '4': 8},
            'C57 VPA EH (1Biol) 28.10.2014': None,
            'BALB CTRL EH (4) 09.09.14': None,
            'BALB CTRL EH (5) 15.09.14': None,
            'FX WT males EH (1) 24.09.14': None,
            'FX KO females EH 2 29.09.14': None,
            'FX KO males EH (2) 19.09.14': None,
            'C57 CTRL EH (5) 20.11.14': None,
            'FX KO males EH (3) 14.01.15': None,
            'FX WT males EH (2) 07.01.15': None,
            'FX WT females EH (1) 16.02.15': None, # 24
            'FX WT females EH (2) 20.02.15': None,
            'FX KO females EH (3) 11.02.15': None,
            'BTBR females (2) EH - 13.03.2015': None,
            'BALB VPA EH Biol (1) 23.03.15': None,
            'C57 VPA EH Biol (2) 27.03.15': None,
            'BTBR males EH (2) 31.03.15': None,
            'C57 NaCl Biol (1) 26.06.15': None,
            'C57 NaCl Biol (2) 01.07.15': None,
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep1': None,
            'C57 CTRL Rep (1) 03.08.15 i 12.08.15/rep2': None,
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep1': None,
            'C57 CTRL Rep (2) 07.08.15 i 16.08.15/rep2': None,
            'BALB NaCl (1) Rep (1) 09.09.15': None,
            'BALB NaCl (2) Rep (1) 02.10.15': None,
            'BALB NaCl (2) Rep (2) 11.10.15': None,
            'FX WT females EH Nen (1) 07.12.15': None,
            'FAMB KO males (1) rep 2 26.02.16': None,
            'FAMB WT males (1) rep 1 22.02.16': None,
            'FAMB WT males (1) rep 2 02.03.16': None,
            'FX WT females EH (B) 05.05.16': None,
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16': None,
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16': None,
            'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16': None,
            'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16': None,
            'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16': None,
            'PV Cre males HT EH (1) 09.06.16': None,
            }

def plot_files(ehd, ax=None):
    """Plot which files were loaded."""
    if ax is None:
        figg = plt.figure()
        ax = plt.subplot(111)
    else:
        figg = ax.get_figure()
    files = ehd._fnames
    days = sorted(set(map(lambda x: x.split('_')[0], files)))
    for didx, day in enumerate(days): 
        ax.text(-5, didx, str(day))
        for time in range(24):
            ax.add_patch(Rectangle((time, didx), 0.5, 0.5, fc='none', ec='k', lw=1))
    for ff in files:
        day = ff.split('_')[0]
        didx = days.index(day)
        time = int(ff.split('.')[0].split('_')[1][0:2])
        ax.add_patch(Rectangle((time, didx), 0.5, 0.5, fc='k', ec='k', lw=1))
    plt.setp(ax, 'xlim', [-5, 24], 'ylim', [0, len(days)])
    # ax.draw()

def plot_antennas_diags(ehd, ax=None, tstart=None, tend=None, binsize=6*3600.):
    """Plot activity level for each antenna"""
    if ax is None:
        figg = plt.figure()
        ax = plt.subplot(111)
    else:
        figg = ax.get_figure()
    tt = ehd.gettimes(ehd.mice)
    an = ehd.getantennas(ehd.mice)
    if tstart is None:
        tstart = min(tt)
    if tend is None:
        tend = max(tt)
    res = {}
    bins = np.arange(tstart, tend+binsize, binsize)
    for aa in range(1, 9):
        taa = map(lambda x: x[0], filter(lambda x: x[1] == aa, zip(tt, an)))
        res[aa] = np.histogram(taa, bins=bins)[0]
        ax.plot((bins[:-1] - tstart)/3600., res[aa], 'k')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Antenna signals')

def plot_antennas(ehd, ax=None, smell_comps=[1, 3], comp_labels={}):
    """Plot how the antennas are numbered."""
    if ax is None:
        figg = plt.figure()
        ax = plt.subplot(111)
    else:
        figg = ax.get_figure()
    ax.axis([-2, 5, -2, 5])
    ant_labels_pos_x = {1: -1, 2: -1, 3: 1, 4: 2, 5: 4, 6: 4, 7: 2, 8: 1}
    ant_labels_pos_y = {1: 1, 2: 2, 3: 4, 4: 4, 5: 2, 6: 1, 7: -1, 8: -1}
    comp_labels_pos_x = {1: -0.5, 2: 3.5, 3: 3.5, 4: -0.5}
    comp_labels_pos_y = {1: 3.5, 2: 3.5, 3: -0.5, 4: -0.5}
    for ant in ehd._ant_pos:
        ax.text(ant_labels_pos_x[ehd._ant_pos[ant]], ant_labels_pos_y[ehd._ant_pos[ant]], ant)
    ax.plot([-0.5, -0.5], [0, 3], '0.5', lw=4)
    ax.plot([3.5, 3.5], [0, 3], '0.5', lw=4)
    ax.plot([0, 3], [-0.5, -0.5], '0.5', lw=4)
    ax.plot([0, 3], [3.5, 3.5], '0.5', lw=4)
    ax.add_patch(Rectangle((-1, 3), 1, 1, fc='none', ec='k', lw=2))
    ax.add_patch(Rectangle((3, 3), 1, 1, fc='none', ec='k', lw=2))
    ax.add_patch(Rectangle((-1, -1), 1, 1, fc='none', ec='k', lw=2))
    ax.add_patch(Rectangle((3, -1), 1, 1, fc='none', ec='k', lw=2))
    if 1 in smell_comps:
        ax.plot([-1, -0.5], [3.5, 4], 'k:', lw=2)
    if 2 in smell_comps:
        ax.plot([3.5, 4], [4, 3.5], 'k:', lw=2)
    if 3 in smell_comps:
        ax.plot([3.5, 4], [-1, -0.5], 'k:', lw=2)
    if 4 in smell_comps:
        ax.plot([-1, -0.5], [-0.5, -1], 'k:', lw=2)
    for comp in comp_labels:
        ax.text(comp_labels_pos_x[comp], comp_labels_pos_y[comp],
                 comp_labels[comp], ha='center', va='center', weight='bold')
    plt.setp(ax, 'xticks', [], 'yticks', [], 'frame_on', False)

def data2csv(name,headers, data, data_type = None):
    f = open(name, 'w')
    header = '\t'.join(['Mouse',] + map(time2title, tts))
    for i in range(len(data)):
        f.write(headers[i] %smells[path]['soc'])
        f.write(header + '\n')
        for mm in ehd.mice: 
            f.write('%s' %(mm))
            for idx in range(len(tts)):
                try:                
                    f.write('\t%s' %locale.format("%f", data_type[mm][idx]))
                except KeyError:
                    f.write('\t')
            f.write('\n')
        f.write('\n\n')
    f.close()
            
if __name__ == '__main__':
    binsizes = [12 * 3600., 2 * 3600., 1 * 3600.]
    bintitles = ['12', '2', '1']


    for path in datasets[datarange]:
        path1 = '../RawData/' + path
        ehd = EcoHab.EcoHabData(path1, _ant_pos=antenna_positions[path])
        ehs = EcoHab.EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path1)
        tstart, tend = cf.gettime('ALL')
        
        plt.figure(figsize=(8, 8))
        plot_antennas(ehd, ax=plt.subplot(311), smell_comps=smells[path].values(), 
                      comp_labels={v: k for k, v in smells[path].items()})
        plt.title(path)
        plot_files(ehd, ax=plt.subplot(312))
        plot_antennas_diags(ehd, ax=plt.subplot(313), tstart=tstart, tend=tend, 
                        binsize=6*3600.)
    
        plt.savefig('../Results/'+'setup_%s.pdf' %path)
    
        def time2title(tt):
            return cf(epoch2num(tt + 300.))

        for binsize, btit in zip(binsizes, bintitles):
            tt = tstart
            tts = []
            results_soc = dict([(mm, []) for mm in ehd.mice])
            results_nsoc = dict([(mm, []) for mm in ehd.mice])
            results_soc_t = dict([(mm, []) for mm in ehd.mice])
            results_nsoc_t = dict([(mm, []) for mm in ehd.mice])
            test_exploration = dict([(mm, []) for mm in ehd.mice])
            while tt < tend:
                ehs.unmask_data()
                ehs.mask_data(tt, tt + binsize)
                for mm in ehd.mice:
                    adds = ehs.getaddresses(mm)
                    durs = ehs.getdurations(mm)
                    test_exploration[mm] += durs
                    results_nsoc[mm].append(adds.count(smells[path]['nsoc']))
                    results_soc[mm].append(adds.count(smells[path]['soc']))
                    results_nsoc_t[mm].append(sum([d for d, ad in zip(durs, adds)
                                  if ad == smells[path]['nsoc']]))
                    results_soc_t[mm].append(sum([d for d, ad in zip(durs, adds)
                                  if ad == smells[path]['soc']]))
                tts.append(tt)
                tt += binsize
                # Save data to cvs file
                name = 'res_%s_%s.csv' %(path, btit)
                data = [results_soc_t,results_soc_t,
                             results_soc,results_nsoc]
                headers = ['Total time with social smell (box %d), seconds\n',
                           'Total time with non-social smell (box %d), seconds\n',
                           'Number of visits to social smell (box %d)\n',
                           'Number of visits to non-social smell (box %d)\n']
                #data2csv(name, headers, data)
            colors = ['black', 'blue','green','cyan','magenta','yellow','red']
            t2 = np.arange(0,4000,0.1)
            tuples = []
            for c in range(len(colors)):
                mm = list(ehd.mice)[c]
                d = test_exploration[mm]
                if d == []:
                    continue
                t = np.copy(d)
                for i in range(1,len(t)):
                    t[i] = t[i]+t[i-1]
                j = 0
                i = 0
                v2 = np.zeros(len(t2))
                while j<len(t2)-1:
                    if t2[j]>t[-1]:
                        break
                    v2[j] = 1.0/np.array(d[i])
                    j+=1
                    if t2[j]>t[i]:
                        i+=1
                plt.plot(t2,v2,'-', color = colors[c])
                if c <2:
                    tuples.append(v2)
            plt.show()
        
            plt.subplot(2,1,1)
            plt.plot(tuples[0],'-', color = colors[2])
            plt.plot(tuples[1],'-', color = colors[-1])
            plt.subplot(2,1,2)
            plt.plot(tuples[0]*tuples[1])
            plt.show()
        ## Totals in dark and light phases
        results_sec = dict([(mm, {}) for mm in ehd.mice])
        for sec in cf.sections()[:-1]:
            ehs.unmask_data()
            ehs.mask_data(*cf.gettime(sec))
            for mm in ehd.mice:
                results_sec[mm][sec] = ehs.getstats(mm)
        #print results_sec
        
        """name = 'activity_%s.csv' %(path)
        headers = ['Total number of visits to boxes\n']
        data = [results_sec[:][:][0]]
        for box in [0, 1, 2, 3]:
            headers.append('Number of visits to box %d\n' %(box + 1))
            #data.append(results_sec[mm][sec][0][box]))"""
        