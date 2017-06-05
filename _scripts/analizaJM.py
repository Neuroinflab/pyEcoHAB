import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
# import IntelliCage_tools as ict
import locale
from ExperimentConfigFile import ExperimentConfigFile
import pickle
import os
from experiments_info import smells, antenna_positions
### How much time mice spend in the 'social' compartment

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

def conv2signal(d,pos, t_max=4000):
    t2 = np.arange(0,t_max,0.1)
    t = np.copy(d)
    for i in range(1,len(t)):
        t[i] = t[i]+t[i-1]
    j = 0
    i = 0
    v2 = np.zeros((len(t2),2))
    while j<len(t2)-1:
        if t2[j]>t[-1]:
            break
        v2[j][0] = (pos[i]-1)/2
        v2[j][1] = (pos[i]-1)%2
        j+=1
        if t2[j]>t[i]:
            i+=1
    return (t2,v2)
    
def norm(x):
    return x/np.linalg.norm(x)
if __name__ == '__main__':
    binsizes = [12 * 3600., 2 * 3600., 1 * 3600.]
    bintitles = ['12', '2', '1']


    for el in next(os.walk('../PreprocessedData/'))[2]:
        with open('../PreprocessedData/'+el, "rb") as input_file:
            ehd = pickle.load(input_file)
            ehs = pickle.load(input_file)
            path = el[:-4]
            path1 = '../RawData/' + path
            cf = ExperimentConfigFile(path1)

            tstart, tend = cf.gettime('ALL')

        for binsize, btit in zip(binsizes, bintitles):
            tt = tstart
            tts = []
            results_soc = dict([(mm, []) for mm in ehd.mice])
            results_nsoc = dict([(mm, []) for mm in ehd.mice])
            results_soc_t = dict([(mm, []) for mm in ehd.mice])
            results_nsoc_t = dict([(mm, []) for mm in ehd.mice])
            test_exploration = dict([(mm, []) for mm in ehd.mice])
            test_position = dict([(mm, []) for mm in ehd.mice])
            while tt < tend:
                ehs.unmask_data()
                ehs.mask_data(tt, tt + binsize)
                for mm in ehd.mice:
                    adds = ehs.getaddresses(mm)
                    durs = ehs.getdurations(mm)
                    test_exploration[mm] += durs
                    test_position[mm] += adds
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
            lt = []
            lv = []
            for c in range(len(colors)):
                mm = list(ehd.mice)[c]
                d = test_exploration[mm]
                pos = test_position[mm]
                if d == []:
                    continue
                t,v = conv2signal(d,pos, t_max=np.max(d))
                lt.append(t)
                lv.append(v)
                plt.plot(t,v[:,0],'-', color = 'r')
                plt.plot(t,v[:,1]-2,'-', color = 'b')
                plt.show()
                print np.dot(norm(t),norm(v))
            plt.show()
            break
        break
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
        