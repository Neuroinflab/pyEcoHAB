from __future__ import division, print_function
import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
# import IntelliCage_tools as ict
import locale
from ExperimentConfigFile import ExperimentConfigFile
import os
import utils
### How much time mice spend in the 'social' compartment

datarange = slice(10, 11, None)
datasets = {
    #'long':'/home/jszmek/EcoHAB_data_November/long_experiment_WT',
    'standard':'/home/jszmek/EcoHAB_data_November/standard_known_stimulus_WT',
}
smells = {
    'long': {'soc': 3, 'nsoc': 1},
    'standard': {'soc': 3, 'nsoc': 1}
}
antenna_positions = {
    'long': None,
    'standard': None,
}
headers = {'soc':['Number of visits to social smell (box %d)\n','Total time with social smell (box %d), seconds\n'],
           'nsoc':['Number of visits to non-social smell (box %d)\n','Total time with non-social smell (box %d), seconds\n',]}

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

def data2csv(name,headers, data,path, data_type=None):
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
                except (TypeError, KeyError):
                    f.write('\t')
            f.write('\n')
        f.write('\n\n')
    f.close()
    
def visits(mouse_address,mouse_durations, key, stim_type):
    result = []
    for address,t in  zip(mouse_address, mouse_durations):
        if address == smells[key][stim_type]:
            result.append(t)
    return result

def save_data_cvs(data,fname,path,key):

    if not os.path.exists(path):
        os.makedirs(path)
    fname = os.path.join(path,fname)
    f = open(fname,'w')
    phases = data['phases']
    mice = data['mice']
    
    header = 'mouse,\"time [h]\"'
    for phase in phases:
        header += ',\"'+phase+"\""
  
    for stim in ['soc','nsoc']:
        
          for j,h in enumerate(headers[stim]):

            f.write(h%smells[key][stim])
            f.write(header+'\n')
            for mouse in mice:
                lines = [mouse for i in data[stim][j][phases[0]][mouse]]
                for phase in phases:
                    if phase != 'ALL':
                        for k,time in enumerate(data['time'][phase]):
                            if phase == phases[0]:
                                lines[k] += ',%3.2f'%(time/3600)
                            lines[k] += ','+str(data[stim][j][phase][mouse][k])
                                                                
                for line in lines:
                    f.write(line)
                    f.write('\n')

def get_time_spent(ehs,cf,key,binsize=12*3600):
  
    tstart, tend = cf.gettime('ALL')
    phases = cf.sections()
    nsoc = dict([(phase, []) for phase in phases])
    soc = dict([(phase, []) for phase in phases])
    t_nsoc = dict([(phase, []) for phase in phases])
    t_soc = dict([(phase, []) for phase in phases])
    tt =  dict([(phase, []) for phase in phases])
    
    for phase in phases:
        nsoc[phase] = dict([(mm, []) for mm in ehs.mice])
        soc[phase] = dict([(mm, []) for mm in ehs.mice])
        t_nsoc[phase] = dict([(mm, []) for mm in ehs.mice])
        t_soc[phase] = dict([(mm, []) for mm in ehs.mice])
        time_start,phase_end = cf.gettime(phase)
        time = time_start
      
        while time < phase_end:
            ehs.unmask_data()
            ehs.mask_data(time, time + binsize)
            for mouse in ehs.mice:
                mouse_address = ehs.getaddresses(mouse)
                mouse_durations = ehs.getdurations(mouse)
 
                nsoc[phase][mouse].append(mouse_address.count(smells[key]['nsoc']))
                soc[phase][mouse].append(mouse_address.count(smells[key]['soc']))
                t_nsoc[phase][mouse].append(visits(mouse_address,mouse_durations, key, 'nsoc'))
                t_soc[phase][mouse].append(visits(mouse_address,mouse_durations, key, 'soc'))
                
                #print  nsoc[phase][mouse],soc[phase][mouse], len(t_nsoc[phase][mouse]),len(t_soc[phase][mouse])
            tt[phase].append((time-time_start)*binsize/3600.)
            time += binsize
 
    for phase in phases:
        for mouse in ehs.mice:
            nsoc[phase][mouse] = np.array(nsoc[phase][mouse])
            soc[phase][mouse] = np.array(soc[phase][mouse])
            t_nsoc[phase][mouse] = np.array(t_nsoc[phase][mouse])
            t_soc[phase][mouse] = np.array(t_soc[phase][mouse])

    return {'nsoc':[nsoc,t_nsoc], 'soc':[soc,t_soc],'time':tt,'mice':ehs.mice,'phases':phases}

def sum_times(data):
    
    phases = data['phases']
    mice = data['mice']
    sum_t_soc = dict([(phase, []) for phase in phases])
    sum_t_nsoc = dict([(phase, []) for phase in phases])

    for phase in phases:

        sum_t_soc[phase] = dict([(mouse, []) for mouse in mice])
        sum_t_nsoc[phase] = dict([(mouse, []) for mouse in mice])
          
        for mouse in mice:

            sum_t_soc[phase][mouse] = np.zeros((len(data['time'][phase])))
            sum_t_nsoc[phase][mouse] = np.zeros((len(data['time'][phase])))
            
            for k,time in enumerate(data['time'][phase]):
                sum_t_soc[phase][mouse][k] = sum(data['soc'][1][phase][mouse][k])
                sum_t_nsoc[phase][mouse][k] = sum(data['nsoc'][1][phase][mouse][k])
 
    return {'phases':phases,'mice':mice,'soc':[data['soc'][0],sum_t_soc],'nsoc':[data['nsoc'][0],sum_t_nsoc],'time':data['time']}
                    
                
def calculate_approach_to_social(data):
    #possibly calculate a square root of sum of squared errors -- if independent errors (Taylor, 1997) probably not
    ats = data[0]/data[1]/(data[2]/data[3])
    d_ats = np.zeros(data[0].shape)
    for i,dat in enumerate(data[4:]):
        d_ats += abs(dat/data[i])
    return ats,d_ats*ats
        
def approach_to_social(data,fname=None,path=None):

    times_nsoc = data['nsoc'][1]
    times_soc = data['soc'][1]
    phases = data['phases']
    last_phase_before_stimulus = [phase for phase in phases if 'EMPTY' in phase and 'dark' in phase][-1]
    first_phase_after_stimulus = [phase for phase in phases if 'SNIFF' in phase and 'dark' in phase][0]

    mice = data['mice']
    mice_no = len(mice)
    shape = (len(times_soc[first_phase_after_stimulus][mice[0]]),len(mice))
    T_S, T_NS, t_s, t_ns = np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(shape)
    
    for i,mouse in enumerate(mice):
        T_S [:,i] = times_soc[first_phase_after_stimulus][mouse]
        T_NS[:,i] = times_nsoc[first_phase_after_stimulus][mouse]
        t_s[:,i] = times_soc[last_phase_before_stimulus][mouse]
        t_ns[:,i] = times_nsoc[last_phase_before_stimulus][mouse]
        
    mean_T_S, mean_T_NS, mean_t_s, mean_t_ns = T_S.mean(axis=1), T_NS.mean(axis=1), t_s.mean(axis=1), t_ns.mean(axis=1)
    std_T_S, std_T_NS, std_t_s, std_t_ns = (T_S.var(axis=1)/mice_no)**.5, (T_NS.var(axis=1)/mice_no)**.5, (t_s.var(axis=1)/mice_no)**.5, (t_ns.var(axis=1)/mice_no)**.5
    ats,d_ats = calculate_approach_to_social([mean_T_S, mean_T_NS, mean_t_s, mean_t_ns,std_T_S, std_T_NS, std_t_s, std_t_ns])
    l = len(data['time'][last_phase_before_stimulus])
    try:
        binsize = data['time'][last_phase_before_stimulus][1]-data['time'][last_phase_before_stimulus][0]
    except IndexError:
        binsize = 12*3600
        
    time = [binsize*i/3600 for i in range(l)]

    if fname:
        if path:
            if not os.path.exists(path):
                os.makedirs(path)
            fname = os.path.join(path,fname)
        np.savetxt(fname, np.array((time,ats,d_ats)).T,header='time, AtS, dAtS',delimiter=',',comments="")
        
    return time, ats, d_ats


if __name__ == '__main__':
    binsizes = [12 * 3600., 2 * 3600., 1 * 3600.,1.5*3600,3600/4]
    bintitles = ['12', '2', '1']

    for binsize in binsizes:
        for key in datasets:
            path1 = datasets[key]
            ehd = EcoHab.EcoHabData(path1)
            ehs = EcoHab.EcoHabSessions(ehd)
            cf = ExperimentConfigFile(path1)
            tstart, tend = cf.gettime('ALL')
            
            print('Binsize ',binsize/3600)
            data = get_time_spent(ehs,cf,key,binsize=binsize)
            data = sum_times(data)
            path = utils.results_path(path1)
            fname = 'collective_results_social_non_social_binsize_%f_h.csv'%(binsize/3600)
            save_data_cvs(data,fname,path,key)
            fname =  'AtS_%f_h.csv'%(binsize/3600)
            approach_to_social(data,fname,path)
        #print data
        
        

