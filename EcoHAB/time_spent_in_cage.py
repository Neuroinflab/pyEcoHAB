import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
# import IntelliCage_tools as ict
import locale
from ExperimentConfigFile import ExperimentConfigFile
import os

### How much time mice spend in the 'social' compartment

datarange = slice(10, 11, None)
datasets = {
    'long':'/home/jszmek/EcoHAB_data_November/long_experiment_WT',
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

def save_data_cvs(data,fname,cf,key):
    
    f = open(fname,'w')
    phases = cf.sections()
    phases = cf.sections()
    
    header = 'mouse,\"time [h]\"'
    for phase in phases:
        header += ',\"'+phase+"\""
  
    for stim in ['soc','nsoc']:
        
          for j,h in enumerate(headers[stim]):

            f.write(h%smells[key][stim])
            f.write(header+'\n')
            for mouse in ehs.mice:
                lines = [mouse for i in range(len(data[stim][j][phases[0]][mouse]))]
                for phase in phases:
                    
                    if phase != 'ALL':
                        for k,time in enumerate(data['time'][phase]):
                            if phase == phases[0]:
                                lines[k] += ',' + str(time)
                            if isinstance(data[stim][j][phase][mouse][k],list):
                                lines[k] += ','+str(sum(data[stim][j][phase][mouse][k]))
                            else:
                                lines[k] += ','+str(data[stim][j][phase][mouse][k])
                                                                
                for line in lines:
                    f.write(line)
                    f.write('\n')
                        
                    
    
        
        

def get_data(ehs,cf,key,binsize=12*3600):
  
    tstart, tend = cf.gettime('ALL')
    phases = cf.sections()
    nsoc = dict([(phase, []) for phase in phases])
    soc = dict([(phase, []) for phase in phases])
    t_nsoc = dict([(phase, []) for phase in phases])
    t_soc = dict([(phase, []) for phase in phases])
    test_exploration = dict([(phase, []) for phase in phases])
    tt =  dict([(phase, []) for phase in phases])
    for phase in phases:
        nsoc[phase] = dict([(mm, []) for mm in ehs.mice])
        soc[phase] = dict([(mm, []) for mm in ehs.mice])
        t_nsoc[phase] = dict([(mm, []) for mm in ehs.mice])
        t_soc[phase] = dict([(mm, []) for mm in ehs.mice])
        test_exploration[phase] = dict([(mm, []) for mm in ehs.mice])
        time_start,phase_end = cf.gettime(phase)
        time = time_start
      
        while time < phase_end:
            ehs.unmask_data()
            ehs.mask_data(time, time + binsize)
            for mouse in ehs.mice:
                mouse_address = ehs.getaddresses(mouse)
                mouse_durations = ehs.getdurations(mouse)
                test_exploration[phase][mouse] += mouse_durations
                nsoc[phase][mouse].append(mouse_address.count(smells[key]['nsoc']))
                soc[phase][mouse].append(mouse_address.count(smells[key]['soc']))
                t_nsoc[phase][mouse].append(visits(mouse_address,mouse_durations, key, 'nsoc'))
                t_soc[phase][mouse].append(visits(mouse_address,mouse_durations, key, 'soc'))
                
                #print  nsoc[phase][mouse],soc[phase][mouse], len(t_nsoc[phase][mouse]),len(t_soc[phase][mouse])
            tt[phase].append((time-time_start)*binsize/3600.)
            time += binsize
            
    
    
    return {'nsoc':[nsoc,t_nsoc], 'soc':[soc,t_soc],'time':tt}

        
                
        
    
  
if __name__ == '__main__':
    binsizes = [12 * 3600., 2 * 3600., 1 * 3600.]
    bintitles = ['12', '2', '1']


    for key in datasets:
        path1 = datasets[key]
        ehd = EcoHab.EcoHabData(path1)#, _ant_pos=antenna_positions[key])
        ehs = EcoHab.EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path1)
        tstart, tend = cf.gettime('ALL')
        
              
        plt.figure(figsize=(8, 8))
        plot_antennas(ehd, ax=plt.subplot(311), smell_comps=smells[key].values(), 
                      comp_labels={v: k for k, v in smells[key].items()})
        plt.title(key)
        plot_files(ehd, ax=plt.subplot(312))
        plot_antennas_diags(ehd, ax=plt.subplot(313), tstart=tstart, tend=tend, 
                        binsize=6*3600.)
    
        plt.savefig(os.path.join(path1,'Results/','setup_%s.pdf' %key))

        data = get_data(ehs,cf,key,binsize=12*3600)
        fname = os.path.join(path1,'Results','collective_results_social_non_social.csv')
        save_data_cvs(data,fname,cf,key)
        
        # def time2title(tt):
        #     return cf(epoch2num(tt + 300.))

        # for binsize, btit in zip(binsizes, bintitles):
        #     print(binsize,btit)
        #     tt = tstart
        #     tts = []
        #     results_soc = dict([(mm, []) for mm in ehd.mice])
        #     results_nsoc = dict([(mm, []) for mm in ehd.mice])
        #     results_soc_t = dict([(mm, []) for mm in ehd.mice])
        #     results_nsoc_t = dict([(mm, []) for mm in ehd.mice])
        #     test_exploration = dict([(mm, []) for mm in ehd.mice])
        #     while tt < tend:
        #         ehs.unmask_data()
        #         ehs.mask_data(tt, tt + binsize)
        #         for mm in ehd.mice:
        #             adds = ehs.getaddresses(mm)
        #             durs = ehs.getdurations(mm)
                    
        #             test_exploration[mm] += durs
                    
        #             results_nsoc[mm].append(adds.count(smells[key]['nsoc']))
                   
        #             results_soc[mm].append(adds.count(smells[key]['soc']))
        #             results_nsoc_t[mm].append(sum([d for d, ad in zip(durs, adds)
        #                           if ad == smells[key]['nsoc']]))
        #             results_soc_t[mm].append(sum([d for d, ad in zip(durs, adds)
        #                           if ad == smells[key]['soc']]))
        #         tts.append(tt)
        #         tt += binsize
        #         # Save data to cvs file
        #     name = 'res_%s_%s.csv' %(key, btit)
        #     data = [results_soc_t,results_soc_t,
        #             results_soc,results_nsoc]
            
        #     headers = ['Total time with social smell (box %d), seconds\n',
        #                'Total time with non-social smell (box %d), seconds\n',
        #                'Number of visits to social smell (box %d)\n',
        #                'Number of visits to non-social smell (box %d)\n']
        #     data2csv(name, headers, data,key)
        #     colors = ['black', 'blue','green','cyan','magenta','yellow','red']
        #     t2 = np.arange(0,4000,0.1)
        #     tuples = []
        #     for c in range(len(colors)):
        #         mm = list(ehd.mice)[c]
        #         d = test_exploration[mm]
        #         if d == []:
        #             continue
        #         t = np.copy(d)
        #         for i in range(1,len(t)):
        #             t[i] = t[i]+t[i-1]
        #         j = 0
        #         i = 0
        #         v2 = np.zeros(len(t2))
        #         while j<len(t2)-1:
        #             if t2[j]>t[-1]:
        #                 break
        #             v2[j] = 1.0/np.array(d[i])
        #             j+=1
        #             if t2[j]>t[i]:
        #                 i+=1
        #         plt.plot(t2,v2,'-', color = colors[c])
        #         if c <2:
        #             tuples.append(v2)
        #     #plt.show()
        
        #     plt.subplot(2,1,1)
        #     plt.plot(tuples[0],'-', color = colors[2])
        #     plt.plot(tuples[1],'-', color = colors[-1])
        #     plt.subplot(2,1,2)
        #     plt.plot(tuples[0]*tuples[1])
        #     #plt.show()
        # ## Totals in dark and light phases
        # results_sec = dict([(mm, {}) for mm in ehd.mice])
        # for sec in cf.sections()[:-1]:
        #     ehs.unmask_data()
        #     ehs.mask_data(*cf.gettime(sec))
        #     for mm in ehd.mice:
        #         results_sec[mm][sec] = ehs.getstats(mm)
        # #print results_sec
   
        # """name = 'activity_%s.csv' %(key)
        # headers = ['Total number of visits to boxes\n']
        # data = [results_sec[:][:][0]]
        # for box in [0, 1, 2, 3]:
        #     headers.append('Number of visits to box %d\n' %(box + 1))
        #     #data.append(results_sec[mm][sec][0][box]))"""
    plt.show()     
