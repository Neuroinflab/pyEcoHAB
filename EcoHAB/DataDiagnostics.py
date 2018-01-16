import sys
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
import numpy as np
import sys
from ExperimentConfigFile import ExperimentConfigFile
import EcoHab
smells = {'soc': 3, 'nsoc': 1}

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
        ax.plot((bins[:-1] - tstart)/3600., res[aa], label=str(aa))
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Antenna signals')
	ax.legend()

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

if __name__ == '__main__':
    try:
        path = sys.argv[1]
    except IndexError:
        sys.exit('No data directory')
        
    ehd = EcoHab.EcoHabData(path)
    ehs = EcoHab.EcoHabSessions(ehd)
    cf = ExperimentConfigFile(path)
    tstart, tend = cf.gettime('ALL')

    plt.figure(figsize=(8, 8))
    plot_antennas(ehd, ax=plt.subplot(311), smell_comps=smells.values(), 
              comp_labels={v: k for k, v in smells.items()})
    plt.title(path)
    plot_files(ehd, ax=plt.subplot(312))
    plot_antennas_diags(ehd, ax=plt.subplot(313), tstart=tstart, tend=tend, 
                    binsize=6*3600.)

    #plt.savefig('../Results/'+'setup_%s.pdf' %path)
    plt.show()
