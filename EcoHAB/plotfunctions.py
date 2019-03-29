# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:38:58 2017

@author: Jan Maka
"""
from __future__ import division, print_function
from matplotlib.dates import epoch2num
import matplotlib.patches as patches
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as st
import utility_functions as utils
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import matplotlib.colors as mcol
import matplotlib.patches as patches
from write_to_file import make_table_of_pairs
nbins = 10
def autolabel(rects,ax):
    """

    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')

def single_barplot(stats,directory, groups,colors,name = "",ylab = "",titles=""):

    path = utils.check_directory(directory,"Results")
        
    N = len(groups)
    ind = np.arange(N)
    width = 0.35
    fig, ax = plt.subplots()
    rects = []
   
    means = [np.mean(stats[group]) for group in groups]
    errs = [st.sem(stats[group]) for group in groups]
    for i, mean in enumerate(means):
        
        if isinstance(titles,list):
            rects.append(ax.bar(ind[i]+width, mean, width, color=colors[i], yerr=errs[i],label=titles[i]))
        else:
            rects.append(ax.bar(ind[i]+width, mean, width, color=colors[i], yerr=errs[i]))

    ax.set_ylabel(ylab)
    ax.set_title('Group comparison')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(groups)
    ax.set_xlabel('Parameter')
    
    if isinstance(titles,list):
        ax.legend()
    
    #autolabel(rects1)
    #autolabel(rects2)

    new_fname = os.path.join(path,(name+'.png'))
    plt.savefig(new_fname)
    #plt.show()

def barplot(stats, names, groups,colors, directory = "Barplots",name = "",ylab = ""):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    N = len(groups)
    ind = np.arange(N)
    width = 0.35
    fig, ax = plt.subplots()
    rects = []
    for i,key in enumerate(names.keys()):
        mean = [np.mean(stats[key][group]) for group in groups]
        err = [st.sem(stats[key][group]) for group in groups]
       
        rects.append(ax.bar(ind+i*width, mean, width, color=colors[key], yerr=err))

    ax.set_ylabel(ylab)
    ax.set_title('Group comparison')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(groups)
    ax.set_xlabel('Parameter')
    
    ax.legend([rects[i][0] for i in range(len(rects))], [key for key in names.keys()])
    
    #autolabel(rects1)
    #autolabel(rects2)
    plt.savefig('../Results/%s/%s.png'%(directory, name))
 

def scaling(p):
    if p<0.05 and p>0:
        return 2.5-p*50
    elif p>-0.05 and p<0:
        return -2.5-p*50
    else:
        return 0

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


    
def make_RasterPlot(main_directory,
                    subdirectory,
                    FAM,
                    phases,
                    name,
                    old_mice,
                    prefix='',
                    to_file=True,
                    vmin=None,
                    vmax=None,
                    title=None):
    
    mice = [mouse[5:] for mouse in old_mice]
    subdirectory = os.path.join(subdirectory, 'raster_plots')
    new_path = utils.check_directory(main_directory, subdirectory)
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, aspect='equal')
    if title:
        plt.suptitle(title, fontsize=14, fontweight='bold')
    assert FAM.shape[0] == len(phases)
    assert FAM.shape[1] == len(mice)
    assert FAM.shape[2] == len(mice)
    
    output, pair_labels = make_table_of_pairs(FAM, phases, mice)
                
    if not vmax and not vmin:
        vmax = FAM.max()
        vmin = FAM.min()
        if vmax*vmin <= 0:
            if abs(vmax) > abs(vmin):
                vmin = -vmax
            else:
                vmax = -vmin
    if vmin*vmax < 0:
        colormap = plt.cm.bwr
    else:
        colormap = plt.cm.Reds
    cax = ax.imshow(output,
                    interpolation='none',
                    origin='lower',
                    aspect='auto',
                    vmin=vmin,
                    vmax=vmax,
                    cmap=colormap)
    cbar = fig.colorbar(cax, ax=ax, ticks=[vmin, 0, vmax])    
    fig.subplots_adjust(left=0.25)
    ax.set_xticks([i for i in range(len(phases))])
    ax.set_yticks([i for i in range(len(pair_labels))])
    assert(output.shape[0]==len(pair_labels))
    ax.set_xticklabels(phases)
    ax.set_yticklabels(pair_labels)
    for tick in ax.get_xticklabels():
            tick.set_rotation(90)
            
    fig.subplots_adjust(left=0.25)
    fig.subplots_adjust(bottom=0.25)
    name = name + prefix
    new_name = os.path.join(new_path, name+'.png')
    plt.savefig(new_name,
                transparent=False,
                bbox_inches=None,
                pad_inches=2,
                frameon=None)

def raster_interactions(directory,
                        FAM,
                        IPP,
                        phases,
                        suffix,
                        mice,
                        scalefactor=0):

    new_path = utils.check_directory(directory, 'interactions/figs/raster_plots')
    fname = os.path.join(new_path, 'Interactions_%s.png' % suffix)
    if not scalefactor:
        scalefactor = IPP.max()

    dt = np.dtype((np.float, (4,)))
    new_shape = (len(mice)*(len(mice)-1)//2, FAM.shape[0])

    followings = np.zeros(new_shape, dtype=dt)
    avoidings = np.zeros(new_shape, dtype=dt)
    pair_labels = utils.list_of_pairs(mice)
    phase_range = FAM.shape[0]
    for i in range(phase_range):
        l = 0
        phase = phases[i]
        for j, mouse in enumerate(mice):
            for k in range(j+1, len(mice)):
                if FAM[i, j, k] < 0.05:
                    if FAM[i, j, k] > 0:
                        followings[l, i] = (1., 0, 0, IPP[i, j, k]/scalefactor/2 )
                    elif FAM[i, j, k] < 0:
                        avoidings [l, i] = (0, 0, 1., IPP[i, j, k]/scalefactor/2)
                if FAM[i, k, j] < 0.05:
                    if FAM[i, k, j] > 0:
                        if np.any(followings[l, i] != 0):
                            followings[l, i][3] += IPP[i, k, j]/scalefactor/2
                        else:
                            followings[l, i] = (1., 0, 0, IPP[i, k, j]/scalefactor/2)
                    elif FAM[i, k, j] < 0:
                        if np.any(avoidings != 0):
                            avoidings[l, i][3] += IPP[i, k, j]*0.5/scalefactor/2
                        else:
                            avoidings[l, i] = (0, 0, 1., IPP[i, k, j]/scalefactor/2)
                l += 1
                
   
    fig = plt.figure(figsize=(phase_range*2, 12))
    ax = fig.add_subplot(111, aspect='equal')
    fig.subplots_adjust(left=0.25)
    plt.imshow(followings, 
               interpolation='nearest',
               origin='upper',
               aspect='auto')
    plt.imshow(avoidings, 
               interpolation='nearest',
               origin='upper',
               aspect='auto')
    ax.set_xticks([i for i in range(len(phases))])
    ax.set_yticks([i for i in range(len(pair_labels))])
    xlabel = phases[FAM.shape[0]-1::-1]
    ax.set_xticklabels(xlabel[::-1])
    ax.set_yticklabels(pair_labels)
    for tick in ax.get_xticklabels():
            tick.set_rotation(90)
    fig.subplots_adjust(left=0.25)
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(fname,
                transparent=False,
                bbox_inches=None,
                pad_inches=.2,
                frameon=None)

def oneRasterPlot(directory,
                  FAM,
                  IPP,
                  phases,
                  name,
                  scalefactor,
                  mice=[],
                  to_file=True):
    
    if not name:
        name = "oneRasterPlot"
    subdirectory = 'RasterPlots'
    n_s,n_l,n_f = FAM.shape
    new_path = utils.check_directory(directory, subdirectory)
    fig = plt.figure(figsize=(n_s*2, 12))
    ax = fig.add_subplot(111, aspect='equal')

    if name:
        plt.suptitle(name, fontsize=14, fontweight='bold')

    xlabels = []
    
    for s in range(n_s):
        xlabels.append(phases[s])
        #phases
        plt.text(0.06+s*0.125,
                 1.025,
                 phases[s],
                 horizontalalignment='center',
                 verticalalignment='center',
                 fontsize=10,
                 transform = ax.transAxes)

        _FAM = FAM[s, :, :]
        _IPP = IPP[s, :, :]
        pair_labels = []
        pos = 0
        for i in range(n_l): #mice
            for j in range(i,n_f): #mice
                if i != j and abs(_FAM[i, j])<0.05 and _FAM[i, j] > 0:
                    ax.add_patch(patches.Rectangle((
                        s,-1*pos),
                        1,
                        1,
                        facecolor=(1, 0, 0,  _IPP[i,j]*0.5/scalefactor)))  
                elif i != j and abs(_FAM[i, j]) < 0.05 and _FAM[i, j] < 0:
                        ax.add_patch(patches.Rectangle((
                            s, -1*pos),
                            1 , 1,facecolor=(0,0,1,_IPP[i,j]*0.5/scalefactor)))
                        
                if i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]>0:
                    ax.add_patch(patches.Rectangle((
                            s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[j,i]*0.5/scalefactor))) 
                elif i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[j,i]*0.5/scalefactor)))
                if i!=j:
     
                    if not mice:
                        pair_labels.append(str(i+1)+'|'+str(j+1))
                    else:
                        pair_labels.append(mice[i]+'|'+mice[j])
                    pos += 1
     
    plt.axis([0.5,n_s,-pos-1,1])
    #plt.tight_layout()
    fig.subplots_adjust(left=0.25)
    fig.subplots_adjust(bottom=0.25)
    ax.set_aspect('auto')
    ax.xaxis.grid()
    ax.get_yaxis().set_ticks([-1*i+0.5 for i in range(pos)])
    ax.set_yticklabels(pair_labels,fontsize=10)
    plt.xlabel("session")
    plt.ylabel("following strength in pair")
    plt.savefig(os.path.join(new_path,name+'.png'),transparent=False, bbox_inches=None, pad_inches=2,frameon=None)

def plot_graph(FAPmatrix, phase, directory, labels=None):
    
    new_path = utils.check_directory(directory, 'interactions/figs/graphs')
    d2,d3 = FAPmatrix.shape
    
    pairs = []       
    for i in range(d2):
        for j in range(d3):
            if i!=j:
                pairs.append([FAPmatrix[i][j],i+1,j+1])
    pairs.sort(key=lambda x: -x[0])
    conn = [(scaling(x[0]),x[1],x[2]) for x in pairs]
    G = nx.MultiDiGraph(multiedges=True, sparse=True)
    for i in range(len(conn)):
        G.add_edges_from([(conn[i][1],conn[i][2])], weight=conn[i][0])
    edge_colors = [conn[i][0] for i in range(len(conn))]
    size = 10
    pos = nx.circular_layout(G)
    headers = ''
    if labels:
        node_labels = {}
        headers = ';'
        for i,node in enumerate(G.nodes()):
            node_labels[node] = labels[i]
            headers += labels[i]+';'
    else:
        
        node_labels = {node:node for node in G.nodes()}     
    fig = plt.figure(figsize=(10*size,10*size)) 
    ax = fig.add_subplot(111, aspect='equal')
 
    nx.draw_networkx_labels(G, pos, labels=node_labels,font_size=120)
    cmap = mcol.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                             colors =[(0, 0, 1), 
                                                      (1, 1., 1), 
                                                      (1, 0, 0)],
                                             N=20-1,
                                             )
    nx.draw(G,pos,ax, node_size=500*size**2,node_color= 'grey',edge_color=edge_colors,edge_cmap=cmap,width = 0,arrows=False, edge_style='dashed',zorder=10)
    for c in conn:
        if abs(c[0]-1)>0.01:
            #print c[0]

            p = patches.FancyArrowPatch(pos[c[2]],pos[c[1]],connectionstyle='arc3, rad=-0.3',arrowstyle="simple",shrinkA=10.2*size, shrinkB=10.2*size,mutation_scale=20*size*abs(c[0]), color = cmap(c[0]+0.5),zorder=1,alpha=0.5)
            ax.add_patch(p)
    new_fname = os.path.join(new_path, 'Interactions_graph_%s.png' % phase)
    plt.savefig(new_fname,
                transparent=False,
                bbox_inches=None,
                pad_inches=2,
                frameon=None)
    plt.close(fig) 


def single_heat_map(result,
                    name,
                    directory,
                    mice,
                    prefix,
                    phase,
                    xlabel='',
                    ylabel='',
                    subdirectory=None,
                    vmax=None,
                    vmin=None,
                    xticks=None,
                    yticks=None):
    name = '%s_%s_%s.png' % (name, prefix, phase)
    fig, ax = plt.subplots()
    if not vmin:
        vmin = result.min()
    if not vmax:
        vmax = result.max()
        
    cax = ax.imshow(result,interpolation='none', aspect='auto', cmap="viridis", origin="lower", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cax)
    turn = False
    if not xticks:
        xticks = [mouse.split('-')[-1] for mouse in mice]
        turn = True
    if not yticks:
        yticks = [mouse.split('-')[-1] for mouse in mice]
        fig.subplots_adjust(left=0.25)
        
    ax.get_yaxis().set_ticks([i for i,x in enumerate(yticks)])
    ax.get_xaxis().set_ticks([i for i,x in enumerate(xticks)])
    ax.set_xticklabels(xticks)
    ax.set_yticklabels(yticks)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    
    if subdirectory:
        subdirectory = os.path.join(subdirectory, 'figs')
    dir_name = utils.check_directory(directory, subdirectory)
    new_name = os.path.join(dir_name, name)
    fig.subplots_adjust(left=0.25)
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(new_name, transparent=False, bbox_inches=None, pad_inches=2, frameon=None)
    
def single_in_cohort_soc_plot(results,
                              results_exp,
                              mice,
                              phase,
                              fname,
                              main_directory,
                              directory,
                              prefix,
                              vmin=0,
                              vmax=1,
                              vmin1=-1,
                              vmax1=1,
                              hist=True,
                              titles=['% time together',
                                      'Expected % time together',
                                      'Excess % time together',
                                      'Histogram of excess % time together'],
                              labels=['', '']):
    new_name = os.path.join(directory, 'figs')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s_%s'% (fname, prefix, phase))
    label_mice = [mouse[-4:] for mouse in mice]
    fig = plt.figure(figsize=(10, 6))
    ax = []
    for i in range(1,5):
        ax.append(fig.add_subplot(2,2,i))
    plt.subplot(221)
    im2 = ax[0].imshow(results, vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
    cbar = fig.colorbar(im2)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_title(titles[0])
    ax[1].imshow(results_exp, vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_title(titles[1])

    deltas = results[results > 0] - results_exp[results > 0]

    plt.subplot(223)
    try:
        if -abs(vmin1) != abs(vmax1):
            maxi = max(abs(vmax1), abs(vmin1))
            mini = -maxi
        else:
            mini = vmin1
            maxi = vmax1
        im = ax[2].imshow((results - results_exp),
                          vmin=mini,
                          vmax=maxi,
                          interpolation='none',
                          cmap='bwr',
                          origin='lower')
            
        ax[2].set_title(titles[2])
        cbar = fig.colorbar(im)
        ax[2].yaxis.tick_left()
        ax[2].get_yaxis().set_ticks([i for i,x in enumerate(mice)])
        ax[2].get_xaxis().set_ticks([i for i,x in enumerate(mice)])
        ax[2].set_xticklabels(label_mice)
        ax[2].set_yticklabels(label_mice)
        ax[2].set_xlabel(labels[0])
        ax[2].set_ylabel(labels[1])
        for label in ax[2].xaxis.get_ticklabels():
            label.set_rotation(90)

    except ValueError:
        pass
    try:
        ax[3].hist(deltas)
    except ValueError:
        pass
    ax[3].set_title(titles[3])
    if hist == True:
        ax[3].set_xlim([-0.1, 0.5])
        ax[3].get_xaxis().set_ticks([-0.1, 0., 0.1, 0.2, 0.3])
        ax[3].set_xticklabels([-0.1, 0., 0.1, 0.2, 0.3])
    fig.suptitle(phase)
    fig.subplots_adjust(left=0.3)
    fig.subplots_adjust(bottom=0.3)
    fig.subplots_adjust(wspace=0.25)
    fig.subplots_adjust(hspace=0.3)
    
    fig.savefig(fname+'.pdf')
    fig.savefig(fname+'.png', dpi=300)


def make_pooled_histograms(res,
                           res_exp,
                           phases,
                           fname,
                           main_directory,
                           directory,
                           prefix,
                           additional_info):

    fig, ax = plt.subplots(1, len(phases), figsize=(len(phases)//2*5, 5))
    max_bins = 0
    min_bins = 100
    max_count = 0
    min_count = 0
    for i, phase in enumerate(phases):
        results = res[i]
        results_exp = res_exp[i]
        deltas = results[results > 0] - results_exp[results > 0]
        n, bins, patches = ax[i].hist(deltas, bins=nbins)
        ax[i].set_title(phases[i], fontsize=14)
        if bins.min() < min_bins:
            min_bins = bins.min()
        if bins.max() > max_bins:
            max_bins = bins.max()
        if max(n) > max_count:
            max_count = max(n)
        if min(n) < min_count:
            min_count = min(n)
        if i:
            ax[i].set_yticklabels([])
        else:
            for tick in ax[i].yaxis.get_major_ticks():
                    tick.label.set_fontsize(14)
        for tick in ax[i].xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 

    for x in ax:
        x.set_xlim([min_bins, max_bins+1])
        x.set_ylim([min_count, max_count+3])

        
    new_name = os.path.join(directory, 'figs')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s_%s'% (fname, prefix, phase))
    print(fname)
    fig.subplots_adjust(wspace=0.15)
    fig.savefig(fname+'.png', dpi=300)
    
def histo():
    figure, ax = plt.subplots(2, 2)
    n, bins1, patches = ax[0][0].hist(out[3], bins=bins)
    ax[0][0].set_title('12')
    ax[0][0].set_xscale("log", nonposx='clip')
    n, bins2, patches = ax[0][1].hist(out[7], bins=bins)
    ax[0][1].set_title('34')
    ax[0][1].set_xscale("log", nonposx='clip')
    n, bins3, patches = ax[1][0].hist(out[11], bins=bins)
    ax[1][0].set_title('56')
    ax[1][0].set_xscale("log", nonposx='clip')
    n, bins4, patches = ax[1][1].hist(out[15], bins=bins)
    ax[1][1].set_title('78')
    ax[1][1].set_xscale("log", nonposx='clip')

    max_lim = max(max(bins1[-1], bins2[-1]), max(bins3[-1], bins4[-1]))
    
    ax[0][0].set_xlim([0, max_lim])
    ax[0][1].set_xlim([0, max_lim])
    ax[1][0].set_xlim([0, max_lim])
    ax[1][1].set_xlim([0, max_lim])

