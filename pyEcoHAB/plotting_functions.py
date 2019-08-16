# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:38:58 2017

@author: Jan Maka
"""
from __future__ import division, print_function, absolute_import
import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcol
import matplotlib.patches as patches
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from . import utility_functions as utils

nbins = 10

    
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
    
    mice = [mouse[-4:] for mouse in old_mice]
    subdirectory = os.path.join(subdirectory, 'raster_plots')
    new_path = utils.check_directory(main_directory, subdirectory)
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, aspect='equal')
    if title:
        plt.suptitle(title, fontsize=14, fontweight='bold')
    assert FAM.shape[0] == len(phases)
    assert FAM.shape[1] == len(mice)
    assert FAM.shape[2] == len(mice)
    
    output, pair_labels = utils.make_table_of_pairs(FAM, phases, mice)
                
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
    assert(output.shape[0] == len(pair_labels))
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
    if len(phases) == 1:
        ax = [ax]
        new_phase = phases[0].replace(' ', '__')
    else:
        new_phase = "%s_%s" %(phases[0].replace(' ', '__'),
                              phases[-1].replace(' ', '__'))
    bins, counts = [], []
    xticks = True
    for i, phase in enumerate(phases):
        results = res[i]
        results_exp = res_exp[i]
        deltas = results[results > 0] - results_exp[results > 0]
        new_title = phases[i]
        yticks = not i
        minb, maxb, minc, maxc = make_single_histogram(ax[i],
                                                       deltas,
                                                       new_title,
                                                       xticks=xticks,
                                                       yticks=yticks,
                                                       xlogscale=False)
        bins.extend([minb, maxb])
        counts.extend([minc, maxc])
    min_bins, max_bins = min(bins), max(bins)
    min_count, max_count = min(counts), max(counts
    for x in ax:
        x.set_xlim([min_bins, max_bins + 1])
        x.set_ylim([min_count, max_count + 3])
    new_name = os.path.join(directory, 'figs')
    directory = utils.check_directory(main_directory, new_name)
    fname =  os.path.join(directory, '%s_%s_%s'% (fname, prefix, new_phase))
    if len(phases) > 1:
        fig.subplots_adjust(wspace=0.15)
    fig.savefig(fname+'.png', dpi=300)


def make_histograms_for_every_mouse(results, fname, mice, main_directory,
                                    directory, prefix, additional_info):
    """
    results should be a dictionary of lists
    """
    fig, ax = plt.subplots(len(mice), len(mice),
                           figsize=(len(mice)//2*5, len(mice)//2*5))
    bins, counts = [], []
    new_name = "all_mice"
    for i, mouse1 in enumerate(mice):
        for j, mouse2 in enumerate(mice):
            if mouse1 == mouse2:
                continue
            if j:
                yticks = False
                ylabel = None
            else:
                yticks = True
                ylabel = mouse1[-4:]
            if i == len(mice) - 1:
                xticks = True
                xlabel = mouse2[-4:]

            else:
                xticks = False
                xlabel = None

            new_title = "%s following %s" % (mouse2[-4:], mouse1[-4:])
            key = "%s_%s" % (mouse1, mouse2)
            intervals = results[key]
            minb, maxb, minc, maxc = make_single_histogram(ax[i, j],
                                                           intervals,
                                                           new_title,
                                                           xticks=xticks,
                                                           yticks=yticks,
                                                           xlabel=xlabel,
                                                           ylabel=ylabel,
                                                           xlogscale=True)
            bins.extend([minb, maxb])
            counts.extend([minc, maxc])
    min_bin, max_bin = min(bins), max(bins)
    min_count, max_count = min(counts), max(counts)
    for s in ax:
        for x in s:
            x.set_xlim([min_bin, max_bin + 1])
            x.set_ylim([min_count, max_count + 3])

    new_dir = utils.check_directory(directory, 'figs')
    dir_name =  utils.check_directory(main_directory, new_dir)

    if prefix != "":
         new_name= '%s_%s_%s.png'% (fname, prefix, new_name)
    else:
        new_name =  '%s_%s.png'% (fname, new_name)

    fname = os.path.join(dir_name, new_name)
    print(fname)
    fig.subplots_adjust(wspace=0.15)
    fig.savefig(fname, dpi=300)


def pool_results_following(res_dict, mice):
    pooled_results = {mouse:[] for mouse in mice}
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            key = "%s_%s" %(mouse1, mouse2)
            pooled_results[mouse2] += res_dict[key]
    return pooled_results


def pool_results_followed(res_dict, mice):
    pooled_results = {mouse:[] for mouse in mice}
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            key = "%s_%s" %(mouse1, mouse2)
            pooled_results[mouse1] += res_dict[key]
    return pooled_results

def make_single_histogram(ax, single_results, title, xticks=False,
                          yticks=False, xlabel=None, ylabel=None
                          xlogscale=False, ylogscale=False):
    hist, bins = np.histogram(single_results, bins=30)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]),len(bins))
    n, bins, patches = ax.hist(single_results, bins=logbins)
    if xlogscale:
        ax.set_xscale('log')
    if ylogscale:
        ax.set_yscale('log')
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=14)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14)
    if yticks:
        for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
    else:
        ax.set_yticklabels([])
    if xticks:
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
    else:
        ax.set_xticklabels([])
    ax.set_title(title, fontsize=14)

    return bins.min(), bins.max(), min(n), max(n)

def make_fig_histogram(results, path, title):
    mice = results.keys()
    fig, ax = plt.subplots(1, len(mice), figsize=(len(mice)//2*5, 5))
    bins = []
    counts = []
    xticks = True
    for i, mouse in enumerate(mice):
        yticks = not i
        new_title = "%s %s" % (mouse[-4:], title)
        intervals = results[mouse]
        minb, maxb, minc, maxc = make_single_histogram(ax,
                                                       intervals,
                                                       new_title,
                                                       xticks=xticks,
                                                       yticks=yticks,
                                                       xlogscale=True)
        bins.extend([minb, maxb])
        counts.extend([minc, maxc])
    min_bin, max_bin = min(bins), max(bins)
    min_count, max_count = min(counts), max(counts)
    for x in ax:
        x.set_xlim([min_bin, max_bin + 1])
        x.set_ylim([min_count, max_count + 3])
    fig.subplots_adjust(wspace=0.15)
    fig.savefig(path, dpi=300)
    print(path)


def make_pooled_histograms_for_every_mouse(results, fname,
                                           mice, main_directory,
                                           directory, prefix,
                                           additional_info):
    results_following = pool_results_following(results, mice)
    results_followed = pool_results_followed(results, mice)
    new_name_following = "pooled_following"
    new_name_followed = "pooled_followed"

    new_dir = utils.check_directory(directory, 'figs')
    dir_name =  utils.check_directory(main_directory, new_dir)

    if prefix != "":
        new_name_following = '%s_%s_%s.png'% (fname, prefix, new_name_following)
        new_name_followed = '%s_%s_%s.png'% (fname, prefix, new_name_followed)
    else:
        new_name_following =  '%s_%s.png'% (fname, new_name_following)
        new_name_followed =  '%s_%s.png'% (fname, new_name_followed)

    fname_following = os.path.join(dir_name, new_name_following)
    fname_followed = os.path.join(dir_name, new_name_followed)
    make_fig_histogram(results_following, fname_following,
                       "following")
    make_fig_histogram(results_followed, fname_followed,
                       "followed")
