# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:38:58 2017

@author: Jan Maka, Joanna Jędrzejewska-Szmek
"""
from __future__ import division, print_function, absolute_import
import os
import numpy as np
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
from . import utility_functions as utils


nbins = 10


def make_labels(my_mice):
    if len(set([mouse[-4:] for mouse in my_mice])) == len(my_mice):
        return [mouse[-4:] for mouse in my_mice]
    return [mouse[:-5] for mouse in my_mice]


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
                    title=None,
                    symmetrical=True,
                    full_dir_tree=True):
    mice = make_labels(old_mice)

    if len(set(old_mice)) != len(set(mice)):
        return
    if full_dir_tree:
        subdirectory = os.path.join(subdirectory, 'raster_plots')
    new_path = utils.check_directory(main_directory, subdirectory)
    if title:
        plt.suptitle(title, fontsize=14, fontweight='bold')
    assert FAM.shape[0] == len(phases)
    assert FAM.shape[1] == len(mice)
    assert FAM.shape[2] == len(mice)

    if symmetrical:
        output, pair_labels = utils.make_table_of_pairs(FAM,
                                                        phases,
                                                        mice)
    else:
        output, pair_labels = utils.make_table_of_all_mouse_pairs(FAM,
                                                                  phases,
                                                                  mice)
    fig = plt.figure(figsize=(len(phases), 0.5*FAM.shape[0]))
    ax = fig.add_subplot(111, aspect='equal')

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
    if prefix:
        name = "%s_%s" % (prefix, name)
    new_name = os.path.join(new_path, name)
    print(new_name)
    plt.savefig(new_name+".png",
                transparent=False,
                bbox_inches=None,
                pad_inches=.5,
                dpi=100)
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
    name = '%s_%s_%s' % (name, prefix, phase)
    fig, ax = plt.subplots()
    if vmin is None:
        vmin = result.min()
    if vmax is None:
        vmax = result.max()

    cax = ax.imshow(result, interpolation='none', aspect='auto',
                    cmap="viridis", origin="lower", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cax)
    turn = False
    if not xticks:
        xticks = [mouse.split('-')[-1] for mouse in mice]
        turn = True
    if not yticks:
        yticks = [mouse.split('-')[-1] for mouse in mice]
        fig.subplots_adjust(left=0.25)

    ax.get_yaxis().set_ticks([i for i, x in enumerate(yticks)])
    ax.get_xaxis().set_ticks([i for i, x in enumerate(xticks)])
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
    fig.savefig(new_name + ".png", transparent=False, bbox_inches=None,
                pad_inches=0.5,  dpi=100)
    plt.close(fig)


def single_timeline_heat_map(result,
                             directory,
                             mice,
                             prefix,
                             phase,
                             binsize,
                             antenna,
                             subdirectory):
    name = 'duration_of_antenna_%s_registration_%s_%4.2f_%s' % (antenna,
                                                                prefix,
                                                                binsize/3600,
                                                                phase)
    fig, ax = plt.subplots()
    vmin = 0
    vmax = 0.5*binsize
    length = len(result[mice[0]])
    out = np.zeros((len(mice), length))
    for j, mouse in enumerate(mice):
        out[j] = np.array(result[mouse])

    cax = ax.imshow(out, interpolation='none', aspect='auto', cmap="viridis",
                    origin="lower", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(cax)
    turn = False
    yticks = [mouse.split('-')[-1] for mouse in mice]

    fig.subplots_adjust(left=0.25)
    ax.get_yaxis().set_ticks([i for i, x in enumerate(yticks)])
    ax.set_yticklabels(yticks)
    ax.set_xlabel("time (h)")
    ax.set_title("In range of antenna %s (in sec)" % antenna)
    if subdirectory:
        subdirectory = os.path.join(subdirectory, 'figs')
    dir_name = utils.check_directory(directory, subdirectory)
    new_name = os.path.join(dir_name, name)
    print(new_name)
    fig.subplots_adjust(left=0.25)
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(new_name + ".png", transparent=False, bbox_inches=None,
                pad_inches=0.5,  dpi=100)
    plt.close(fig)


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
                              labels=['', ''],
                              full_dir_tree=True):
    if full_dir_tree:
        new_name = os.path.join(directory, 'figs')
    else:
        new_name = directory
    directory = utils.check_directory(main_directory, new_name)
    fname = os.path.join(directory, '%s_%s_%s' % (fname, prefix, phase))
    print(fname)
    label_mice = make_labels(mice)
    fig = plt.figure(figsize=(10, 6))
    ax = []
    for i in range(1, 5):
        ax.append(fig.add_subplot(2, 2, i))
    im2 = ax[0].imshow(results, vmin=vmin, vmax=vmax,
                       interpolation='none', origin='lower')
    cbar = fig.colorbar(im2)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_title(titles[0])
    ax[1].imshow(results_exp, vmin=vmin, vmax=vmax,
                 interpolation='none', origin='lower')
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_title(titles[1])
    deltas = results[results > 0] - results_exp[results > 0]

    try:
        if -abs(vmin1) != abs(vmax1):
            maxi = max(abs(vmax1), abs(vmin1))
            mini = -maxi
        else:
            mini = vmin1
            maxi = vmax1
        im = ax[2].imshow(results - results_exp,
                          vmin=mini,
                          vmax=maxi,
                          interpolation='none',
                          cmap='bwr',
                          origin='lower')
        ax[2].set_title(titles[2])
        cbar = fig.colorbar(im)
        ax[2].yaxis.tick_left()
        ax[2].get_yaxis().set_ticks([i for i, x in enumerate(mice)])
        ax[2].get_xaxis().set_ticks([i for i, x in enumerate(mice)])
        ax[2].set_xticklabels(label_mice)
        ax[2].set_yticklabels(label_mice)
        ax[2].set_xlabel(labels[0])
        ax[2].set_ylabel(labels[1])
        for label in ax[2].xaxis.get_ticklabels():
            label.set_rotation(90)

    except ValueError:
        pass
    try:
        ax[3].hist(deltas, density=1)
    except ValueError:
        pass
    ax[3].set_title(titles[3])
    if hist is True:
        ax[3].set_xlim([-0.1, 0.5])
        ax[3].get_xaxis().set_ticks([-0.1, 0., 0.1, 0.2, 0.3])
        ax[3].set_xticklabels([-0.1, 0., 0.1, 0.2, 0.3])
    else:
        ax[3].set_xlim([vmin1, vmax1])
        ticks = np.linspace(vmin1, vmax1, 4)
        ax[3].get_xaxis().set_ticks(ticks)
        ax[3].set_ylim([0, 1])
        if vmax1 > 20:
            ax[3].set_xticklabels([np.round(tick) for tick in ticks])
        elif vmax1 < 0.01:
            ax[3].set_xticklabels([np.round(tick, 3) for tick in ticks])

    fig.suptitle(phase)
    fig.subplots_adjust(bottom=0.1)
    fig.subplots_adjust(wspace=0.25)
    fig.subplots_adjust(hspace=0.3)
    fig.savefig(fname+'.png', dpi=100,  bbox_inches=None,
                pad_inches=2)
    plt.close(fig)
    print(fname+'.png')


def pooled_hists(res, res_exp, phases, fname, main_directory, directory,
                 prefix, additional_info, full_dir_tree=True):

    fig, ax = plt.subplots(1, len(phases),
                           figsize=(2.5 + (len(phases) - 1)//2*5, 5))
    if len(phases) == 1:
        ax = [ax]
        new_phase = phases[0].replace(' ', '__')
    else:
        new_phase = "%s_%s" % (phases[0].replace(' ', '__'),
                               phases[-1].replace(' ', '__'))
    bins, counts = [], []
    xticks = True
    for i, phase in enumerate(phases):
        results = utils.dict_to_array_2D(res[phase][0],
                                         list(res[phase][0].keys()),
                                         list(res[phase][0].keys()))
        results_exp = utils.dict_to_array_2D(res_exp[phase][0],
                                             list(res[phase][0].keys()),
                                             list(res[phase][0].keys()))

        deltas = results[results > 0] - results_exp[results > 0]
        new_title = phases[i]
        yticks = not i
        minb, maxb, minc, maxc = make_single_histogram(ax[i],
                                                       deltas,
                                                       10,
                                                       title=new_title,
                                                       xticks=xticks,
                                                       yticks=yticks,
                                                       xlogscale=False)
        bins.extend([minb, maxb])
        counts.extend([minc, maxc])
    min_bins, max_bins = min(bins), max(bins)
    min_count, max_count = min(counts), max(counts)
    for x in ax:
        x.set_xlim([min_bins, max_bins + 1])
        x.set_ylim([min_count, max_count + 3])
    if full_dir_tree:
        directory = os.path.join(directory, 'figs')
    directory = utils.check_directory(main_directory, directory)
    fname = os.path.join(directory, '%s_%s_%s' % (fname, prefix,
                                                  new_phase))
    if len(phases) > 1:
        fig.subplots_adjust(wspace=0.15)
    fig.savefig(fname + '.png', dpi=100, bbox_inches=None,
                pad_inches=0.5)
    plt.close(fig)


def make_histograms_for_every_mouse(results, fname, mice, main_directory,
                                    directory, prefix, additional_info,
                                    full_dir_tree=True):
    """
    results should be a dictionary of lists
    """
    fig, ax = plt.subplots(len(mice), len(mice),
                           figsize=(len(mice)//2*5, len(mice)//2*5))
    bins, counts = [], []
    new_name = "all_mice"
    xlabel = None
    for i, mouse1 in enumerate(mice):
        for j, mouse2 in enumerate(mice):
            if mouse1 == mouse2:
                if i == 0:
                    ax[i, j].set_ylabel(mouse1, fontsize=14)
                    ax[i, j].set_title(mouse2, fontsize=14)
                ax[i, j].set_yticklabels([])
                ax[i, j].set_xticklabels([])
                continue
            if j:
                yticks = False
                ylabel = None
            else:
                yticks = True
                ylabel = mouse1
            if i == len(mice) - 1:
                xticks = True
            else:
                xticks = False
            if i == 0:
                new_title = mouse2
            else:
                new_title = ""

            key = "%s|%s" % (mouse1, mouse2)
            intervals = results[key]
            minb, maxb, minc, maxc = make_single_histogram(ax[i, j],
                                                           intervals,
                                                           30,
                                                           title=new_title,
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
    if full_dir_tree:
        directory = os.path.join(directory, 'figs')
    dir_name = utils.check_directory(main_directory, directory)
    plt.gcf().text(0.02, 0.5, "Followed mouse", fontsize=28, rotation=90)
    plt.gcf().text(0.5, 0.02, "Following mouse", fontsize=28)

    if prefix != "":
        new_name = '%s_%s_%s' % (fname, prefix, new_name)
    else:
        new_name = '%s_%s' % (fname, new_name)

    fname = os.path.join(dir_name, new_name)
    print(fname)
    fig.subplots_adjust(wspace=0.15)
    fig.savefig(fname + ".png",
                bbox_inches=None,
                pad_inches=.5,
                dpi=100)
    plt.close(fig)


def pool_results_following(res_dict, mice):
    pooled_results = {mouse: [] for mouse in mice}
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            key = "%s|%s" % (mouse1, mouse2)
            pooled_results[mouse2] += res_dict[key]
    return pooled_results


def pool_results_followed(res_dict, mice):
    pooled_results = {mouse: [] for mouse in mice}
    for mouse1 in mice:
        for mouse2 in mice:
            if mouse1 == mouse2:
                continue
            key = "%s|%s" % (mouse1, mouse2)
            pooled_results[mouse1] += res_dict[key]
    return pooled_results


def single_histogram_figures(single_results, fname, main_directory,
                             path, title, nbins=10,
                             xlabel=None, ylabel=None,
                             fontsize=14, xlogscale=False,
                             ylogscale=False, xmin=None, xmax=None,
                             ymin=None, ymax=None,
                             median_mean=False, add_text="",
                             full_dir_tree=True):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    if full_dir_tree:
        new_dir = os.path.join(path, 'figs')
    else:
        new_dir = path
    dir_name = utils.check_directory(main_directory, new_dir)
    new_fname = os.path.join(dir_name, fname)
    if nbins is False:
        nbins = int(max(single_results))
    if add_text:
        title += add_text

    make_single_histogram(ax, single_results, nbins, title=title,
                          xticks=True,
                          yticks=True, xlabel=xlabel, ylabel=ylabel,
                          xlogscale=xlogscale, ylogscale=ylogscale,
                          fontsize=fontsize, xmin=xmin, xmax=xmax,
                          ymin=ymin, ymax=ymax,
                          median_mean=median_mean)
    fig.savefig(new_fname + ".png",
                bbox_inches=None,
                pad_inches=.5,
                dpi=100)
    print(new_fname)
    plt.close(fig)


def make_single_histogram(ax, single_results, nbins, title="", xticks=False,
                          yticks=False, xlabel=None, ylabel=None,
                          xlogscale=False, ylogscale=False, fontsize=14,
                          xmin=None, xmax=None, ymin=None, ymax=None,
                          median_mean=False):
    if len(single_results) == 0:
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        return 0, 0, 0, 0
    if xlogscale:
        hist, bins = np.histogram(single_results, nbins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
        n, bins, patches = ax.hist(single_results, bins=logbins)
        ax.set_xscale("log")
    else:
        n, bins, patches = ax.hist(single_results, bins=nbins)

    if ylogscale:
        ax.set_yscale('log')
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    if yticks:
        for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize)
    else:
        ax.set_yticklabels([])
    if xticks:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
    else:
        ax.set_xticklabels([])
    ax.set_title(title, fontsize=fontsize)
    if xmin is None:
        xmin = min(bins) - 0.5
    if xmax is None:
        xmax = max(bins) + 0.5
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = max(n)+1
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    if median_mean is True:
        mean = np.mean(single_results)
        median = np.median(single_results)
        ax.axvline(mean, color='k', linestyle='dashed', linewidth=1)
        ax.axvline(median, color='r', linestyle='dashed', linewidth=1)
        ylims = ax.get_ylim()
        ax.text(mean*0.9, (ylims[1]-ylims[0])*0.8, "mean = %4.1f" % mean)
        ax.text(median*0.9, (ylims[1]-ylims[0])*0.9, "median = %4.1f" % median)
    return bins.min(), bins.max(), min(n), max(n)


def make_fig_histogram(results, path, title):
    mice = list(results.keys())
    fig, ax = plt.subplots(1, len(mice), figsize=(len(mice)//2*5, 5))
    bins = []
    counts = []
    xticks = True
    for i, mouse in enumerate(mice):
        yticks = not i
        new_title = "%s %s" % (title, mouse)
        intervals = results[mouse]
        minb, maxb, minc, maxc = make_single_histogram(ax[i],
                                                       intervals,
                                                       30,
                                                       title=new_title,
                                                       xticks=xticks,
                                                       yticks=yticks,
                                                       xlogscale=True)
        bins.extend([minb, maxb])
        counts.extend([minc, maxc])
    min_bin, max_bin = min(bins), max(bins)
    min_count, max_count = min(counts), max(counts)
    for x in ax:
        if min_bin == max_bin == max_count == min_count:
            continue
    fig.subplots_adjust(wspace=0.15)
    fig.savefig(path + ".png", dpi=100,
                bbox_inches=None,
                pad_inches=.5)
    print(path)
    plt.close(fig)


def pooled_hists_for_every_mouse(results, fname, mice, main_directory,
                                 directory, prefix, additional_info,
                                 full_dir_tree=True):
    results_following = pool_results_following(results, mice)
    results_followed = pool_results_followed(results, mice)
    new_name_following = "pooled_following"
    new_name_followed = "pooled_followed"
    if full_dir_tree:
        directory = os.path.join(directory, 'figs')
    dir_name = utils.check_directory(main_directory, directory)

    if prefix != "":
        new_name_following = '%s_%s_%s' % (fname, prefix, new_name_following)
        new_name_followed = '%s_%s_%s' % (fname, prefix, new_name_followed)
    else:
        new_name_following = '%s_%s' % (fname, new_name_following)
        new_name_followed = '%s_%s' % (fname, new_name_followed)

    fname_following = os.path.join(dir_name, new_name_following)
    fname_followed = os.path.join(dir_name, new_name_followed)
    make_fig_histogram(results_following, fname_following,
                       "following")
    make_fig_histogram(results_followed, fname_followed,
                       "followed")


def make_visit_duration_histogram(results, time, phase, mice,
                                  fname, main_directory,
                                  directory, prefix, additional_info):

    dir_name = os.path.join(main_directory, directory)
    dir_name = utils.check_directory(dir_name, "figs")

    for mouse in mice:
        if prefix != "":
            new_name = '%s_%s_%s_%s' % (fname, mouse, prefix, phase)
        else:
            new_name = '%s_%s_%s' % (fname, mouse, phase)
        new_name = os.path.join(dir_name, new_name)
        ncols = len(results.keys())
        nrows = len(results[list(results.keys())[0]][mouse])
        fig, ax = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*2.5))
        fontsize = 14 + nrows
        if nrows == 1:
            ax = np.expand_dims(ax, 0)
        bins, counts = [], []
        for i, key in enumerate(sorted(results.keys())):
            for j, out in enumerate(results[key][mouse]):
                if i:
                    yticks = False
                    ylabel = None
                else:
                    yticks = True
                    ylabel = "# at %2.0f h" % time[j]
                if j == nrows - 1:
                    xticks = True
                    xlabel = "visit duration"
                else:
                    xticks = False
                    xlabel = None

                minb, maxb, minc,\
                    maxc = make_single_histogram(ax[j, i], out,
                                                 10,
                                                 title=key,
                                                 xticks=xticks,
                                                 yticks=yticks,
                                                 xlabel=xlabel,
                                                 ylabel=ylabel,
                                                 xlogscale=True,
                                                 fontsize=fontsize)
                bins.extend([minb, maxb])
                counts.extend([minc, maxc])
        min_bin, max_bin = min(bins), max(bins)
        min_count, max_count = min(counts), max(counts)
        for s in ax:
            for x in s:
                x.set_xlim([min_bin, max_bin + 1])
                x.set_ylim([min_count, max_count + 3])
        fig.suptitle("%s visit durations in %s" % (mouse, phase),
                     fontsize=fontsize+2)
        fig.subplots_adjust(wspace=0.15)
        fig.savefig(new_name + ".png", dpi=100,
                    bbox_inches=None,
                    pad_inches=.5)
        print(new_name)
        plt.close(fig)


def histograms_antenna_transitions(t_times, setup_config,
                                   main_directory, directory, fname, prefix,
                                   xmin=None, xmax=None, ymin=None, ymax=None,
                                   additional_info=""):
    dir_correct_same = os.path.join(directory, "correct_antenna_transition",
                                    "same_antenna")
    dir_correct_cage = os.path.join(directory,
                                    "correct_antenna_transition",
                                    "cages")
    dir_correct_tunnel = os.path.join(directory,
                                      "correct_antenna_transition",
                                      "tunnels")
    dir_incorrect = os.path.join(directory,
                                 "incorrect_antenna_transition")
    cage_pd = setup_config.cage_pair_dict()
    tunnel_pd = setup_config.tunnel_pair_dict()
    title_double_a = "consecutive reading by antenna %s (%s) %s at %s"
    title_other = "passing through %s phase %s at %s"
    fname = "%s_%s_%s" % (prefix, fname, additional_info)
    max_count = 0
    nbins = {}
    title = {}
    dest_fname = {}
    dir_name = {}
    xlogscale = {}
    new_xmin = 1000
    new_xmax = 0
    incorrect_transitions = []
    allowed_keys = {}

    for ph in t_times.keys():
        allowed_keys[ph] = []
        xlogscale[ph] = {}
        nbins[ph] = {}
        dest_fname[ph] = {}
        title[ph] = {}
        for lab in t_times[ph].keys():
            nbins[ph][lab] = {}
            dest_fname[ph][lab] = {}
            title[ph][lab] = {}
            xlogscale[ph][lab] = {}
            for key in t_times[ph][lab].keys():
                if not len(t_times[ph][lab][key]):
                    continue
                try:
                    first, last = key.split(" ")
                    gen_key = "%s %s" % (min(first, last), max(first, last))
                    if gen_key in setup_config.mismatched_pairs:
                        incorrect_transitions.extend(t_times[ph][lab][key])
                        continue
                    else:
                        allowed_keys[ph].append(key)
                    if first == last:
                        dir_name[key] = dir_correct_same
                        title[ph][lab][key] = title_double_a %\
                            (first,
                             setup_config.address[first],
                             ph,
                             lab)
                        dest_fname[ph][lab][key] =\
                            "%s_%s_%s_%s_start_at_%d" %\
                            (fname,
                             setup_config.address[first].replace(" ", "_"),
                             key.replace(" ", "_"),
                             ph.replace(" ", "_"),
                             lab)
                    else:
                        if key in setup_config.tunnel_pairs():
                            dir_name[key] = dir_correct_tunnel
                            title[ph][lab][key] = title_other %\
                                (tunnel_pd[key],
                                 ph,
                                 lab)
                            dest_fname[ph][lab][key] =\
                                "%s_%s_%s_%s_start_at_%d" %\
                                (fname,
                                 tunnel_pd[key].replace(" ", "_"),
                                 key.replace(" ", "_"),
                                 ph.replace(" ", "_"),
                                 lab)

                        else:
                            dir_name[key] = dir_correct_cage
                            title[ph][lab][key] = title_other %\
                                (cage_pd[key],
                                 ph,
                                 lab)
                            dest_fname[ph][lab][key] =\
                                "%s_%s_%s_%s_start_at_%s" %\
                                (fname,
                                 cage_pd[key].replace(" ", "_"),
                                 key.replace(" ", "_"),
                                 ph.replace(" ", "_"),
                                 lab)

                except ValueError:
                    allowed_keys[ph].append(key)
                    if key == "tunnels":
                        dir_name[key] = dir_correct_tunnel
                    elif key == "cages":
                        dir_name[key] = dir_correct_cage
                    title[ph][lab][key] = "%s phase %s at %s" %\
                        (key, ph, lab)
                    dest_fname[ph][lab][key] = "%s_%s_%s_start_at_%s" %\
                        (fname,
                         key.replace(" ", "_"),
                         ph.replace(" ", "_"),
                         lab)
                xlogscale[ph][lab][key] = True
                nbins[ph][lab][key] = 10
                if len(t_times[ph][lab][key]) > 1000:
                    nbins[ph][lab][key] = 40
                while True:
                    if 0 in t_times[ph][lab][key]:
                        t_times[ph][lab][key].remove(0)
                    else:
                        break
                hist, bins = np.histogram(t_times[ph][lab][key],
                                          nbins[ph][lab][key])
                if xlogscale[ph][lab][key]:
                    logbins = np.logspace(np.log10(bins[0]),
                                          np.log10(bins[-1]),
                                          len(bins))
                    hist, bins = np.histogram(t_times[ph][lab][key],
                                              bins=logbins)
                if max(hist) > max_count:
                    max_count = max(hist)
                if new_xmin > min(t_times[ph][lab][key]):
                    new_xmin = min(t_times[ph][lab][key])
                if new_xmax < max(t_times[ph][lab][key]):
                    new_xmax = max(t_times[ph][lab][key])
    if xmin is None:
        xmin = new_xmin
    if xmax is None:
        xmax = new_xmax
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = max_count
    for ph in t_times.keys():
        for lab in t_times[ph].keys():
            for key in allowed_keys[ph]:
                if not len(t_times[ph][lab][key]):
                    continue
                single_histogram_figures(t_times[ph][lab][key],
                                         dest_fname[ph][lab][key],
                                         main_directory, dir_name[key],
                                         title[ph][lab][key],
                                         nbins=nbins[ph][lab][key],
                                         xlogscale=xlogscale[ph][lab][key],
                                         xlabel="Transition times (s)",
                                         ylabel="count", xmin=xmin, xmax=xmax,
                                         ymin=0, ymax=ymax,
                                         fontsize=14, median_mean=True)
    single_histogram_figures(incorrect_transitions, "incorrect_transitions",
                             main_directory, dir_incorrect,
                             "Incorrect antenna transitions",
                             nbins=30,
                             xlogscale=True,
                             xlabel="Transition times (s)",
                             ylabel="count",
                             fontsize=14, median_mean=True)
    return incorrect_transitions
