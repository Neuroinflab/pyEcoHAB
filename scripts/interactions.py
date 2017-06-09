#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:48:42 2016

@author: Jan Maka
"""
# EcoHAB libraries
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from experiments_info import smells, antenna_positions
# from plotfunctions import *
import plotfunctions
from data_managment import *
# Third party libraries
import networkx as nx
import numpy as np
import os
import pickle
# Debug libraries
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcol
import matplotlib.patches as patches
from time import sleep
from datetime import datetime

import threading

threads = []
started = 0
finished = 0


class Experiment(object):
    """Class, which represent data from one of the experiments"""

    def __init__(self, path, exp_name=None):
        if not exp_name:
            exp_name = path
        self.path = path
        self.directory = os.path.join('..', 'Results/' + path)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.ehd = EcoHab.EcoHabData(os.path.join('..', 'RawData', exp_name), _ant_pos=antenna_positions[exp_name])
        self.ehs = EcoHab.EcoHabSessions9states(self.ehd, shortest_session_threshold=0)
        self.fs = self.ehs.fs
        self.sd = self.ehs.signal_data
        self.cf = ExperimentConfigFile(os.path.join('..', 'RawData', exp_name))
        mice = list(self.ehd.mice)
        self.mice = filter(lambda x: len(self.ehs.getstarttimes(x)) > 30, mice)
        self.lm = len(self.mice)

    def calculate_fvalue(self, window='default', treshold=2, min_interactions=1, force=False, fols=None, ops=None):
        self.fols = fols
        self.ops = ops
        self.treshold = treshold
        self.min_interactions = min_interactions
        self.tstart, self.tend = self.cf.gettime('ALL')
        self.fpatterns = []
        self.opatterns = []
        if window == 'default':
            sessions = filter(lambda x: x.endswith('dark') or x.endswith('light'), self.cf.sections())
            self.phases = [(self.cf.gettime(sec)[0] - self.tstart, self.cf.gettime(sec)[1] - self.tstart) for sec in
                           sessions]
        else:
            if isinstance(window, float) or isinstance(window, int):
                self.phases = [(i * window * 3600, np.min([(i + 1) * window * 3600, self.sd.shape[0]])) for i in
                               range(int(np.ceil(self.sd.shape[0] * 1.0 / (window * 3600 * self.fs))))]
            elif isinstance(window, list):
                self.phases = [(st * window[0] * 3600, (st + 1) * window[0] * 3600) for st in window[1]]
            else:
                raise TypeError
        self.f = np.zeros((len(self.mice), len(self.mice), len(self.phases)))
        self.interactions = np.zeros((len(self.phases), len(self.mice), len(self.mice), 8, 3))
        self.f_sum = np.zeros((len(self.phases), self.lm))
        if 'following.npy' in list(os.walk(os.path.join('..', 'Results', self.path)))[0][2] and not force:
            self.f = np.load(os.path.join('..', 'Results/' + self.path + '/') + 'following.npy')
            return self.f
        else:
            for s in range(len(self.phases)):
                ts, te = self.phases[s]
                print 'Phase %s. from %sh, to %sh' % (s + 1, np.round(ts / 3600., 2), np.round(te / 3600., 2))
                self.interactions[s, :, :, :, :] = self.interaction_matrix(ts, te)
            if not os.path.exists(os.path.join('..', 'PreprocessedData/InteractionsData/')):
                os.makedirs(os.path.join('..', 'PreprocessedData/InteractionsData/'))
            np.save(os.path.join('..', 'PreprocessedData/InteractionsData/') + '%s.npy' % self.path, self.interactions)
            return self.f

    def interaction_matrix(self, ts, te):
        """Calculates fvalue matrix for the given period from ts to te
        """
        imatrix = np.zeros((self.lm, self.lm, 8, 3))
        for ii in range(len(self.mice)):
            for jj in range(len(self.mice)):
                if ii < jj:
                    imatrix[ii, jj, :, :], patterns = self.findpatterns((ii, jj), ts, te)
                    self.fpatterns += patterns[0]
                    self.opatterns += patterns[1]
                    imatrix[jj, ii, :, :], patterns = self.findpatterns((jj, ii), ts, te)
                    self.fpatterns += patterns[0]
                    self.opatterns += patterns[1]
        return imatrix

    def validatePatterns(self, plots_nr=9, trange=[-3, 3]):
        size = np.ceil(np.sqrt(plots_nr))
        t = np.arange(trange[0], trange[1], 1.0 / self.fs)
        frandom_idx = [np.random.randint(0, len(self.fpatterns) - 1) for i in range(plots_nr)]
        orandom_idx = [np.random.randint(0, len(self.opatterns) - 1) for i in range(plots_nr)]
        plt.suptitle("Random following patterns", fontsize=14, fontweight='bold')
        for i, idx in enumerate(frandom_idx):
            ax = plt.subplot(size, size, i + 1)
            print self.fpatterns[idx]
            ii, jj, s = self.fpatterns[idx]
            ax.set_title("%s|%s|t=%s" % (ii, jj, s * 1. / self.fs))
            plt.plot(t, self.sd[s - 3 * self.fs:s + 3 * self.fs, ii] - 0.05, 'ro', label="leader")
            plt.plot(t, self.sd[s - 3 * self.fs:s + 3 * self.fs, jj] + 0.05, 'bo', label="follower")
            plt.axis([-3.1, 3.1, -0.5, 9.5])
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()
        plt.suptitle("Random avoiding patterns", fontsize=14, fontweight='bold')
        for i, idx in enumerate(orandom_idx):
            ax = plt.subplot(size, size, i + 1)
            print self.opatterns[idx]
            ii, jj, s = self.opatterns[idx]
            ax.set_title("%s|%s|t=%s" % (ii, jj, s * 1. / self.fs))
            plt.plot(t, self.sd[s - 3 * self.fs:s + 3 * self.fs, ii] - 0.05, 'ro', label="leader")
            plt.plot(t, self.sd[s - 3 * self.fs:s + 3 * self.fs, jj] + 0.05, 'bo', label="follower")
            plt.axis([-3.1, 3.1, -0.5, 9.5])
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

    def findpatterns(self, (m1, m2), t1, t2):
        fs = self.fs
        sd = self.sd
        detected_idx = [[], []]
        follow_stat = {}
        m1_idx = np.where(((np.roll(sd[:, m1], 1) - sd[:, m1]) != 0))[0]
        m2_stats = self.ehs.statistics[self.mice[m2]]["preference"]
        moves = ['24', '42', '46', '64', '68', '86', '82', '28']
        for m in moves:
            follow_stat[m] = np.zeros(3)
            start_st, end_st = int(m[0]), int(m[1])
            if end_st == 2 and start_st == 8:
                index = 0
            elif end_st == 8 and start_st == 2:
                index = 1
            elif start_st < end_st:
                index = 0
            else:
                index = 1
            follow_stat[m][2] = m2_stats[int(m[0])][index] / np.sum(m2_stats[int(m[0])])
        # print m2_stats
        # print 'ruchliwosc', len(m1_idx)
        for i in range(2, len(m1_idx)):
            if m1_idx[i] > t1 * fs and m1_idx[i] < t2 * fs:
                s = m1_idx[i]
                start_st = int(sd[m1_idx[i - 2], m1])
                end_st = int(sd[s, m1])
                e = s + self.treshold * fs
                try:
                    period1 = list(sd[s:e, m2])
                    period2 = list(sd[s - 2 * fs:s, m2])
                    period3 = list(sd[s - int(0.1 * fs):s, m2])
                    # define conditions
                    unknown_state = sd[s, m1] == 0 or int(sd[m1_idx[i - 1], m1]) == 0
                    unknown_previous_states = (sd[m1_idx[i - 1], m1] == 0 and (sd[m1_idx[i - 2], m1] == 0)) or (
                    sd[m1_idx[i - 1], m1] != 0 and (sd[m1_idx[i - 2], m1] == 0))
                    in_pipe = sd[s, m1] % 2 == 1
                    if unknown_state or unknown_previous_states or in_pipe or start_st == end_st:
                        continue

                    op_idx = (2 * sd[m1_idx[i - 2], m1] - sd[s, m1] - 1) % 8 + 1
                    same_start = sd[m1_idx[i - 2], m1] in period2  # POPRAWIC!!!!!!!
                    first_m1 = sd[s, m1] not in period3 and op_idx not in period3
                    followed = period1.count(sd[s, m1]) > 0
                    go_oposite = op_idx in period1
                    if followed and go_oposite:
                        followed = period1.index(sd[s, m1]) < period1.index(op_idx)
                        go_oposite = not followed
                    if same_start and first_m1:
                        # print start_st,end_st, op_idx
                        # print start_st,end_st,index
                        # print sd[m1_idx[i-2],m1],sd[s,m1], (2*sd[m1_idx[i-2],m1]-sd[s,m1]-1)%8+1,followed,go_oposite
                        if followed:
                            #####POPRAWIC
                            # print np.ceil((period1.index(sd[s,m1]))/self.fs)
                            if self.fols != None:
                                self.fols[np.ceil((period1.index(sd[s, m1])) / self.fs)] += 1
                            follow_stat[str(start_st) + str(end_st)][0] += 1
                            # print p, sd[m1_idx[i-2],m1],[index]
                            detected_idx[0].append((m1, m2, s))
                        elif go_oposite:
                            # print np.ceil((period1.index(op_idx))/self.fs)
                            if self.ops != None:
                                self.ops[np.ceil((period1.index(op_idx)) / self.fs)] += 1
                            follow_stat[str(start_st) + str(end_st)][1] += 1
                            detected_idx[1].append((m1, m2, s))
                except IndexError:
                    print 'Err'
                    continue
        return np.array(follow_stat.values()), detected_idx


def createRandomExpepiments(paths, core='Random_test'):
    if not os.path.exists('../PreprocessedData/RandomData/'):
        os.makedirs('../PreprocessedData/RandomData/')
    statistics_rd = {}
    for j in range(len(paths)):
        rdE = Experiment('Random_test%s' % (j + 1), exp_name=core)
        rdname = 'Random_test%s.pkl' % (j + 1)
        print rdE.sd.shape
        rd = np.zeros(rdE.sd.shape)
        print '####%s#####' % j
        ######Create test file#########
        for i in range(j, len(paths) + j):
            E = Experiment(paths[i % len(paths)])
            print i, i % E.lm, E.sd.shape
            end = np.min([len(rd[:, i % rdE.lm]), len(E.sd[:, i % E.lm])])
            rd[:end, i % rdE.lm] = E.sd[:end, i % E.lm]
            statistics_rd[rdE.mice[i % rdE.lm]] = E.ehs.statistics[E.mice[i % E.lm]]
            # print statistics_rd[E_t.mice[i%E_t.lm]]["preference"]
        rdE.mice = rdE.mice[:len(paths)]
        rdE.lm = len(paths)
        rdE.sd = rd[:, :end]
        rdE.ehs.statistics = statistics_rd
        with open('../PreprocessedData/RandomData/' + rdname, 'wb') as output:
            pickle.dump(rdE, output, pickle.HIGHEST_PROTOCOL)
            del rdE


def singleThreadProcessData(key, path, window, ts):
    global started, finished
    started += 1
    print key, window
    print path
    if key == "RD":
        with open('../PreprocessedData/RandomData/' + path + '.pkl', "rb") as input_file:
            E = pickle.load(input_file)
    else:
        E = Experiment(path)
    E.calculate_fvalue(window=window, treshold=ts, force=True)
    # E.validatePatterns()
    finished += 1


def preprocessData(names, window, ts=3):
    global threads
    for key in names.keys():
        for path in names[key]:
            threads.append(threading.Thread(target=singleThreadProcessData, args=[key, path, window, ts]))
            print ("yay, new thread")
            threads[-1].daemon = True
            threads[-1].start()
    while started != finished:
        sleep(1)
        """
        for path in names[key]:
            print path
            if key=="RD":
                with open('../PreprocessedData/RandomData/'+path+'.pkl', "rb") as input_file:
                    E = pickle.load(input_file)
            else:
                E = Experiment(path)
            E.calculate_fvalue(window = window, treshold =ts, force=True)
            E.validatePatterns()
        """


def Interpersec(names, ts=3, directory='InterPerSec'):
    if not os.path.exists('../Results/%s/' % directory):
        os.makedirs('../Results/%s/' % directory)
    for key in names.keys():
        fols = np.zeros(ts - 1)
        ops = np.zeros(ts - 1)
        for path in names[key]:
            print path
            if key == "RD":
                with open('../PreprocessedData/RandomData/' + path + '.pkl', "rb") as input_file:
                    E = pickle.load(input_file)
            else:
                E = Experiment(path)
            E.calculate_fvalue(window=24, treshold=ts, force=True, fols=fols, ops=ops)
            fols = E.fols
            ops = E.ops
        print '###########%s###########' % key
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1, ts, 1), fols * 1. / np.sum(fols))
        plt.ylim(0, 0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of followings")
        plt.savefig('../Results/' + directory + '/' + key + 'folinsec.png')
        plt.show()
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1, ts, 1), ops * 1. / np.sum(ops))
        plt.ylim(0, 0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of avoidance")
        plt.savefig('../Results/' + directory + '/' + key + 'opsinsec.png')
        plt.show()
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1, ts, 1), (fols + ops) * 1. / np.sum(fols + ops))
        plt.ylim(0, 0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of interactions")
        plt.savefig('../Results/' + directory + '/' + key + 'interinsec.png')
        plt.show()


def InteractionsPerPair(names):
    IPP = {}
    for key in names.keys():
        IPP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..', 'PreprocessedData/InteractionsData/') + '%s.npy' % path)
            i8states = np.sum(patterns[:, :, :, :, :2], axis=4)
            interactions = np.sum(i8states, axis=3)
            IPP[key].append(interactions)
    return IPP


def FollowingPerPair(names):
    FPP = {}
    for key in names.keys():
        FPP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..', 'PreprocessedData/InteractionsData/') + '%s.npy' % path)
            i8states = np.sum(patterns[:, :, :, :, :1], axis=4)
            interactions = np.sum(i8states, axis=3)
            FPP[key].append(interactions)
    return FPP


def AvoidingPerPair(names):
    APP = {}
    for key in names.keys():
        APP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..', 'PreprocessedData/InteractionsData/') + '%s.npy' % path)
            i8states = np.sum(patterns[:, :, :, :, 1:2], axis=4)
            interactions = np.sum(i8states, axis=3)
            APP[key].append(interactions)
    return APP


def factorial(n):
    if n < 2: return 1
    return reduce(lambda x, y: x * y, xrange(2, int(n) + 1))


def prob(s, p, n):
    x = 1.0 - p
    a = n - s
    b = s + 1
    c = a + b - 1
    prob = 0.0
    for j in xrange(a, c + 1):
        prob += factorial(c) / (factorial(j) * factorial(c - j)) \
                * x ** j * (1 - x) ** (c - j)
    return prob


def FAprobablity(a):
    pf, pa = prob(int(a[1]), 0.5, int(a[0] + a[1])), prob(int(a[0]), 0.5, int(a[0] + a[1]))
    if pf < pa and pf < 0.05:
        v = round(pf, 3)  # 0.5-pf
    elif pa < pf and pa < 0.05:
        print pa
        v = round(-pa, 3)  # pa-0.5
    else:
        v = 0.
    return v


def easyFAP(patterns, p=0.5):
    rawFA = np.sum(patterns[:, :, :, :, :2], axis=3, dtype='float32')
    FAPmatrix = np.apply_along_axis(FAprobablity, 3, rawFA)
    return FAPmatrix


def FollowingAvoidingMatrix(names):
    FAP = {}
    for key in names.keys():
        print key
        FAP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..', 'PreprocessedData/InteractionsData/') + '%s.npy' % path)
            FAP[key].append(easyFAP(patterns, p=0.5))
    return FAP


def longest_sequence(FAM, n_s=6):
    DARK = []
    ALL = []
    n_s, n_l, n_f = FAM.shape
    for i in range(n_l):
        for j in range(i, n_f):
            t_sec_ALL = 0
            sec_ALL = 0
            t_sec_DARK = 0
            sec_DARK = 0
            for s in range(n_s):
                if FAM[s, i, j] > 0 or FAM[s, j, i] > 0:
                    t_sec_ALL += 1
                else:
                    if sec_ALL < t_sec_ALL:
                        sec_ALL = t_sec_ALL
                    t_sec_ALL = 0
                if s % 2 == 0:
                    if FAM[s, i, j] > 0 or FAM[s, j, i] > 0:
                        t_sec_DARK += 1
                    else:
                        if sec_DARK < t_sec_DARK:
                            sec_DARK = t_sec_DARK
                        t_sec_DARK = 0
            ALL.append(sec_ALL)
            DARK.append(sec_DARK)
    return DARK, ALL


def follsactive(FAM, n_s=6):
    DARK = []
    ALL = []
    n_s, n_l, n_f = FAM.shape
    for i in range(n_l):
        for j in range(i, n_f):
            sec_ALL = 0.5
            sec_DARK = 0.5
            for s in range(n_s):
                if FAM[s, i, j] > 0 or FAM[s, j, i] > 0:
                    sec_ALL += 1
                if s % 2 == 0:
                    if FAM[s, i, j] > 0 or FAM[s, j, i] > 0:
                        sec_DARK += 1
            ALL.append(sec_ALL)
            DARK.append(sec_DARK)
    return DARK, ALL


if __name__ == "__main__":
    ts = 3
    experiments = load_experiments_info("experiments_desc.csv")
    comparisons = load_comparisons_info("comparisons.csv")
    # print comparisons.keys()
    ##for i in comparisons:
    ##    names, colors = group_data(i,comparisons,experiments, color_lst = ["red","green", "blue"])
    #    print nameswin
    #    preprocessData(names,window = 12,ts=3)
    names, colors = group_data('KO-WT|mouse|FX|females|1', comparisons, experiments, color_lst=["red", "green", "blue"])
    # names, colors = group_data('VPA-CTRL-NaCl|mouse|C57|males|1',comparisons,experiments, color_lst = ["red","green", "blue"])

    # createRandomExpepiments(exp_paths)
    # Interpersec(names,ts=200)
    startingpoint = datetime.now()
    #preprocessData(names, window=6, ts=3)
    print 'processing took '+ str(-(startingpoint-datetime.now()).total_seconds()) + ' seconds'
    IPP = InteractionsPerPair(names)
    scalefactor = np.max([np.max(i) for i in IPP["KO"]] + [np.max(i) for i in IPP["WT"]])
    # scalefactor = np.max([np.max(i) for i in IPP["BALB_VPA"]]+[np.max(i) for i in IPP["BALB_CTRL"]]+[np.max(i) for i in IPP["BALB_NaCl"]])
    FPP = FollowingPerPair(names)
    APP = AvoidingPerPair(names)
    FAM = FollowingAvoidingMatrix(names)
    result = []
    pair_inc = []
    pair_long = []

    stats = {}
    stats["KO"] = {}
    stats["KO"]["SLD"] = []
    stats["KO"]["NFD"] = []
    stats["KO"]["SLA"] = []
    stats["KO"]["NFA"] = []
    """
    for i in range(3):
        _FAM = FAM["KO"][i]
        LSD, LSA = longest_sequence(_FAM, n_s=6)
        FSD, FSA = follsactive(_FAM, n_s=6)
        stats["KO"]["SLA"] += LSA
        stats["KO"]["NFA"] += FSA
        stats["KO"]["NFD"] += FSD
        stats["KO"]["SLD"] += LSD
        print np.mean(LSA), np.mean(FSA), FAM["KO"][i].shape[2]
        plt.hist(LSD)
        plt.hist(FSD)
        plt.show()
    """
    print '#####################WT#######################'
    stats["WT"] = {}
    stats["WT"]["SLD"] = []
    stats["WT"]["NFD"] = []
    stats["WT"]["SLA"] = []
    stats["WT"]["NFA"] = []
    for i in range(4):
        _FAM = FAM["WT"][i]
        LSD, LSA = longest_sequence(_FAM, n_s=6)
        FSD, FSA = follsactive(_FAM, n_s=6)
        stats["WT"]["SLA"] += LSA
        stats["WT"]["NFA"] += FSA
        stats["WT"]["NFD"] += FSD
        stats["WT"]["SLD"] += LSD
        print np.mean(LSA), np.mean(FSA), FAM["WT"][i].shape[2]
    # plt.hist(LSA)
    #    plt.hist(FSA)
    #    plt.show()


    plotfunctions.createRasterPlots(FAM,IPP,names,scalefactor)
    plotfunctions.createRasterPlotsSUM(FAM,IPP,names,scalefactor)
    #plotfunctions.CreateRelationGraphs(FAM,IPP,names,scalefactor/50)

    """
    plotfunctions.barplot(stats,names,["SLA","SLD",], colors, name="AverageLengthofSec",ylab = "Average length of sequence per pair")
    plotfunctions.barplot(stats,names,["NFA","NFD" ], colors, name="AverageNumberofFols", ylab="Average number of followings per pair")
    #createRasterPlots(FAM,IPP,names,scalefactor)           
    #createRasterPlotsSUM(FAM,IPP,names,scalefactor)
    #CreateRelationGraphs(FAM,IPP,names,scalefactor/50)
    statsIPP = plotfunctions.plotphist(IPP,names,colors,to_file = True,directory = 'Interactions',vrange = [0,120], prange = [0,0.11])
    statsFPP = plotfunctions.plotphist(FPP,names,colors,to_file = True,directory = 'Followings',vrange = [0,120], prange = [0,0.11])
    statsAPP = plotfunctions.plotphist(APP,names,colors,to_file = True,directory = 'Avoidings',vrange = [0,120], prange = [0,0.11])
    print statsIPP
    plotfunctions.barplot(statsIPP,names,["Interactions"], colors, name="InteractionPerPairplotfunctions.barplot",ylab="Average number of interactions per pair")
    plotfunctions.barplot(statsFPP,names,["Followings"], colors,name="FollowingsPerPairplotfunctions.barplot",ylab="Average number of followings per pair")
    plotfunctions.barplot(statsAPP,names,["Avoidings"], colors,name="AvoidingsPerPairplotfunctions.barplot", ylab="Average number of avoidings per pair")
    """
    #print  st.mannwhitneyu(statsIPP["KO"]["mean"], statsIPP["WT"]["mean"], use_continuity=True)
    #print  st.mannwhitneyu(statsFPP["KO"]["mean"], statsFPP["WT"]["mean"], use_continuity=True)
    #print  st.mannwhitneyu(statsAPP["KO"]["mean"], statsAPP["WT"]["mean"], use_continuity=True)
    #print  st.mannwhitneyu(stats["KO"]["median"], stats["WT"]["median"], use_continuity=True)
    #FAP = FollowingAvoidingMatrix(names)
    #stats = plotfunctions.plotphist(FAP,names,colors,to_file = True,directory = 'FA',vrange = [-1,1], prange = [0,1])
    #print  st.mannwhitneyu(stats["KO"]["mean"], stats["WT"]["mean"], use_continuity=True)
    #print  st.mannwhitneyu(stats["KO"]["median"], stats["WT"]["median"], use_continuity=True)
    #m = np.round(FAP["KO"][0][0],3)
    #print m


