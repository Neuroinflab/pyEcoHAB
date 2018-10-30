#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:48:42 2016

@author: Jan Maka, Asia Jędrzejewska-Szmek
"""
#EcoHAB libraries
from __future__ import division,print_function

import EcoHab
import sys
if sys.version_info < (3, 0):
    import ConfigParser
else:
    import configparser as ConfigParser
from ExperimentConfigFile import ExperimentConfigFile
from experiments_info import smells, antenna_positions
from data_managment import *
#Third party libraries
import networkx as nx 
import numpy as np
import scipy
import os
#Debug libraries
import matplotlib.pyplot as plt   
import matplotlib.patches as patches
import plotfunctions
import write_to_file as wtf
import utils
from numba import jit

class Experiment(object):
    """Class, which represent data from one experiment"""
    def _fix_config(self):
        tstarts = []
        tends = []
        phases = self.cf.sections()
        for phase in phases:
            if 'dark' in phase or 'light' in phase or 'DARK' in phase or 'LIGHT' in phase:
                st, en = self.cf.gettime(phase)
                tstarts.append(st)
                tends.append(en)
        t_min = np.argmin(tstarts)
        t_max = np.argmax(tends)
        if 'ALL' not in phases :
            self.cf.add_section('ALL')
            
        startdate = self.cf.get(phases[t_min],'startdate')
        starttime = self.cf.get(phases[t_min],'starttime')
        enddate = self.cf.get(phases[t_max],'enddate')
        endtime = self.cf.get(phases[t_max],'endtime')
        
        self.cf.set('ALL','startdate',startdate)
        self.cf.set('ALL','starttime',starttime)
        self.cf.set('ALL','enddate',enddate)
        self.cf.set('ALL','endtime',endtime)
        
    def _merge_config(self,other_experiment):
        for phase in other_experiment.cf.sections():
            if phase != 'ALL':
                self.cf.add_section(phase)
                startdate = other_experiment.cf.get(phase,'startdate')
                starttime = other_experiment.cf.get(phase,'starttime')
                enddate = other_experiment.cf.get(phase,'enddate')
                endtime = other_experiment.cf.get(phase,'endtime')
                
                self.cf.set(phase,'startdate',startdate)
                self.cf.set(phase,'starttime',starttime)
                self.cf.set(phase,'enddate',enddate)
                self.cf.set(phase,'endtime',endtime)
        self._fix_config()
            
    def _remove_phases(self,mask):
        if mask:
            for phase in self.cf.sections():
                st, en = self.cf.gettime(phase)
                if en <= mask[0]:
                    self.cf.remove_section(phase)
                elif st >= mask[1]:
                    self.cf.remove_section(phase)
 
            self._fix_config()
    
    def __init__(self, path,**kwargs):
        self.path = path
        _ant_pos = kwargs.pop('_ant_pos', None)
        which_phase = kwargs.pop('which_phase', 'ALL')
        mask = kwargs.pop('mask', None)
        from_file = kwargs.pop('from_file', False)
        h_m_a = kwargs.pop('how_many_appearances', 1000)
        factor = kwargs.pop('factor', 2)
        tags = kwargs.pop('remove_mice', [])
        self.compensate_for_lost_antenna = kwargs.pop('compensate_for_lost_antenna', False)
        if kwargs:
            raise Exception('Too many arguments for class Experiment')
        self.directory = utils.results_path(path)
        self.from_file = from_file
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.cf = ExperimentConfigFile(self.path)
        if which_phase and mask == None:
            try:
                mask = self.cf.gettime(which_phase)
            except ConfigParser.NoSectionError:
                mask = None
        self.ehs = EcoHab.EcoHabSessions9states(path=self.path,
                                                _ant_pos=_ant_pos,
                                                mask=mask,
                                                shortest_session_threshold=0,
                                                how_many_appearances=h_m_a,
                                                factor=factor,
                                                remove_mice=tags,
                                                compensate_for_lost_antenna=self.compensate_for_lost_antenna)
        self._remove_phases(mask)
        self.fs = self.ehs.fs
        mice = list(self.ehs.mice)
        self.mice = [mouse for mouse in mice if len(self.ehs.getstarttimes(mouse)) > 30]
        self.lm = len(self.mice)
        self.t_start = self.ehs.data['AbsStartTimecode'][0]
        self.t_end = self.ehs.data['AbsStartTimecode'][-1]
        
        self.fname_ending = which_phase
        self.phases = None
        self.key_list = ['Address',
                         'Tag',
                         'AbsStartTimecode',
                         'AbsEndTimecode',
                         'VisitDuration',
                         'ValidVisitSolution']
        self.prefix = utils.make_prefix(self.path)
        
    def add_mouse_data(self, other_experiment):
        other_sd = other_experiment.ehs.signal_data
        for mouse in other_experiment.mice:
            if mouse in self.ehs.signal_data:
                np.concatenate((self.ehs.signal_data[mouse],
                                other_sd[mouse]))
            else:
                self.ehs.signal_data[mouse] = other_sd[mouse]
            
        np.concatenate((self.ehs.signal_data['time'],
                        other_sd['time']))
        
        for key in self.key_list:
            self.ehs.data[key].extend(other_experiment.ehs.data[key])
                                                   
    def merge_experiment(self, other_experiment, same_mice=True):
        """This is only supposed to be used with one experiment (or an
        experiment with the same mice), when e.g. the antennae setup has
        changed.
        
        """
        other_cf = other_experiment.cf
        assert self.ehs.fs == other_experiment.ehs.fs
        if same_mice:
            assert self.mice==other_experiment.mice
        else:
            self.mice.append(other_experiment.mice)
        self.add_mouse_data(other_experiment)
        self.ehs.mice = self.mice
        self.lm = len(self.mice)
        self.ehs.t_start_exp = np.min(self.ehs.data['AbsStartTimecode'])
        self.ehs.t_end_exp = np.max(self.ehs.data['AbsStartTimecode'])    
        self._merge_config(other_experiment)
        self.t_start = min(self.ehs.data['AbsStartTimecode'])
        self.t_end = max(self.ehs.data['AbsStartTimecode'])
        add = ""#"_"+other_experiment.fname_ending+"_merged_"
        self.fname_ending = self.fname_ending + add

    def calculate_phases(self, window='default', which_phase='ALL'):
        self.tstart, self.tend = self.cf.gettime('ALL')
        sd = self.ehs.signal_data
        if window == 'default':
            sessions = filter(lambda x: x.endswith('dark')
                              or x.endswith('light'),
                              self.cf.sections())
            self.phases = [(self.cf.gettime(sec)[0]-self.tstart,
                            self.cf.gettime(sec)[1]-self.tstart)
                           for sec in sessions]
        elif window == "ALL":
            self.phases = [[0, int(np.ceil(self.tend-self.tstart))]]
        else:
            if isinstance(window,float) or isinstance(window,int):
                phase_nb = int(np.ceil((self.t_end - self.t_start)/(window*3600)))
                self.phases = [(i*window*3600,
                                np.min([(i+1)*window*3600,
                                        len(sd[self.mice[0]])]))
                               for i in range(phase_nb)]
            elif isinstance(window, list):
                self.phases = [(st*window[0]*3600,
                                (st+1)*window[0]*3600) for st in window[1]]
            else:
                raise TypeError

    def calculate_antenna_errors_phase(self, t1, t2):
        outs = np.zeros(len(self.mice))
        for i, mouse in enumerate(self.mice):
            outs[i] = self.calculate_antenna_errors_mouse_phase(mouse, t1, t2)
        return outs
    
    def calculate_antenna_errors(self):    
        self.calculate_phases(window='default', which_phase='ALL')
        print(self.phases)
        errors = np.zeros((len(self.phases), len(self.mice)))
        for s, (ts, te) in enumerate(self.phases):
            errors[s] = self.calculate_antenna_errors_phase(ts, te)
            
        self.write_to_file(errors, 'Antennae_errors.csv', subdirectory="antennae_errors")
        return errors
    def calculate_antenna_errors_mouse_phase(self, mouse, t1, t2):
        start = int(self.fs*t1)
        end =  int(self.fs*t2)
        mouse_sig = self.ehs.signal_data[mouse][start:end+1]
        if mouse_sig[0]:
            errors = 0
        else:
            errors = 1
        for i, s in enumerate(mouse_sig[1:]):
            if s == 0:
                errors += 1
            elif abs(s - mouse_sig[i]) == 1:
                continue
            elif abs(s - mouse_sig[i]) == 7:
                continue
            elif abs(s - mouse_sig[i]) == 0:
                continue            
            else:
                errors += 1
        return errors/(self.fs*(t2-t1))

    def calculate_fvalue(self,
                         window='default',
                         threshold=2,
                         force=False,
                         fols=None,
                         ops=None,
                         which_phase='ALL'):
        self.calculate_phases(window=window, which_phase=which_phase)
        self.fols = fols
        self.ops = ops
        self.threshold = threshold
        self.fpatterns = []
        self.opatterns = []
        l = len(self.phases)
        print(self.phases)
        self.f = np.zeros((self.lm, self.lm, l))
        self.interactions = np.zeros((l, self.lm, self.lm, 8, 3))
        self.f_sum = np.zeros((l,self.lm))
        new_path = os.path.join(self.directory,
                                'PreprocessedData/IteractionsData/')
        new_fname_patterns = os.path.join(new_path,
                                          'Patterns_cut_%s_%s.npy' % (self.fname_ending,
                                                                      which_phase))
        new_fname_fpatterns = os.path.join(new_path,
                                           'fpatterns_cut_%s_%s.npy' % (self.fname_ending,
                                                                        which_phase))
        new_fname_opatterns = os.path.join(new_path,
                                           'opatterns_cut_%s_%s.npy' % (self.fname_ending,
                                                                        which_phase))
        new_fname_mice = os.path.join(new_path,
                                      'mice_cut_%s_%s.npy' % (self.fname_ending,
                                                              which_phase))
        
        if not os.path.exists(os.path.dirname(new_path)):
            os.makedirs(os.path.dirname(new_path))
        if self.from_file:
            self.interactions = np.load(new_fname_patterns)
            self.fpatterns = np.load(new_fname_fpatterns)
            self.opatterns = np.load(new_fname_opatterns)
            self.mice_list =  np.load(new_fname_mice)
        else:
            for s in range(len(self.phases)):
                ts, te = self.phases[s]
                print('Phase %s. from %sh, to %sh'%(s + 1,
                                                    np.round(ts/3600., 2),
                                                    np.round(te/3600., 2)))
                self.interactions[s, :, :, :, :] = self.interaction_matrix(ts,
                                                                           te)
                np.save(new_fname_patterns, self.interactions)
                np.save(new_fname_fpatterns, self.fpatterns)
                np.save(new_fname_opatterns, self.opatterns)
                np.save(new_fname_mice, self.mice)
        return self.f

    def tube_dominance_test(self,
                            window='default',
                            which_phase='ALL'):
        self.calculate_phases(window=window,
                              which_phase=
                              which_phase)
        self.tube_dominance_matrix = np.zeros((len(self.phases),
                                               self.lm,
                                               self.lm))
        new_path = utils.check_directory(self.directory,
                                         'PreprocessedData/IteractionsData/')
        new_fname_ = os.path.join(new_path,
                                  'Patterns_%s.npy' % self.fname_ending)
        new_fname_mice = os.path.join(new_path,
                                      'mice_cut_%s_%s.npy' % (self.fname_ending,
                                                              which_phase))
        for s in range(len(self.phases)):
            ts, te = self.phases[s]
            print('Phase %s. from %sh, to %sh'%(s + 1,
                                                np.round(ts/3600., 2),
                                                np.round(te/3600., 2)))
            out = self.calculate_tube_dominance_matrix(ts, te, normalize=False)
            self.tube_dominance_matrix[s, :, :] = out
            
        np.save(new_fname_mice, self.mice)
        np.save(new_fname_, self.tube_dominance_matrix)

    def calculate_tube_dominance_matrix(self,ts, te, normalize):
        """
        Calculates tube dominance matrix for the given period from ts to te
        check_tube_dominance(mouse1,mouse2,t1,t2) checks if mouse2
        dominates over mouse 1. That is why mice are switched in the
        loop.

        """
        imatrix = np.zeros((self.lm, self.lm))
        for ii, mouse1 in enumerate(self.mice):
            for jj,mouse2 in enumerate(self.mice):
                if ii < jj:
                    imatrix[ii, jj] = self.check_tube_dominance(mouse2,
                                                                mouse1,
                                                                ts,
                                                                te)
                    imatrix[jj, ii] = self.check_tube_dominance(mouse1,
                                                                mouse2,
                                                                ts,
                                                                te)
                    total = imatrix[ii, jj] + imatrix[jj, ii]
                    if normalize:
                        if total:
                            imatrix[ii, jj] = imatrix[ii, jj]/total
                            imatrix[jj, ii] = imatrix[jj, ii]/total
        return imatrix
        
    def check_tube_dominance(self, mouse1, mouse2, t1, t2):
        """
        This is an implementation of the tube dominance test,
        where two mice are placed on opposite sides of a tube
        and one of the mice forces its opponent out of the tube.
        The test checks if mouse1 entered a tube and backed
        out and if mouse2 remains in the tube and entered the tube
        from the other side. Effectively this function checks
        how many times mouse2 domineered over mouse1.
        """
        dominance_stats = 0
        sd = self.ehs.signal_data
        try:
            [mouse1_idx_start,
             mouse1_idx_stop,
             mouse1_indices] = self.mouse_indices(mouse1, t1, t2)
        except (ValueError, TypeError):
            return dominance_stats
        i = mouse1_idx_start - 1
        for end in mouse1_indices[mouse1_idx_start:mouse1_idx_stop+1]:
            i = i + 1
            if i < 2:
                continue
            start = mouse1_indices[i-2]
            start_state = int(sd[mouse1][start]) #from start to mouse1_indices[i-1]-1
            middle_state = int(sd[mouse1][mouse1_indices[i-1]])#from mouse1_indices[i-1] to end-1
            end_state = int(sd[mouse1][end])#from end
            if start_state != end_state:
                continue
            if not middle_state %2:
                continue
            if not start_state*middle_state*end_state:
                continue
            try:
                end_state_mouse_2 = sd[mouse2][end-1]
                post_end_state_mouse_2 = sd[mouse2][end]
            except IndexError:
                print("Index error")
                return dominance_stats
            if end_state_mouse_2 != middle_state or post_end_state_mouse_2 != middle_state:
                continue
            #find previous state
            idx = end - 1
            while True:
                idx += -1
                if end_state_mouse_2 != sd[mouse2][idx]:
                    t_middle_end = idx
                    middle_state_mouse2 = sd[mouse2][idx]
                    break
                if idx < self.fs*t1:
                    t_middle_end = 0
                    middle_state_mouse2 = 0
                    break

            if not middle_state_mouse2 or middle_state_mouse2 == start_state:
                continue
            dominance_stats +=1
        return dominance_stats

    def interaction_matrix(self, ts, te):
        """Calculates fvalue matrix for the given period from ts to te
        """
        imatrix = np.zeros((self.lm,self.lm,8,3))
        for ii, mouse1 in enumerate(self.mice):
            for jj,mouse2 in enumerate(self.mice):
                if ii < jj:
                    imatrix[ii,jj,:,:],patterns = self.findpatterns(ii,jj,ts, te)
                    self.fpatterns+=patterns[0]
                    self.opatterns+=patterns[1]
                    imatrix[jj,ii,:,:],patterns = self.findpatterns(jj,ii,ts, te)
                    self.fpatterns+=patterns[0]
                    self.opatterns+=patterns[1]
        return imatrix
    
    @jit
    def convolution(self, ts, te, tau=60): # tau equal to 60 s
        sd = self.ehs.signal_data
        shift = tau*self.fs
        new_ts, new_te = ts*self.fs, te*self.fs
        length = new_te - new_ts
        
        fig, ax = plt.subplots(self.lm, self.lm)
        print(ax.shape)
        time = np.linspace(-tau, tau, 2*shift+1)
        for ii, mouse1 in enumerate(self.mice):
            for jj, mouse2 in enumerate(self.mice):
                if ii != jj:
                    s2 = sd[mouse2][new_ts:new_te]
                    s1 = np.zeros(2*shift+length)
                    if ts < tau:
                        if te - self.t_end > tau:
                            s1[shift:length+2*shift] = sd[mouse1][new_ts:new_te+shift]
                        else:
                            s1[shift:length+shift] = sd[mouse1][new_ts:new_te]
                    else:
                        if te - self.t_end > tau:
                            s1 = sd[mouse1][new_ts-shift:new_te+shift]
                        else:
                            s1[0:length+shift] = sd[mouse1][new_ts-shift:new_te]
                            
                    out = np.zeros(2*shift+1)
                    m1 = sd[mouse1][new_ts:new_te].mean()
                    std1 = sd[mouse1][new_ts:new_te].var()**.5
                    m2 = sd[mouse2][new_ts:new_te].mean()
                    std2 = sd[mouse2][new_ts:new_te].var()**.5
                    
                    for i, new_tau in enumerate(range(-shift, shift+1)):
                        out[i] = sum((s1[shift-new_tau:length+shift-new_tau]-m1)*(s2-m2))/(std1*std2)/length
                    print(mouse1, mouse2, out.max(), time[out.argmax()])
                    ax[ii, jj].plot(time, out)
                    ax[ii, jj].set_ylim([-.5, .5])
                    if not ii:
                        ax[ii, jj].set_title(mouse2)
                    if not jj:
                        ax[ii, jj].set_title(mouse1)
                        

    plt.show()
    def validatePatterns(self,plots_nr = 9, trange = [-3,3]):
        sd = self.ehs.signal_data
        size = np.ceil(np.sqrt(plots_nr))
        t = np.arange(trange[0],trange[1],1.0/self.fs)
        frandom_idx = [np.random.randint(0,len(self.fpatterns)-1) for i in range(plots_nr)]
        orandom_idx = [np.random.randint(0,len(self.opatterns)-1) for i in range(plots_nr)]
        plt.figure()
        plt.suptitle("Random following patterns",
                     fontsize=14,
                     fontweight='bold')
        for i,idx in enumerate(frandom_idx):
            ax = plt.subplot(size, size, i+1)

            ii,jj,s = self.fpatterns[idx]
            ax.set_title("%s|%s|t=%s"%(ii,jj,s*1./self.fs))
            mouse1 = self.mice[ii]
            mouse2 = self.mice[jj]
            plt.plot(t,sd[mouse1][s-3*self.fs:s+3*self.fs]-0.05,'ro',label="leader")
            plt.plot(t,sd[mouse2][s-3*self.fs:s+3*self.fs]+0.05,'bo',label="follower")
            plt.axis([-3.1,3.1,-0.5,9.5])
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()
        plt.figure()
        plt.suptitle("Random avoiding patterns", fontsize=14, fontweight='bold')
        for i,idx in enumerate(orandom_idx):
            ax = plt.subplot(size, size,i+1)

            ii,jj,s = self.opatterns[idx]
            mouse1 = self.mice[ii]
            mouse2 = self.mice[jj]
            ax.set_title("%s|%s|t=%s"%(ii,jj,s*1./self.fs))
            plt.plot(t,sd[mouse1][s-3*self.fs:s+3*self.fs]-0.05,'ro',label="leader")
            plt.plot(t,sd[mouse2][s-3*self.fs:s+3*self.fs]+0.05,'bo',label="follower")
            plt.axis([-3.1,3.1,-0.5,9.5])
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

    @jit    
    def calculate_mouse_stats(self,mouse):
        
        follow_stat = {}
        mouse_stats = self.ehs.statistics[mouse]["preference"]
        moves = ['24','42','46','64','68','86','82','28']
        for move in moves:
            follow_stat[move] = np.zeros(3)
            start_st, end_st = int(move[0]),int(move[1])
            if end_st == 2 and start_st==8:
                index = 0
            elif end_st == 8 and start_st==2:
                index = 1
            elif start_st <end_st:
                index = 0
            else:
                index = 1
            if np.sum(mouse_stats[int(move[0])]) > 0:
                follow_stat[move][2] = mouse_stats[int(move[0])][index]/np.sum(mouse_stats[int(move[0])])
            else:
                follow_stat[move][2] = 0
        return follow_stat

    def mouse_indices(self, mouse1, t1, t2):
        sd = self.ehs.signal_data
        mouse1_indices = np.where(((np.roll(sd[mouse1], 1) - sd[mouse1]) != 0))[0]
        try:
            mouse1_idx_start = np.where(mouse1_indices > t1*self.fs)[0][0]
        except IndexError:
            return None
        
        try:    
            mouse1_idx_stop =  np.where(mouse1_indices < t2*self.fs)[0][-1]
        except IndexError:
            return None
        return [mouse1_idx_start, mouse1_idx_stop, mouse1_indices]

    @jit
    def findpatterns(self,ii,jj,t1, t2):
    
        sd = self.ehs.signal_data
        detected_idx = [[],[]]
        mouse1 = self.mice[ii]
        mouse2 = self.mice[jj]
    
        follow_stat = self.calculate_mouse_stats(mouse2)
        #print(follow_stat)
        out = self.mouse_indices(mouse1, t1,t2)
        if out is None:
            return np.array(follow_stat.values()), detected_idx
        else:
            [mouse1_idx_start,mouse1_idx_stop,mouse1_indices] = self.mouse_indices(mouse1, t1,t2)
            
        i = mouse1_idx_start - 1
       
        for start in mouse1_indices[mouse1_idx_start:mouse1_idx_stop+1]:
            
            i = i + 1
           
            if i < 2:
                continue

            start_state = int(sd[mouse1][mouse1_indices[i-2]])
            middle_state = int(sd[mouse1][mouse1_indices[i-1]])
            end_state = int(sd[mouse1][start])
           
            end = start + self.threshold*self.fs

            unknown_state = end_state == 0 or middle_state == 0 or start_state == 0
            
            in_pipe = end_state%2==1
            # define conditions
            if unknown_state or start_state == end_state or in_pipe:
                continue
            
            period1 = list(sd[mouse2][start:end])
            
            period2 = list(sd[mouse2][start-2*self.fs:start])

            period3 = list(sd[mouse2][start-int(0.1*self.fs):start])
              
            op_idx = (2*start_state-end_state-1)%8+1 #opposite direction
                
            same_start = start_state in period2 #POPRAWIC!!!!!!!

            first_mouse1 = end_state not in period3 and op_idx not in period3
                
            followed = period1.count(end_state) > 0
                
            go_oposite = op_idx in period1
                
            if followed and go_oposite:
                followed = period1.index(end_state) < period1.index(op_idx)
                go_oposite = not followed
                    
                    
            if same_start and first_mouse1:

                if followed:

                    if self.fols!=None:
                        self.fols[np.ceil((period1.index(end_state))/self.fs)]+=1
                    if str(start_state)+str(end_state) in follow_stat:
                        follow_stat[str(start_state)+str(end_state)][0] += 1
                       
                    detected_idx[0].append((ii,jj,start))

                elif go_oposite:

                    if self.ops!=None:
                        self.ops[np.ceil((period1.index(op_idx))/self.fs)]+=1

                    if str(start_state)+str(end_state) in follow_stat:
                        follow_stat[str(start_state)+str(end_state)][1] += 1
                    detected_idx[1].append((ii,jj,start))
                                                
      
        return np.array(list(follow_stat.values())), detected_idx
    
    def InteractionsPerPair(self,start,end):

        i8states = np.sum(self.interactions[:,:,:,:,start:end],axis=4) # po kierunku -- 0 following, 1 -- avoiding
        return np.sum(i8states ,axis=3) #po wszystkich stanach (24,42, etc.)

    def FollowingAvoidingMatrix(self):
        
        patterns = self.interactions
        rawFA = np.sum(patterns[:,:,:,:,:2] ,axis=3,dtype = 'float32')
        FAPmatrix = np.apply_along_axis(FAprobablity, 3, rawFA)
        return FAPmatrix

    
    def normalize(self,matrix):
        nx,ny = matrix.shape
        for i in range(nx):
            for j in range(i+1,ny):
                total = matrix[i,j]+matrix[j,i]
                if total:
                    matrix[i,j] /= total
                    matrix[j,i] /= total
        return matrix

    def write_to_file(self, output, fname, subdirectory=""):
        if self.compensate_for_lost_antenna:
            fname = '%s_comp_for_lost_antenna' % fname
        if subdirectory:
            new_path = utils.check_directory(self.directory, subdirectory)
        else:
            new_path = self.directory
        
        new_fname = os.path.join(new_path, fname)
        
        f = open(new_fname, 'w')
        
        header = 'mouse'
        for phase in self.cf.sections():
            if 'dark' in phase:
                header += ',' + phase
            elif 'DARK' in phase:
                header += ',' + phase
            elif 'light' in phase:
                header += ',' + phase
            elif 'LIGHT' in phase:
                header += ',' + phase
        header = header + '\n'
        f.write(header)
        for i, mouse in enumerate(self.mice):
            f.write(mouse)
            for j in range(len(self.phases)):
                f.write(',')
                f.write(str(output[j, i]))
            f.write('\n')
        f.close()

    def get_interaction(self, what):
        if what == "Following" or what == "following":
            what = "following"
            return self.InteractionsPerPair(0,1), what
        if what == "Avoiding" or what == "avoiding":
            what = "avoiding"
            return self.InteractionsPerPair(1,2), what
        if what == "FAM":
            what = "interactions"
            return self.FollowingAvoidingMatrix(), what
       
        print("Unknown value %s, available values are following or avoiding" % what)
        return
    
    def get_phases(self, phases):
        if phases:
            if isinstance(phases, list):
                phases_comment = 'from_%s_to_%s' % (phases[0].replace(' ',''), phases[-1].replace(' ',''))
                return phases, phases_comment
            if isinstance(phases, str):
                phases_comment = phases
                phases = [phases]
                return phases, phases_comment
            
            print("Unknown phases type, leaving")
            return
    
        phases_comment = "all_phases"
        phases = []
        for phase in self.cf.sections():
            if 'dark' in phase or 'DARK' in phase:
                phases.append(phase)
            elif 'light' in phase or 'LIGHT' in phase:
                    phases.append(phase)

        return phases, phases_comment

    def write_tables_to_file(self, what, phases=None, write_all=False):
        try:
            result, what = self.get_interaction(what)
        except ValueError:
            return
        try:
            phases, phases_comment = self.get_phases(phases)
        except ValueError:
            return
        print('Write %s to'%what)

        fname = '%s_%s' % (what, phases_comment)
        if self.compensate_for_lost_antenna:
            fname = '%s_comp_for_lost_antenna' %fname
        wtf.write_csv_tables(result, phases, self.mice, self.directory, what, fname, self.prefix)
       
        if write_all:
            new_directory = utils.check_directory(self.directory, 'single')
            for i, small_matrix in enumerate(result):
                print(i, small_matrix)
                new_fname = '%s_%s' % (what, phases[i])
                wtf.write_csv_tables([result[i]],
                                     [phases[i]],
                                     self.mice,
                                     new_directory,
                                     what,
                                     new_fname,
                                     self.prefix)
                
                
        
    def plot_fam(self, phases=None, scalefactor=0):
        FAM = self.FollowingAvoidingMatrix()
        IPP = self.InteractionsPerPair(0,2)

        try:
            phases, phases_comment = self.get_phases(phases)
        except ValueError:
            return
        suffix = '%s_%s' % (self.prefix, phases_comment)
        
        plotfunctions.raster_interactions(self.directory,
                                          FAM,
                                          IPP,
                                          phases,
                                          suffix,
                                          self.mice,
                                          scalefactor=scalefactor)
        
        self.write_tables_to_file("FAM", phases=phases)
        mice = [mouse.split('-')[-1] for mouse in self.mice]
        for i, fam in enumerate(FAM):
            plotfunctions.plot_graph(fam, phases[i], self.directory, labels=mice)
            
    def generate_heatmaps(self, what, phases=None):
        
        try:
            result, what = self.get_interaction(what)
        except ValueError:
            return
        try:
            phases, phases_comment = self.get_phases(phases)
        except ValueError:
            return
        vmin, vmax = result.min(), result.max()
        
        for i, phase in enumerate(phases):
            suffix = '%s_%s' % (self.prefix, phase.replace(' ', '_'))
            plotfunctions.single_heat_map(result[i],
                                          what,
                                          self.directory,
                                          self.mice,
                                          suffix,
                                          subdirectory=what,
                                          vmin=vmin,
                                          vmax=vmax)
            
    def plotTubeDominanceRasters(self, name=None, mice=[]):
        subdirectory = 'tube_dominance/figs'
        n_s = self.tube_dominance_matrix.shape[0]
        phases = self.cf.sections()[:n_s]
        plotfunctions.make_RasterPlot(self.directory,
                                      subdirectory,
                                      self.tube_dominance_matrix,
                                      phases,
                                      'tube_dominance',
                                      self.mice,
                                      self.prefix,
                                      to_file=True,
                                      vmin=None,
                                      vmax=None,
                                      title=None)
        wtf.write_csv_tables(self.tube_dominance_matrix,
                             phases,
                             self.mice,
                             self.directory,
                             'tube_dominance',
                             'tube_dominance',
                             self.prefix)
        
    def plot_all(self,t1,t2):
        sd = self.ehs.signal_data
        fig = plt.figure()
        ax = []
        lm = self.lm
        if lm < 4:
            how_many = lm
            nrows = 1
            ncolumns = how_many
        else:
            if lm % 3:
                how_many = ((lm + 1) % 3 == 0)*(lm + 1) + ((lm + 2) % 3 == 0)*(lm+2)
            else:
                how_many = lm
                    
            nrows = how_many//3
            ncolumns = 3
            k = 0
    
            for i in range(nrows):
                for j in range(ncolumns):
                    if k < lm:
                        ax.append(fig.add_subplot(nrows,ncolumns,k+1))
                        mouse = self.mice[k]
                        sig = sd[mouse][t1*self.fs:t2*self.fs]
                        time = np.linspace(t1,t2,len(sig))
                        ax[k].plot(time,sig,'dk')
                        ax[k].set_title(mouse)
                        ax[k].set_xlabel('time [s]')
                        ax[k].set_ylabel('State')
                        k = k+1                        
    

def binomial_probability(s, p, n):
    prob = 0.0

    for j in range(s, n + 1):
        prob += scipy.special.binom(n, j) * p**j * (1 - p)**(n-j)
   
    return prob


def FAprobablity(a, p=0.5):
    
    number_of_following = int(a[0])
    number_of_avoiding = int(a[1])
    total = number_of_following + number_of_avoiding
    probability_of_following = 0.5
    probability_of_avoiding = 0.5
    
    pf = binomial_probability(number_of_following, probability_of_following, total)
    pa = binomial_probability(number_of_avoiding, probability_of_avoiding, total)

    if pf < pa and pf < 0.05:
        v = round(pf,6)#0.5-pf
    elif pa < pf and pa < 0.05:
        v = round(-pa,6)#pa-0.5
    else:
        v = 0.
    
    return v

def easyFAP(patterns):
    rawFA = np.sum(patterns[:,:,:,:,:2] ,axis=3,dtype = 'float32')
    FAPmatrix = np.apply_along_axis(FAprobablity, 3, rawFA)
    return FAPmatrix
    
def longest_sequence(FAM):
    DARK = []
    ALL = []
    n_s,n_l,n_f = FAM.shape
 
    for i in range(n_l):
        for j in range(i,n_f):
            t_sec_ALL = 0
            sec_ALL = 0
            t_sec_DARK = 0
            sec_DARK = 0
            for s in range(n_s):
                #print(s,i,j,FAM[s,i,j],FAM[s,j,i],t_sec_ALL,sec_ALL)
                if FAM[s,i,j]>0 or FAM[s,j,i] >0:
                    t_sec_ALL+=1
                else:
                    if sec_ALL<t_sec_ALL:
                        sec_ALL = t_sec_ALL
                    t_sec_ALL = 0
                if s%2==0:
                    if FAM[s,i,j]>0 or FAM[s,j,i] >0:
                        t_sec_DARK+=1
                    else:
                        if sec_DARK<t_sec_DARK:
                            sec_DARK = t_sec_DARK
                        t_sec_DARK = 0
            ALL.append(sec_ALL)
            DARK.append(sec_DARK)
    return DARK, ALL

def follsactive(FAM):
    DARK = []
    ALL = []
    n_s,n_l,n_f = FAM.shape
    for i in range(n_l):
        for j in range(i,n_f):
            sec_ALL = 0.5
            sec_DARK = 0.5
            for s in range(n_s):
                if FAM[s,i,j]>0 or FAM[s,j,i] >0:
                    sec_ALL+=1
                if s%2==0:
                    if FAM[s,i,j]>0 or FAM[s,j,i] >0:
                        sec_DARK+=1
            ALL.append(sec_ALL)
            DARK.append(sec_DARK)
    return DARK, ALL



if __name__ == "__main__":
    antenna_pos = {
        "/home/jszmek/EcoHAB_data_November/long_experiment_KO":{'1':1,'2':5,'3':3,'4':6,'5':4,'6':2,'7':7,'8':8},
        "/home/jszmek/EcoHAB_data_November/long_experiment_WT":{'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8}
}
    a_dirs  = [
        "/home/jszmek/EcoHAB_data_November/long_experiment_KO",
        #"/home/jszmek/EcoHAB_data_November/long_experiment_WT"
    ]

    phases = {"/home/jszmek/EcoHAB_data_November/long_experiment_KO":["BEGINNING","MIDDLE"],
         "/home/jszmek/EcoHAB_data_November/long_experiment_WT":["WRONG_ANTENNAS", "CORRECT_ANTENNAS"]}
    #compensate_for_lost_antenna=True
    ts = 3   
    window=12
    IPP = {}
    FAM = {}
    sections = {}
    directories = {}
    endings = {}
    mice = {}
    for a_dir in a_dirs:
        IPP[a_dir] = []
        FAM[a_dir] = []
        sections[a_dir] = []
        directories[a_dir] = []
        endings[a_dir] = []
        mice[a_dir] = []
        if a_dir == "/home/jszmek/EcoHAB_data_November/long_experiment_KO":
            E1 = Experiment(path=a_dir,_ant_pos=antenna_pos[a_dir],which_phase="WRONG_ANTENNAS", compensate_for_lost_antenna=True)
            E2 = Experiment(path=a_dir,_ant_pos=None, which_phase="CORRECT_ANTENNAS", compensate_for_lost_antenna=True)
            E1.merge_experiment(E2)
        else:
             E1 = Experiment(path=a_dir,_ant_pos=antenna_pos[a_dir], compensate_for_lost_antenna=False)
        window = 12
        E1.calculate_fvalue(window=window, threshold=ts)
        E1.write_tables_to_file("following")
        E1.write_tables_to_file("avoiding")
        E1.write_tables_to_file("FAM")
        E1.generate_heatmaps("following")
        E1.generate_heatmaps("avoiding")
        E1.plot_fam()
        

        window = "ALL"
        E1.calculate_fvalue(window=window, threshold=ts)
        
        E1.write_tables_to_file("following", phases="ALL")
        E1.write_tables_to_file("avoiding", phases="ALL")
        E1.write_tables_to_file("FAM", phases="ALL")
        E1.generate_heatmaps("following", phases="ALL")
        E1.generate_heatmaps("avoiding", phases="ALL")
        E1.plot_fam(phases="ALL")
                
