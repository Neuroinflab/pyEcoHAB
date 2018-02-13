#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:48:42 2016

@author: Jan Maka
"""
#EcoHAB libraries
from __future__ import division,print_function

import EcoHab
import ConfigParser
from ExperimentConfigFile import ExperimentConfigFile
from experiments_info import smells, antenna_positions
from plotfunctions import *
from data_managment import *
#Third party libraries
import networkx as nx 
import numpy as np
import scipy
import os
import pickle
#Debug libraries
import time
import matplotlib.pyplot as plt   
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcol
import matplotlib.patches as patches
import numpy.random

class Experiment(object):
    """Class, which represent data from one of the experiments"""
    def _fix_config(self):
        tstarts = []
        tends = []
        phases = self.cf.sections()
        for phase in phases:
            if 'dark' in phase or 'light' in phase:
                st, en = self.cf.gettime(phase)
                tstarts.append(st)
                tends.append(en)
                
        t_min = np.argmin(tstarts)
        t_max = np.argmax(tends)
        if 'ALL' in phases :
            
            startdate = self.cf.get(phases[t_min],'startdate')
            starttime = self.cf.get(phases[t_min],'starttime')
            enddate = self.cf.get(phases[t_max],'enddate')
            endtime = self.cf.get(phases[t_max],'endtime')
            
            self.cf.set('ALL','startdate',startdate)
            self.cf.set('ALL','starttime',starttime)
            self.cf.set('ALL','enddate',enddate)
            self.cf.set('ALL','endtime',endtime)
            
    def _remove_phases(self,mask):
        if mask:
            for phase in self.cf.sections():
                st, en = self.cf.gettime(phase)
                if en <= mask[0]:
                    self.cf.remove_section(phase)
                elif st >= mask[1]:
                    self.cf.remove_section(phase)
 
            self._fix_config()

    def __init__(self, path,**kwargs):#_ant_pos=None,which_phase='ALL',mask=None,from_file=False):
        
        self.path = path
        _ant_pos = kwargs.pop('_ant_pos',None)
        which_phase = kwargs.pop('which_phase','ALL')
        mask = kwargs.pop('mask',None)
        from_file = kwargs.pop('from_file',False)
        how_many_appearances = kwargs.pop('how_many_appearances',1000)
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])
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
        
        self.ehs = EcoHab.EcoHabSessions9states(path=self.path, _ant_pos=_ant_pos,mask=mask,shortest_session_threshold=0,how_many_appearances=how_many_appearances,factor=factor,remove_mice=tags)
        self._remove_phases(mask) 
        self.fs = self.ehs.fs

       
        mice = list(self.ehs.mice)
        
        self.mice = filter(lambda x: len(self.ehs.getstarttimes(x)) > 30, mice)
        
        self.sd = {}
        # for mouse in self.mice:
        #     self.sd[mouse] = np.array(self.ehs.signal_data[mouse],dtype=np.int)
        # self.sd['time'] = self.ehs.signal_data['time']
        self.sd = self.ehs.signal_data.copy()
        self.lm = len(self.mice)

        if mask:
            self.t_start = self.ehs.data['AbsStartTimecode'][0]
            self.t_end = self.ehs.data['AbsStartTimecode'][-1]
        else:
            self.t_start = 0
            self.t_end = len(self.ehs.signal_data[self.mice[0]])/self.fs

        self.fname_ending = which_phase
        
        self.phases = None
        self.remove_zeros()
        self.state_probabilities = {}
        
        
    def calculate_phases(self,window='default',which_phase='ALL'):
        
        self.tstart, self.tend = self.cf.gettime('ALL')
        sd = self.ehs.signal_data
        if window=='default':
            sessions = filter(lambda x: x.endswith('dark') or x.endswith('light'), self.cf.sections())
            self.phases = [(self.cf.gettime(sec)[0]-self.tstart ,self.cf.gettime(sec)[1]-self.tstart) for sec in sessions]
        else:
            if isinstance(window,float) or isinstance(window,int):
                phase_nb = int(np.ceil((self.t_end-self.t_start)/(window*3600)))
                self.phases = [(i*window*3600,np.min([(i+1)*window*3600,len(sd[self.mice[0]])])) for i in range(phase_nb)]
            elif isinstance(window, list):
                self.phases = [(st*window[0]*3600,(st+1)*window[0]*3600) for st in window[1]]
            else:
                raise TypeError
        
       
    def calculate_fvalue(self,window='default',treshold = 2, force=False,fols=None,ops=None,which_phase='ALL'):
        
        self.calculate_phases(window=window,which_phase=which_phase)
        
        self.fols = fols
        self.ops = ops
        self.treshold = treshold
 
        
        self.fpatterns = []
        self.opatterns = []
            
        self.f = np.zeros((len(self.mice),len(self.mice),len(self.phases)))
        self.interactions = np.zeros((len(self.phases),len(self.mice),len(self.mice),8,3))
        self.f_sum = np.zeros((len(self.phases),self.lm))
                              
        new_path = os.path.join(self.directory,'PreprocessedData/IteractionsData/')
    
        new_fname_patterns = os.path.join(new_path, 'Patterns_cut_%s_%s.npy'%(self.fname_ending,which_phase))
        new_fname_fpatterns = os.path.join(new_path, 'fpatterns_cut_%s_%s.npy'%(self.fname_ending,which_phase))
        new_fname_opatterns = os.path.join(new_path, 'opatterns_cut_%s_%s.npy'%(self.fname_ending,which_phase))
        new_fname_mice = os.path.join(new_path, 'mice_cut_%s_%s.npy'%(self.fname_ending,which_phase))
        
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
                print('Phase %s. from %sh, to %sh'%(s+1,np.round(ts/3600.,2), np.round(te/3600.,2)))
                self.interactions[s,:,:,:,:] = self.interaction_matrix(ts, te)
                self.independent_data_comparison(ts+self.t_start,te+self.t_start)
                np.save(new_fname_patterns,self.interactions)
                np.save(new_fname_fpatterns,self.fpatterns)
                np.save(new_fname_opatterns,self.opatterns)
                np.save(new_fname_mice,self.mice)
        return self.f

    def tube_dominance_test(self,window='default',which_phase='ALL'):
        
        self.calculate_phases(window=window,which_phase=which_phase)
        self.tube_dominance_matrix = np.zeros((len(self.phases),self.lm,self.lm))
        new_path = utils.check_directory(self.directory,'PreprocessedData/IteractionsData/')
        
        new_fname_ = os.path.join(new_path, 'Patterns_%s.npy'%self.fname_ending)
        new_fname_mice = os.path.join(new_path, 'mice_cut_%s_%s.npy'%(self.fname_ending,which_phase))
  
        for s in range(len(self.phases)):
            ts, te = self.phases[s]
            print('Phase %s. from %sh, to %sh'%(s+1,np.round(ts/3600.,2), np.round(te/3600.,2)))
            self.tube_dominance_matrix[s,:,:] = self.calculate_tube_dominance_matrix(ts, te)
            

        np.save(new_fname_mice,self.mice)
        np.save(new_fname_,self.tube_dominance_matrix)
            
    def interaction_matrix(self,ts, te):
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
    
    def validatePatterns(self,plots_nr = 9, trange = [-3,3]):
        sd = self.ehs.signal_data
        size = np.ceil(np.sqrt(plots_nr))
        t = np.arange(trange[0],trange[1],1.0/self.fs)
        frandom_idx = [np.random.randint(0,len(self.fpatterns)-1) for i in range(plots_nr)]
        orandom_idx = [np.random.randint(0,len(self.opatterns)-1) for i in range(plots_nr)]
        plt.figure()
        plt.suptitle("Random following patterns", fontsize=14, fontweight='bold')
        for i,idx in enumerate(frandom_idx):
            ax = plt.subplot(size, size,i+1)

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

    def calculate_mouse_stats(self,mouse):
        
        follow_stat = {}
        mouse_stats = self.ehs.statistics[mouse]["preference"]
        #print(mouse2_stats)
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
            follow_stat[move][2] = mouse_stats[int(move[0])][index]/np.sum(mouse_stats[int(move[0])])
        return follow_stat


    def mouse_indices(self, mouse1, t1,t2):
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
        return [mouse1_idx_start,mouse1_idx_stop,mouse1_indices]


    def plot_all(self,t1,t2):
        sd = self.ehs.signal_data
        fig = plt.figure()
        ax = []
        lm = len(self.mice)
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
                        print(k)
                        mouse = self.mice[k]
                        sig = sd[mouse][t1*self.fs:t2*self.fs]
                        time = np.linspace(t1,t2,len(sig))
                        ax[k].plot(time,sig,'dk')
                        ax[k].set_title(mouse)
                        ax[k].set_xlabel('time [s]')
                        ax[k].set_ylabel('State')
                        k = k+1                        
                        print(k)

    def calculate_tube_dominance_matrix(self,ts, te):
        """Calculates fvalue matrix for the given period from ts to te
        """
        self.plot_all(ts,te)
        imatrix = np.zeros((self.lm,self.lm))
        for ii, mouse1 in enumerate(self.mice):
            for jj,mouse2 in enumerate(self.mice):
                if ii < jj:
                    imatrix[ii,jj] = self.check_tube_dominance((mouse1,mouse2),ts, te)
                    imatrix[jj,ii] = self.check_tube_dominance((mouse2,mouse1),ts, te)
                    total = imatrix[ii,jj] + imatrix[jj,ii]
                    if total:
                        imatrix[ii,jj] = imatrix[ii,jj]/total
                        imatrix[jj,ii] = imatrix[jj,ii]/total
        return imatrix
    
    def check_tube_dominance(self, (mouse1,mouse2),t1, t2):
        """
        This is an implementation of the tube dominance test, where two mice are placed on opposite sides of a tube
        and one of the mice forces its opponent out of the tube. The test checks if mouse1 entered a tube and backed
        out and if mouse2 remains in the tube and entered the tube from the other side.
        Effectively this function checks how many times mouse2 domineered over mouse1.
        """
        dominance_stats = 0
        sd = self.ehs.signal_data
        try:
            [mouse1_idx_start,mouse1_idx_stop,mouse1_indices] = self.mouse_indices(mouse1, t1,t2)
        except ValueError:
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
                
            idx = end-1
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
    def findpatterns(self,ii,jj,t1, t2):
    
        sd = self.ehs.signal_data
        detected_idx = [[],[]]
        mouse1 = self.mice[ii]
        mouse2 = self.mice[jj]
    
        follow_stat = self.calculate_mouse_stats(mouse2)
        #print(follow_stat)
        try:
            [mouse1_idx_start,mouse1_idx_stop,mouse1_indices] = self.mouse_indices(mouse1, t1,t2)
        except ValueError:
            return np.array(follow_stat.values()), detected_idx
            
        i = mouse1_idx_start - 1
       
        for start in mouse1_indices[mouse1_idx_start:mouse1_idx_stop+1]:
            
            i = i + 1
           
            if i < 2:
                continue

            start_state = int(sd[mouse1][mouse1_indices[i-2]])
            middle_state = int(sd[mouse1][mouse1_indices[i-1]])
            end_state = int(sd[mouse1][start])
            # if abs(start_state-middle_state) !=1 or abs(middle_state-end_state)!=1:
            #     continue
            #print(start_state,middle_state,end_state)
            end = start + self.treshold*self.fs

            unknown_state = end_state == 0 or middle_state == 0 or start_state == 0
            
            in_pipe = end_state%2==1
            # define conditions
            if unknown_state or start_state == end_state or in_pipe:
                continue
            
            try:
                
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
                        #####POPRAWIC
                        #print(np.ceil((period1.index(sd[mouse1][start]))/self.fs))
                        if self.fols!=None:
                            self.fols[np.ceil((period1.index(end_state))/self.fs)]+=1
                        follow_stat[str(start_state)+str(end_state)][0] += 1 
                        #print p, sd[mouse1_indices[i-2],mouse1],[index]
                        detected_idx[0].append((ii,jj,start))

                    elif go_oposite:
                        #print np.ceil((period1.index(op_idx))/self.fs)
                        if self.ops!=None:
                            self.ops[np.ceil((period1.index(op_idx))/self.fs)]+=1
                        follow_stat[str(start_state)+str(end_state)][1] += 1
                        detected_idx[1].append((ii,jj,start))
                                                
            except IndexError:
                print('Err')
                continue
       
        return np.array(follow_stat.values()), detected_idx
    
    def InteractionsPerPair(self,start,end):

        i8states = np.sum(self.interactions[:,:,:,:,start:end],axis=4) # po kierunku -- 0 following, 1 -- avoiding
        return np.sum(i8states ,axis=3) #po wszystkich stanach (24,42, etc.)

    def FollowingAvoidingMatrix(self):
        
        patterns = self.interactions
        rawFA = np.sum(patterns[:,:,:,:,:2] ,axis=3,dtype = 'float32')
        FAPmatrix = np.apply_along_axis(FAprobablity, 3, rawFA)
        return FAPmatrix
    
    def plotTubeDominanceRasters(self,name=None,mice=[]):
        
        if not name:
            name = "TubeDominance_"+self.fname_ending
        subdirectory = 'RasterPlots'
    
        new_path = utils.check_directory(self.directory,subdirectory)
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, aspect='equal')
        if name:
            plt.suptitle(name, fontsize=14, fontweight='bold')

        n_s,n_l,n_f = self.tube_dominance_matrix.shape
        
        for s in range(n_s): #phases
            
            plt.text(0.06+s*0.125, 1.025,self.cf.sections()[s], horizontalalignment='center', verticalalignment='center', fontsize=10,  transform = ax.transAxes)
            # MakeRelationGraph(FAM[s,:,:],IPP[s,:,:],exp,s,key,directory,scalefactor)
            _TDT = self.tube_dominance_matrix[s,:,:]
            pair_labels = []
            pos = 0
            for i in range(n_l): #mice
                
                for j in range(i,n_f): #mice
                
                    if i!=j and  np.isclose(_TDT[i,j],1-_TDT[j,i]):
                        
                        if _TDT[i,j] > 0.5:
                            ax.add_patch(patches.Rectangle((
                                s, -1*pos),1 , 1,facecolor=(1,0,0,0.5*np.round(_TDT[i,j],2))))
                        elif _TDT[i,j] < 0.5:
                            ax.add_patch(patches.Rectangle((
                                s, -1*pos),1 , 1,facecolor=(0,0,1,0.5*np.round(1-_TDT[i,j],2))))
                        else:
                            ax.add_patch(patches.Rectangle((
                                s, -1*pos),1 , 1,facecolor=(0.5,.5,.5,.5)))#grey patch if both equal and 0.5

                    
                            
                    if i!=j:
                        if not mice:
                            pair_labels.append(str(j+1)+'|'+str(i+1))
                        else:
                            pair_labels.append(mice[j]+'|'+mice[i])
                        pos+=1
            for i in range(8-n_s):
                ax.add_patch(patches.Rectangle((
                                        n_s+i, -pos+1),1, pos,facecolor="lightgrey"))
        plt.axis([0,8,-pos+1,1])
        fig.subplots_adjust(left=0.25)

        ax.set_aspect('auto')
        ax.xaxis.grid()
        ax.xaxis.set_ticklabels([])
        ax.get_yaxis().set_ticks([-1*i+0.5 for i in range(pos)])
        ax.set_yticklabels(pair_labels,fontsize=10)
        plt.xlabel("session")
        plt.ylabel("following strength in pair")
        plt.savefig(os.path.join(new_path,name+'.png'),transparent=False, bbox_inches=None, pad_inches=3,frameon=None)
        # plt.show()
 
    def remove_zeros(self):
        indices = set() #set of indices, where position is zero
        for mouse in self.mice:
            zeros = np.where(self.sd[mouse] == 0)[0]
            if len(zeros):
                indices.update(set(zeros))
        indices = sorted(list(indices))
        mask = np.ones(len(self.sd['time']), dtype=bool)
        mask[indices] = False
        for key in self.sd.keys():
            self.sd[key] = self.sd[key][mask]

    def get_mask(self, starttime,endtime):
        arr = self.sd['time']
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            return (idcs[0], idcs[-1] + 1)
        if len(idcs) == 1:
            return (idcs[0], idcs[0] + 1)
        return (0,0)
         
    def single_state_probability(self,mouse,starttime,endtime):
        result = np.zeros((8,))
        _mask_slice = self.get_mask(starttime,endtime)
        for i in range(1,9):
            result[i-1] += len(np.where(self.sd[mouse][_mask_slice[0]:_mask_slice[1]] == i)[0])/len(self.sd[mouse][_mask_slice[0]:_mask_slice[1]])
        return result

    def vector_state_probability(self,signal,starttime,endtime):
        shape = tuple([8 for i in self.mice])
        print(shape)
        result = np.empty(shape)
        _mask_slice = self.get_mask(starttime,endtime)
        length = len(signal['time'][_mask_slice[0]:_mask_slice[1]])
        print(length)
        new_signal = np.empty((length,len(self.mice)),dtype=np.int)
        for i,mouse in enumerate(self.mice):
            new_signal[:,i] = np.array(signal[mouse][_mask_slice[0]:_mask_slice[1]],dtype=int)
        for row in new_signal:
            print(row)
            result[row-1] += 1
        return result/length
    
    def generate_independent_data_sequences(self,starttime,endtime):
        probabilities = {}
        for mouse in self.mice:
            probabilities[mouse] = self.single_state_probability(mouse,starttime,endtime)
        mask = self.get_mask(starttime,endtime)
        length = len(self.sd['time'][mask[0]:mask[1]])
        surrogate_data = {}
        positions = [int(x) for x in range(1,9)]
        for mouse in self.mice:
            surrogate_data[mouse] = numpy.random.choice(positions,length,p=probabilities[mouse])

        print('generated surrogate data')
        surrogate_data['time'] = np.linspace(starttime,endtime,length)
        return surrogate_data


    def independent_data_comparison(self,starttime,endtime):
        sur_data = self.generate_independent_data_sequences(starttime,endtime)
        independent_states = self.vector_state_probability(sur_data,starttime,endtime)
        sur_data_vec = independent_states[np.where(independent_states !=0)]
        independent_states = sur_data_vec.sort()
        print(starttime,endtime,(endtime-starttime)/3600)
        states = self.vector_state_probability(self.sd,starttime,endtime)
        print('where')
        non_zeros = states[np.where(states !=0)]
        states = non_zeros.sort()
        print(states)
        

        plt.figure()
        plt.semilogy(independent_states,label='Independent model')
        plt.semilogy(states, label='Empirical')
        plt.xlabel('Location pattern rank')
        plt.ylabel('Location pattern probability')
        plt.show()
        
        
def createRandomExpepiments(paths,core='Random_test' ):
    if not os.path.exists('../PreprocessedData/RandomData/'):
        os.makedirs('../PreprocessedData/RandomData/')
    statistics_rd = {}
    for j in range(len(paths)):
        rdE = Experiment('Random_test%s'%(j+1),exp_name=core)
        rdname = 'Random_test%s.pkl'%(j+1)
        print(rdE.sd.shape)
        rd = np.zeros(rdE.sd.shape)
        print('####%s#####'%j)
        ######Create test file#########
        for i in range(j,len(paths)+j):
            E = Experiment(paths[i%len(paths)])
            print(i, i%E.lm,E.sd.shape)
            end = np.min([len(rd[:,i%rdE.lm]),len(E.sd[:,i%E.lm])])
            rd[:end,i%rdE.lm] = E.sd[:end,i%E.lm]
            statistics_rd[rdE.mice[i%rdE.lm]] =  E.ehs.statistics[E.mice[i%E.lm]]
            #print statistics_rd[E_t.mice[i%E_t.lm]]["preference"]
        rdE.mice = rdE.mice[:len(paths)]
        rdE.lm = len(paths)
        rdE.sd = rd[:,:end]
        rdE.ehs.statistics = statistics_rd
        with open('../PreprocessedData/RandomData/'+rdname, 'wb') as output:
           pickle.dump(rdE,output,pickle.HIGHEST_PROTOCOL)
           del rdE

def preprocessData(names,window,ts=3):
    IPP = {}
    APP = {}
    FPP = {}
    FAM = {}
    directory = {}
    phases = {}
    for key in names.keys():
        IPP[key] = {}
        APP[key] = {}
        FPP[key] = {}
        FAM[key] = {}
        directory[key] = {}
        phases[key] = {}
        for path in names[key]:
            print(path)
            # if key=="RD":
            #     with open('../PreprocessedData/RandomData/'+path+'.pkl', "rb") as input_file:
            #         E = pickle.load(input_file)
            # else:
            E = Experiment(path, _ant_pos=antenna_positions[path])
            E.calculate_fvalue(window=window,treshold=ts,force=True)
            IPP[key][path] = E.InteractionsPerPair(0,2)
            APP[key][path] = E.InteractionsPerPair(1,2)
            FPP[key][path] = E.InteractionsPerPair(0,1)
            FAM[key][path] = E.FollowingAvoidingMatrix()
            directory[key][path] = E.directory
            phases[key][path] = E.cf.sections()
            #E.validatePatterns()
    return IPP,FPP,APP,FAM,directory, phases

def Interpersec(names,ts=3, directory='InterPerSec'):
    if not os.path.exists('../Results/%s/'%directory):
        os.makedirs('../Results/%s/'%directory)
    for key in names.keys():
        fols = np.zeros(ts-1)
        ops  = np.zeros(ts-1)
        for path in names[key]:
            print(path)
            if key=="RD":
                with open('../PreprocessedData/RandomData/'+path+'.pkl', "rb") as input_file:
                    E = pickle.load(input_file)
            else:
                E = Experiment(path)
            E.calculate_fvalue(window = 24, treshold =ts, force=True, fols=fols,ops=ops)
            fols = E.fols
            ops = E.ops
        print('###########%s###########'%key)
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1,ts,1),fols*1./np.sum(fols))
        plt.ylim(0,0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of followings")
        plt.savefig('../Results/'+directory+'/'+key+'folinsec.png')
        #plt.show()
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1,ts,1),ops*1./np.sum(ops))
        plt.ylim(0,0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of avoidance")
        plt.savefig('../Results/'+directory+'/'+key+'opsinsec.png')
        #plt.show()
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1,ts,1),(fols+ops)*1./np.sum(fols+ops))
        plt.ylim(0,0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of interactions")
        plt.savefig('../Results/'+directory+'/'+key+'interinsec.png')
        #plt.show()

def binomial_probability(s, p, n):
    prob = 0.0

    for j in range(s, n + 1):
        prob += scipy.special.binom(n, j) * p**j * (1 - p)**(n-j)
   
    return prob


def FAprobablity(a,p=0.5):
    
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
    
def FollowingAvoidingMatrix(names):
    FAP = {}
    for key in names.keys(): 
        print(key)
        FAP[key] = []
        for path in names[key]:
            new_path = utils.check_directory(path,'PreprocessedData/IteractionsData/')
            fname = os.path.join(new_path,'interactions.npy')
            patterns = np.load(fname)
            FAP[key].append(easyFAP(patterns,p=0.5))
    return FAP

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
    ts = 3
    experiments = load_experiments_info("experiments_desc.csv")
    comparisons = load_comparisons_info("comparisons.csv")
    #print(experiments)
    #print comparisons.keys()
    ##for i in comparisons:
    ##    names, colors = group_data(i,comparisons,experiments, color_lst = ["red","green", "blue"])
    #    print names
    #    preprocessData(names,window = 12,ts=3)
    names, colors = group_data('KO-WT|mouse|FX|females|1',comparisons,experiments, color_lst = ["red","green", "blue"])
    #names, colors = group_data('VPA-CTRL-NaCl|mouse|C57|males|1',comparisons,experiments, color_lst = ["red","green", "blue"])
    
    #createRandomExpepiments(exp_paths)
    #Interpersec(names,ts=200)
    IPP, FPP, APP, FAM, directories, phases = preprocessData(names,window=12,ts=ts)
    print(IPP.keys())
    
    
    scalefactor = np.max([np.max(IPP["KO"][path]) for path in IPP["KO"]]+[np.max(IPP["WT"][path]) for path in IPP["WT"]])
    for key in IPP:
        for path in IPP[key]:
            oneRasterPlot(directories[key][path],FAM[key][path],IPP[key][path],phases[key][path],'Interactions',scalefactor)

        #oneRasterPlot(directories[key],FAM[key],FPP[key],phases[key],'Following',scalefactor)

        #oneRasterPlot(directories[key],FAM[key],APP[key],phases[key],'Avoiding',scalefactor)

        #plot_graph(FAM[key],1, base_name = 'following_graph', nr='')

    result = []
    pair_inc = []
    pair_long = []
    
    # stats = {}
    # stats["KO"] = {}
    # stats["KO"]["SLD"] = []
    # stats["KO"]["NFD"] = []
    # stats["KO"]["SLA"] = []
    # stats["KO"]["NFA"] = []
    # for i in range(3):
    #     _FAM=FAM["KO"][i]
    #     LSD, LSA = longest_sequence(_FAM)
    #     FSD, FSA = follsactive(_FAM)
    #     stats["KO"]["SLA"]+=LSA
    #     stats["KO"]["NFA"]+=FSA
    #     stats["KO"]["NFD"]+=FSD
    #     stats["KO"]["SLD"]+=LSD
    #     print(np.mean(LSA), np.mean(FSA), FAM["KO"][i].shape[2])
    #     #plt.hist(LSD)
    #     #plt.hist(FSD)
    #     ###plt.show()
    # print('#####################WT#######################')
    # stats["WT"] = {}
    # stats["WT"]["SLD"] = []
    # stats["WT"]["NFD"] = []
    # stats["WT"]["SLA"] = []
    # stats["WT"]["NFA"] = []
    # for i in range(4):
    #     _FAM=FAM["WT"][i]
    #     LSD, LSA = longest_sequence(_FAM)
    #     FSD, FSA = follsactive(_FAM)
    #     stats["WT"]["SLA"]+=LSA
    #     stats["WT"]["NFA"]+=FSA
    #     stats["WT"]["NFD"]+=FSD
    #     stats["WT"]["SLD"]+=LSD
    
    #    plt.hist(LSA)
    #    plt.hist(FSA)
    #    ##plt.show()
    
    
    
    
