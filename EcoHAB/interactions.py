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

    def __init__(self, path,_ant_pos=None,which_phase='ALL',mask=None):
        
        self.path = path
        self.directory = utils.results_path(path)
        print(self.directory)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        self.cf = ExperimentConfigFile(self.path)
        
        if which_phase != 'ALL' and mask == None:
            try:
                mask = self.cf.gettime(which_phase)
            except ConfigParser.NoSectionError:
                mask = None
        
        self._remove_phases(mask)            

        
        self.ehs = EcoHab.EcoHabSessions9states(self.path, _ant_pos=_ant_pos,mask=mask,shortest_session_threshold=0)

        self.fs = self.ehs.fs
        self.sd =  self.ehs.signal_data
       
        mice = list(self.ehs.mice)
        
        self.mice = filter(lambda x: len(self.ehs.getstarttimes(x)) > 30, mice)
        
        if mice != self.mice:
            print("False")
            for mouse in mice:
                if mouse not in self.mice:
                    del self.sd[mouse]
                    del self.ehs.signal_data[mouse]
                    
        self.lm = len(self.mice)
        
    def calculate_fvalue(self,window='default',treshold = 2, force=False,fols=None,ops=None,which_phase='ALL'):
        self.fols = fols
        self.ops = ops
        self.treshold = treshold
 
        self.tstart, self.tend = self.cf.gettime('ALL')
        self.fpatterns = []
        self.opatterns = []
        


        if window=='default':
            sessions = filter(lambda x: x.endswith('dark') or x.endswith('light'), self.cf.sections())
            self.phases = [(self.cf.gettime(sec)[0]-self.tstart ,self.cf.gettime(sec)[1]-self.tstart) for sec in sessions]
        else:
            if isinstance(window,float) or isinstance(window,int):
                self.phases = [(i*window*3600,np.min([(i+1)*window*3600,len(self.sd[self.mice[0]])])) for i in range(int(np.ceil(len(self.sd[self.mice[0]])*1.0/(window*3600*self.fs))))]
            elif isinstance(window, list):
                self.phases = [(st*window[0]*3600,(st+1)*window[0]*3600) for st in window[1]]
            else:
                raise TypeError
            
        self.f = np.zeros((len(self.mice),len(self.mice),len(self.phases)))
        self.interactions = np.zeros((len(self.phases),len(self.mice),len(self.mice),8,3))
        self.f_sum = np.zeros((len(self.phases),self.lm))
                              
        new_path = os.path.join(self.directory,'PreprocessedData/IteractionsData/')
        new_fname_patterns = os.path.join(new_path, 'Patterns_%s.npy'%which_phase)
        new_fname_fpatterns = os.path.join(new_path, 'fpatterns_%s.npy'%which_phase)
        new_fname_opatterns = os.path.join(new_path, 'opatterns_%s.npy'%which_phase)
        
        if not os.path.exists(os.path.dirname(new_path)):
            os.makedirs(os.path.dirname(new_path))
    
        # try:
        #     self.interactions = np.load(new_fname_patterns)
        #     self.fpatterns = np.load(new_fname_fpatterns)
        #     self.opatterns = np.load(new_fname_opatterns)
            
        # except IOError:
        
        for s in range(len(self.phases)):
            ts, te = self.phases[s]
            #print('Phase %s. from %sh, to %sh'%(s+1,np.round(ts/3600.,2), np.round(te/3600.,2)))
            self.interactions[s,:,:,:,:] = self.interaction_matrix(ts, te)
            np.save(new_fname_patterns,self.interactions)
            np.save(new_fname_fpatterns,self.fpatterns)
            np.save(new_fname_opatterns,self.opatterns)

        return self.f
     
    
    def interaction_matrix(self,ts, te):
        """Calculates fvalue matrix for the given period from ts to te
        """
        imatrix = np.zeros((self.lm,self.lm,8,3))
        for ii, mouse1 in enumerate(self.mice):
            for jj,mouse2 in enumerate(self.mice):
                if ii < jj:
                    imatrix[ii,jj,:,:],patterns = self.findpatterns((ii,jj),ts, te)
                    self.fpatterns+=patterns[0]
                    self.opatterns+=patterns[1]
                    imatrix[jj,ii,:,:],patterns = self.findpatterns((jj,ii),ts, te)
                    self.fpatterns+=patterns[0]
                    self.opatterns+=patterns[1]
        return imatrix
    
    def validatePatterns(self,plots_nr = 9, trange = [-3,3]):
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
            plt.plot(t,self.sd[mouse1][s-3*self.fs:s+3*self.fs]-0.05,'ro',label="leader")
            plt.plot(t,self.sd[mouse2][s-3*self.fs:s+3*self.fs]+0.05,'bo',label="follower")
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
            plt.plot(t,self.sd[mouse1][s-3*self.fs:s+3*self.fs]-0.05,'ro',label="leader")
            plt.plot(t,self.sd[mouse2][s-3*self.fs:s+3*self.fs]+0.05,'bo',label="follower")
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
    
    def findpatterns(self, (ii,jj),t1, t2):
        
        fs = self.fs

        detected_idx = [[],[]]
        mouse1 = self.mice[ii]
        mouse2 = self.mice[jj]

        follow_stat = self.calculate_mouse_stats(mouse2)
        mouse1_indices = np.where(((np.roll(self.sd[mouse1], 1) - self.sd[mouse1]) != 0))[0]
        try:
            mouse1_idx_start = np.where(mouse1_indices > t1*fs)[0][0]
        except IndexError:
            return follow_stat.values(),detected_idx

        try:    
            mouse1_idx_stop =  np.where(mouse1_indices < t2*fs)[0][-1]
        except IndexError:
            return follow_stat.values(),detected_idx
            
        i = mouse1_idx_start - 1
       
        for start in mouse1_indices[mouse1_idx_start:mouse1_idx_stop+1]:
            
            i = i + 1
            start_state = int(self.sd[mouse1][mouse1_indices[i-2]])
            middle_state = int(self.sd[mouse1][mouse1_indices[i-1]])
            end_state = int(self.sd[mouse1][start])
 
            end = start + self.treshold*fs

            unknown_state = end_state == 0 or middle_state == 0 or start_state == 0
            
            in_pipe = end_state%2==1
            # define conditions
            if unknown_state or start_state == end_state or in_pipe:
                continue
            
            try:
                
                period1 = list(self.sd[mouse2][start:end])
                
                period2 = list(self.sd[mouse2][start-2*fs:start])

                period3 = list(self.sd[mouse2][start-int(0.1*fs):start])
              
                op_idx = (2*start_state-end_state-1)%8+1 #opposite direction
                
                same_start = start_state in period2 #POPRAWIC!!!!!!!

                first_mouse1 = end_state not in period3 and op_idx not in period3
                
                followed = period1.count(end_state) > 0
                
                go_oposite = op_idx in period1
                
                if followed and go_oposite:
                    followed = period1.index(end_state) < period1.index(op_idx)
                    go_oposite = not followed
                    
                    
                if same_start and first_mouse1:
                    #print start_state,end_state, op_idx
                    #print start_state,end_state,index
                    #print self.sd[mouse1_indices[i-2],mouse1],self.sd[s,mouse1], (2*self.sd[mouse1_indices[i-2],mouse1]-self.sd[s,mouse1]-1)%8+1,followed,go_oposite
                    if followed:
                        #####POPRAWIC
                        #print(np.ceil((period1.index(self.sd[mouse1][start]))/self.fs))
                        if self.fols!=None:
                            self.fols[np.ceil((period1.index(end_state))/self.fs)]+=1
                        follow_stat[str(start_state)+str(end_state)][0] += 1 
                        #print p, self.sd[mouse1_indices[i-2],mouse1],[index]
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
        v = round(pf,3)#0.5-pf
    elif pa < pf and pa < 0.05:
        v = round(-pa,3)#pa-0.5
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
    
    
    
    
