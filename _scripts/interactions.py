#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:48:42 2016

@author: jmaka
"""
import random
import time
import EcoHab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
import locale
from ExperimentConfigFile import ExperimentConfigFile
from experiments_info import smells, antenna_positions
import os
import pickle
from mpl_toolkits.mplot3d import Axes3D
from experiments_info import smells, antenna_positions
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import matplotlib.colors as mcol
import matplotlib.patches as patches
import scipy.stats as st
from plotfunctions import barplot


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def plotfsec(self,fols, ops, to_file = True):
    x = np.arange(1,ts+1)
    plt.stem(x,self.fols/np.sum(self.fols)*100)
    plt.axis([0,max_ts+1,0,15])
    plt.xlabel('1 s time bin')
    plt.ylabel('Percetage of positive f')
    plt.savefig('../Results/rdn_fols%s'%max_ts+'.png')
    plt.show()
    plt.stem(x,self.ops/np.sum(self.ops)*100)
    plt.axis([0,max_ts+1,0,15])
    plt.xlabel('1 s time bin')
    plt.ylabel('Percetage of negative f')
    plt.savefig('../Results/rdn_ops%s'%max_ts+'.png')
    plt.show()
    
    
    
class Experiment(object):
    def __init__(self, path, exp_name=None):
        if not exp_name:
            exp_name = path
        self.path = path
        self.directory = os.path.join('..','Results/'+path)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.ehd = EcoHab.EcoHabData(os.path.join('..','RawData',exp_name), _ant_pos=antenna_positions[exp_name])
        self.ehs = EcoHab.EcoHabSessions9states(self.ehd,shortest_session_threshold=0)
        self.fs = self.ehs.fs
        self.sd =  self.ehs.signal_data
        self.cf = ExperimentConfigFile(os.path.join('..','RawData',exp_name))
        mice = list(self.ehd.mice)
        self.mice = filter(lambda x: len(self.ehs.getstarttimes(x)) > 30, mice)
        self.lm = len(self.mice)
        
    def calculate_fvalue(self,window='default',treshold = 2,min_interactions = 1, force=False,fols=None,ops=None):
        self.fols = fols
        self.ops = ops
        self.treshold = treshold
        self.min_interactions = min_interactions
        self.tstart, self.tend = self.cf.gettime('ALL')
        if window=='default':
            sessions = filter(lambda x: x.endswith('dark') or x.endswith('light'), self.cf.sections())
            self.phases = [(self.cf.gettime(sec)[0]-self.tstart ,self.cf.gettime(sec)[1]-self.tstart) for sec in sessions]
        else:
            if isinstance(window,float) or isinstance(window,int):
                self.phases = [(i*window*3600,np.min([(i+1)*window*3600,self.sd.shape[0]])) for i in range(int(np.ceil(self.sd.shape[0]*1.0/(window*3600*self.fs))))]
            elif isinstance(window, list):
                self.phases = [(st*window[0]*3600,(st+1)*window[0]*3600) for st in window[1]]
            else:
                raise TypeError
        self.f = np.zeros((len(self.mice),len(self.mice),len(self.phases)))
        self.interactions = np.zeros((len(self.phases),len(self.mice),len(self.mice),8,3))
        self.f_sum = np.zeros((len(self.phases),self.lm))
        if 'following.npy' in list(os.walk(os.path.join('..','Results',self.path)))[0][2] and not force:
            self.f = np.load(os.path.join('..','Results/'+self.path+'/')+'following.npy')
            return self.f
        else:
            for s in range(len(self.phases)):
                ts, te = self.phases[s]
                print 'Phase %s. from %sh, to %sh'%(s+1,np.round(ts/3600.,2), np.round(te/3600.,2))
                self.interactions[s,:,:,:,:] = self.interaction_matrix(ts, te)
            if not os.path.exists(os.path.join('..','PreprocessedData/IteractionsData/')):
                os.makedirs(os.path.join('..','PreprocessedData/IteractionsData/'))
            np.save(os.path.join('..','PreprocessedData/IteractionsData/')+'%s.npy'%self.path,self.interactions)
            return self.f
     
    
    def interaction_matrix(self,ts, te):
        """Calculates fvalue matrix for the given period from ts to te
        """
        imatrix = np.zeros((self.lm,self.lm,8,3))
        for ii in range(len(self.mice)):
            for jj in range(len(self.mice)):
                if ii < jj:
                    imatrix[ii,jj,:,:],_ = self.findpatterns((ii,jj),ts, te)
                    imatrix[jj,ii,:,:],_ = self.findpatterns((jj,ii),ts, te)   
#                    t = np.arange(-3,3,1.0/self.fs)
#                    for i in range(len(_[0])):
#                        s = _[0][i]
#                        plt.plot(t,self.sd[s-3*self.fs:s+3*self.fs,ii]-0.05,'ro',label="leader")
#                        plt.plot(t,self.sd[s-3*self.fs:s+3*self.fs,jj]+0.05,'bo',label="follower")
#                        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#
#                        plt.axis([-3.1,3.1,-0.5,9.5])
#                        plt.show()
                    
        return imatrix
    
    def findpatterns(self, (m1,m2),t1, t2):
        fs = self.fs
        sd = self.sd
        detected_idx = [[],[]]
        follow_stat = {}
        m1_idx = np.where(((np.roll(sd[:,m1], 1) - sd[:,m1]) != 0))[0]
        m2_stats = self.ehs.statistics[self.mice[m2]]["preference"]
        moves = ['24','42','46','64','68','86','82','28']
        for m in moves:
            follow_stat[m] = np.zeros(3)
            start_st, end_st = int(m[0]),int(m[1])
            if end_st == 2 and start_st==8:
                index = 0
            elif end_st == 8 and start_st==2:
                index = 1
            elif start_st <end_st:
                index = 0
            else:
                index = 1
            follow_stat[m][2] = m2_stats[int(m[0])][index]/np.sum(m2_stats[int(m[0])])
        #print m2_stats
        #print 'ruchliwosc', len(m1_idx)
        for i in range(2,len(m1_idx)):
            if m1_idx[i] > t1*fs and m1_idx[i]<t2*fs:
                s = m1_idx[i]
                start_st = int(sd[m1_idx[i-2],m1])
                end_st = int(sd[s,m1])
                e =s+self.treshold*fs
                try:
                    period1 = list(sd[s:e,m2])
                    period2 = list(sd[s-2*fs:s,m2])
                    period3 = list(sd[s-int(0.1*fs):s,m2])
                    # define conditions
                    unknown_state = sd[s,m1]==0 or int(sd[m1_idx[i-1],m1])==0
                    unknown_previous_states = (sd[m1_idx[i-1],m1] ==0 and (sd[m1_idx[i-2],m1] ==0)) or (sd[m1_idx[i-1],m1] !=0 and (sd[m1_idx[i-2],m1] ==0))
                    in_pipe =  sd[s,m1]%2==1
                    if unknown_state or unknown_previous_states or in_pipe or start_st==end_st:
                        continue
                    
                    op_idx = (2*sd[m1_idx[i-2],m1]-sd[s,m1]-1)%8+1
                    same_start = sd[m1_idx[i-2],m1] in period2 #POPRAWIC!!!!!!!
                    first_m1 = sd[s,m1] not in period3 and op_idx not in period3
                    followed = period1.count(sd[s,m1])>0
                    go_oposite = op_idx in period1
                    if followed and go_oposite:
                        followed = period1.index(sd[s,m1])<period1.index(op_idx)
                        go_oposite = not followed
                    if same_start and first_m1:
                        #print start_st,end_st, op_idx
                        #print start_st,end_st,index
                        #print sd[m1_idx[i-2],m1],sd[s,m1], (2*sd[m1_idx[i-2],m1]-sd[s,m1]-1)%8+1,followed,go_oposite
                        if followed:
                            #####POPRAWIC
                            #print np.ceil((period1.index(sd[s,m1]))/self.fs)
                            if self.fols!=None:
                                self.fols[np.ceil((period1.index(sd[s,m1]))/self.fs)]+=1
                            follow_stat[str(start_st)+str(end_st)][0] += 1 
                            #print p, sd[m1_idx[i-2],m1],[index]
                            detected_idx[0].append(s)
                        elif go_oposite:
                            #print np.ceil((period1.index(op_idx))/self.fs)
                            if self.ops!=None:
                                self.ops[np.ceil((period1.index(op_idx))/self.fs)]+=1
                            follow_stat[str(start_st)+str(end_st)][1] += 1
                            detected_idx[1].append(s)
                except IndexError:
                    print 'Err'
                    continue
        return np.array(follow_stat.values()), detected_idx
        


def plotphist(data,names,colors,to_file = True,directory = 'Interactions',vrange = [], prange = [],nbins = 24):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    stats = {}
    for key in data.keys():
        stats[key] = {}
        stats[key][directory] = []
        stats[key]["mean"] = []
        stats[key]["median"] = []
        for i in range(len(data[key])):
            x = data[key][i].flatten()
            stats[key][directory]+=list(x[x!=0])
            stats[key]["mean"].append(np.mean(x[x!=0]))
            stats[key]["median"].append(np.median(x[x!=0]))
            plt.suptitle(names[key][i], fontsize=14, fontweight='bold')
            plt.hist(x, bins = np.arange(vrange[0],vrange[1],(vrange[1]-vrange[0])*1./nbins),normed=True)
            plt.axvline(np.mean(x), color=colors[key], linestyle='solid', linewidth=2)
            plt.axvline(np.median(x), color=colors[key], linestyle='dashed', linewidth=2)
            plt.ylim(prange)
            plt.savefig('../Results/'+directory+'/'+names[key][i]+'.png')
            plt.show()
    for key in data.keys():
        for i in range(len(data[key])):
            x = data[key][i].flatten()
            plt.subplot(2,1,1)
            plt.axvline(np.mean(x[x!=0]), color=colors[key], linestyle='solid', linewidth=2)
            plt.subplot(2,1,2)
            plt.axvline(np.median(x[x!=0]), color=colors[key], linestyle='dashed', linewidth=2)
    plt.savefig('../Results/'+directory+'/'+"WT_KO_RD"+'.png')
    plt.show()
    plt.show()
    return stats
        
def createRandomExpepiments(paths,core='Random_test' ):
    if not os.path.exists('../PreprocessedData/RandomData/'):
        os.makedirs('../PreprocessedData/RandomData/')
    statistics_rd = {}
    for j in range(len(paths)):
        rdE = Experiment('Random_test%s'%(j+1),exp_name=core)
        rdname = 'Random_test%s.pkl'%(j+1)
        print rdE.sd.shape
        rd = np.zeros(rdE.sd.shape)
        print '####%s#####'%j
        ######Create test file#########
        for i in range(j,len(paths)+j):
            E = Experiment(paths[i%len(paths)])
            print i, i%E.lm,E.sd.shape
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
    for key in names.keys():
        for path in names[key]:
            print path
            if key=="RD":
                with open('../PreprocessedData/RandomData/'+path+'.pkl', "rb") as input_file:
                    E = pickle.load(input_file)
            else:
                E = Experiment(path)
            E.calculate_fvalue(window = window, treshold =ts, force=True)
            
def Interpersec(names,ts=3, directory='InterPerSec'):
    if not os.path.exists('../Results/%s/'%directory):
        os.makedirs('../Results/%s/'%directory)
    for key in names.keys():
        fols = np.zeros(ts-1)
        ops  = np.zeros(ts-1)
        for path in names[key]:
            print path
            if key=="RD":
                with open('../PreprocessedData/RandomData/'+path+'.pkl', "rb") as input_file:
                    E = pickle.load(input_file)
            else:
                E = Experiment(path)
            E.calculate_fvalue(window = 24, treshold =ts, force=True, fols=fols,ops=ops)
            fols = E.fols
            ops = E.ops
        print '###########%s###########'%key
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1,ts,1),fols*1./np.sum(fols))
        plt.ylim(0,0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of followings")
        plt.savefig('../Results/'+directory+'/'+key+'folinsec.png')
        plt.show()
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1,ts,1),ops*1./np.sum(ops))
        plt.ylim(0,0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of avoidance")
        plt.savefig('../Results/'+directory+'/'+key+'opsinsec.png')
        plt.show()
        plt.suptitle(key, fontsize=14, fontweight='bold')
        plt.stem(np.arange(1,ts,1),(fols+ops)*1./np.sum(fols+ops))
        plt.ylim(0,0.5)
        plt.xlabel("time bin [s]")
        plt.ylabel("avarage number of interactions")
        plt.savefig('../Results/'+directory+'/'+key+'interinsec.png')
        plt.show()
          
def InteractionsPerPair(names):
    IPP = {}
    for key in names.keys():
        IPP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..','PreprocessedData/IteractionsData/')+'%s.npy'%path)
            i8states = np.sum(patterns[:,:,:,:,:2],axis=4)
            interactions = np.sum(i8states ,axis=3)
            IPP[key].append(interactions)
    return IPP

def FollowingPerPair(names):
    FPP = {}
    for key in names.keys():
        FPP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..','PreprocessedData/IteractionsData/')+'%s.npy'%path)
            i8states = np.sum(patterns[:,:,:,:,:1],axis=4)
            interactions = np.sum(i8states ,axis=3)
            FPP[key].append(interactions)
    return FPP    

def AvoidingPerPair(names):
    APP = {}
    for key in names.keys():
        APP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..','PreprocessedData/IteractionsData/')+'%s.npy'%path)
            i8states = np.sum(patterns[:,:,:,:,1:2],axis=4)
            interactions = np.sum(i8states ,axis=3)
            APP[key].append(interactions)
    return APP 
    
def factorial(n): 
    if n < 2: return 1
    return reduce(lambda x, y: x*y, xrange(2, int(n)+1))

def prob(s, p, n):
    x = 1.0 - p
    a = n - s
    b = s + 1
    c = a + b - 1
    prob = 0.0
    for j in xrange(a, c + 1):
        prob += factorial(c) / (factorial(j)*factorial(c-j)) \
                * x**j * (1 - x)**(c-j)
    return prob
def FAprobablity(a):
    pf,pa=prob(int(a[1]), 0.5, int(a[0]+a[1])), prob(int(a[0]), 0.5, int(a[0]+a[1]))
    if pf<pa and pf<0.05:
        v = round(pf,3)#0.5-pf
    elif pa<pf and pa<0.05:
        print pa
        v = round(-pa,3)#pa-0.5
    else:
        v = 0.
    return v
def easyFAP(patterns,p=0.5):
    rawFA = np.sum(patterns[:,:,:,:,:2] ,axis=3,dtype = 'float32')
    FAPmatrix = np.apply_along_axis(FAprobablity, 3, rawFA)
    return FAPmatrix
    
def FollowingAvoidingMatrix(names):
    FAP = {}
    for key in names.keys(): 
        print key
        FAP[key] = []
        for path in names[key]:
            print path
            patterns = np.load(os.path.join('..','PreprocessedData/IteractionsData/')+'%s.npy'%path)
            FAP[key].append(easyFAP(patterns,p=0.5))
    return FAP
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

def createRasterPlots(FAM,IPP,names,scalefactor,to_file = True,directory = 'RasterPlots'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    for key in names.keys():
            for exp in range(len(names[key])):
                fig = plt.figure(figsize=(12,12))
                ax = fig.add_subplot(111, aspect='equal')
                plt.suptitle('%s'%(names[key][exp]), fontsize=14, fontweight='bold')
                n_s,n_l,n_f = FAM[key][exp].shape
                for s in range(n_s):
                    plt.text(0.06+s*0.125, 1.025,ExperimentConfigFile('../RawData/'+names[key][exp]).sections()[s],
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=10,
                     transform = ax.transAxes)
                    #MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor)
                    _FAM = FAM[key][exp][s,:,:]
                    _IPP = IPP[key][exp][s,:,:]
                    n_l,n_f = _FAM.shape
                    print n_l, n_f
                    pair_labels = []
                    pos = 0
                    for i in range(n_l):
                        for j in range(i,n_f):
                            if i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[i,j]*0.5/scalefactor)))  
                            elif i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[i,j]*0.5/scalefactor)))
                            if i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[j,i]*0.5/scalefactor))) 
                            elif i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[j,i]*0.5/scalefactor)))
                            if i!=j:
                                #ax.add_patch(patches.Rectangle((0, -1*pos+1),8,1,facecolor="black",fill=False))
                                pair_labels.append(str(i+1)+'|'+str(j+1))
                                #pair_labels.append(str(j+1)+'|'+str(i+1))
                                pos+=1
                for i in range(8-n_s):
                    ax.add_patch(patches.Rectangle((
                                        n_s+i, -pos+1),1, pos,facecolor="lightgrey"))
                plt.axis([0,8,-pos+1,1])
                ax.set_aspect('auto')
                ax.xaxis.grid()
                ax.xaxis.set_ticklabels([])
                ax.get_yaxis().set_ticks([-1*i+0.5 for i in range(pos)])
                ax.set_yticklabels(pair_labels)
                plt.xlabel("session")
                plt.ylabel("following strength in pair")
                plt.savefig('../Results/%s/%s.png'%(directory,names[key][exp]))
                #plt.show()
                plt.close(fig)   

def createRasterPlotsSUM(FAM,IPP,names,scalefactor,to_file = True,directory = 'RasterPlotsSUM'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    for key in names.keys():
            fig = plt.figure(figsize=(12,12))
            ax = fig.add_subplot(111, aspect='equal')
            plt.suptitle(key, fontsize=14, fontweight='bold')
            exp_pos = 0
            for exp in range(len(names[key])):
                n_s,n_l,n_f = FAM[key][exp].shape
                for s in range(n_s):
                    plt.text(0.06+s*0.125, 1.025,ExperimentConfigFile('../RawData/'+names[key][exp]).sections()[s],
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=10,
                     transform = ax.transAxes)
                    #MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor)
                    _FAM = FAM[key][exp][s,:,:]
                    _IPP = IPP[key][exp][s,:,:]
                    n_l,n_f = _FAM.shape
                    print n_l, n_f
                    pair_labels = []
                    pos = exp_pos+0
                    for i in range(n_l):
                        for j in range(i,n_f):
                            if i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[i,j]*0.5/scalefactor)))  
                            elif i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[i,j]*0.5/scalefactor)))
                            if i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[j,i]*0.5/scalefactor))) 
                            elif i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[j,i]*0.5/scalefactor)))
                            if i!=j:
                                #ax.add_patch(patches.Rectangle((0, -1*pos+1),8,1,facecolor="black",fill=False))
                                pair_labels.append(str(i+1)+'|'+str(j+1))
                                #pair_labels.append(str(j+1)+'|'+str(i+1))
                                pos+=1
                for i in range(8-n_s):
                    ax.add_patch(patches.Rectangle((
                                        n_s+i, -pos+1),1, pos,facecolor="lightgrey"))
                pos -=exp_pos
                exp_pos += pos
            plt.axis([0,8,-exp_pos+1,1])
            ax.set_aspect('auto')
            ax.xaxis.grid()
            ax.xaxis.set_ticklabels([])
            #ax.get_yaxis().set_ticks([-1*i+0.5 for i in range(pos)])
            ax.set_yticklabels([])
            plt.xlabel("session")
            #plt.ylabel("following strength in pair")
            plt.savefig('../Results/%s/%s.png'%(directory,key))
            #plt.show()
            plt.close(fig) 


def CreateRelationGraphs(FAM,IPP,names,scalefactor,to_file = True,directory = 'InteractionsGraphs'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    for key in names.keys():
        for exp in range(len(names[key])):
            n_s,n_l,n_f = FAM[key][exp].shape
            for s in range(n_s):
                MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor,names)
    
def MakeRelationGraph(FAM,IPP,exp,s,key,directory,scalefactor,names, fig = None, xy0 = 0):
    print exp
    n_l,n_f = FAM.shape
    power = []
    conn = []
    for i in range(n_l):
        for j in range(n_f):
            if i!=j and abs(FAM[i,j])<0.05:
                power.append(IPP[i,j])
                conn.append([IPP[i,j]*FAM[i,j]/abs(FAM[i,j]),i+1,j+1])
    #scalefactor = np.max(power)/50
    G = nx.MultiDiGraph(multiedges=True, sparse=True)
    for i in range(len(conn)):
        G.add_edges_from([(conn[i][1],conn[i][2])], weight=conn[i][0])
    edge_colors = [conn[i][0] for i in range(len(conn))]
    size = 10
    pos=nx.circular_layout(G)
    for key in pos.keys():
        pos[key]+=np.array([2,0])
    node_labels = {node:node for node in G.nodes()}  
    if not fig:
        fig = plt.figure(figsize=(10*size,10*size))
    plt.suptitle('%s'%(names[key][exp]), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111, aspect='equal')
    #fig.suptitle(self.path, fontsize=14*size, fontweight='bold')
    nx.draw_networkx_labels(G, pos, labels=node_labels,font_size=120)
    cmap = mcol.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                             colors =[(0, 0, 1), 
                                                      (1, 1., 1), 
                                                      (1, 0, 0)],
                                             N=20-1,
                                             )
    nx.draw(G,pos,ax, node_size=500*size**2,node_color= 'grey',edge_color=edge_colors,edge_cmap=cmap,width = 0,arrows=False, edge_style='dashed',zorder=10)
    for c in conn:
        p = patches.FancyArrowPatch(pos[c[2]],pos[c[1]],connectionstyle='arc3, rad=-0.3',arrowstyle="simple",shrinkA=12*size, shrinkB=12*size,mutation_scale=size*c[0]/scalefactor, color = cmap(c[0]+0.5),zorder=1,alpha=0.5)
        ax.add_patch(p)
    plt.savefig('../Results/%s/%s_exp%s_s%s_graph.png'%(directory,key,exp,s))
    #plt.show()
    plt.close(fig)    


def create_group_graph(FAM,IPP,names,scalefactor,to_file = True,directory = 'InteractionsGraphs'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    max_s = 0
    for key in names.keys():
        for exp in range(len(names[key])):
            n_s,n_l,n_f = FAM[key][exp].shape
            if n_s>max_s:
                max_s = n_s
    figs = []
    size = 10
    for s in range(max_s):
        fig = plt.figure(figsize=(10*size,10*size)) 
        figs.append(fig)       
    for key in names.keys():
        for exp in range(len(names[key])):
            n_s,n_l,n_f = FAM[key][exp].shape
            for s in range(n_s):
                MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor)


ts = 3    
paths_WT = [#'FX WT females EH Nen (1) 07.12.15',
            'FX WT females EH (1) 16.02.15',
            'FX WT females EH (B) 05.05.16',
            'FX WT females EH (B) - powtorzenie 2. - 18.05.16',
            'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16',
            #'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16',
            ]
paths_KO =  ['FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16',
            #'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16',
            #'FX KO F (3) EH (1) 15.05.14',
            'FX KO females EH 2 29.09.14',
            'FX KO females EH (3) 11.02.15'
            ]
exp_paths = paths_WT+paths_KO
paths_RD = ['Random_test%s'%(j+1) for j in range(len(exp_paths))]    
colors = {
          "RD":'grey',
          "WT":'green',
          "KO":'red',
          }
names = {}
#names["RD"] = paths_RD
names["WT"] = paths_WT
names["KO"] = paths_KO            
            
#paths_BALB_VPA = [
#                    'BALB VPA EH (1) 05.05.14', # ***
#                    'BALB VPA EH (2) 09.05.14', # 5
#                    'BALB VPA EH (4) 16.06.14 - 19.06.14',
#                    'BALB VPA EH Biol (1) 23.03.15'
#                  ]
#paths_BALB_CTRL = [
#                   'BALB CTRL EH (3) 12.06.14 - 15.06.14',
#                   'BALB CTRL EH (4) 09.09.14',
#                   'BALB CTRL EH (5) 15.09.14'
#                   
#                   ]
#paths_BALB_NaCl = [
#                  'BALB NaCl (1) Rep (1) 09.09.15',
#                  'BALB NaCl (2) Rep (1) 02.10.15'
#                   ]
#            
#names = {}
#names["BALB_VPA"] = paths_BALB_VPA
#names["BALB_CTRL"] = paths_BALB_CTRL
#names["BALB_NaCl"] = paths_BALB_NaCl
#colors = {
#          "BALB_VPA":'blue',
#          "BALB_CTRL":'green',
#          "BALB_NaCl":'red',
#          }            
#createRandomExpepiments(exp_paths)
#Interpersec(names,ts=200)
#preprocessData(names,window = 12,ts=3)
IPP = InteractionsPerPair(names)
scalefactor = np.max([np.max(i) for i in IPP["KO"]]+[np.max(i) for i in IPP["WT"]])
#scalefactor = np.max([np.max(i) for i in IPP["BALB_VPA"]]+[np.max(i) for i in IPP["BALB_CTRL"]]+[np.max(i) for i in IPP["BALB_NaCl"]])
FPP = FollowingPerPair(names)
APP = AvoidingPerPair(names)
FAM = FollowingAvoidingMatrix(names)
result = []
pair_inc = []
pair_long = []
def longest_sequence(FAM,n_s=6):
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

def follsactive(FAM,n_s=6):
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
stats = {}
stats["KO"] = {}
stats["KO"]["SLD"] = []
stats["KO"]["NFD"] = []
stats["KO"]["SLA"] = []
stats["KO"]["NFA"] = []
for i in range(3):
    _FAM=FAM["KO"][i]
    LSD, LSA = longest_sequence(_FAM,n_s=6)
    FSD, FSA = follsactive(_FAM,n_s=6)
    stats["KO"]["SLA"]+=LSA
    stats["KO"]["NFA"]+=FSA
    stats["KO"]["NFD"]+=FSD
    stats["KO"]["SLD"]+=LSD
    print np.mean(LSA), np.mean(FSA), FAM["KO"][i].shape[2]
    #plt.hist(LSD)
    #plt.hist(FSD)
    #plt.show()
print '#####################WT#######################'
stats["WT"] = {}
stats["WT"]["SLD"] = []
stats["WT"]["NFD"] = []
stats["WT"]["SLA"] = []
stats["WT"]["NFA"] = []
for i in range(4):
    _FAM=FAM["WT"][i]
    LSD, LSA = longest_sequence(_FAM,n_s=6)
    FSD, FSA = follsactive(_FAM,n_s=6)
    stats["WT"]["SLA"]+=LSA
    stats["WT"]["NFA"]+=FSA
    stats["WT"]["NFD"]+=FSD
    stats["WT"]["SLD"]+=LSD
    print np.mean(LSA), np.mean(FSA), FAM["WT"][i].shape[2]
#    plt.hist(LSA)
#    plt.hist(FSA)
#    plt.show()
barplot(stats,names,["SLA","SLD",], colors, name="AverageLengthofSec",ylab = "Average length of sequence per pair")
barplot(stats,names,["NFA","NFD" ], colors, name="AverageNumberofFols", ylab="Average number of followings per pair")
#createRasterPlots(FAM,IPP,names,scalefactor)           
#createRasterPlotsSUM(FAM,IPP,names,scalefactor)
#CreateRelationGraphs(FAM,IPP,names,scalefactor/50)
statsIPP = plotphist(IPP,names,colors,to_file = True,directory = 'Interactions',vrange = [0,120], prange = [0,0.11])
statsFPP = plotphist(FPP,names,colors,to_file = True,directory = 'Followings',vrange = [0,120], prange = [0,0.11])
statsAPP = plotphist(APP,names,colors,to_file = True,directory = 'Avoidings',vrange = [0,120], prange = [0,0.11])
print statsIPP
barplot(statsIPP,names,["Interactions"], colors, name="InteractionPerPairBarplot",ylab="Average number of interactions per pair")
barplot(statsFPP,names,["Followings"], colors,name="FollowingsPerPairBarplot",ylab="Average number of followings per pair")
barplot(statsAPP,names,["Avoidings"], colors,name="AvoidingsPerPairBarplot", ylab="Average number of avoidings per pair")

#print  st.mannwhitneyu(statsIPP["KO"]["mean"], statsIPP["WT"]["mean"], use_continuity=True)
#print  st.mannwhitneyu(statsFPP["KO"]["mean"], statsFPP["WT"]["mean"], use_continuity=True)
#print  st.mannwhitneyu(statsAPP["KO"]["mean"], statsAPP["WT"]["mean"], use_continuity=True)
#print  st.mannwhitneyu(stats["KO"]["median"], stats["WT"]["median"], use_continuity=True)
#FAP = FollowingAvoidingMatrix(names)
#stats = plotphist(FAP,names,colors,to_file = True,directory = 'FA',vrange = [-1,1], prange = [0,1])
#print  st.mannwhitneyu(stats["KO"]["mean"], stats["WT"]["mean"], use_continuity=True)
#print  st.mannwhitneyu(stats["KO"]["median"], stats["WT"]["median"], use_continuity=True)
#m = np.round(FAP["KO"][0][0],3)
#print m
