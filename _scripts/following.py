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


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

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
        self.interactions = np.zeros((len(self.mice),len(self.mice),len(self.phases)))
        self.f_sum = np.zeros((len(self.phases),self.lm))
        if 'following.npy' in list(os.walk(os.path.join('..','Results',self.path)))[0][2] and not force:
            self.f = np.load(os.path.join('..','Results/'+self.path+'/')+'following.npy')
            return self.f
        else:
            for s in range(len(self.phases)):
                ts, te = self.phases[s]
                print 'Phase %s. from %sh, to %sh'%(s+1,np.round(ts/3600.,2), np.round(te/3600.,2))
                self.f[:,:,s], self.interactions[:,:,s] = self.fvalue_matrix(ts, te,s)
            np.save(os.path.join('..','Results/'+self.path+'/')+'following.npy',self.f)
            return self.f
     
    
    def fvalue_matrix(self,ts, te,s):
        """Calculates fvalue matrix for the given period from ts to te
        """
        fmatrix = np.zeros((self.lm,self.lm))
        imatrix = np.zeros((self.lm,self.lm))
        for ii in range(len(self.mice)):
            for jj in range(len(self.mice)):
                if ii < jj:
                    fmatrix[ii,jj],imatrix[ii,jj] = self.fvalue((ii,jj),ts, te,s)
                    fmatrix[jj,ii],imatrix[jj,ii] = self.fvalue((jj,ii),ts, te,s)
                    
        return fmatrix,imatrix  
    
    def fvalue(self, (m1,m2),t1, t2,s ):
        follow_stat, detected_idx = self.fpatterns((m1,m2),t1, t2)
        #print follow_stat
        #print follow_stat
        '''t = np.arange(-3,3,1.0/self.fs)
        for i in range(len(detected_idx[0])):
            s = detected_idx[1][i]
            plt.plot(t,self.sd[s-3*self.fs:s+3*self.fs,m1]-0.05,'ro')
            plt.plot(t,self.sd[s-3*self.fs:s+3*self.fs,m2]+0.05,'bo')
            plt.axis([-3.1,3.1,-0.5,9.5])
            plt.show()
        #print follow_stat'''
        print '#################%s%s##################'%(m1,m2)
        gpos = 0
        gneg = 0
        gp = []
        pp = [0]
        for key in follow_stat.keys():
            pos,neg,p = follow_stat[key]
            gpos+=pos
            gneg+=neg
            gp.append(p)
            pf,pa =round(prob(int(neg), 1-p, int(pos+neg)),3), round(prob(int(pos), p, int(pos+neg)),3)
#            if pf <0.05:
#                pp.append(2.5-pf*50)
#                print '\x1b[6;30;41m', key, follow_stat[key], pf,pa,'\x1b[0m'
#            elif pa <0.05:
#                pp.append(-2.5+pa*50)
#                print '\x1b[6;30;44m', key, follow_stat[key], pf,pa,'\x1b[0m'
#            else:
#                print key, follow_stat[key], pf,pa
        gpmean,gpstd = round(np.mean(gp),3), round(np.std(gp),3)
        pf,pa =round(prob(int(gneg), 1-gpmean, int(gpos+gneg)),3), round(prob(int(gpos), gpmean, int(gpos+gneg)),3)
        if pf <0.05:
            pp = 2.5-pf*50
            print '\x1b[6;30;41m',"global", [gpos,gneg,gpmean], pf,pa,'\x1b[0m'
        elif pa <0.05:
            pp = -2.5+pa*50
            print '\x1b[6;30;44m', "global", [gpos,gneg,gpmean], pf,pa,'\x1b[0m'
        else:
            pp = 0
            print "global", [gpos,gneg,gpmean], pf,pa
        print gpstd
        return np.sum(pp), gpos+gneg

    
    def fpatterns(self, (m1,m2),t1, t2):
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
                    

                    same_start = sd[m1_idx[i-2],m1] in period2 #POPRAWIC!!!!!!!
                    first_m1 = sd[s,m1] not in period3
                    op_idx = (2*sd[m1_idx[i-2],m1]-sd[s,m1]-1)%8+1
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
                            self.fols[np.ceil((period1.index(sd[s,m1]))/self.fs)]+=1
                            follow_stat[str(start_st)+str(end_st)][0] += 1 
                            #print p, sd[m1_idx[i-2],m1],[index]
                            detected_idx[0].append(s)
                        elif go_oposite:
                            #print np.ceil((period1.index(op_idx))/self.fs)
                            self.ops[np.ceil((period1.index(op_idx))/self.fs)]+=1
                            follow_stat[str(start_st)+str(end_st)][1] += 1
                            detected_idx[1].append(s)
                except IndexError:
                    print 'Err'
                    continue
        return follow_stat, detected_idx
        
    def plot_graphs(self, nr=''):
        for s in range(len(self.phases)):
            self.plot_graph(s, nr=nr)
            
            
    def plot_fsums(self,base_name = 'following_sum',nr=''):
        plt.figure(figsize=(15,15))
        plt.suptitle(self.path, fontsize=14, fontweight='bold')
        for s in range(len(self.phases)):
            plt.subplot(len(self.phases),1,s+1)
            plt.title('Session %s'%(s+1))
            self.plot_fsum(s)
        plt.savefig(self.directory+'/'+'following_sum_%s_%s%s.png'%(s,self.treshold,nr))
        plt.show()
        
        plt.figure(figsize=(15,15))
        plt.suptitle(self.path, fontsize=14, fontweight='bold')
        for s in range(len(self.phases)):
            plt.subplot(len(self.phases),1,s+1)
            plt.title('Session %s'%(s+1))
            self.plot_fsum2(s)
        plt.savefig(self.directory+'/'+'all_sum_%s_%s%s.png'%(s,self.treshold,nr))
        plt.show()
        
    def plot_fsum(self,s):
        for i in range(self.lm):
            data = []
            mouse_nr = []
            for j in range(self.lm):
                if i!=j:
                    if self.f[i,j,s] > 1:
                        data.append(self.f[i,j,s]-1)
                        mouse_nr.append(j)
            self.f_sum[s,i] = np.sum(data)
        nr = np.arange(1,self.lm+1)
        #print nr, f_sum
        plt.stem(nr, self.f_sum[s,:])
        plt.axis([0,self.lm+1.5,0,5])
        plt.xlabel('Mouse number')
        plt.ylabel('Sum of f parameter')
        plt.xticks(nr)
        
    def plot_fsum2(self,s):
        for i in range(self.lm):
            data = []
            mouse_nr = []
            for j in range(self.lm):
                if i!=j:
                    if self.f[i,j,s] > 0:
                        data.append(self.f[i,j,s]-1)
                        mouse_nr.append(j)
            self.f_sum[s,i] = np.sum(data)
        nr = np.arange(1,self.lm+1)
        #print nr, f_sum
        plt.stem(nr, self.f_sum[s,:])
        plt.axis([0,self.lm+1.5,-5,5])
        plt.xlabel('Mouse number')
        plt.ylabel('Sum of f parameter')
        plt.xticks(nr)
        
    
    def plot_graph(self,s, base_name = 'following_graph', nr=''):
        pairs = []       
        for i in range(self.lm):
            for j in range(self.lm):
                if i!=j:
                    pairs.append([self.f[i][j][s],i+1,j+1])
        pairs.sort(key=lambda x: -x[0])
        conn = [(np.round((x[0]),2),x[1],x[2]) for x in pairs]
        G = nx.MultiDiGraph(multiedges=True, sparse=True)
        for i in range(len(conn)):
            G.add_edges_from([(conn[i][1],conn[i][2])], weight=conn[i][0])
        edge_colors = [conn[i][0] for i in range(len(conn))]
        size = 10
        pos=nx.circular_layout(G)
        node_labels = {node:node for node in G.nodes()}     
        fig = plt.figure(figsize=(10*size,10*size)) 
        ax = fig.add_subplot(111, aspect='equal')
        fig.suptitle(self.path, fontsize=14*size, fontweight='bold')
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
                p = patches.FancyArrowPatch(pos[c[2]],pos[c[1]],connectionstyle='arc3, rad=-0.3',arrowstyle="simple",shrinkA=12*size, shrinkB=12*size,mutation_scale=20*size*abs(c[0]), color = cmap(c[0]+0.5),zorder=1,alpha=0.5)
                ax.add_patch(p)
        plt.savefig(self.directory+'/'+base_name+'_%s_ts%s%s.png'%(s,self.treshold,nr))
        plt.close(fig)
        
    def plotfhist(self,ts,to_file = True,nr = ''):
        plt.suptitle(self.path, fontsize=14, fontweight='bold')
        x = self.f[:,:,:].flatten()
        print 'Neutral%s'%len(x[x==0])
        plt.hist(x[x!=0])#,bins=np.arange(-5,5,0.1))
        plt.savefig(self.directory+'/fhist_%s%s.png'%(ts,nr))
        plt.show()
    
    def plotihist(data,names,colors,to_file = True,directory = 'Interactions'):
        if not os.path.exists('../Results/'+directory):
            os.makedirs('../Results/'+directory)
        stats = {}
        for key in data.keys():
            stats[key] = {}
            stats[key]["mean"] = []
            stats[key]["median"] = []
            for i in range(len(data[key])):
                x = data[key][i].flatten()
                stats[key]["mean"].append(np.mean(x))
                stats[key]["median"].append(np.median(x))
                plt.suptitle(names[key][i], fontsize=14, fontweight='bold')
                plt.hist(x, bins = np.arange(0,120,5),normed=True)#,bins=np.arange(-5,5,0.1))
                plt.axvline(np.mean(x), color=colors[key], linestyle='solid', linewidth=2)
                plt.axvline(np.median(x), color=colors[key], linestyle='dashed', linewidth=2)
                plt.ylim((0,0.11))
                plt.savefig('../Results/'+directory+'/'+names[key][i]+'.png')
                plt.show()
        for key in data.keys():
            for i in range(len(data[key])):
                x = data[key][i].flatten()
                plt.subplot(2,1,1)
                plt.axvline(np.mean(x), color=colors[key], linestyle='solid', linewidth=2)
                plt.subplot(2,1,2)
                plt.axvline(np.median(x), color=colors[key], linestyle='dashed', linewidth=2)
        plt.show()
        return None
            
def createRandomExpepiments(paths,core='Random_test' ):
    statistics_rd = {}
    for j in range(len(paths)):
        rdE = Experiment('Random_test%s'%j,exp_name=core)
        rdname = 'Random_test%s.pkl'%j
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
        rdE.sd = rd
        rdE.ehs.statistics = statistics_rd
        with open('../PreprocessedData/'+rdname, 'wb') as output:
           pickle.dump(rdE,output,pickle.HIGHEST_PROTOCOL)
           rd_paths.append(rdname)
           del rdE
    return rd_paths

def calculatestats():
    f_arr = []
    for i in range(len(fexps)):
        f_arr+=fexps[i]
        el = np.array(fexps[i])
        x = el[el!=0]
        fstats.append([len(el[el==0])]+[round(f(x),5) for f in [np.mean, np.std, st.skew, st.kurtosis]])
        istats.append([round(f(iexps[i]),5) for f in [np.mean, np.std, st.skew, st.kurtosis]])
        print fstats[-1]
           
def analyzeData(paths,ts=3,rd = False):
    fols = np.zeros(ts)
    ops  = np.zeros(ts)
    fexps = []
    iexps = []
    for path in paths:
        print path
        if rd:
            with open('../PreprocessedData/'+path+'.pkl', "rb") as input_file:
                E = pickle.load(input_file)
        else:
            E = Experiment(path)
        E.calculate_fvalue(window = 24, treshold =ts, force=True, fols=fols,ops=ops)
        fols = E.fols
        ops = E.ops
        E.plot_fsums()
        #E.plot_graphs()
        E.plotfhist(ts)
        #E.plotihist(ts)
        fexps.append(E.f)
        print E.f[:,:,0]
        iexps.append(E.interactions)
        break
    return fexps, iexps
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
paths_RD = ['Random_test%s'%j for j in range(len(exp_paths))]    
data = {}
colors = {
          "RD":'grey',
          "WT":'green',
          "KO":'red',
          }
names = {}
names["RD"] = paths_RD
names["WT"] = paths_WT
names["KO"] = paths_KO
            
#createRandomExpepiments(exp_paths)
#fRD,data["RD"] = analyzeData(paths_RD,rd = True)
#fstatsWT,data["WT"] = analyzeData(paths_WT)
fstatsKO,data["KO"] = analyzeData(paths_KO)
plotihist(data,names,colors,to_file = True,directory = 'Interactions')


