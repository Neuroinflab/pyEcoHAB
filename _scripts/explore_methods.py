#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:30:26 2016

@author: jmaka
"""

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

def prepare_data(ehs, mice, times, margin=12*3600.):
    """Prepare masked data."""
    t1, t2 = times
    ehs.mask_data(t1 - margin, t2)
    data = {}
    for mm in mice:
        ads = ehs.getaddresses(mm)
        sts = ehs.getstarttimes(mm)
        ens = ehs.getendtimes(mm)
        data_mm = []
        for ad, st, en in zip(ads, sts, ens):
            if en > t1:
                data_mm.append((ad, max(st, t1), min(en, t2)))
        data[mm] = data_mm
    return data


def conv2dicrete(ehs, mice, times, margin=12*3600.,t_max=43200):
    t = np.arange(0,t_max,0.1)
    v = np.zeros((len(t),len(mice)))
    t1, t2 = times
    ehs.mask_data(t1 - margin, t2)
    data = {}
    for mm in mice:
        ads = ehs.getaddresses(mm)
        sts = ehs.getstarttimes(mm)
        ens = ehs.getendtimes(mm)
        data_mm = []
        for ad, st, en in zip(ads, sts, ens):
            if en > t1:
                data_mm.append((ad, max(st, t1), min(en, t2)))
        data[mm] = data_mm
    i =0
    for mm in mice:
        for rec in data[mm]:
            s = int((rec[1]-t1)*10)
            e = int((rec[2]-t1)*10)
            v[s:e,i] = rec[0]
        i+=1
    return t,v

def following(signal_data,(m1,m2),t1, t2 ,fs,treshold = 2, min_overlap = 0.0):
    sd = signal_data
    detected_idx = [[],[]]
    follow_stat = np.zeros((8,3))
    m1_idx = np.where(((np.roll(sd[:,m1], 1) - sd[:,m1]) != 0))[0]

    print 'ruchliwosc', len(m1_idx)
    for i in range(2,len(m1_idx)):
        if m1_idx[i] > t1*fs and m1_idx[i]<t2*fs:
            s = m1_idx[i]
            #print s
            e =s+treshold*fs
            try:
                period1 = list(sd[s:e,m2])
                period2 = list(sd[s-2*fs:e-2*fs,m2])
                period3 = list(sd[s-int(0.1*fs):s,m2])
                # define conditions
                unknown_state = sd[s,m1]==0
                unknown_previous_states = (sd[m1_idx[i-1],m1] ==0 and (sd[m1_idx[i-2],m1] ==0)) or (sd[m1_idx[i-1],m1] !=0 and (sd[m1_idx[i-2],m1] ==0))
                in_pipe =  sd[s,m1]%2==1
                same_start = sd[m1_idx[i-2],m1] in period2 #POPRAWIC!!!!!!!
                first_m1 = sd[s,m1] not in period3
                followed = period1.count(sd[s,m1])/fs>min_overlap
                go_oposite = (2*sd[m1_idx[i-2],m1]-sd[s,m1]-1)%8+1 in period1
        
                if unknown_state or unknown_previous_states or in_pipe:
                    continue
                elif same_start and first_m1:
                    if followed:
                        follow_stat[sd[m1_idx[i-2],m1]-1,0] += 1 
                        index = 0 if sd[m1_idx[i-2],m1] <sd[s,m1] else 1
                        follow_stat[sd[m1_idx[i-2],m1]-1,2] += 1./ehs.statistics[list(ehd.mice)[m2]]["preference"][sd[m1_idx[i-2],m1]][index]
                        detected_idx[0].append(s)
                    elif go_oposite:
                        follow_stat[sd[m1_idx[i-2],m1]-1,1] += 1
                        detected_idx[1].append(s)
            except IndexError:
                break
    print follow_stat
    #print follow_stat
    #print np.sum(follow_stat[:,0])
    if np.sum(follow_stat[:,0:2])>1:
        #print m1, m2
        for state in range(8):
            follow_stat[state,2]/=np.sum(follow_stat[state,:-1])
            follow_stat[state,2]-=1
            #print follow_stat[state,2],
        #print '####' , np.nansum(follow_stat[:,2])
        return np.nanmean(follow_stat[:,2])*10000/len(m1_idx),detected_idx
    else:
        return 10, detected_idx


def random_data(signal_data):
    sd = signal_data
    rd = np.zeros(sd.shape)
    for m in range(sd.shape[1]):
        print m
        mm = list(ehd.mice)[m]
        st = int(sd[0,m])
        ts = 0
        t = np.int(np.random.choice(ehs.statistics[mm]["state_time"][st],size=1)*ehs.fs)
        te = ts+t
        while te<sd.shape[0]:
            rd[ts:te,m] = st
            r = np.random.random()
            if r< ehs.statistics[mm]["preference"][st][1]:
                st = (st+1)%9
            else:
                st = (st-1)%9
            if st ==0:
                st = 1
            ts = te
            t = np.int(np.random.choice(ehs.statistics[mm]["state_time"][st],size=1)*ehs.fs)
            te = ts+t 
    return rd    

path = 'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16'
ehd = EcoHab.EcoHabData(os.path.join('..','RawData',path), _ant_pos=antenna_positions[path])
ehs = EcoHab.EcoHabSessions9states(ehd,shortest_session_threshold=0)
sd =  ehs.signal_data
print sd.shape
cf = ExperimentConfigFile(os.path.join('..','RawData',path))
tstart, tend = cf.gettime('ALL')
diff = ehs.t_start_exp-tstart
mice = list(ehd.mice)
print ehs.statistics[mice[0]]["preference"]