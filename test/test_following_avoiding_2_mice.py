from __future__ import division, print_function

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import make_circular
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from EcoHAB import interactions
from EcoHAB import plotfunctions



a_dirs  = ["/home/jszmek/EcoHAB/test/Circular_data_following"]
masks = {"/home/jszmek/EcoHAB/test/Circular_data_following":[]}
phases = {"/home/jszmek/EcoHAB/test/Circular_data_following":['ALL']}
antenna_pos = {"/home/jszmek/EcoHAB/test/Circular_data_following":{'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8}
}

data_phases = ['FOLLOWING','AVOIDING']
ending = ['','']
m2_function = [make_circular.new_state_clockwise,make_circular.new_state_counter_clockwise]
phase_duration = 12*60*60
exp_duration = [4,4]

ts = 3   
window=12
IPP = {}
FAM = {}
sections = {}
directories = {}
endings = {}

for a_dir in a_dirs:
    
    make_circular.make_data(a_dir,data_phases,ending,m2_function,phase_duration,exp_duration,state_m1=1,state_m2=1,previous_m1=8,previous_m2=8,mouse1 = "0065-0136661698",mouse2 = "0065-0136656570",followings=[2,10])
    
    IPP[a_dir] = []
    FAM[a_dir] = []
    sections[a_dir] = []
    directories[a_dir] = []
    endings[a_dir] = []

    if masks[a_dir] == []:
        for phase in phases[a_dir]:
            E = interactions.Experiment(a_dir,_ant_pos=antenna_pos[a_dir],which_phase=phase,how_many_appearances=10)#,mask=m1)
            mouse_positions = E.sd
            fname = 'Interactions_' + E.fname_ending
            E.calculate_fvalue(window=window, treshold=ts, force=True)
            IPP[a_dir].append(E.InteractionsPerPair(0,2))
            FAM[a_dir].append(E.FollowingAvoidingMatrix())
            sections[a_dir].append(E.cf.sections())
            directories[a_dir].append(E.directory)
            endings[a_dir].append(E.fname_ending)
    else:
        for mask in masks[a_dir]:
            E = interactions.Experiment(a_dir,_ant_pos=antenna_pos[a_dir],mask=mask,how_many_appearances=10)
            mouse_positions = E.sd
            fname = 'Interactions_' + E.fname_ending

            E.calculate_fvalue(window=window, treshold=ts, force=True)
            IPP[a_dir].append(E.InteractionsPerPair(0,2))
            FAM[a_dir].append(E.FollowingAvoidingMatrix())
            sections[a_dir].append(E.cf.sections())
            directories[a_dir].append(E.directory)
            endings[a_dir].append(E.fname_ending)
                        

scalefactor = 0

for a_dir in a_dirs:
    maxi = 0

    for ipp in IPP[a_dir]:
        
        if np.max(ipp) > maxi:
            maxi = np.max(ipp)
    scalefactor += maxi
print(directories)
for a_dir in a_dirs:
    for i, ipp in enumerate(IPP[a_dir]):
        interactions.oneRasterPlot(directories[a_dir][i],FAM[a_dir][i],ipp,sections[a_dir][i],'Interactions_ts_'+str(ts)+'_s_'+endings[a_dir][i],scalefactor)

        for k,l in enumerate(FAM[a_dir][i]):
            plotfunctions.plot_graph(FAM[a_dir][i],k,sections[a_dir][i],directories[a_dir][i])
            

