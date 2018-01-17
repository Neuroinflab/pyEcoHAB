from __future__ import division, print_function
import interactions
import sys
import plotfunctions
import numpy as np
import matplotlib.pyplot as plt

antenna_pos = {"/home/jszmek/EcoHAB_data_November/long_experiment_KO":{'1':1,'2':5,'3':3,'4':6,'5':4,'6':2,'7':7,'8':8},
               "/home/jszmek/EcoHAB_data_November/long_experiment_WT":{'1':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8}
}
# if len(sys.argv) < 2:
#     sys.exit("No data directory given")
a_dirs  = ["/home/jszmek/EcoHAB_data_November/long_experiment_KO","/home/jszmek/EcoHAB_data_November/long_experiment_WT"]
#[0,1507905985.0]
masks = {"/home/jszmek/EcoHAB_data_November/long_experiment_KO":[],#[0,1508410232],[1508410232.,1508740930.0]],
         "/home/jszmek/EcoHAB_data_November/long_experiment_WT":[]}
phases = {"/home/jszmek/EcoHAB_data_November/long_experiment_KO":["BEGINNING","MIDDLE"],
         "/home/jszmek/EcoHAB_data_November/long_experiment_WT":["BEGINNING","MIDDLE"]}
ts = 3   
window=12
IPP = {}
FAM = {}
sections = {}
directories = {}
endings = {}
for a_dir in a_dirs:
    IPP[a_dir] = []
    FAM[a_dir] = []
    sections[a_dir] = []
    directories[a_dir] = []
    endings[a_dir] = []

    if masks[a_dir] == []:
        for phase in phases[a_dir]:
            E = interactions.Experiment(a_dir,_ant_pos=antenna_pos[a_dir],which_phase=phase)#,mask=m1)
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
            E = interactions.Experiment(a_dir,_ant_pos=antenna_pos[a_dir],mask=mask)
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

for a_dir in a_dirs:
    for i, ipp in enumerate(IPP[a_dir]):
        interactions.oneRasterPlot(directories[a_dir][i],FAM[a_dir][i],ipp,sections[a_dir][i],'Interactions_ts_'+str(ts)+'_s_'+endings[a_dir][i],scalefactor)
        for k,l in enumerate(FAM[a_dir][i]):
            plotfunctions.plot_graph(FAM[a_dir][i],k,sections[a_dir][i],directories[a_dir][i])
            

