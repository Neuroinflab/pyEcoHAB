# -*- coding: utf-8 -*-
from __future__ import division, print_function
import EcoHab
import numpy as np
import os
import utility_functions as utils
from write_to_file import save_data_cvs
from data_info import *
from ExperimentConfigFile import ExperimentConfigFile

### How much time mice spend in the 'social' compartment
    
def get_visits(mouse_address,mouse_durations,stim_address):
    result = []
    for address, t in  zip(mouse_address, mouse_durations):
        if address == stim_address:
            result.append(t)
    return result

def initialize_data(cf,chamber_dict):
    phases = find_phases(cf)
    visits = {}
    times = {}
    tt = phase_dictionary(cf)
    for key in chamber_dict:
        visits[key] = phase_dictionary(cf)
        times[key] = phase_dictionary(cf)
 
    return visits, times, tt

def find_phases(cf):
    return utils.filter_dark_light(cf.sections())

def phase_dictionary(cf):
    return dict([(phase, []) for phase in find_phases(cf)])

def phase_dictionary_phases(phases):
    return dict([(phase, []) for phase in phases])

def mouse_dict(mice):
    return dict([(mm, []) for mm in mice])

def loop_through_mice(v, t, address, ehs):
    for mouse in ehs.mice:
        mouse_address = ehs.getaddresses(mouse)
        mouse_durations = ehs.getdurations(mouse)
        v[mouse].append(mouse_address.count(address))
        t[mouse].append(get_visits(mouse_address, mouse_durations, address))

def loop_through(visits,times,tt,cf,ehs,chamber_dict,binsize):
    phases = find_phases(cf)
    i = 0

    for key in chamber_dict:  
        for phase in phases:
            visits[key][phase] = mouse_dict(ehs.mice)
            times[key][phase] = mouse_dict(ehs.mice)
    

    for phase in phases:
        time_start, phase_end = cf.gettime(phase)
        time = time_start
            
        while time < phase_end:
            ehs.unmask_data()
            ehs.mask_data(time, time + binsize)
            for key in chamber_dict:  
                loop_through_mice(visits[key][phase], times[key][phase], chamber_dict[key], ehs)

            tt[phase].append((time-time_start)*binsize/3600.)
            time += binsize
                
def results_to_array(visits):
    for key in visits:
        for phase in visits[key]:
            for mouse in visits[key][phase]:
                visits[key][phase][mouse] = np.array(visits[key][phase][mouse])
                 

def get_time_spent_in_each_chamber(ehs, cf, chamber_dict, binsize):
    
    visits, times, tt = initialize_data(cf,chamber_dict)
    loop_through(visits,times,tt,cf,ehs,chamber_dict,binsize)   
    results_to_array(visits)
    results_to_array(times)
    out = {}
    for key in visits:
        out[key] = [visits[key], times[key]]
    out['time'] = tt
    out['mice'] = ehs.mice
    out['phases'] = find_phases(cf)
    return out

def sum_data(data):
    phases = data.pop('phases')
    mice = data.pop('mice')
    time = data.pop('time')
    st = {}
    new_data = {'phases':phases, 'mice':mice, 'time':time}
    for key in data:
        st[key] = phase_dictionary_phases(phases)
        new_data[key] = []
        new_data[key].append(data[key][0])
        new_data[key].append(st[key])
        for phase in phases:
            st[key][phase] = mouse_dict(mice)
            for mouse in mice:
                st[key][phase][mouse] = np.zeros((len(time[phase])))
                for k, t in enumerate(time[phase]):
                    st[key][phase][mouse][k] = sum(data[key][1][phase][mouse][k])
    return new_data
        
def get_time_spent(ehs,cf,key,binsize=12*3600):
  
    tstart, tend = cf.gettime('ALL')
    phases = find_phases(cf)
    nsoc = phase_dictionary(cf)
    soc = phase_dictionary(cf)
    t_nsoc = phase_dictionary(cf)
    t_soc = phase_dictionary(cf)
    tt =  phase_dictionary(cf)
    
    for phase in phases:
        nsoc[phase] = dict([(mm, []) for mm in ehs.mice])
        soc[phase] = dict([(mm, []) for mm in ehs.mice])
        t_nsoc[phase] = dict([(mm, []) for mm in ehs.mice])
        t_soc[phase] = dict([(mm, []) for mm in ehs.mice])
        time_start,phase_end = cf.gettime(phase)
        time = time_start
      
        while time < phase_end:
            ehs.unmask_data()
            ehs.mask_data(time, time + binsize)
            for mouse in ehs.mice:
                mouse_address = ehs.getaddresses(mouse)
                mouse_durations = ehs.getdurations(mouse)
 
                nsoc[phase][mouse].append(mouse_address.count(smells[key]['nsoc']))
                soc[phase][mouse].append(mouse_address.count(smells[key]['soc']))
                t_nsoc[phase][mouse].append(get_visits(mouse_address,mouse_durations, smells[key]['nsoc']))
                t_soc[phase][mouse].append(get_visits(mouse_address,mouse_durations, smells[key][ 'soc']))
                
            tt[phase].append((time-time_start)*binsize/3600.)
            time += binsize
 
    for phase in phases:
        for mouse in ehs.mice:
            nsoc[phase][mouse] = np.array(nsoc[phase][mouse])
            soc[phase][mouse] = np.array(soc[phase][mouse])
            t_nsoc[phase][mouse] = np.array(t_nsoc[phase][mouse])
            t_soc[phase][mouse] = np.array(t_soc[phase][mouse])

    return {'nsoc':[nsoc,t_nsoc], 'soc':[soc,t_soc],'time':tt,'mice':ehs.mice,'phases':phases}

def sum_times(data):
    phases = data['phases']
    mice = data['mice']
    sum_t_soc = phase_dictionary_phases(phases)
    sum_t_nsoc = phase_dictionary_phases(phases)

    for phase in phases:

        sum_t_soc[phase] = mouse_dict(mice)
        sum_t_nsoc[phase] = mouse_dict(mice)
          
        for mouse in mice:

            sum_t_soc[phase][mouse] = np.zeros((len(data['time'][phase])))
            sum_t_nsoc[phase][mouse] = np.zeros((len(data['time'][phase])))
            
            for k,time in enumerate(data['time'][phase]):
                sum_t_soc[phase][mouse][k] = sum(data['soc'][1][phase][mouse][k])
                sum_t_nsoc[phase][mouse][k] = sum(data['nsoc'][1][phase][mouse][k])
 
    return {'phases':phases,'mice':mice,'soc':[data['soc'][0],sum_t_soc],'nsoc':[data['nsoc'][0],sum_t_nsoc],'time':data['time']}
                    
                
def calculate_approach_to_social(data):
    #possibly calculate a square root of sum of squared errors -- if independent errors (Taylor, 1997) probably not
    ats = data[0]/data[1]/(data[2]/data[3])
    d_ats = np.zeros(data[0].shape)
    for i,dat in enumerate(data[4:]):
        d_ats += abs(dat/data[i])
    return ats,d_ats*ats
        
def approach_to_social(data,fname=None,path=None):

    times_nsoc = data['nsoc'][1]
    times_soc = data['soc'][1]
    phases = data['phases']
    last_phase_before_stimulus = [phase for phase in phases if 'EMPTY' in phase and 'dark' in phase or 'EMPTY' in phase and 'DARK' in phase][-1]
    first_phase_after_stimulus = [phase for phase in phases if 'SNIFF' in phase and 'dark' in phase or 'SNIFF' in phase and 'DARK' in phase][0]

    mice = data['mice']
    mice_no = len(mice)
    shape = (len(times_soc[first_phase_after_stimulus][mice[0]]),len(mice))
    T_S, T_NS, t_s, t_ns = np.zeros(shape), np.zeros(shape), np.zeros(shape), np.zeros(shape)
    
    for i,mouse in enumerate(mice):
        T_S [:,i] = times_soc[first_phase_after_stimulus][mouse]
        T_NS[:,i] = times_nsoc[first_phase_after_stimulus][mouse]
        t_s[:,i] = times_soc[last_phase_before_stimulus][mouse]
        t_ns[:,i] = times_nsoc[last_phase_before_stimulus][mouse]
        
    mean_T_S, mean_T_NS, mean_t_s, mean_t_ns = T_S.mean(axis=1), T_NS.mean(axis=1), t_s.mean(axis=1), t_ns.mean(axis=1)
    std_T_S, std_T_NS, std_t_s, std_t_ns = (T_S.var(axis=1)/mice_no)**.5, (T_NS.var(axis=1)/mice_no)**.5, (t_s.var(axis=1)/mice_no)**.5, (t_ns.var(axis=1)/mice_no)**.5
    ats,d_ats = calculate_approach_to_social([mean_T_S, mean_T_NS, mean_t_s, mean_t_ns,std_T_S, std_T_NS, std_t_s, std_t_ns])
    l = len(data['time'][last_phase_before_stimulus])
    try:
        binsize = data['time'][last_phase_before_stimulus][1]-data['time'][last_phase_before_stimulus][0]
    except IndexError:
        binsize = 12*3600
        
    time = [binsize*i/3600 for i in range(l)]

    if fname:
        if path:
            if not os.path.exists(path):
                os.makedirs(path)
            fname = os.path.join(path,fname)
        np.savetxt(fname, np.array((time,ats,d_ats)).T,header='time; AtS; dAtS',delimiter=';',comments="")
        
    return time, ats, d_ats
homepath = os.path.expanduser("~/")
threshold = 3

if __name__ == '__main__':
    datapaths = ['EcoHAB_data_November/FxKO_Bodziec_1st_hour',
                'EcoHAB_data_November/FxWT_bodziec_1st_hour_13.07.18',
    ]
    for new_path in datapaths:
        path = os.path.join(homepath, new_path)
        prefix = utils.make_prefix(path)
        print(prefix)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
        else:
            remove_mouse = None
        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 10
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        ehs = EcoHab.EcoHabSessions(ehd)
        cf = ExperimentConfigFile(path)
        for binsize in [3600]:
            print('Binsize ',binsize/3600)
            results_path = utils.results_path(path)
            fname_all_chambers = '%scollective_results_all_chambers_binsize_%f_h.csv'%(prefix, binsize//3600)
            try:
                cages = non_standard_cages[path]
            except KeyError:
                cages = standard_cages
            try:
                headers = non_standard_headers[path]
            except KeyError:
                headers = standard_headers
            data = get_time_spent_in_each_chamber(ehs, cf, cages, binsize=binsize)
            data = sum_data(data)
            save_data_cvs(data, fname_all_chambers, results_path, cages, headers)
