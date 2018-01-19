from __future__ import division, print_function
import time
import numpy as np
import os

def new_state_clockwise(state):
    if state < 8:
        return state+1
    return 1

def new_state_counter_clockwise(state):
    if state > 1:
        return state-1
    return 8

def increment(previous):
    if previous%2:
        return 0.5
    return 1

def make_str(k):
    if k < 10:
        return "0"+str(k)
    return str(k)

def make_date(t):
    st = time.localtime(t)
    date = '%d.%s.%s'%(st.tm_year,make_str(st.tm_mon),make_str(st.tm_mday))
    return date

def make_date_inverse(t):
    st = time.localtime(t)
    date = '%s.%s.%d'%(make_str(st.tm_mday),make_str(st.tm_mon),st.tm_year)
    return date

def make_hour(t):
    stru = time.localtime(t)
    po_kropce = str(int(1000*(t - np.floor(t))))
    return make_str(stru.tm_hour)+':'+make_str(stru.tm_min)+':'+make_str(stru.tm_sec)+'.'+po_kropce

def write_config(config,phase,t,end):
    date = make_date_inverse(t)
    hour = make_hour(t)[:5]
    if end:
         config.write('enddate = %s\n'%date)
         config.write('endtime = %s\n'%hour)
         config.write('\n')
    else:
        phase_name = '['+phase+']\n'
        config.write(phase_name)
        config.write('startdate = %s\n'%date)
        config.write('starttime = %s\n'%hour)

def open_new_file(t,dir_name):
    st = time.localtime(t)
    date = str(st.tm_year)+'.'+make_str(st.tm_mon)+'.'+make_str(st.tm_mday)
    new_fname = str(st.tm_year)+make_str(st.tm_mon)+make_str(st.tm_mday)+"_"+make_str(st.tm_hour)+"0000.txt"
    location = os.path.join(dir_name,new_fname)
    return open(location,'w')

def make_data(dir_name,phases,endings,m2_function,phase_duration,exp_duration,state_m1=1,state_m2=1,previous_m1=8,previous_m2=8,mouse1 = "0065-0136661698",mouse2 = "0065-0136656570"):
    
    counter = 1
    t0 = 1516030860.-41*60-4*60*60
    time_m1 = t0
    time_m2 = t0 + 0.21
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    fname_config = os.path.join(dir_name,'config.txt')
    config = open(fname_config,'w')

    t = t0
    t_start = t0
    for i,phase in enumerate(phases):
        for j in range(exp_duration[i]):

            phase_name = phase+' '+str(j+1)+' '+endings[j%2]
            write_config(config,phase_name,t,False)
            
            while (t < t_start + phase_duration):
                date = make_date(t)
                f = open_new_file(t,dir_name)
                time_counter = t
                while t < time_counter+60*60:
                    f.write(str(counter)+"\t"+date+"\t")
                    if time_m1 < time_m2:
                        f.write(make_hour(time_m1)+"\t"+str(state_m1)+"\t200\t"+mouse1+"\n")
                        previous_m1 = state_m1
                        state_m1 = new_state_clockwise(state_m1)
                        time_m1 += increment(previous_m1)
                        t = time_m1
                        counter += 1
                    else:
                        f.write(make_hour(time_m2)+"\t"+str(state_m2)+"\t200\t"+mouse2+"\n")
                        previous_m2 = state_m2
                        state_m2 = m2_function[i](state_m2)
                        time_m2 += increment(previous_m2)
                        t = time_m2
                        counter +=1
            t_start = t
            write_config(config,phase_name,t,True)

    write_config(config,'ALL',t0,False)
    write_config(config,'ALL',t,True)
    
