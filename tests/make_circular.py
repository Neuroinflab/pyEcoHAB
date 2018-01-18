import time
import numpy as np
import os

def new_state(state):
    if state < 8:
        return state+1
    return 1

def increment(previous):
    if previous%2:
        return 0.5
    return 1

def make_str(k):
    if k < 10:
        return "0"+str(k)
    return str(k)

def make_hour(t):
    stru = time.localtime(t)
    po_kropce = str(int(1000*(t - np.floor(t))))
    return make_str(stru.tm_hour)+':'+make_str(stru.tm_min)+':'+make_str(stru.tm_sec)+'.'+po_kropce


if __name__ == '__main__':
    t0 = 1516030871.-41*60

    #time.localtime(t0) = time.struct_time(tm_year=2018, tm_mon=1, tm_mday=15, tm_hour=16, tm_min=41, tm_sec=11, tm_wday=0, tm_yday=15, tm_isdst=0)
    
    st = time.localtime(t0)
    
    fname = str(st.tm_year)+make_str(st.tm_mon)+make_str(st.tm_mday)+"_"+make_str(st.tm_hour)+"0000.txt"
    
    date = str(st.tm_year)+'.'+make_str(st.tm_mon)+'.'+make_str(st.tm_mday)
    
    f = open(fname,"w")

    state_m1 = 1
    state_m2 = 1
    previous_m1 = 8
    previous_m2 = 8
    mouse1 = "0065-0136661698"
    mouse2 = "0065-0136656570"

    counter = 1
    time_m1 = t0
    time_m2 = t0+0.21

    while True:
        f.write(str(counter)+"\t"+date+"\t")
        if time_m1 < time_m2:
            f.write(make_hour(time_m1)+"\t"+str(state_m1)+"\t200\t"+mouse1+"\n")
            previous_m1 = state_m1
            state_m1 = new_state(state_m1)
            time_m1 += increment(previous_m1)
        else:
            f.write(make_hour(time_m2)+"\t"+str(state_m2)+"\t200\t"+mouse2+"\n")
            previous_m2 = state_m2
            state_m2 = new_state(state_m2)
            time_m2 += increment(previous_m2)
            counter += 1
            
        if counter > 3600:
            break
    

            
