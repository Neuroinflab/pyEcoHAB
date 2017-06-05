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
### How much time mice spend with each other


datarange = slice(0, 10, None)
datasets = [
    'BALB CTRL EH (3) 12.06.14 - 15.06.14', #XXX
    'BALB VPA EH (4) 16.06.14 - 19.06.14', #XXX
    'C57 VPA EH (1) 20.06.14',  #XXX
    'C57 CTRL EH (3) 25.07.14',
    'C57 CTRL EH (4) 30.07.14',  #XXX
    'C57 VPA EH (1Biol) 28.10.2014',  #XXX
    'BALB CTRL EH (4) 09.09.14', #XXX
    'C57 CTRL EH (5) 20.11.14',  #XXX
    'BALB VPA EH Biol (1) 23.03.15', #XXX
    'C57 VPA EH Biol (2) 27.03.15',  #XXX
    'FX KO females EH 2 29.09.14', #XXX
    'FX KO females EH (3) 11.02.15', #XXX
    'FX WT females EH (1) 16.02.15', #XXX
    'FX WT females EH (2) 20.02.15',
    'FX WT females EH (B) 05.05.16',
    'FX WT females EH (B) - powtorzenie 2. - 18.05.16', # 15
    'FX WT  females EH Lab 2 (1) - powtorzenie 1. - 14.05.16',
    'FX WT  females EH Lab 2 (1) - powtorzenie 2. - 26.05.16',
    'FX KO  females EH Lab 2 (1) - powtorzenie 1. - 10.05.16',
    'FX KO  females EH Lab 2 (1) - powtorzenie 2. - 22.05.16',
    'PV Cre males HT EH (1) 09.06.16',
    ]

# # address = {1: 4, 2: 1, 3: 1, 4: 2, 5: 2, 6: 3, 7: 3, 8: 4}
# smells = {  'C57 CTRL EH (3) 25.07.14': {'soc': 1, 'nsoc': 3},
#             'C57 CTRL EH (4) 30.07.14': {'soc': 3, 'nsoc': 1},
#             'C57 CTRL EH (5) 20.11.14': {'soc': 3, 'nsoc': 1},
#             }
# antenna_positions = {'C57 CTRL EH (3) 25.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
#             'C57 CTRL EH (4) 30.07.14': {'7': 1, '8': 2, '4': 3, '3': 4, '1': 5, '2': 6, '6': 7, '5': 8},
#             'C57 CTRL EH (5) 20.11.14': None,
#             }

def interval_overlap(int1, int2):
    """Return overlap between two intervals."""
    ints = sorted([int1, int2], key=lambda x: x[0])
    if ints[1][0] > ints[0][1]:
        return 0
    else:
        return min(ints[0][1], ints[1][1]) - ints[1][0]

def mice_overlap(data_mice, m1, m2, address):
    """Return time overlap of mice m1 and m2 in cage <address>."""
    ints1 = [(s, e) for a, s, e in data_mice[m1] if a == address]
    ints2 = [(s, e) for a, s, e in data_mice[m2] if a == address]
    durs1 = [x[1] - x[0] for x in ints1]
    durs2 = [x[1] - x[0] for x in ints2]
    total = 0.
    for int1 in ints1:
        for int2 in ints2:
            total += interval_overlap(int1, int2)
    return total, sum(durs1), sum(durs2)

def mice_together(data_mice, m1, m2, total_time=43200.):
    """Return the time spent together by two mice and expected time
    assuming independence."""
    result = np.zeros((4, 3))
    for address in [1, 2, 3, 4]:
        result[address-1] = mice_overlap(data_mice, m1, m2, address)
    fracs = result[:, 1:] / total_time
    time_together = result[:, 0].sum()
    exp_time_together = (fracs[:, 0] * fracs[:, 1]).sum()
    return time_together/total_time, exp_time_together

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

def following(data_mice,times,t_max=43200):
    t1, t2 = times
    t = np.arange(0,t_max,0.1)
    v = np.zeros((len(t),len(mice)))
    i =0
    for mm in mice:
        for rec in data[mm]:
            s = int((rec[1]-t1)*10)
            e = int((rec[2]-t1)*10)
            v[s:e,i] = rec[0]
        i+=1
    return t,v
def norm(x):
    return x/np.linalg.norm(x)

def best_val(x,y):
    maxi = 0
    for i in range(1,11):
        new = np.dot(x[i:],y[:-i])
        if new>maxi:
            maxi = new
        new = np.dot(y[i:],x[:-i])
        if new>maxi:
            maxi = new

    return maxi

def following2(data,v,(idx1,idx2),times,delay = 2, power_step = 3.33):
    foll_signalls = []
    t1, t2 = times
    score_val = 0
    power = 1.
    for i in range(1,len(data)):
        s = int((data[i][1]-t1)*10)
        e = int((data[i][1]-t1+delay)*10)
        period1 = list(v[s:e,idx2])
        period2 = list(v[s-20:e-20,idx2])
        period3 = list(v[s-5:s,idx2])
        if data[i-1][0]==data[i][0]:
            continue
        if (data[i][0] in period1) and (data[i-1][0] in period2) and (data[i][0] not in period3):
            score_val +=power
            power*=power_step
            foll_signalls.append((s,e))
        else:
            power = 1
    try:
        score_val/=len(data)
    except ZeroDivisionError:
        score_val=0
    return score_val, foll_signalls

if __name__ == '__main__':


    alldeltas = {}
    allratios = {}
    allresults = {}


    for el in next(os.walk('../PreprocessedData/'))[2]:
        with open('../PreprocessedData/'+el, "rb") as input_file:
            ehd = pickle.load(input_file)
            ehs = pickle.load(input_file)
            path = el[:-4]
            path1 = '../RawData/' + path
            cf = ExperimentConfigFile(path1)

        tstart, tend = cf.gettime('ALL')
        alldeltas[path] = {}
        allratios[path] = {}
        allresults[path] = {}
        mice = list(ehd.mice)
        mice = filter(lambda x: len(ehs.getstarttimes(x)) > 30, mice)



        phases = filter(lambda x: x.endswith('dark'), cf.sections())
        print el
        #print phases
        # phases = ['SNIFF 1 dark']
        # phases = filter(lambda x: x.endswith('dark') or x.endswith('light'), cf.sections())
        incohort = {}
        scalar = {}
        s_i =0
        following_array = np.zeros((len(mice),len(mice),len(phases)))
        for sec in phases:
            print sec
            #print sec
            data = prepare_data(ehs, mice, cf.gettime(sec))
            t,v = following(data,cf.gettime(sec),t_max=43200)
            results = np.zeros((len(mice), len(mice)))
            results_exp = np.zeros((len(mice), len(mice)))
            for ii in range(len(mice)):
                for jj in range(len(mice)):
                    if ii < jj:
                        key =str(ii)+':'+str(jj)
                        following_array[ii,jj,s_i],f_idx = following2(data[mice[ii]],v,(ii,jj),cf.gettime(sec),delay = 2, power_step = 0.1)
                        following_array[jj,ii,s_i],f_idx = following2(data[mice[jj]],v, (jj,ii),cf.gettime(sec),delay = 2, power_step = 0.1)
                        if ii == 4 and  jj ==8:
                            t = np.arange(-2,2,0.1)
                            for i in range(len(f_idx)):
                                plt.plot(t,v[f_idx[i][0]-20:f_idx[i][1],jj]-0.05,'ro')
                                plt.plot(t,v[f_idx[i][0]-20:f_idx[i][1],ii]+0.05,'bo')
                                plt.axis([-2.1,2.1,-0.5,5])
                                plt.show()
                        """res = mice_together(data, mice[ii], mice[jj])
                        results[ii, jj] = res[0]
                        results_exp[ii, jj] = res[1]
                        try:
                            incohort[key].append(res[0]-res[1])
                        except:
                            incohort[key] = [res[0]-res[1]]"""
                        """sca = best_val(norm(v[:,ii]),norm(v[:,jj]))
                        if sca<0.5:
                            print sec, key, sca
                        try:
                            scalar[key].append(sca)
                        except:
                            scalar[key] = [sca]"""
            """np.savetxt('results_%s_%s.csv' %(path, sec), results, fmt='%.6f', delimiter=';')
            np.savetxt('results_exp_%s_%s.csv' %(path, sec), results_exp,
                    fmt='%.6f', delimiter=';')
            np.savetxt('results_final_%s_%s.csv' %(path, sec), results-results_exp,
                    fmt='%.6f', delimiter=';')"""
            #plt.imshow(following_array[:,:,s_i],interpolation='none')
            #plt.savefig('Results/'+path+'_'+str(s_i)+'.png')
            #s_i+=1
        #fig = plt.figure()
        """experiment_fol = np.zeros((len(mice)**2-len(mice),6))
        #ax = fig.add_subplot(111, projection='3d')
        for ii in range(len(mice)):
            for jj in range(len(mice)):
                if ii !=jj:
                    x = np.ones(ii*len(mice)+jj)
                    y = np.arange(len(phases))
                    z = following_array[ii,jj,:]
                    experiment_fol[ii*(len(mice)-1)+jj,:] = following_array[ii,jj,:]
            #ax.scatter(x,y,z)
        plt.imshow(experiment_fol,interpolation='none')
        plt.show()"""
        """fig = plt.figure()
        for key in scalar.keys():
            x = np.arange(len(phases))
            plt.plot(x,scalar[key],'o-')
        my_xticks = phases
        plt.xlim([-0.2,len(phases)-0.8])
        plt.xticks(x, my_xticks)
        fig.set_size_inches(10,6)
        fig.savefig('Results/'+path+'_scalar.png')"""
        """fig = plt.figure()
        for key in incohort.keys():
            x = np.arange(len(phases))
            plt.plot(x,incohort[key],'o-')
        my_xticks = phases
        plt.xlim([-0.2,len(phases)-0.8])
        plt.xticks(x, my_xticks)
        fig.set_size_inches(10,6)
        fig.savefig('Results/'+path+'_incohort.png')"""



        """ plt.figure(figsize=(10, 6))
            plt.subplot(221)
            plt.imshow(results, vmin=0, vmax=0.5)
            plt.xticks([])
            plt.yticks([])
            plt.title('% time together')
            plt.subplot(222)
            plt.imshow(results_exp, vmin=0, vmax=0.5)
            plt.xticks([])
            plt.yticks([])
            plt.title('Expected % time together')
            deltas = results[results > 0] - results_exp[results > 0]
            alldeltas[path][sec] = deltas
            ratios = results[results > 0] / results_exp[results > 0]
            allratios[path][sec] = ratios
            allresults[path][sec] = results[results > 0]
            plt.subplot(223)
            try:
                # plt.imshow(results - results_exp, vmin=-max(np.abs(deltas)),
                #  vmax=max(np.abs(deltas)))
                plt.imshow(results - results_exp, vmin=-0.25, vmax=0.25)
                plt.xticks([])
                plt.yticks([])
            except ValueError:
                pass
            plt.title('Excess % time together')
            plt.colorbar()
            plt.subplot(224)
            try:
                plt.hist(deltas)
            except ValueError:
                pass
            plt.title('Histogram of excess % time together')
            plt.xlim([-0.1, 0.3])
            plt.suptitle(path)
            plt.savefig('friends_%s_%s.pdf' %(path.translate(None, '/'), sec))
            plt.savefig('friends_%s_%s.png' %(path.translate(None, '/'), sec), dpi=300)"""
        # np.save('alldeltas_EH_new_%s.npy' %path, alldeltas[path])
        # np.save('allratios_EH_new_%s.npy' %path, allratios[path])
        # np.save('allesults_EH_new_%s.npy' %path, allresults[path])


    # np.save('alldeltas.npy', alldeltas)
    # np.save('allratios.npy', allratios)
    # np.save('allesults.npy', allresults)
