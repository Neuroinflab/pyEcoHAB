from __future__ import print_function, division
import EcoHab
from ExperimentConfigFile import ExperimentConfigFile
from data_info import *
import os
import utils
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
bins = 2000
homepath = os.path.expanduser("~/")
threshold = 2

def in_tube(antenna, next_antenna):
    if antenna % 2:
        if next_antenna  == antenna  + 1:
            return True
    else:
        if next_antenna == antenna - 1:
            return True
    
    return False

if __name__ == '__main__':

    for new_path in datasets:
       
        path = os.path.join(homepath, new_path)
        prefix = utils.make_prefix(path)
        if new_path in remove_tags:
            remove_mouse = remove_tags[new_path]
        else:
            remove_mouse = None
        if new_path not in antenna_positions:
            antenna_positions[new_path] = None
        if new_path not in how_many_appearances:
            how_many_appearances[new_path] = 500
        if remove_mouse:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    remove_mice=remove_mouse,
                                    how_many_appearances=how_many_appearances[new_path])
        else:
            ehd = EcoHab.EcoHabData(path=path,
                                    _ant_pos=antenna_positions[new_path],
                                    how_many_appearances=how_many_appearances[new_path])

        new_data = pd.DataFrame.from_dict(ehd.data)
        
        titles = ['12', '34', '56', '78']
        mice = ehd.mice
        tag_data = np.array(ehd.data['Tag'])
        outs = {}
        cf = ExperimentConfigFile(path)
        phases = cf.sections()
        followings = {}
        for phase in phases:
            st, en = cf.gettime(phase)
            ehd.mask_data(st, en)
            followings[phase] = {}
            for mouse1 in mice:
                antennas1 = ehd.getantennas(mouse1)
                times1 = ehd.gettimes(mouse1)
                for mouse2 in mice:
                    if mouse2 != mouse1:
                        antennas2 = ehd.getantennas(mouse2)
                        times2 = ehd.gettimes(mouse2)
                        key = mouse1 + '_' + mouse2
                        followings[phase][key] = {'antennas':[],
                                                  'mouse1_time':[],
                                                  'mouse2_time':[],
                                                  'mouse1_delta':[],
                                                  'mouse2_delta':[]}
               

                        for i, antenna1 in enumerate(antennas1[:-1]):
                            next_antenna1 = antennas1[i+1]
                            delta_t1 = times1[i+1] - times1[i]
                            key1 = str(antenna1) + str(next_antenna1)
                            which = in_tube(antenna1, next_antenna1)                
                            if which:
                                
                                 idxs1 = np.where(np.array(times2) >= times1[i])[0]
                                 idxs2 = np.where(np.array(times2) <= times1[i]+threshold)[0]
                                 common_idxs = list(set(idxs1)&set(idxs2))
                                 if len(common_idxs):
                                     for ci in common_idxs:
                                         antenna2 = antennas2[ci]
                                         if ci + 1 < len(antennas2):
                                             next_antenna2 = antennas2[ci+1]
                                             delta_t2 = times2[ci+1] - times2[ci]
                                             if antenna2 == antenna1:
                                                 if in_tube(antenna2, next_antenna2):
                                                     followings[phase][key]['antennas'].append(key1)
                                                     followings[phase][key]['mouse1_time'].append(times1[i])
                                                     followings[phase][key]['mouse2_time'].append(times2[ci])
                                                     followings[phase][key]['mouse1_delta'].append(delta_t1)
                                                     followings[phase][key]['mouse2_delta'].append(delta_t2)
                                                     
            for key in followings[phase].keys():                                 
                print(followings[phase][key])
                                     
                                
    data = pd.DataFrame.from_dict(followings)
    data.to_csv('followings_wyniki')

            
            # fig, ax = plt.subplots(1, 1)
            # fig.suptitle(mouse)
            # n, bins, patches = ax.hist(outs[mouse], bins=bins)
            # #ax.set_xlim([0, 10])
            # ax.set_xscale("log", nonposx='clip')
            # plt.show()
                
                                      
 
            
