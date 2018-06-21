from __future__ import print_function, division, absolute_import

import os
import glob
import time
import numpy as np
import sys
from DataBase import Data, AntennaReadOut, Animal
import utils
from collections import Container, Counter
from operator import methodcaller, attrgetter


    


class EcoHabData(Data):
    """Reads in a folder with data from Eco-HAB"""
    standard_antenna_positions = {'1': 1, '2': 2, '3': 3, '4': 4,
                                  '5': 5, '6': 6, '7': 7, '8': 8}
    pattern1 = '*0000.txt'
    pattern2 = '*0000_???.txt'
    read_in_dict  = {
        5: utils.process_line_5,
        6: utils.process_line_6,
        7: utils.process_line_7,             
    }
    
    def __init__(self, path, antenna_positions=None, verbose=False, **kwargs):
        super(EcoHabData,self).__init__()       
        self.path = path
        self.verbose = verbose
        
        if isinstance(antenna_positions, dict):
            self.antenna_positions = antenna_positions
            
        elif antenna_positions is None:
            self.antenna_positions = self.standard_antenna_positions
            
        how_many_appearances = kwargs.pop('how_many_appearances',100) # to exclude fake tags aka ghost mice
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])
        max_break = kwargs.pop('max_break', 60*60)
        raw_data = self.read_in_data(self.path, how_many_appearances, factor, tags)
        readouts = self.extract_readouts(raw_data)
        self.add_readouts(readouts)

        # antenna_breaks = self.check_antenna_presence(max_break)
        # if antenna_breaks:
        #     print('Antenna not working')
        #     for antenna in antenna_breaks:
        #         print(antenna,antenna_breaks[antenna])
        # self.antenna_mismatch()
        
        # if mask:
        #     self._cut_out_data(mask)


    @staticmethod    
    def remove_ghost_tags(data, how_many_appearances,factor,tags=[]):
        new_data = []
        ghost_mice = set()
        counters = Counter()
        dates = {}
        
        if isinstance(tags, list):
            for tag in tags:
                ghost_mice.add(tag)
        elif isinstance(tags, str):
            ghost_mice.add(tags)
        elif tags in None:
            pass
        else:
            print('Unknown tag format for removal of mice', tags) 
            
        for d in data:
            mouse = d[4]
            print(mouse)
            if mouse not in dates:
                dates[mouse] = Counter()
            counters[mouse] += 1
            date = d[1].split()[0]
 
            dates[mouse][date] += 1

        how_many = utils.how_many_days(dates, factor)
        
        for mouse in counters:
            if counters[mouse] < how_many_appearances or len(dates[mouse]) <= how_many:
                    ghost_mice.add(mouse)
   
        for d in data:
            mouse = d[4]
            if mouse not in ghost_mice:
                new_data.append(d)

        return new_data[:]
    
    def find_files(self, path):
        """Finds files in provided directory (path) and reads in data"""
        fname1 = os.path.join(path, self.pattern1)
        fname2 = os.path.join(path, self.pattern2)
        return glob.glob(fname1) + glob.glob(fname2)
 
    def append_data(self, file_list):
        """Reads in data file"""
        results = []
        for fname in file_list:
            hour, date, datenext = self.parse_fname(fname)
            single_hour_data = self.read_in_single_file(fname, (hour, date, datenext))
            results += single_hour_data
        return results
    
 
    def read_in_data(self, path, how_many_appearances, factor, tags):
        files = self.find_files(path)
        raw_data = self.append_data(files)
        return self.remove_ghost_tags(raw_data, how_many_appearances, factor, tags)

        
    @staticmethod
    def extract_readouts(data):
        readouts = np.zeros_like(data, dtype=AntennaReadOut)
        for i, data_line in enumerate(data):
            readout[i] = AntennaReadOut(EventId=data_line[0],
                                         Start=data_line[1],
                                         AntennaId=data_line[2],
                                         Duration=data_line[3],
                                         Animal=Animal(Tag=data_line[4],))
        return readouts
        

  

    @staticmethod
    def parse_fname(fname):
        """"Extracts time and date from data's filename"""
        fname = os.path.split(fname)[-1]
            
        hour = fname[9:11]
        date = fname[:8]
        datenext = time.strftime('%Y%m%d',
                                 time.localtime(time.mktime(time.strptime(date, '%Y%m%d')) + 24 * 3600.))

        return hour, date, datenext

    def read_in_single_file(self, fname, date_hour_obj):
        file_obj = open(fname)
        out = []
        for line in file_obj:
            elements = str.strip(line).split()
            length = len(elements)
            try:
                out += [self.read_in_dict[length](elements, date_hour_obj)]
            except KeyError:
                raise(IOError('Unknown data format in file %s' %f.name))
        return out
    
    
    
        
    def check_antenna_presence(self):
        t_start = self.get_start()
        all_times = self.readouts.getAttributes('Start')
        breaks = {}
        
        for antenna in range(1, 9):
            
            antenna_idx = np.where(np.array(self.data['Antenna']) == antenna)[0]
            times = all_times[antenna_idx] 
            breaks[antenna] = []
            if len(times):
                if times[0] - t_start > self.max_break:
                    breaks[antenna].append([0, np.round(times[0])])
        
                intervals = times[1:]-times[0:-1]

                where_breaks = np.where(intervals > self.max_break)[0]

                if len(where_breaks):
                    for i in where_breaks:
                        breaks[antenna].append([np.round(times[i]), np.round(times[i+1])])
            else:
                breaks[antenna].append([np.round(t_start),self.data['Time'][-1]])         
                    
        # for antenna in breaks:
        #     print(antenna, breaks[antenna])
        return breaks
    
    def antenna_mismatch(self):
        
        t_start = self.data['Time'][0]
        all_times = np.array(self.data['Time'])
        weird_transit = [[], []]
        for mouse in self.mice:
        
            mouse_idx = np.where(np.array(self.data['Tag']) == mouse)[0]
            times = all_times[mouse_idx]
            antennas = np.array(self.data['Antenna'])[mouse_idx]

            for i, a in enumerate(antennas[:-1]):
                if abs(a - antennas[i+1]) not in [0,1,7]:
                    weird_transit[0].append(times[i])
                    if a < antennas[i+1]:
                        weird_transit[1].append(str(a)+' '+str(antennas[i+1]))
                    else:
                        weird_transit[1].append(str(antennas[i+1])+' '+str(a))
                    
        pairs = list(set(weird_transit[1]))
        
        mismatches = {}

        for pair in pairs:
            mismatches[pair] = weird_transit[1].count(pair)
            print(pair, mismatches[pair],np.round(100*mismatches[pair]/len(self.data['Antenna'])),'per 100')

        return weird_transit

    
                        
    def __repr__ (self):
        """Nice string representation for printing this class."""
        mystring = 'Eco-HAB data loaded from:\n%s\nin the folder%s\n' %(
                   self._fnames.__str__(), self.path) 
        return mystring
                            
    @staticmethod
    def convert_time(s): 
        """Convert date and time to seconds since epoch"""
        actual_date, millisec = s.split('.')
        sec_to_epoch = time.mktime(time.strptime(actual_date, '%Y%m%d %H:%M:%S'))
        return sec_to_epoch + float(millisec)/1000
    
    def checkData_one_mouse(self,mouse):
        antennas = self.getantennas(mouse)
        times  = self.gettimes(mouse)
        wrong_times = []
        for i,next_antenna in enumerate(antennas[1:]):
            if abs(next_antenna - antennas[i]) not in [0,1]:
                if next_antenna !=8 and antennas[i] != 8:
                    wrong_times.append(times[i])
                    #print('mouse',mouse,'wrong antenna',antennas[i],'time',times[i],'antenna',next_antenna,'next time',times[i+1],times[i+1]-times[i])
        
        return wrong_times
    
class IEcoHabSession(object):
    def getstarttimes(self, *arg, **kwarg):
        raise NotImplementedError("Virtual method called")

    def getendtimes(self, *arg, **kwarg):
        raise NotImplementedError("Virtual method called")

    def getdurations(self, *arg, **kwarg):
        raise NotImplementedError("Virtual method called")

    def getaddresses(self, *arg, **kwarg):
        raise NotImplementedError("Virtual method called")

    def getproperty(self, *arg, **kwarg):
        raise NotImplementedError("Virtual method called")

class EcoHabSessions(IEcoHabSession):
    """Calculates 'visits' to Eco-HAB compartments."""
    def __init__(self, ehd, **kwargs):
        
        self._ehd = ehd
        self.mask = None
        self._mask_slice = None
        # No longer used after 14 May 2014
        # self.same_antenna_threshold = kwargs.pop('same_antenna_threshold', 2.)
        self.shortest_session_threshold = kwargs.pop('shortest_session_threshold', 2.)
        
        tempdata = []
        same_pipe = {1: [1, 2], 2: [1, 2], 3: [3, 4], 4: [3, 4],
                     5: [5, 6], 6: [5, 6], 7: [7, 8], 8: [7, 8]}
        opposite_pipe = {1: [5, 6], 2: [5, 6], 3: [7, 8], 4: [7, 8],
                     5: [1, 2], 6: [1, 2], 7: [3, 4], 8: [3, 4]}
        address = {1: 4, 2: 1, 3: 1, 4: 2, 5: 2, 6: 3, 7: 3, 8: 4}
        address_not_adjacent = {1: 1, 2: 4, 3: 2, 4: 1, 5: 3, 6: 2, 7: 4, 8: 3}
        # Surrounding: difference between antennas only 2 or 6 
        surrounding = {(1, 3): 1, (1, 7): 4, (2, 4): 1, (2, 8): 4,
                       (3, 5): 2, (4, 6): 2, (5, 7): 3, (6, 8): 3}
        
        for mm in ehd.mice:
            tt = self._ehd.gettimes(mm)
            an = self._ehd.getantennas(mm)

            for tstart, tend, anstart, anend in zip(tt[:-1], tt[1:], an[:-1], an[1:]):
                # # Old code - until 14 May 2014
                # if tend - tstart < self.shortest_session_threshold:
                #     continue
                # if anend in same_pipe[anstart]:
                #     if (anend != anstart) or (tend - tstart < self.same_antenna_threshold):
                #         continue
                #         # Jesli odczyty sa w tej samej rurze, ale w wiekszym odstepie niz
                #         # 2 sekundy, to raczej przyjmujemy ze sesja byla - poprawic.
                # valid_solution = np.abs(anstart - anend) in [0, 1, 7]
                # tempdata.append((address[anstart], mm, tstart, tend, tend-tstart,
                #              valid_solution))

                
                # Workflow as agreed on 14 May 2014
                if tend - tstart < self.shortest_session_threshold:
                    continue
                diff = np.abs(anstart - anend)
                if diff == 0:
                    tempdata.append((address[anstart], mm, tstart, tend, tend-tstart,
                                 True))
                elif diff in [1, 7]:
                    if anend in same_pipe[anstart]:
                        continue
                    else:
                        tempdata.append((address[anstart], mm, tstart, tend, tend-tstart,
                                     True))
                elif diff in [2, 6]:
                    tempdata.append((surrounding[(min(anstart, anend), max(anstart, anend))], 
                                mm, tstart, tend, tend-tstart,
                                False))
                elif diff in [3, 4, 5]:
                    if anend in opposite_pipe[anstart]:
                        continue
                    else:
                        tempdata.append((address_not_adjacent[anstart], 
                                    mm, tstart, tend, tend-tstart,
                                    False))
                
                            
        tempdata.sort(key=lambda x: x[2])

        self.data = {'Tag': [],
             'Address': [],
             'AbsStartTimecode': [],
             'AbsEndTimecode': [],
             'VisitDuration': [],
             'ValidVisitSolution': [],}
        self.data['Address'] = [x[0] for x in tempdata]
        self.data['Tag'] = [x[1] for x in tempdata]
        self.data['AbsStartTimecode'] = [x[2] for x in tempdata]
        self.data['AbsEndTimecode'] = [x[3] for x in tempdata]
        self.data['VisitDuration'] = [x[4] for x in tempdata]
        self.data['ValidVisitSolution'] = [x[5] for x in tempdata]
        mice = self._ehd.mice
        self.mice = [mouse for mouse in mice if len(self.getstarttimes(mouse)) > 30]
        
    def unmask_data(self):
        """Remove the mask - future queries will not be clipped"""
        self.mask = None
        self._mask_slice = None

    def mask_data(self, *args):
        """mask_data(endtime) or mask_data(starttime, endtime)
        All future queries will be clipped to the visits starting between
        starttime and endtime."""
        try:
            starttime = args[0]
            endtime = args[1]
        except IndexError:   
            starttime = min(self.getstarttimes(self._ehd.mice))
            endtime = args[0]
        self.mask = (starttime, endtime) 
        arr = np.array(self.data['AbsStartTimecode'])
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            self._mask_slice = (idcs[0], idcs[-1] + 1)
        elif len(idcs) == 1:
            self._mask_slice = (idcs[0], idcs[0] + 1)
        else:
            self._mask_slice = (0, 0)

    def getproperty(self, mice, propname, astype=None):
        if sys.version_info < (3, 0):
            if isinstance(mice, (str, unicode)):
                mice = [mice]
        else:
            if isinstance(mice, str):
                mice = [mice]

        if self.mask is None:
            if astype is None:
                return [x[0] for x in zip(self.data[propname], 
                        self.data['Tag']) if x[1] in mice]
            elif astype == 'float':                          
                return [float(x[0]) for x in zip(self.data[propname], 
                        self.data['Tag']) if x[1] in mice]
        else:
            if astype is None:
                return [x[0] for x in zip(
                        self.data[propname][self._mask_slice[0]:self._mask_slice[1]], 
                        self.data['Tag'][self._mask_slice[0]:self._mask_slice[1]]) 
                        if x[1] in mice] 
            elif astype == 'float':
                return [float(x[0]) for x in zip(
                        self.data[propname][self._mask_slice[0]:self._mask_slice[1]], 
                        self.data['Tag'][self._mask_slice[0]:self._mask_slice[1]]) 
                        if x[1] in mice]
                    
    def getstarttimes(self, mice): 
        return self.getproperty(mice, 'AbsStartTimecode', 'float')
                    
    def getendtimes(self, mice):
        return self.getproperty(mice, 'AbsEndTimecode', 'float')
                    
    def getdurations(self, mice): 
        return self.getproperty(mice, 'VisitDuration', 'float')
    
    def getaddresses(self, mice): 
        return self.getproperty(mice, 'Address')
    
    def getstats(self, mm):
        """Return total number of visits 
        and total time spent in compartments."""
        durations = self.getdurations(mm)
        adds = self.getaddresses(mm)
        totv = [0, 0, 0, 0]
        tott = [0., 0., 0., 0.]
        for idx, ad in enumerate([1, 2, 3, 4]):
            durs = [x for x, y in zip(durations, adds) if y == ad]
            totv[idx] = len(durs)
            tott[idx] = sum(durs)
        return totv, tott
        
        
        
class EcoHabSessions9states(Data,IEcoHabSession):
    """Calculates 'visits' to Eco-HAB compartments."""
    def plot_histograms_time_spent_each_state(self,mouse):
        fig = plt.figure()
        ax = []
        fig.suptitle(mouse)
        for i in range(8):
            ax.append(fig.add_subplot(2,4,i+1))
       
        for i in range(4):
            chamber = np.array(self.statistics[mouse]["state_time"][2*(i+1)])
            corridor  = np.array(self.statistics[mouse]["state_time"][2*i+1])
            ax[i].hist(corridor,bins=100)
            ax[i].set_title(str(2*i+1))
            ax[i+4].hist(chamber,bins=100)
            ax[i+4].set_title(str(2*(i+1)))
            
        plt.show()
        
    def _initialize_statistics(self,statistics,mouse):
        
        statistics[mouse] = {}
        statistics[mouse]["state_freq"]= np.zeros(9)
        statistics[mouse]["state_time"]= [[] for i in range(9)]
        statistics[mouse]["preference"]={}
        for i in range(9):
            statistics[mouse]["preference"][i] = np.zeros(2)
            
    def _calculate_visitis(self, too_fast=2,compensate_for_lost_antenna=False,in_pipe=0.7):
        
        statistics = {}
        tempdata = []
 
        
        for mm in self.mice:
            
            tt = self._ehd.gettimes(mm)
            an = self._ehd.getantennas(mm)
            
            self._initialize_statistics(statistics,mm)
                 
            previous = 0
            
            for tstart, tend, anstart, anend in zip(tt[:-1], tt[1:], an[:-1], an[1:]):
                
                t_diff = tend - tstart

                if  t_diff < self.shortest_session_threshold:
                    state = 0
                    previous = 0
                    statistics[mm]["state_freq"][state]+=1
                    statistics[mm]["state_time"][state].append(tend - tstart)
                    continue
                
                diff = np.abs(anstart - anend)
                s = int((tstart - self.t_start_exp)*self.fs)
                e = int((tend- self.t_start_exp)*self.fs)

                if diff in [0,1,7]:
                    ############Most obvious reading##########
                    if diff in [1, 7]:
                        #printprevious,)
                        state = self._find_state(anstart,anend)
                        tempdata.append((state, mm, tstart, tend, t_diff, True))
                        self.signal_data[mm][int(s):int(e)] = state
                        previous = self._update_statistics(statistics,mm,state,anend,t_diff)
                        
                    elif diff == 0 and previous != 0:
                       
                        if t_diff < too_fast:
                            continue
                        diff2 = anstart - previous
                        state  = 0
                        if diff2 >=0:
                            if diff2 == 0 and previous != 1:
                                state = previous-1
                            elif diff2 == 0 and previous == 1:
                                state = 8
                            elif diff2 == 1:
                                state = previous+1
                            elif diff2 == 7:
                                state = 1
                            
                        tempdata.append((state, mm, tstart, tend, t_diff, True))
                        self.signal_data[mm][int(s):int(e)] = state
                        previous = self._update_statistics(statistics,mm,state,anend,t_diff)
                    
                elif compensate_for_lost_antenna and diff in [2,6]:
                    
                        if t_diff < too_fast:
                            continue
                        
                        if anstart == 8 and anend == 2 :
                            middle_antenna = 1
                        elif anstart == 2 and anend == 8:
                            middle_antenna = 1
                        elif anstart == 7 and anend == 1 :
                            middle_antenna = 8
                        elif anstart == 1 and anend == 7:
                            middle_antenna = 8
                        else:
                            middle_antenna = anstart + (anend-anstart)//2

                        previous = self._find_state(anstart,middle_antenna)
                        state = self._find_state(middle_antenna,anend)
                    
                        duration = e-s
                   
                        if anstart % 2:
                            middle = s + int(in_pipe*self.fs)
                            t_diff_prev = in_pipe
                            t_diff_post = t_diff - in_pipe
                        else:
                            middle = s + duration - int(in_pipe*self.fs)
                            t_diff_prev = t_diff - in_pipe
                            t_diff_post = in_pipe
                            
                        self.signal_data[mm][s:middle] = previous
                        self.signal_data[mm][middle:e] = state
                    
                        statistics[mm]["state_freq"][previous]+=1
                        statistics[mm]["state_freq"][state]+=1
                    
                        statistics[mm]["state_time"][previous].append(t_diff_prev)
                        statistics[mm]["state_time"][state].append(t_diff_post)
                        tempdata.append((previous,mm, tstart,tend,t_diff_prev ,False))
                        tempdata.append((previous,mm, tstart,tend,t_diff_post ,False))
                        previous = state
            

        tempdata.sort(key=lambda x: x[2])
        self.data = {'Tag': [],
             'Address': [],
             'AbsStartTimecode': [],
             'AbsEndTimecode': [],
             'VisitDuration': [],
             'ValidVisitSolution': [],}
        self.data['Address'] = [x[0] for x in tempdata]
        self.data['Tag'] = [x[1] for x in tempdata]
        self.data['AbsStartTimecode'] = [x[2] for x in tempdata]
        self.data['AbsEndTimecode'] = [x[3] for x in tempdata]
        self.data['VisitDuration'] = [x[4] for x in tempdata]
        self.data['ValidVisitSolution'] = [x[5] for x in tempdata]

        return statistics

    def _find_state(self,anstart,anend):
        if abs(anstart-anend) == 1:
            return int((anstart + anend)/2.-0.5)
        if abs(anstart-anend) == 7:
            return 8
        return 0

    def _update_statistics(self,statistics,mm,state,anend,t_diff):
        statistics[mm]["state_freq"][state]+=1
        statistics[mm]["state_time"][state].append(t_diff)
        if anend == state:
            statistics[mm]["preference"][state][1]+=1
        else:
            statistics[mm]["preference"][state][0]+=1  

        return state
        
    def __init__(self, **kwargs):
        path = kwargs.pop("path")
        _ant_pos = kwargs.pop('_ant_pos', None)
        _mask = kwargs.pop('mask', None)
        how_many_appearances = kwargs.pop('how_many_appearances',1000)
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])
        compensate_for_lost_antenna = kwargs.pop('compensate_for_lost_antenna',False)
        ehd = EcoHabData(path=path,_ant_pos=_ant_pos,mask=_mask,how_many_appearances=how_many_appearances, factor=factor,remove_mice=tags)
        self._ehd = ehd
        data = {}
        mask = None
        self.mice = self._ehd.mice
        super(EcoHabSessions9states,self).__init__(data=data,mask=mask)
        self.shortest_session_threshold = kwargs.pop('shortest_session_threshold', 2)
        self.fs = 10
        self.t_start_exp = np.min(self._ehd.data['Time'])
        self.t_end_exp = np.max(self._ehd.data['Time'])
        t = np.arange(self.t_start_exp,self.t_end_exp,1/self.fs)
        self.signal_data = {}
        self.signal_data['time'] = t
        for mouse in self.mice:
            self.signal_data[mouse] = np.zeros(len(t),dtype=np.int8)          
        
        self.statistics = self._calculate_visitis(compensate_for_lost_antenna=compensate_for_lost_antenna)


            
            
    def mask_data(self, *args):
        """mask_data(endtime) or mask_data(starttime, endtime)
        All future queries will be clipped to the visits starting between
        starttime and endtime."""
        try:
            starttime = args[0]
            endtime = args[1]
        except IndexError:   
            starttime = min(self.getstarttimes(self.mice))
            endtime = args[0]
        self.mask = (starttime, endtime) 
        arr = np.array(self.data['AbsStartTimecode'])
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            self._mask_slice = (idcs[0], idcs[-1] + 1)
        elif len(idcs) == 1:
            self._mask_slice = (idcs[0], idcs[0] + 1)
        else:
            self._mask_slice = (0, 0)

                     
    def getstarttimes(self, mice): 
        return self.getproperty(mice, 'AbsStartTimecode', 'float')
                    
    def getendtimes(self, mice):
        return self.getproperty(mice, 'AbsEndTimecode', 'float')
                    
    def getdurations(self, mice): 
        return self.getproperty(mice, 'VisitDuration', 'float')
    
    def getaddresses(self, mice): 
        return self.getproperty(mice, 'Address')
    
    def getstats(self, mm):
        """Return total number of visits 
        and total time spent in compartments."""
        durations = self.getdurations(mm)
        adds = self.getaddresses(mm)
        totv = [0, 0, 0, 0] #total number of visits
        tott = [0., 0., 0., 0.] #total number of time
        for idx, ad in enumerate([1, 2, 3, 4]):
            durs = [x for x, y in zip(durations, adds) if y == ad]
            totv[idx] = len(durs)
            tott[idx] = sum(durs)
        return totv, tott
    def getproperty(self, mice, propname, astype=None):
        if sys.version_info < (3, 0):
            if isinstance(mice, (str, unicode)):
                mice = [mice]
        else:
            if isinstance(mice, str):
                mice = [mice]

 
        if self.mask is None:
            if astype is None:
                return [x[0] for x in zip(self.data[propname], 
                        self.data['Tag']) if x[1] in mice]
            elif astype == 'float':                          
                return [float(x[0]) for x in zip(self.data[propname], 
                        self.data['Tag']) if x[1] in mice]
        else:
            if astype is None:
                return [x[0] for x in zip(
                        self.data[propname][self._mask_slice[0]:self._mask_slice[1]], 
                        self.data['Tag'][self._mask_slice[0]:self._mask_slice[1]]) 
                        if x[1] in mice] 
            elif astype == 'float':
                return [float(x[0]) for x in zip(
                        self.data[propname][self._mask_slice[0]:self._mask_slice[1]], 
                        self.data['Tag'][self._mask_slice[0]:self._mask_slice[1]]) 
                        if x[1] in mice]
        
if __name__ == '__main__':
    eco = EcoHabData('/home/jszmek/EcoHAB_data_November/Social structure males 02.03/')
   
