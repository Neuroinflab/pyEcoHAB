from __future__ import print_function,division

import os
import csv
import time
import numpy as np
import pylab as plt


def parse_fname(fname):
    """"Extracts time and date from data's filename"""

    hour = fname[9:11]
    date = fname[:8]
    datenext = time.strftime('%Y%m%d',
                             time.localtime(time.mktime(time.strptime(date, '%Y%m%d')) + 24 * 3600.))

    return hour, date, datenext


class EcoHabData(object):
    """Reads in a folder with data from Eco-HAB"""

    def process_line_6(self,elements):
        """remove point from 2nd column of new data files"""
        return [elements[0],' '.join([elements[1].replace('.', ''), elements[2]]), elements[3], elements[4], elements[5]]
    
    def process_line_7(self,elements):
        """remove point from 2nd column of new data files"""
        return [elements[0],' '.join([elements[1].replace('.', ''), elements[2]]), elements[3], elements[4], elements[5], elements[6]]
    def process_line_5(self, elements, date):
        """Add date to data (old data files)"""
        elements[1] = ' '.join([date, elements[1]])
        return elements

    def read_file(self,fname):
        """Reads in data file"""
        hour, date, datenext = parse_fname(fname)

        f = open(os.path.join(self.path, fname),'r')

        for line in f:
            elements = line.split()

            if len(elements) == 6:
                self.rawdata.append(self.process_line_6(elements))
            elif len(elements) == 7:
                self.rawdata.append(self.process_line_7(elements))
            elif len(elements) == 5:
                if hour == '23' and elements[1][:2] == '00':
                    self.rawdata.append(self.process_line_5(elements,datenext))
                else:
                    self.rawdata.append(self.process_line_5(elements,date))
            else:
                    raise(IOError('Unknown data format in file %s' %f))

    def get_data(self):
        """Finds files in provided directory (path) and reads in data"""

        self._fnames =  filter(lambda x: x.endswith('0000.txt') or
                             (x[-12:-7] == '0000_' and x.endswith('.txt')),
                             os.listdir(self.path))
        
        for f_name in self._fnames:
            self.read_file(f_name)
            self.days.add(f_name.split('_')[0])

    def __init__(self, path, _ant_pos=None):
        self.days = set()
        self.path = path
        self.rawdata = []
        self.get_data()

        self.rawdata = self.remove_ghost_tags()
        
        self.mice = set([d[4] for d in self.rawdata])        
        self.rawdata.sort(key=lambda x: self.convert_time(x[1]))
  
        self.data = {}
        self.data['Id'] = [d[0] for d in self.rawdata]
        self.data['Time'] = [self.convert_time(d[1]) for d in self.rawdata]
        if _ant_pos is None:
            self._ant_pos = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8}
        else:
            self._ant_pos = _ant_pos
        self.data['Antenna'] = [self._ant_pos[d[2]] for d in self.rawdata]
        self.data['Tag'] = [d[4] for d in self.rawdata]
        # Masking - by default do not mask
        self.mask = None 
        self._mask_slice = None
        
        
    def remove_ghost_tags(self, how_many_appearances=20,factor=3):
        new_data = []
        ghost_mice = []
        counters = {}
        dates = {}
        for d in self.rawdata:
            mouse = d[4]
            if mouse not in counters:
                counters[mouse] = 0
            if mouse not in dates:
                dates[mouse] = set()
                
            counters[mouse] += 1
            dates[mouse].add(d[1].split()[0])
        how_many_days = len(self.days)/factor
        for mouse in counters:
            if counters[mouse] < how_many_appearances or dates[mouse] <= how_many_days:
                ghost_mice.append(mouse)
   
        for d in self.rawdata:
            mouse = d[4]
            if mouse not in ghost_mice:
                new_data.append(d)

        return new_data[:]
                        
    def __repr__ (self):
        """Nice string representation for prtinting this class."""
        mystring = 'Eco-HAB data loaded from:\n%s\nin the folder%s\n' %(
                   self._fnames.__str__(), self.path) 
        return mystring

    def mask_data(self, starttime, endtime):
        """mask_data(starttime, endtime)
        All future queries will be clipped to the visits starting between
        starttime and endtime."""
        self.mask = (starttime, endtime) 
        arr = np.array(self.data['Time'])
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            self._mask_slice = (idcs[0], idcs[-1] + 1)
        elif len(idcs) == 1:
            self._mask_slice = (idcs[0], idcs[0] + 1)
        else:
            self._mask_slice = (0, 0)

    def unmask_data(self):
        """Remove the mask - future queries will not be clipped"""
        self.mask = None
        self._mask_slice = None

    def getproperty(self, mice, propname, astype=None):
        if isinstance(mice, (str, unicode)):
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

    def getantennas(self, mice):
        return self.getproperty(mice, 'Antenna')
                                                  
    def gettimes(self, mice): 
        return self.getproperty(mice, 'Time', 'float')
                            
    @staticmethod
    def convert_time(s): 
        """Convert date and time to seconds since epoch""" 
        return (time.mktime(time.strptime(s[:-4], '%Y%m%d %H:%M:%S'))
                    + float(s[-3:])/1000)

    def checkData_one_mouse(self,mouse):
        antennas = self.getantennas(mouse)
        times  = self.gettimes(mouse)
        wrong_times = []
        for i,next_antenna in enumerate(antennas[1:]):
            if abs(next_antenna - antennas[i]) not in [0,1]:
                if next_antenna !=8 and antennas[i] != 8:
                    wrong_times.append(times[i])
        
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
        if isinstance(mice, (str, unicode)):
            mice = [mice]
        # if not isinstance(mice, collections.Container):
        #     mice = [mice]
 
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
        
        
        
class EcoHabSessions9states(IEcoHabSession):
    """Calculates 'visits' to Eco-HAB compartments.
    state: 0 -- visit too fast
    state: 1 -- pipe between antennae 1 and 2, chamber 1 and 2
    state: 2 -- chamber 2 (antennae 2 and 3)
    state: 3 -- pipe between antennae 3 and 4, chamber 2 and 3
    state: 4 -- chamber 3 (antennae 4 and 5)
    state: 5 -- pipe between antennae 5 and 6, chamber 3 and 4
    state: 6 -- chamber 4 (antennae 6 and 7)
    state: 7 -- pipe between antennae 7 and 8, chamber 4 and 1
    state: 8 -- chamber 1 (antennae 8 and 1)

    """
    

    def __init__(self, ehd, **kwargs):
        self._ehd = ehd
        self.mask = None
        self._mask_slice = None
        self.shortest_session_threshold = kwargs.get('shortest_session_threshold', 2)
        self.fs = kwargs.get('fs',10)
        
        self.t_start_exp = np.min(self._ehd.data['Time'])
        self.t_end_exp = np.max(self._ehd.data['Time'])
        
        t = np.arange(self.t_start_exp,self.t_end_exp,1/self.fs)
        self.signal_data = np.zeros((len(t),len(ehd.mice)),dtype =np.int8)      
        tempdata, self.statistics = self.get_statistics()

            
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
        
    def time_to_experiment_beginning(self,t):
        return (t - self.t_start_exp)*self.fs
    
    def get_statistics(self):
        
        tempdata = []
        statistics = {}
        
        for n, mm in enumerate(self._ehd.mice):
            tt = self._ehd.gettimes(mm)
            an = self._ehd.getantennas(mm)

            statistics[mm] = {}
            statistics[mm]["state_freq"]= np.zeros(9)
            statistics[mm]["state_time"]= [[] for i in range(9)]
            statistics[mm]["preference"]={}
            
            for i in range(9):
                 statistics[mm]["preference"][i] = np.zeros(2)
                 
            previous = 0
            for tstart, tend, anstart, anend in zip(tt[:-1], tt[1:], an[:-1], an[1:]):
                
                t_diff = tend - tstart

                if t_diff < self.shortest_session_threshold:
                    state = 0
                    previous = 0
                    statistics[mm]["state_freq"][state]+=1
                    statistics[mm]["state_time"][state].append(t_diff)
                    continue

                antennae_diff = np.abs(anstart - anend)
            
                if antennae_diff in [0,1,7]:
                    
                    s = self.time_to_experiment_beginning(tstart)
                    e = self.time_to_experiment_beginning(tend)

                    ############Most obvious reading##########
                    if antennae_diff in [1, 7]:
                        if antennae_diff == 1:
                            state = int((anstart + anend)/2-0.5)
                        else:
                            state = 8
                        #Save state to stadard and signal data
                        tempdata.append((state, mm, tstart, tend, t_diff,
                                         True))
                        self.signal_data[int(s):int(e),n] = state
                                         
                        statistics[mm]["state_freq"][state]+=1
                        statistics[mm]["state_time"][state].append(t_diff)
                        
                        if anend == state:
                            statistics[mm]["preference"][state][1]+=1 #counterclockwise
                        else:
                            statistics[mm]["preference"][state][0]+=1  #clockwise
                        previous = state
                        
                    elif antennae_diff == 0 and previous != 0:

                        if t_diff < 2:
                            continue

                        antennae_diff2 = abs(anstart - previous)
          
                        if antennae_diff2 == 0 and previous != 1:
                            state = previous-1
                        elif antennae_diff2 == 0 and previous == 1:
                            state = 8
                        elif antennae_diff2 == 1:
                            state = previous+1
                            if state > 8:
                                state = 1
                        elif antennae_diff2 == 7:
                            state = 1
                            
                            
                        #print previous,state,anend
                        #Save state to stadard and signal data
                        tempdata.append((state, mm, tstart, tend, t_diff,
                                         True))
                        self.signal_data[int(s):int(e),n] = state
                                         
                        statistics[mm]["state_freq"][state]+=1
                        statistics[mm]["state_time"][state].append(tend - tstart)
                        if anend == state:
                            statistics[mm]["preference"][state][1]+=1
                        else:
                            statistics[mm]["preference"][state][0]+=1  
                        previous = state       

        return tempdata, statistics         
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
        if isinstance(mice, (str, unicode)):
            mice = [mice]
        # if not isinstance(mice, collections.Container):
        #     mice = [mice]
 
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
        totv = [0, 0, 0, 0] #total number of visits
        tott = [0., 0., 0., 0.] #total number of time
        for idx, ad in enumerate([1, 2, 3, 4]):
            durs = [x for x, y in zip(durations, adds) if y == ad]
            totv[idx] = len(durs)
            tott[idx] = sum(durs)
        return totv, tott
        
if __name__ == '__main__':
    eco = EcoHabData('/home/jszmek/EcoHAB/eh')
    sessions = EcoHabSessions9states(eco)
    for mouse in eco.mice:
        print(sessions.statistics[mouse]['state_freq'])
   
