from __future__ import print_function, division

import os
import time
import numpy as np
import sys
max_break = 60*60

class Data(object):
    
    def __init__(self,data,mask):
        self.mask = None
        self._mask_slice = None
        self.data = data
        if self.mask:
            self._cut_out_data(mask)
            
    def _cut_out_data(self,new_mask):
        mask = self._find_mask_indices(new_mask)
        self.data['Id'] = self.data['Id'][mask[0]:mask[1]]
        self.data['Time'] = self.data['Time'][mask[0]:mask[1]]
        
        self.data['Antenna'] = self.data['Antenna'][mask[0]:mask[1]]
        self.data['Tag'] = self.data['Tag'][mask[0]:mask[1]]
        
    def _find_mask_indices(self,mask):
        
        starttime, endtime = mask
        arr = np.array(self.data['Time'])
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]

        if len(idcs) >= 2:
            return (idcs[0], idcs[-1] + 1)
        if len(idcs) == 1:
            return (idcs[0], idcs[0] + 1)
        
        return (0,0)
    
    def mask_data(self, starttime, endtime):
        """mask_data(starttime, endtime)
        All future queries will be clipped to the visits starting between
        starttime and endtime."""
        self.mask = (starttime, endtime) 
        self._mask_slice = self._find_mask_indices(self.mask)

    def unmask_data(self):
        """Remove the mask - future queries will not be clipped"""
        self.mask = None
        self._mask_slice = None

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

    def getantennas(self, mice):
        return self.getproperty(mice, 'Antenna')
                                                  
    def gettimes(self, mice): 
        return self.getproperty(mice, 'Time', 'float')
    
def parse_fname(fname):
    """"Extracts time and date from data's filename"""

    hour = fname[9:11]
    date = fname[:8]
    datenext = time.strftime('%Y%m%d',
                             time.localtime(time.mktime(time.strptime(date, '%Y%m%d')) + 24 * 3600.))

    return hour, date, datenext

class EcoHabData(Data):
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
                    self.rawdata.append(self.process_line_5(elements, datenext))
                else:
                    self.rawdata.append(self.process_line_5(elements, date))
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
        
    def check_antenna_presence(self):
        t_start = self.data['Time'][0]
        all_times = np.array(self.data['Time'])
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

    def __init__(self,**kwargs):# path, _ant_pos=None,mask=None):
        self.days = set()
        self.path = kwargs.pop("path")
        self.rawdata = []
        self.get_data()
        self.max_break = max_break
        how_many_appearances = kwargs.pop('how_many_appearances',1000)
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])
        self.rawdata = self.remove_ghost_tags(how_many_appearances,factor,tags=tags)
        self.mice = list(set([d[4] for d in self.rawdata]))
        self.rawdata.sort(key=lambda x: self.convert_time(x[1]))
        _ant_pos = kwargs.pop('_ant_pos',None)
        mask = kwargs.pop('mask',None)
        
        if _ant_pos is None:
            self._ant_pos = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8}
        else:
            self._ant_pos = _ant_pos
        
        #self.mask = None 
        #self._mask_slice = None

        data = {}
        data['Time'] = [self.convert_time(d[1]) for d in self.rawdata]
        data['Id'] = [d[0] for d in self.rawdata]
        
        data['Antenna'] = [self._ant_pos[d[2]] for d in self.rawdata]
        data['Tag'] = [d[4] for d in self.rawdata]
        super(EcoHabData,self).__init__(data, mask)
        antenna_breaks = self.check_antenna_presence()
        if antenna_breaks:
            print('Antenna not working')
            for antenna in antenna_breaks:
                print(antenna,antenna_breaks[antenna])
        self.antenna_mismatch()
        
        if mask:
            self._cut_out_data(mask)
          
    def merge_experiment(self, other):
        assert self.mice == other.mice
        self.mask = None
        self._mask_slice = None
        t_mine = min(self.data['Time'])
        t_other = min(other.data['Time'])
        
        if t_other < t_mine:
            assert max(other.data['Time']) < t_mine
            new_data = {}
            for key in other.data:
                new_data[key] = other_data[key][:]
                new_data[key].append(self.data['key'])
            self.data = new_data
        elif t_mine < t_other:
            assert max(self.data['Time']) < t_other
            for key in self.data:
                self.data[key].append(other.data[key])
                
    def remove_ghost_tags(self, how_many_appearances,factor,tags=[]):
        new_data = []
        ghost_mice = []
        counters = {}
        dates = {}
        if len(tags):
            for tag in tags:
                ghost_mice.append(tag)
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
            if counters[mouse] < how_many_appearances or len(dates[mouse]) <= how_many_days:
                if mouse not in ghost_mice:
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

    @staticmethod
    def convert_time(s): 
        """Convert date and time to seconds since epoch"""
        return (time.mktime(time.strptime(s.split('.')[0], '%Y%m%d %H:%M:%S'))  + float(s.split('.')[-1])/1000.)
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
    eco = EcoHabData('/home/jszmek/EcoHAB/eh')
    sessions = EcoHabSessions(eco)
    i = 0
    for mouse in eco.mice:
        totv,tott = sessions.getstats(mouse)
        if totv == [0,0,0,0] or tott == [0,0,0,0]:
            
            i = i+1
        print(mouse,totv, tott, sum(totv))

