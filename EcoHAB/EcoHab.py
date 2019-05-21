from __future__ import print_function, division

import os
import time
import numpy as np
import sys
import utility_functions as utils
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
        self.data['Duration'] = self.data['Duration'][mask[0]:mask[1]]

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

    def getdurations(self, mice):
        return self.getproperty(mice, 'Duration')

def parse_fname(fname):
    """"Extracts time and date from data's filename"""

    hour = fname[9:11]
    date = fname[:8]
    datenext = time.strftime('%Y%m%d',
                             time.localtime(time.mktime(time.strptime(date, '%Y%m%d')) + 24 * 3600.))

    return hour, date, datenext

class EcoHabData(Data):
    """Reads in a folder with data from Eco-HAB"""
    def get_home_cage_antenna(self):
        antennas = []
        for mouse in self.mice:
            antennas.append(self.getantennas(mouse)[0])
        return max(set(antennas), key=antennas.count)

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
        #     printantenna, breaks[antenna])
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

    def get_mice(self):
        mouse_list = list(set([d[4] for d in self.rawdata]))
        mouse_dict = {mouse[-4:]:mouse for mouse in mouse_list}
        mouse_keys = sorted(mouse_dict.keys())
        return [mouse_dict[mouse] for mouse in mouse_keys]

    def __init__(self, **kwargs):# path, _ant_pos=None,mask=None):
        self.days = set()
        self.path = kwargs.pop("path")
        self.rawdata = []
        self.get_data()
        self.max_break = max_break
        how_many_appearances = kwargs.pop('how_many_appearances', 50)
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])
        self.rawdata = self.remove_ghost_tags(how_many_appearances,
                                              factor,
                                              tags=tags)
         
        self.mice = self.get_mice()
        self.rawdata.sort(key=lambda x: self.convert_time(x[1]))
        _ant_pos = kwargs.pop('_ant_pos',None)
        mask = kwargs.pop('mask',None)
        
        if _ant_pos is None:
            self._ant_pos = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8}
        else:
            self._ant_pos = _ant_pos
 
        data = {}
        data['Time'] = [self.convert_time(d[1]) for d in self.rawdata]
        data['Id'] = [d[0] for d in self.rawdata]
        data['Antenna'] = [self._ant_pos[d[2]] for d in self.rawdata]
        data['Tag'] = [d[4] for d in self.rawdata]
        data['Duration'] = [d[3] for d in self.rawdata]
        super(EcoHabData,self).__init__(data, mask)
        antenna_breaks = self.check_antenna_presence()
        if antenna_breaks:
            print('Antenna not working')
            for antenna in antenna_breaks:
                print(antenna, ':')
                for breaks in antenna_breaks[antenna]:
                    print(self.print_time_human(breaks[0]),  self.print_time_human(breaks[1]))
                    print((breaks[1] - breaks[0])/3600, 'h')
        self.antenna_mismatch()
        
        if mask:
            self._cut_out_data(mask)
        self.prefix = utils.make_prefix(self.path)
        self.res_dir = utils.results_path(self.path)
          
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
                
    def remove_ghost_tags(self, how_many_appearances, factor, tags=[]):
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
    def print_time_human(tt):
        """convert seconds to date and time since epoch """
        st = time.localtime(tt)
        return time.asctime(st)
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
    """Calculates 'visits' to Eco-HAB chambers."""

    def _calculate_states(self):
        tempdata = []
        for mouse in self._ehd.mice:
            times, antennas = utils.get_times_antennas(self._ehd,
                                                       mouse,
                                                       0,
                                                       -1)
            tempdata.extend(utils.get_states_for_ehs(times, antennas,
                                                     mouse,
                                                     self.threshold))
        tempdata.sort(key=lambda x: x[2])
        return tempdata

    def __init__(self, ehd, **kwargs):
        
        self._ehd = ehd
        self.mask = None
        self._mask_slice = None
        self.threshold = kwargs.pop('shortest_session_threshold', 2.)
        temp_data = self._calculate_states()
        self.data = {}
        self.data['Address'] = [x[0] for x in temp_data]
        self.data['Tag'] = [x[1] for x in temp_data]
        self.data['AbsStartTimecode'] = [x[2] for x in temp_data]
        self.data['AbsEndTimecode'] = [x[3] for x in temp_data]
        self.data['VisitDuration'] = [x[4] for x in temp_data]
        self.data['ValidVisitSolution'] = [x[5] for x in temp_data]
        self.mice = self._ehd.mice
        
    def unmask_data(self):
        """Remove the mask - future queries will not be clipped"""
        self.mask = None
        self._mask_slice = None

    def mask_data(self, *args):
        """mask_data(endtime) or mask_data(starttime, endtime)
        All future queries will be clipped to the visits starting between
        starttime and endtime."""
        try:
            start = args[0]
            end = args[1]
        except IndexError:   
            start = min(self.getstarttimes(self._ehd.mice))
            end = args[0]
        self.mask = (start, end)
        idcs = utils.get_idx_between(start, end, self.data['AbsStartTimecode'])
        if len(idcs) >= 2:
            self._mask_slice = (idcs[0], idcs[-1])
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
        # starts = self.getproperty(mice, 'AbsStartTimecode', 'float')
        # ends = self.getproperty(mice, 'AbsEndTimecode', 'float')
        # return [abs(x[1] - x[0]) for x in zip(starts, ends)]
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
    

