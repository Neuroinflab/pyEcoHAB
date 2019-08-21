from __future__ import print_function, division, absolute_import
import os
import sys
import numpy as np
from . import utility_functions as utils
max_break = 60*60

class DataBase(object):
    
    def __init__(self, data, mask):
        self.mask = None
        self._mask_slice = None
        self.data = data
        if self.mask:
            self._cut_out_data(mask)
            
    def _cut_out_data(self, new_mask):
        mask = self._find_mask_indices(new_mask)
        for key in self.data.keys():
            self.data[key] = self.data[key][mask[0]:mask[1]]
    
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


class Data(DataBase):
    def __init__(self, data, mask):
        super(Data, self).__init__(data, mask)

    def _find_mask_indices(self, mask):
        starttime, endtime = mask
        arr = np.array(self.data['Time'])
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            return (idcs[0], idcs[-1] + 1)
        if len(idcs) == 1:
            return (idcs[0], idcs[0] + 1)
        return (0, 0)

    def getantennas(self, mice):
        return self.getproperty(mice, 'Antenna')
                                                  
    def gettimes(self, mice): 
        return self.getproperty(mice, 'Time', 'float')

    def getdurations(self, mice):
        return self.getproperty(mice, 'Duration')

class Visits(DataBase):

    def __init__(self, data, mask):
        super(Visits, self).__init__(data, mask)

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
        arr = np.array(self.data['AbsStartTimecode'])
        idcs = np.where((arr >= start) & (arr < end))[0]
        if len(idcs) >= 2:
            self._mask_slice = (idcs[0], idcs[-1] +1)
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


class EcoHabData(Data):
    """Reads in EcoHAB data"""

    def __init__(self, **kwargs):# path, _ant_pos=None,mask=None):
        self.days = set()
        self.path = kwargs.pop("path")
        self.rawdata = []
        self.read_in_data()
        self.max_break = max_break
        how_many_appearances = kwargs.pop('how_many_appearances', 50)
        remove_antennas = kwargs.pop('remove_antennas', [])
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])
        self.rawdata = self.remove_ghost_tags(how_many_appearances,
                                              factor,
                                              tags=tags)
         
        self.mice = self.get_mice()
        self.rawdata.sort(key=lambda x: utils.time_to_sec(x[1]))
        _ant_pos = kwargs.pop('_ant_pos',None)
        mask = kwargs.pop('mask',None)
        
        if _ant_pos is None:
            self._ant_pos = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8}
        else:
            self._ant_pos = _ant_pos
 
        data = {}
        data['Time'] = [utils.time_to_sec(d[1]) for d in self.rawdata]
        data['Id'] = [d[0] for d in self.rawdata]
        data['Antenna'] = [self._ant_pos[d[2]] for d in self.rawdata]
        data['Tag'] = [d[4] for d in self.rawdata]
        data['Duration'] = [d[3] for d in self.rawdata]
        data = self._remove_antennas(data, remove_antennas)
        super(EcoHabData,self).__init__(data, mask)
        antenna_breaks = self.check_antenna_presence()
        if antenna_breaks:
            print('Antenna not working')
            for antenna in antenna_breaks:
                print(antenna, ':')
                for breaks in antenna_breaks[antenna]:
                    print(utils.print_human_time(breaks[0]),
                          utils.print_human_time(breaks[1]))
                    print((breaks[1] - breaks[0])/3600, 'h')
        self.antenna_mismatch()
        
        if mask:
            self._cut_out_data(mask)
        self.prefix = utils.make_prefix(self.path)
        self.res_dir = utils.results_path(self.path)

    def get_home_cage_antenna(self):
        """
        Finds home antenna. This is a function used to calculate one
        of the measures of dominance in two cage experiments, when mice
        are living in an environment 
        """
        antennas = []
        for mouse in self.mice:
            antennas.append(self.getantennas(mouse)[0])
        return max(set(antennas), key=antennas.count)

    @staticmethod
    def _remove_antennas(data, antennas=[]):
        keys = data.keys()
        if isinstance(antennas, int):
            antennas = [antennas]
        if isinstance(antennas, list):
            new_data = {key:[] for key in keys}
            for i, antenna in enumerate(data['Antenna']):
                if antenna not in antennas:
                    for key in keys:
                        new_data[key].append(data[key][i])
                else:
                    print([data[key][i] for key in keys])
            return new_data
        return data

    def read_single_file(self, fname):
        """Reads in single data file"""
        hour, date, datenext = utils.parse_fname(fname)
        f = open(os.path.join(self.path, fname),'r')
        for line in f:
            elements = line.split()
            if len(elements) == 5:
                if hour == '23' and elements[1][:2] == '00':
                    line = utils.process_line_5_elements(elements,
                                                         datenext)
                else:
                    line = utils.process_line_5_elements(elements,
                                                         date)
            elif len(elements) > 5:
                line = utils.process_line_more_elements(elements)
            else:
                raise(IOError('Unknown data format in file %s' %f))
            self.rawdata += [line]

    def read_in_data(self):
        """Finds files in provided directory (path) and reads in data"""
        self._fnames = utils.get_filenames(self.path)
        for f_name in self._fnames:
            self.read_single_file(f_name)
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

    def check_single_mouse_data(self,mouse):
        antennas = self.getantennas(mouse)
        times  = self.gettimes(mouse)
        error_crossing_times = []
        for i,next_antenna in enumerate(antennas[1:]):
            if abs(next_antenna - antennas[i]) not in [0, 1]:
                if next_antenna !=8 and antennas[i] != 8:
                    error_crossing_times.append(times[i])
        return error_crossing_times
    


class EcoHabSessions(Visits):
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
        self.threshold = kwargs.pop('shortest_session_threshold', 2.)
        temp_data = self._calculate_states()
        data = {}
        data['Address'] = [x[0] for x in temp_data]
        data['Tag'] = [x[1] for x in temp_data]
        data['AbsStartTimecode'] = [x[2] for x in temp_data]
        data['AbsEndTimecode'] = [x[3] for x in temp_data]
        data['VisitDuration'] = [x[4] for x in temp_data]
        data['ValidVisitSolution'] = [x[5] for x in temp_data]
        super(EcoHabSessions,self).__init__(data, None)
        self.mice = self._ehd.mice
        self.prefix = self._ehd.prefix
        self.res_dir = self._ehd.res_dir
        


    
