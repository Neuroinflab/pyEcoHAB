from __future__ import print_function, division, absolute_import

import os
import glob
import time
import numpy as np
import sys
from DataBase import Data
from Nodes import AntennaReadOut, Animal
import utils
from collections import Container, Counter
from operator import methodcaller, attrgetter

class EcoHabData(Data):
    """Reads in data from a directory with EcoHAB data"""
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

        # to exclude fake tags aka ghost mice
        how_many_appearances = kwargs.pop('how_many_appearances',100) 
        factor = kwargs.pop('factor',2)
        #if you want to remove mice from analysis -- posssibly remove
        tags = kwargs.pop('remove_mice',[])
        max_break = kwargs.pop('max_break', 60*60)
        raw_data = self.read_in_data(self.path, how_many_appearances, factor, tags)

        self.add_mice(raw_data)
        readouts = self.extract_readouts(raw_data, source=self.path)
        self.add_readouts(readouts)
       
        antenna_breaks = self.check_antenna_presence(max_break)
        # if antenna_breaks:
        #     print('Antenna not working')
        #     for antenna in antenna_breaks:
        #         print(antenna,antenna_breaks[antenna])
        # self.antenna_mismatch()
        
        # if mask:
        #     self._cut_out_data(mask)


    @staticmethod    
    def remove_ghost_tags(data, how_many_appearances, factor, tags=[]):
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

    def add_mice(self, raw_data):
        mice = self.get_all_mice(raw_data)
        for mouse in mice:
            mouse_obj = Animal(mouse)
            self.add_animal(mouse_obj)
            
    def get_all_mice(self, raw_data):
        return set(np.array(raw_data)[:,-1])
    
    def read_in_data(self, path, how_many_appearances, factor, tags):
        files = self.find_files(path)
        raw_data = self.append_data(files)
        return self.remove_ghost_tags(raw_data, how_many_appearances, factor, tags)
        
    @staticmethod
    def extract_readouts(data, source=None):
        readouts = []
        for i, data_line in enumerate(data):
            time_sec = utils.convert_time(data_line[1])
            readouts += [AntennaReadOut(EventId=int(data_line[0]),
                                        Start=time_sec,
                                        Date=data_line[1],
                                        AntennaId=int(data_line[2]),
                                        Duration=float(data_line[3]),
                                        Animal=Animal(data_line[4]))]
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
        
    def check_antenna_presence(self, max_break):
        t_start = self.get_start()
        all_times = self.readouts.get_attributes('Start')
        breaks = {}
        
        for antenna in range(1, 9):
            
            antenna_idx = np.where(np.array(self.data['Antenna']) == antenna)[0]
            times = all_times[antenna_idx] 
            breaks[antenna] = []
            if len(times):
                if times[0] - t_start > max_break:
                    breaks[antenna].append([0, np.round(times[0])])
        
                intervals = times[1:]-times[0:-1]

                where_breaks = np.where(intervals > max_break)[0]

                if len(where_breaks):
                    for i in where_breaks:
                        breaks[antenna].append([np.round(times[i]), np.round(times[i+1])])
            else:
                breaks[antenna].append([np.round(t_start), self.data['Time'][-1]])         
                    
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
