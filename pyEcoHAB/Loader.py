from __future__ import print_function, division, absolute_import
import os
import numpy as np
from . import BaseFunctions
from . import utility_functions as utils


class Loader(object):
    """Reads in EcoHAB data"""
    STANDARD_ANTENNAS = {'1': 1, '2': 2,
                         '3': 3, '4': 4,
                         '5': 5, '6': 6,
                         '7': 7, '8': 8}
    MAX_BREAK = 3600

    def __init__(self, path, **kwargs):
        #Read in parameters
        self.path = path
        _ant_pos = kwargs.pop('antenna_positions', None)

        if _ant_pos is None:
            self._ant_pos = self.STANDARD_ANTENNAS
        else:
            self._ant_pos = _ant_pos

        mask = kwargs.pop('mask', None)
        self.threshold = kwargs.pop('antenna_threshold', 2.)
        self.res_dir = kwargs.pop("results_path",
                                  utils.results_path(self.path))
        self.prefix = utils.make_prefix(self.path)
        self.max_break = kwargs.pop("max_break", self.MAX_BREAK)
        how_many_appearances = kwargs.pop('how_many_appearances', 50)
        remove_antennas = kwargs.pop('remove_antennas', [])
        factor = kwargs.pop('factor',2)
        tags = kwargs.pop('remove_mice',[])

        #Read in data
        rawdata = self._read_in_raw_data(factor,
                                        how_many_appearances,
                                        tags)
        data = self._from_raw_data(rawdata,
                                 self._ant_pos,
                                 remove_antennas)
        #As in antenna readings
        self.readings = BaseFunctions.Data(data, mask)
        self.mice = self.get_mice()
        self.run_diagnostics()
        self.visits = self._calculate_visits()

    def get_home_cage_antenna(self):
        """
        Finds home antenna. This is a function used to calculate one
        of the measures of dominance in two cage experiments, when mice
        are living in an environment 
        """
        antennas = []
        for mouse in self.mice:
            antennas.append(self.readings.data.getantennas(mouse)[0])
        return max(set(antennas), key=antennas.count)

    @staticmethod
    def _remove_antennas(data, antennas):
        if isinstance(antennas, int):
            antennas = [antennas]
        if isinstance(antennas, list):
            keys = data.keys()
            new_data = {key:[] for key in keys}
            for i, antenna in enumerate(data['Antenna']):
                if antenna not in antennas:
                    for key in keys:
                        new_data[key].append(data[key][i])
                else:
                    print("Removing record",
                          [data[key][i] for key in keys])
            return new_data
        return data

    def _read_single_file(self, fname):
        """Reads in single data file"""
        hour, date, datenext = utils.parse_fname(fname)
        raw_data = []
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
            raw_data += [line]
        return raw_data

    @staticmethod
    def _remove_ghost_tags(raw_data,
                           how_many_appearances,
                           how_many_days,
                           tags=[]):
        new_data = []
        ghost_mice = []
        counters = {}
        dates = {}
        if len(tags):
            for tag in tags:
                ghost_mice.append(tag)
        print(len(raw_data))
        for d in raw_data:
            mouse = d[4]
            if mouse not in counters:
                counters[mouse] = 0
            if mouse not in dates:
                dates[mouse] = set()
            counters[mouse] += 1
            dates[mouse].add(d[1].split()[0])

        for mouse in counters:
            if counters[mouse] < how_many_appearances or len(dates[mouse]) <= how_many_days:
                if mouse not in ghost_mice:
                    ghost_mice.append(mouse)
        for d in raw_data:
            mouse = d[4]
            if mouse not in ghost_mice:
                new_data.append(d)

        return new_data[:]

    def _read_in_raw_data(self, factor,
                          how_many_appearances, tags):
        """Reads in data from files in self.path.
        Removes ghost tags from data"""
        raw_data = []
        days = set()
        self._fnames = utils.get_filenames(self.path)
        for f_name in self._fnames:
            raw_data += self._read_single_file(f_name)
            days.add(f_name.split('_')[0])
        how_many_days = len(days)/factor
        data = self._remove_ghost_tags(raw_data,
                                       how_many_appearances,
                                       how_many_days,
                                       tags=tags)
        return data

    def _from_raw_data(self, raw_data, antenna_pos,
                       remove_antennas=[]):
        data = {}
        data['Time'] = [utils.time_to_sec(d[1]) for d in raw_data]
        data['Id'] = [d[0] for d in raw_data]
        data['Antenna'] = [antenna_pos[d[2]] for d in raw_data]
        data['Tag'] = [d[4] for d in raw_data]
        data['Duration'] = [d[3] for d in raw_data]

        return self._remove_antennas(data, remove_antennas)

    def run_diagnostics(self):
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

    def check_antenna_presence(self):
        t_start = self.readings.data['Time'][0]
        all_times = np.array(self.readings.data['Time'])
        breaks = {}
        
        for antenna in range(1, 9):
            antenna_idx = np.where(np.array(self.readings.data['Antenna']) == antenna)[0]
            times = all_times[antenna_idx] 
            breaks[antenna] = []
            if len(times):
                if times[0] - t_start > self.max_break:
                    breaks[antenna].append([0, np.round(times[0])])
                intervals = times[1:] - times[0:-1]
                where_breaks = np.where(intervals > self.max_break)[0]
                if len(where_breaks):
                    for i in where_breaks:
                        breaks[antenna].append([np.round(times[i]), np.round(times[i+1])])
            else:
                breaks[antenna].append([np.round(t_start),
                                        self.readings.data['Time'][-1]])
        return breaks
    
    def antenna_mismatch(self):
        
        t_start = self.readings.data['Time'][0]
        all_times = np.array(self.readings.data['Time'])
        weird_transit = [[], []]
        for mouse in self.mice:
            mouse_idx = np.where(np.array(self.readings.data['Tag']) == mouse)[0]
            times = all_times[mouse_idx]
            antennas = np.array(self.readings.data['Antenna'])[mouse_idx]
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
            print(pair, mismatches[pair],
                  np.round(100*mismatches[pair]/len(self.readings.data['Antenna'])),
                  'per 100')
        return weird_transit

    def get_mice(self):
        mouse_list = list(set(self.readings.data["Tag"]))
        mouse_dict = {mouse[-4:]:mouse for mouse in mouse_list}
        mouse_keys = sorted(mouse_dict.keys())
        return [mouse_dict[mouse] for mouse in mouse_keys]
                        
    def __repr__ (self):
        """Nice string representation for prtinting this class."""
        mystring = 'Eco-HAB data loaded from:\n%s\nin the folder%s\n' %(
                   self._fnames.__str__(), self.path) 
        return mystring

    def check_single_mouse_data(self,mouse):
        antennas = self.data.getantennas(mouse)
        times  = self.data.gettimes(mouse)
        error_crossing_times = []
        for i, next_antenna in enumerate(antennas[1:]):
            if abs(next_antenna - antennas[i]) not in [0, 1]:
                if next_antenna != 8 and antennas[i] != 8:
                    error_crossing_times.append(times[i])
        return error_crossing_times

    def _calculate_states(self):
        tempdata = []
        for mouse in self.mice:
            times, antennas = utils.get_times_antennas(self.readings,
                                                       mouse,
                                                       0, -1)
            tempdata.extend(utils.get_states_for_ehs(times, antennas,
                                                     mouse,
                                                     self.threshold))
        tempdata.sort(key=lambda x: x[2])
        return tempdata

    def _calculate_visits(self):
        temp_data = self._calculate_states()
        data = {}
        data['Address'] = [x[0] for x in temp_data]
        data['Tag'] = [x[1] for x in temp_data]
        data['AbsStartTimecode'] = [x[2] for x in temp_data]
        data['AbsEndTimecode'] = [x[3] for x in temp_data]
        data['VisitDuration'] = [x[4] for x in temp_data]
        data['ValidVisitSolution'] = [x[5] for x in temp_data]
        return BaseFunctions.Visits(data, None)

    #add get_visits, get_readings
