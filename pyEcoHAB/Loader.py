from __future__ import print_function, division, absolute_import
import os
import sys
from collections import OrderedDict

try:
    basestring
except NameError:
    basestring = str

import numpy as np

from . import BaseFunctions
from . import utility_functions as utils
from .utils import for_loading as ufl

class EcoHabDataBase(object):

    def __init__(self, data, mask, threshold):
        self.readings = BaseFunctions.Data(data, mask)
        self.threshold = threshold
        self.mice = self.get_mice()
        self.visits = self._calculate_visits()
        self.session_start = sorted(self.get_times(self.mice))[0]
        self.session_end = sorted(self.get_times(self.mice))[-1]

    def _calculate_animal_positions(self):
        """
        Calculate timings of animal visits to Eco-HAB compartments.
        """
        tempdata = []
        for mouse in self.mice:
            times, antennas = utils.get_times_antennas(self.readings,
                                                       mouse,
                                                       0, -1)
            tempdata.extend(utils.get_animal_position(times, antennas,
                                                      mouse,
                                                      self.threshold))
        tempdata.sort(key=lambda x: x[2])
        return tempdata

    def _calculate_visits(self):
        temp_data = self._calculate_animal_positions()
        data = {}
        data['Address'] = [x[0] for x in temp_data]
        data['Tag'] = [x[1] for x in temp_data]
        data['AbsStartTimecode'] = [x[2] for x in temp_data]
        data['AbsEndTimecode'] = [x[3] for x in temp_data]
        data['VisitDuration'] = [x[4] for x in temp_data]
        data['ValidVisitSolution'] = [x[5] for x in temp_data]
        return BaseFunctions.Visits(data, None)

    def mask_data(self, starttime, endtime):
        self.mask = (starttime, endtime)
        self.readings.mask_data(self.mask)
        self.visits.mask_data(self.mask)

    def unmask_data(self):
        """Remove the mask - future queries will not be clipped"""
        self.mask = None
        self.readings.unmask_data()
        self.visits.unmask_data()

    def get_antennas(self, mice):
        return self.readings.getproperty(mice,
                                         'Antenna')

    def get_times(self, mice):
        return self.readings.getproperty(mice,
                                         'Time',
                                         'float')
    def get_durations(self, mice):
        """Return duration of registration
        by antenna"""
        return self.readings.getproperty(mice,
                                         'Duration',
                                         'float')
    #add get_visits, get_readings
    def get_visit_addresses(self, mice):
        return self.visits.getproperty(mice,
                                       'Address')
    def get_starttimes(self, mice):
        return self.visits.getproperty(mice,
                                       'AbsStartTimecode',
                                       'float')

    def get_endtimes(self, mice):
        return self.visits.getproperty(mice,
                                       'AbsEndTimecode',
                                       'float')

    def get_visit_durations(self, mice):
        return self.visits.getproperty(mice,
                                       'VisitDuration',
                                       'float')
    def how_many_antennas(self):
        all_antennas = set(self.get_antennas(self.mice))
        return len(all_antennas)

    def get_mice(self):
        mouse_list = list(set(self.readings.data["Tag"]))
        #new EcoHAB has a different tag nameing convention
        #last five digits are the same whereas in previous version
        #there was a prefix and first digits where the same
        if len(set([mouse[-4:] for mouse in mouse_list])) == len(mouse_list):
            mouse_dict = {mouse[-6:]:mouse for mouse in mouse_list}
            mouse_keys = sorted(mouse_dict.keys())
            return [mouse_dict[mouse] for mouse in mouse_keys]
        return sorted(mouse_list)

    def get_home_cage_antenna(self):
        """
        Finds home antenna. This is a function used to calculate one
        of the measures of dominance in two cage experiments. 
        """
        antennas = []
        for mouse in self.mice:
            antennas.append(self.get_antennas(mouse)[0])
        return max(set(antennas), key=antennas.count)

    def get_visits(self, mice=None, t_start=None, t_end=None):
        """
        Return a list of visits to Eco-HAB compartments. Each visit is 
        a named dictionary with following fields: t_start, t_end, tag
        """
        if isinstance(mice, str):
            if mice in self.mice:
                mice = [mice]
            else:
                print("Could not find animal %s" % mice)
                return []
        if mice is None:
            mice = self.get_mice()
        if t_start is None:
            t_start = self.session_start
        if t_end is None:
            t_end = self.session_end

        self.visits.mask_data([t_start, t_end])
        out = []
        for mouse in mice:
            addresses = self.get_visit_addresses(mouse)
            start_times = self.get_starttimes(mouse)
            end_times = self.get_endtimes(mouse)
            durations = self.get_durations(mouse)
            for i, a in enumerate(addresses):
                visit = ufl.NamedDict("Visit_%s_%d" % (mouse, i),
                                      tag=mouse, address=a,
                                      t_start=start_times[i],
                                      t_end=end_times[i],
                                      duration=durations[i])
                out.append(visit)
        return sorted(out, key = lambda o: o["t_start"])


class Loader(EcoHabDataBase):
    """Reads in Eco-HAB data files that are located in path.

    This class reads in data collected by the Eco-HAB system, parses them
    and removes in-correct registrations. After loading the data Loader triggers
    calculation of timings of animal visits to Eco-HAB compartments. Currently
    Loader assumes that there are 4 Eco-HAB compartments denoted by A, B, C, D).

    Loader converts date and time of registration to float using
    time.localtime()

    Args:
        path: string
           directory containing Eco-HAB data

    Keyword Args:
        antenna_positions: dictionary
           a dictionary specifing conversion of registered antenna ids
           to integers
           int(antenna_id) is the default conversion
        mask: list or tuple of floats
           Loader will read in data registed between mask[0] and mask[1].
           mask[0] and mask[1] need to be expressed seconds from the epoch,
           since Loader converts animal tag registration times to seconds
           since the epoch. By default the whole data is saved by Loader.
        visit_threshold: float
           visits shorter than visit_threshold will be rejected
           Default value is 2 s (parameter based on mouse behavior)
        res_dir: string
           results path directory.
           By default results will be saved in path/Results
        prefix: string
           a prefix (string) added to all generated result files.
           By default an info.txt file in path directory is parsed and added to
           all filenames of results files. If no prefix is provided and path
           directory does not contain an info.txt file, no prefix is added.
        max_break: float
           breaks in antenna registrations longer than max_break
           will be reported, while loading Eco-HAB data.
        how_many_appearances: int
           Animal tags that are registered less times than how_many_appearances
           will be removed from loaded data. By default an animal tag must be
           registered at least 200 times not to be removed from loaded data
        min_appearance_factor: float of value less than 1
           Animal tags that are registered in fraction of the experiment
           duration lower than min_appearance_factor will be removed
           from loaded data. By default an animal tag must be registered
           during at least half of the experiment duration.
        remove_antennas: list
           Registrations by antenna ids in remove_antennas will be removed from
           loaded data. By default Loader keeps all the registrations.
        remove_mice: list
           Animal tag registrations to be removed from loaded data. By default
           no registrations are removed.
    """
    STANDARD_ANTENNAS = {'1': 1, '2': 2,
                         '3': 3, '4': 4,
                         '5': 5, '6': 6,
                         '7': 7, '8': 8}
    MAX_BREAK = 3600

    def __init__(self, path, **kwargs):
        #Read in parameters
        self.path = path
        antenna_positions = kwargs.pop('antenna_positions', None)

        if antenna_positions is None:
            self.antenna_positions = self.STANDARD_ANTENNAS
        else:
            self.antenna_positions = antenna_positions

        self.mask = kwargs.pop('mask', None)
        self.visit_threshold = kwargs.pop('visit_threshold', 2.)
        self.res_dir = kwargs.pop("res_dir",
                                  ufl.results_path(self.path))
        self.prefix = ufl.make_prefix(self.path)
        self.max_break = kwargs.pop("max_break", self.MAX_BREAK)
        how_many_appearances = kwargs.pop('how_many_appearances', 50)
        min_appearance_factor = kwargs.pop('min_appearance_factor', 0.5)
        remove_antennas = kwargs.pop('remove_antennas', [])
        factor = 1/min_appearance_factor
        tags = kwargs.pop('remove_mice',[])

        #Read in data
        rawdata = self._read_in_raw_data(factor,
                                         how_many_appearances,
                                         tags)
        data = ufl.from_raw_data(rawdata,
                                 self.antenna_positions)
        data = ufl.remove_antennas(data, remove_antennas)
        #As in antenna readings
        
        ufl.run_diagnostics(data, self.max_break)
        super(Loader, self).__init__(data, self.mask,
                                     self.visit_threshold)
        self.cages = self.get_cages()

    def get_cages(self):
        return sorted(list(set(self.get_visit_addresses(self.mice))))

    def _read_in_raw_data(self, factor, how_many_appearances, tags):
        """Reads in data from files in self.path.
        Removes ghost tags from data"""
        raw_data = []
        days = set()
        self._fnames = ufl.get_filenames(self.path)
        if not len(self._fnames ):
            sys.exit("%s is empty"% self.path)
        for f_name in self._fnames:
            raw_data += ufl.read_single_file(self.path, f_name)
            days.add(f_name.split('_')[0])
        how_many_days = len(days)/factor
        data = ufl.remove_ghost_tags(raw_data,
                                       how_many_appearances,
                                       how_many_days,
                                       tags=tags)
        data.sort(key=lambda x: ufl.time_to_sec(x[1]))
        return data
                         
    def __repr__ (self):
        """Nice string representation for prtinting this class."""
        mystring = 'Eco-HAB data loaded from:\n%s\nin the folder%s\n' %(
                   self._fnames.__str__(), self.path) 
        return mystring

    def check_single_mouse_data(self, mouse):
        antennas = self.data.get_antennas(mouse)
        times  = self.data.get_times(mouse)
        error_crossing_times = []
        for i, next_antenna in enumerate(antennas[1:]):
            if abs(next_antenna - antennas[i]) not in [0, 1]:
                if next_antenna != 8 and antennas[i] != 8:
                    error_crossing_times.append(times[i])
        return error_crossing_times



class Merger(EcoHabDataBase):
    def __init__(self, *data_sources, **kwargs):
        self._source_list = map(str, data_sources)
        self._ignore_mice_diff = kwargs.pop('ignore_mice_differences',
                                            False)
        self.session_start = None
        self.session_end = None

        grouped_sources = self._group_data_sources(data_sources)
        for key in sorted(grouped_sources.keys()):
            data_source = grouped_sources[key]
            try:
                self._append_data_source(data_source)

            except:
                print("ERROR processing {}".format(data_source))
                raise
        self.mice = self.get_mice()
        self.session_start = sorted(self.get_times(self.mice))[0]
        self.session_end = sorted(self.get_times(self.mice))[-1]
        self.res_dir = kwargs.pop("results_path",
                                  data_sources[0].res_dir  + "_merged")
        self.prefix = kwargs.pop("results_path",
                                 data_sources[0].prefix + "_merged")

    def _group_data_sources(self, data_sources):
        grouped = {}
        for data_source in data_sources:
            key = data_source.session_start
            grouped[key] = data_source
        return grouped


    def _append_data_source(self, new_data):
        if self.session_start is None and self.session_end is None:
            super(Merger, self).__init__(new_data.readings.data, None,
                                         new_data.visit_threshold)

        if new_data.session_start != None:
            if self.session_start != None and self.session_end != None:
                if self.session_start < new_data.session_start\
                   and self.session_end > new_data.session_start:
                    print('Overlap of EcoHAB sessions detected')

        if new_data.session_end != None:
            if self.session_start != None and self.session_end != None:
                if self.session_start < new_data.session_end\
                   and self.session_end > new_data.session_end:
                    print('Overlap of EcoHAN sessions detected')

        for key in new_data.readings.data.keys():
            self.readings.data[key].extend(new_data.readings.data[key])
        for key in new_data.visits.data.keys():
            self.visits.data[key].extend(new_data.visits.data[key])
