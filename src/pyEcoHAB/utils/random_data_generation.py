# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from . import general as utils
from . import BaseFunctions

class PseudoLoader(object):
    def __init__(self, data, setup_config):
        self.registrations = BaseFunctions.Data(data, None)
        self.mice = self.get_mice()
        self.setup_config = setup_config
        self.directions = setup_config.directions
        self.session_start = sorted(self.get_times(self.mice))[0]
        self.session_end = sorted(self.get_times(self.mice))[-1]

    def get_mice(self):
        mouse_list = list(set(self.registrations.data["Tag"]))
        # new Eco-HAB has a different mouse tag naming convention
        # last five digits are the same whereas in previous version
        # there was a prefix and first digits where the same
        if len(set([mouse[-4:] for mouse in mouse_list])) == len(mouse_list):
            mouse_dict = {mouse[-6:]: mouse for mouse in mouse_list}
            mouse_keys = sorted(mouse_dict.keys())
            return [mouse_dict[mouse] for mouse in mouse_keys]
        return sorted(mouse_list)

    def mask_data(self, start_time, end_time):
        """
        Hide registrations and visits in ranges (self.session_start, start_time)
        and (end_time, self.session_end).

        Args:
           start_time: float
           end_time: float
        """
        self.mask = (start_time, end_time)
        self.registrations.mask_data(self.mask)

    def unmask_data(self):
        """Remove the mask - future registrations and visits 
        queries will not be clipped"""
        self.mask = None
        self.registrations.unmask_data()

    def get_antennas(self, mice):
        return self.registrations.getproperty(mice,
                                         'Antenna')

    def get_times(self, mice):
        return self.registrations.getproperty(mice,
                                         'Time',
                                         'float')

def get_shifts(mice_list):
    shift_dict = {}
    for mouse in mice_list:
        shift_dict[mouse] = np.random.uniform(-1800., 1800.)
    return shift_dict


def randomly_shift_data(data): 
    mice = sorted(set(data[:]["Tag"]))
    new_data = data.copy()
    shift_dict = get_shifts(mice)
    for i, line in enumerate(data):
        key = line[-1]
        new_data[i]["Time"] = line["Time"] + shift_dict[key]
    return new_data


def generate_surrogate_data(e_data, timeline, binsize, mice, N, func):
    out_data = []
    for i in range(N):
        new_data = randomly_shift_data(e_data.registrations.data)
        dataE = PseudoLoader(new_data, e_data.setup_config)
        phases,\
            total_time,\
            data,\
            data_keys = utils.get_registrations_bins(dataE, timeline,
                                                     binsize, mice,
                                                     func)
        out_data.append(data)
    return out_data

def reshape_surrogate_data(data):
    # data is a list of dictionaries, we need a dictionary of lists

    out_data = {}
    for x in data:
        # first dict
        for key1 in x.keys():
            if key1 not in out_data:
                out_data[key1] = {}
            for key2 in x[key1].keys():
                if key2 not in out_data[key1]:
                    out_data[key1][key2] = []
                out_data[key1][key2].append(x[key1][key2])
    return out_data
