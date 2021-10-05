# SPDX-License-Identifier: LGPL-2.1-or-later
from __future__ import print_function, division, absolute_import
import sys
import numpy as np


class DataBase(object):

    def __init__(self, data, mask):
        self.mask = None
        self._mask_slice = None
        self.data = data
        if mask:
            self._cut_out_data(mask)

    def _find_mask_indices(self, mask, column_name):
        arr = np.array(self.data[column_name])

        if len(mask) >= 2:
            starttime = mask[0]
            endtime = mask[-1]
        elif len(args) == 1:
            starttime = min(self.data[column_name])
            endtime = mask[0]
        else:
            return (0, len(arr) - 1)
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            return (idcs[0], idcs[-1] + 1)
        if len(idcs) == 1:
            return (idcs[0], idcs[0] + 1)
        return (0, 0)

    def mask_data(self, args, column_name):
        """mask_data(endtime) or mask_data(starttime, endtime)
        All future queries will be clipped to the visits starting between
        starttime and endtime."""
        arr = np.array(self.data[column_name])
        if isinstance(args, int) or isinstance(args, float):
            start = min(arr)
            end = args[0]
        elif len(args) >= 2:
            start = args[0]
            end = args[-1]
        elif len(args) == 1:
            start = min(arr)
            end = args[0]
        else:
            start = min(arr)
            end = max(arr)

        self.mask = (start, end)
        self._mask_slice = self. _find_mask_indices(self.mask,
                                                    column_name)

    def unmask_data(self):
        """Remove the mask - future queries will not be clipped"""
        self.mask = None
        self._mask_slice = None

    def _cut_out_data(self, new_mask):
        mask = self._find_mask_indices(new_mask)
        for key in self.data.keys():
            self.data[key] = self.data[key][mask[0]: mask[1]]

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
            mask_0, mask_1 = self._mask_slice[0], self._mask_slice[1]
            if astype is None:
                return [x[0] for x in zip(self.data[propname][mask_0:mask_1],
                                          self.data['Tag'][mask_0:mask_1])
                        if x[1] in mice]
            elif astype == 'float':
                return [float(x[0]) for x in zip(
                    self.data[propname][mask_0:mask_1],
                    self.data['Tag'][mask_0:mask_1])
                        if x[1] in mice]


class Data(DataBase):
    def __init__(self, data, mask):
        super(Data, self).__init__(data, mask)

    def get_antennas(self, mice):
        return self.getproperty(mice, 'Antenna')

    def get_times(self, mice):
        return self.getproperty(mice, 'Time', 'float')

    def get_durations(self, mice):
        return self.getproperty(mice, 'Duration')

    def mask_data(self, mask):
        super(Data, self).mask_data(mask, column_name="Time")


class Visits(DataBase):
    def __init__(self, data, mask):
        super(Visits, self).__init__(data, mask)

    def get_starttimes(self, mice):
        return self.getproperty(mice, 'AbsStartTimecode', 'float')

    def get_endtimes(self, mice):
        return self.getproperty(mice, 'AbsEndTimecode', 'float')

    def get_durations(self, mice):
        return self.getproperty(mice, 'VisitDuration', 'float')

    def get_visit_addresses(self, mice):
        return self.getproperty(mice, 'Address')

    def mask_data(self, mask):
        super(Visits, self).mask_data(mask, column_name="AbsStartTimecode")
