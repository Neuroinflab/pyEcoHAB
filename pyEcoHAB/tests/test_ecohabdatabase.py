from __future__ import print_function, division, absolute_import
import os
import unittest

import pyEcoHAB.utils.for_loading as uf
import pyEcoHAB.utility_functions as utils
from pyEcoHAB import data_path
from pyEcoHAB.Loader import EcoHabDataBase, Loader
from pyEcoHAB import Timeline
from pyEcoHAB.SetupConfig import SetupConfig


class TestSingleAntennaStats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        setup_config = SetupConfig()
        cls.data = EcoHabDataBase(data, None, 2, setup_config)
        cls.config = Timeline(path)

    def test_antenna1(self):
        times = self.config.get_time_from_epoch("ALL")
        result = self.data.get_registration_stats("mouse_1", times[0],
                                                  times[1], "1", 3600)
        self.assertEqual(result, ([3], [613/1000]))

    def test_antenna2(self):
        times = self.config.get_time_from_epoch("ALL")
        result = self.data.get_registration_stats("mouse_1", times[0],
                                                  times[1], "1", 1800)
        self.assertEqual(result, ([2, 1], [460/1000, 153/1000]))

    def test_antenna8(self):
        times = self.config.get_time_from_epoch("ALL")
        result = self.data.get_registration_stats("mouse_1", times[0],
                                                  times[1], "8", 900)
        self.assertEqual(result, ([0, 0, 0, 1], [0, 0, 0, 1026/1000]))


if __name__ == '__main__':
    unittest.main()
