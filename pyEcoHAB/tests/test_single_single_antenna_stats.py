from __future__ import print_function, division, absolute_import

import unittest
import os
from pyEcoHAB import single_antenna_registrations as sar
from pyEcoHAB import utility_functions as uf
from pyEcoHAB import Loader
from pyEcoHAB import Timeline
from pyEcoHAB import data_path

class TestExecution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sample_data = os.path.join(data_path, "weird_short")
        cls.data = Loader(sample_data)
        cls.config = Timeline(sample_data)

    def test1(self):
        sar.get_single_antenna_stats(self.data, self.config, 3600, 1)

    def test2(self):
        sar.get_single_antenna_stats(self.data, self.config, 3600)

    def test3(self):
        self.assertRaises(Exception, sar.get_single_antenna_stats,
                          self.data, self.config, 3600, "gugu")

    def test4(self):
        sar.get_single_antenna_stats(self.data, self.config, 900, 1)

    def test5(self):
        sar.get_single_antenna_stats(self.data, self.config, 900)

    def test6(self):
        sar.get_single_antenna_stats(self.data, self.config, 1800, 1)

    def test7(self):
        sar.get_single_antenna_stats(self.data, self.config, 1800)
    
if __name__ == '__main__':
    unittest.main()
