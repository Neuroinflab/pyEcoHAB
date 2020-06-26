from __future__ import print_function, division, absolute_import
import os
import unittest
import numpy as np
import pyEcoHAB.utils.for_loading as uf
import pyEcoHAB.utility_functions as utils
from pyEcoHAB import data_path
from pyEcoHAB.SetupConfig import SetupConfig
from pyEcoHAB import Loader, Merger

class TestLoader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path1 = os.path.join(data_path, "modular_1",
                                 "data_setup_additional")
        cls.dataset1 = Loader(cls.path1, visit_threshold=1.5, prefix="gugu")
        cls.setup1 = SetupConfig(cls.path1)
        cls.dataset1_standard = Loader(cls.path1)

    def test_path(self):
        self.assertEqual(self.path1, self.dataset1.path)

    def test_visit_threshold(self):
        self.assertEqual(self.dataset1.visit_threshold, 1.5)
    
    def test_prefix(self):
        self.assertEqual(self.dataset1_standard.prefix, "WT_F_test_1_")

    def test_prefix_2(self):
        self.assertEqual(self.dataset1.prefix, "gugu")

if __name__ == '__main__':
    unittest.main()

