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
        cls.dataset1 = Loader(cls.path1)

    def test_path(self):
        self.assertEqual(self.path1, self.dataset1.path)


if __name__ == '__main__':
    unittest.main()

