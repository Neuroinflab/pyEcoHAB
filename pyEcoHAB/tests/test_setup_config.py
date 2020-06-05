from __future__ import print_function, division, absolute_import
import os
import unittest
from collections import OrderedDict

from pyEcoHAB import SetupConfig 
from pyEcoHAB import data_path

class TestReadingIn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = SetupConfig()
        path1 = os.path.join(data_path, "test_setups")
        cls.d_path = SetupConfig(path=path1)
        path2 = os.path.join(data_path, "test_setups_2")
        cls.d_2 = SetupConfig(path=path2)
        cls.c = SetupConfig(path=path2, fname="setup2.txt")

    def test_read_in_default_path(self):
        self.assertEqual(self.d.path, data_path)

    def test_read_in_default_fname(self):
        self.assertEqual(self.d.fname, "standard_setup.txt")

    def test_read_in_path(self):
        path =  os.path.join(data_path, "test_setups")
        self.assertEqual(self.d_path.path, path)

    def test_read_in_path_fname(self):
        self.assertEqual(self.d_path.fname, "setup.txt")

    def test_raise_no_setup(self):
         path1 = os.path.join(data_path, "test_setup")
         self.assertRaises(Exception, SetupConfig, path=path1)

    def test_read_in_path2(self):
        path =  os.path.join(data_path, "test_setups_2")
        self.assertEqual(self.d_2.path, path)

    def test_read_in_path_fname2(self):
        self.assertEqual(self.d_2.fname, "setup2.txt")

    def test_read_in_custom_path(self):
        path =  os.path.join(data_path, "test_setups_2")
        self.assertEqual(self.c.path, path)

    def test_read_in_custom_fname(self):
        self.assertEqual(self.c.fname, "setup2.txt")


class TestGetCagesTunnels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.default = SetupConfig()

    def test_standard_cages(self):
        out = sorted(["cage A", "cage B", "cage C", "cage D"])
        res = sorted(self.default.get_cages())
        self.assertEqual(out, res)
   
    def test_standard_tunnels(self):
        out = sorted(["tunnel 1", "tunnel 2", "tunnel 3", "tunnel 4"])
        res = sorted(self.default.get_tunnels())
        self.assertEqual(out, res)


class TestGetDicts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.default = SetupConfig()
        path = os.path.join(data_path, "test_setups")
        cls.custom = SetupConfig(path=path, fname="setup_internal.txt")

    def test_default_cages(self):
        out = self.default.get_cages_dict()
        correct = OrderedDict()
        correct["cage A"] = [2, 3]
        correct["cage B"] = [4, 5]
        correct["cage C"] = [6, 7]
        correct["cage D"] = [8, 1]
        self.assertEqual(out, correct)


    def test_default_tunnels(self):
        out = self.default.get_tunnels_dict()
        correct = OrderedDict()
        correct["tunnel 1"] = [1, 2]
        correct["tunnel 2"] = [3, 4]
        correct["tunnel 3"] = [5, 6]
        correct["tunnel 4"] = [7, 8]
        self.assertEqual(out, correct)

    def test_default_internal(self):
        out = self.default.get_compartments_with_additional_antennas()
        self.assertEqual(out, [])

    def test_custom_internal(self):
        out = self.custom.get_compartments_with_additional_antennas()
        self.assertEqual(out, ["cage B"])
        
if __name__ == '__main__':
    unittest.main()
