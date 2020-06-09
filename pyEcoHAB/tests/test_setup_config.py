from __future__ import print_function, division, absolute_import
import os
import unittest

from pyEcoHAB import SetupConfig 
from pyEcoHAB import data_path

SAME_PIPE = { "1": ["1", "2"],
             "2": ["1", "2"],
             "3": ["3", "4"],
             "4": ["3", "4"],
             "5": ["5", "6"],
             "6": ["5", "6"],
             "7": ["7", "8"],
             "8": ["7", "8"]}

SAME_ADDRESS = {
    "1": ["1", "8"],
    "2": ["2", "3"],
    "3": ["2", "3"],
    "4": ["4", "5"],
    "5": ["4", "5"],
    "6": ["6", "7"],
    "7": ["6", "7"],
    "8": ["1", "8"],
}

OPPOSITE_PIPE = {"1": ["5", "6"],
                 "2": ["5", "6"],
                 "3": ["7", "8"],
                 "4": ["7", "8"],
                 "5": ["1", "2"],
                 "6": ["1", "2"],
                 "7": ["3", "4"],
                 "8": ["3", "4"]}

ADDRESS = {"1": "cage A", #"4"
           "2": "cage B", #"1",
           "3": "cage B", #"1",
           "4": "cage C", #"2",
           "5": "cage C", #"2",
           "6": "cage D", #"3",
           "7": "cage D", #"3",
           "8": "cage A", #"4"
}

ADDRESS_NON_ADJACENT = {"1": "cage B", #"1",
                        "2": "cage A", #"4",
                        "3": "cage C", #"2",
                        "4": "cage B", #"1",
                        "5": "cage D", #"3",
                        "6": "cage C", #"2",
                        "7": "cage A", #"4",
                        "8": "cage D", #"3"
}
# Surrounding: difference between antennas only "2" or "6" -- skipped one antenna
SURROUNDING = {("1", "3"): "cage B", #"1",
               ("1", "7"): "cage A", #"4",
               ("2", "4"): "cage B", #"1",
               ("2", "8"): "cage A", #"4",
               ("3", "5"): "cage C", #"2",
               ("4", "6"): "cage C", #"2",
               ("5", "7"): "cage D", #"3",
               ("6", "8"): "cage D", #"3"
}


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
        correct = {}
        correct["cage B"] = ['2', '3']
        correct["cage C"] = ['4', '5']
        correct["cage D"] = ['6', '7']
        correct["cage A"] = ['8', '1']
        self.assertEqual(out, correct)


    def test_default_tunnels(self):
        out = self.default.get_tunnels_dict()
        correct = {}
        correct["tunnel 1"] = ['1', '2']
        correct["tunnel 2"] = ['3', '4']
        correct["tunnel 3"] = ['5', '6']
        correct["tunnel 4"] = ['7', '8']
        self.assertEqual(out, correct)

    def test_default_internal(self):
        out = self.default.internal_antennas
        self.assertEqual(out, [])

    def test_custom_internal(self):
        out = self.custom.internal_antennas
        self.assertEqual(out, ["cage B"])

    def test_same_tunnel_default(self):
        out = {}
        out["1"] = ["1", "2"]
        out["2"] = ["1", "2"]
        out["3"] = ["3", "4"]
        out["4"] = ["3", "4"]
        out["5"] = ["5", "6"]
        out["6"] = ["5", "6"]
        out["7"] = ["7", "8"]
        out["8"] = ["7", "8"]
        self.assertEqual(out, self.default.same_tunnel)

    def test_same_address_default(self):
        out = {}
        out["1"] = ["8", "1"]
        out["2"] = ["2", "3"]
        out["3"] = ["2", "3"]
        out["4"] = ["4", "5"]
        out["5"] = ["4", "5"]
        out["6"] = ["6", "7"]
        out["7"] = ["6", "7"]
        out["8"] = ["8", "1"]
        self.assertEqual(out, self.default.same_address)


class TestOppositePipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.default = SetupConfig()
        path = os.path.join(data_path, "test_setups")
        cls.custom = SetupConfig(path=path, fname="setup_internal.txt")

    def test(self):
        self.assertEqual(self.default.opposite_tunnel, OPPOSITE_PIPE)

    def test_empty(self):
        self.assertEqual(self.custom.opposite_tunnel, {})

    def test_entrance_antennas1(self):
        out = ["1", "2", "3", "4", "5", "6", "7", "8"]
        self.assertEqual(sorted(out), sorted(self.default.entrance_antennas))

    def test_entrance_antennas2(self):
        out = ["1", "2"]
        self.assertEqual(sorted(out), sorted(self.custom.entrance_antennas))

    def test_other_tunnel_antenna_d1(self):
        self.assertEqual(["2"], self.default.other_tunnel_antenna("1"))

    def test_other_tunnel_antenna_d2(self):
        self.assertEqual(["1"], self.default.other_tunnel_antenna("2"))

    def test_other_tunnel_antenna_d3(self):
        self.assertEqual(["4"], self.default.other_tunnel_antenna("3"))

    def test_other_tunnel_antenna_d4(self):
        self.assertEqual(["3"], self.default.other_tunnel_antenna("4"))

    def test_other_tunnel_antenna_d5(self):
        self.assertEqual(["6"], self.default.other_tunnel_antenna("5"))

    def test_other_tunnel_antenna_d6(self):
        self.assertEqual(["5"], self.default.other_tunnel_antenna("6"))

    def test_other_tunnel_antenna_d7(self):
        self.assertEqual(["8"], self.default.other_tunnel_antenna("7"))

    def test_other_tunnel_antenna_d8(self):
        self.assertEqual(["7"], self.default.other_tunnel_antenna("8"))

    def test_other_tunnel_antenna_c1(self):
        self.assertEqual(["2"], self.custom.other_tunnel_antenna("1"))

    def test_other_tunnel_antenna_c2(self):
        self.assertEqual(["1"], self.custom.other_tunnel_antenna("2"))

    def test_other_cage_antenna_d1(self):
        self.assertEqual(["8"], self.default.other_cage_antenna("1"))

    def test_other_cage_antenna_d2(self):
        self.assertEqual(["3"], self.default.other_cage_antenna("2"))

    def test_other_cage_antenna_d3(self):
        self.assertEqual(["2"], self.default.other_cage_antenna("3"))

    def test_other_cage_antenna_d4(self):
        self.assertEqual(["5"], self.default.other_cage_antenna("4"))

    def test_other_cage_antenna_d5(self):
        self.assertEqual(["4"], self.default.other_cage_antenna("5"))

    def test_other_cage_antenna_d6(self):
        self.assertEqual(["7"], self.default.other_cage_antenna("6"))

    def test_other_cage_antenna_d7(self):
        self.assertEqual(["6"], self.default.other_cage_antenna("7"))

    def test_other_cage_antenna_d8(self):
        self.assertEqual(["1"], self.default.other_cage_antenna("8"))

    def test_other_cage_antenna_c1(self):
        self.assertEqual([], self.custom.other_cage_antenna("1"))

    def test_other_cage_antenna_c2(self):
        self.assertEqual([], self.custom.other_cage_antenna("2"))

    def test_next_tunnel_antennas_c1(self):
        self.assertEqual([], self.custom.next_tunnel_antennas("1"))

    def test_next_tunnel_antennsa_c2(self):
        self.assertEqual([], self.custom.next_tunnel_antennas("2"))

    def test_next_tunnel_antennas_d1(self):
        self.assertEqual(sorted(["4", "3", "7", "8"]),
                         self.default.next_tunnel_antennas("1"))

    def test_next_tunnel_antennas_d2(self):
        self.assertEqual(sorted(["4", "3", "7", "8"]),
                         self.default.next_tunnel_antennas("2"))

    def test_next_tunnel_antennas_d3(self):
        self.assertEqual(sorted(["1", "2", "5", "6"]),
                         self.default.next_tunnel_antennas("3"))

    def test_next_tunnel_antennas_d4(self):
        self.assertEqual(sorted(["1", "2", "5", "6"]),
                         self.default.next_tunnel_antennas("4"))

    def test_next_tunnel_antennas_d5(self):
        self.assertEqual(sorted(["3", "4", "7", "8"]),
                         self.default.next_tunnel_antennas("5"))

    def test_next_tunnel_antennas_d6(self):
        self.assertEqual(sorted(["3", "4", "7", "8"]),
                         self.default.next_tunnel_antennas("6"))

    def test_next_tunnel_antennas_d7(self):
        self.assertEqual(sorted(["1", "2", "5", "6"]),
                         self.default.next_tunnel_antennas("7"))

    def test_next_tunnel_antennas_d8(self):
        self.assertEqual(sorted(["1", "2", "5", "6"]),
                         self.default.next_tunnel_antennas("8"))

    def test_cage_address_dict_default(self):
        self.assertEqual(ADDRESS, self.default.get_cage_address_dict())

    def test_cage_address_dict_custom(self):
        out = {"1": "cage A",
               "2": "cage B"}
        self.assertEqual(out, self.custom.get_cage_address_dict())

    def test_cage_adjacent_default(self):
        self.assertEqual(ADDRESS_NON_ADJACENT,
                         self.default.get_address_non_adjacent_dict())

    def test_cage_adjacent_custom(self):
        out = {"1": "cage B",
               "2": "cage A"}
        self.assertEqual(out,
                         self.custom.get_address_non_adjacent_dict())


if __name__ == '__main__':
    unittest.main()
