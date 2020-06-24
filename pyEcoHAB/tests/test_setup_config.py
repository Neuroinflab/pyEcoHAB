from __future__ import print_function, division, absolute_import
import os
import unittest

from pyEcoHAB import SetupConfig, ExperimentSetupConfig
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

ADDRESS = {"1": "cage A", #4
           "2": "cage B", #1,
           "3": "cage B", #1,
           "4": "cage C", #2,
           "5": "cage C", #2,
           "6": "cage D", #3,
           "7": "cage D", #3,
           "8": "cage A", #4
}

ADDRESS_NON_ADJACENT = {"1": "cage B", #1,
                        "2": "cage A", #4,
                        "3": "cage C", #2,
                        "4": "cage B", #1,
                        "5": "cage D", #3,
                        "6": "cage C", #2,
                        "7": "cage A", #4,
                        "8": "cage D", #3
}
# Surrounding: difference between antennas only 2 or 6 -- skipped one antenna
SURROUNDING = {("1", "3"): "cage B", #1,
               ("1", "7"): "cage A", #4,
               ("2", "4"): "cage B", #1,
               ("2", "8"): "cage A", #4,
               ("3", "5"): "cage C", #2,
               ("4", "6"): "cage C", #2,
               ("5", "7"): "cage D", #3,
               ("6", "8"): "cage D", #3
}

KEYS = ['12', '21', '34', '43', '56', '65', '78', '87']
PAIRS = ["1 3", "1 4", "1 5", "1 6", "1 7", "2 4", "2 5", "2 6", "2 7", "2 8",
         "3 5", "3 6", "3 7", "3 8", "4 6", "4 7", "4 8", "5 7", "5 8", "6 8"]


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
        correct["cage B"] = ["2", "3"]
        correct["cage C"] = ["4", "5"]
        correct["cage D"] = ["6", "7"]
        correct["cage A"] = ["8", "1"]
        self.assertEqual(out, correct)


    def test_default_tunnels(self):
        out = self.default.get_tunnels_dict()
        correct = {}
        correct["tunnel 1"] = ["1", "2"]
        correct["tunnel 2"] = ["3", "4"]
        correct["tunnel 3"] = ["5", "6"]
        correct["tunnel 4"] = ["7", "8"]
        self.assertEqual(out, correct)

    def test_default_internal(self):
        out = self.default.internal_antennas
        self.assertEqual(out, [])

    def test_custom_internal(self):
        out = self.custom.internal_antennas
        self.assertEqual(out, ["8"])

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
        self.assertEqual(["8"], self.custom.other_cage_antenna("2"))

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

    def test_cage_address_default(self):
        self.assertEqual(ADDRESS, self.default.get_cage_address_dict())

    def test_cage_address_custom(self):
        out = {"1": "cage A",
               "2": "cage B",
               "8": "cage B"}
        self.assertEqual(out, self.custom.get_cage_address_dict())

    def test_cage_adjacent_default(self):
        self.assertEqual(ADDRESS_NON_ADJACENT,
                         self.default.get_address_non_adjacent_dict())

    def test_cage_adjacent_custom(self):
        out = {"1": "cage B",
               "2": "cage A"}
        self.assertEqual(out,
                         self.custom.get_address_non_adjacent_dict())

    def test_cage_surrounding_default(self):
        self.assertEqual(SURROUNDING,
                         self.default.get_surrounding_dict())

    def test_cage_surrounding_custom(self):
        self.assertEqual({("1", "8"): 'cage B'},
                         self.custom.get_surrounding_dict())

    def test_directions_custom(self):
        self.assertEqual(["12", "21"],
                         self.custom.get_directions_dict())

    def test_directions_default(self):
        self.assertEqual(KEYS,
                         self.default.get_directions_dict())

    def test_unused_antennas_d(self):
        self.assertEqual([],
                         self.default.find_unused_antennas())

    def test_unused_antennas_custom(self):
        self.assertEqual(['3', '4', '5', '6', '7'],
                         self.custom.find_unused_antennas())

    def test_mismatch_pairs_default(self):
        self.assertEqual(PAIRS,
                         self.default.get_mismatched_pairs())

    def test_mismatch_pairs_custom(self):
        self.assertEqual(["1 8"],
                         self.custom.get_mismatched_pairs())

    def test_all_antennas_custom(self):
        self.assertEqual(self.custom.all_antennas, ["1", "2", "8"])

    def test_all_antennas_default(self):
        self.assertEqual(self.default.all_antennas, ["1", "2", "3", "4",
                                                     "5", "6", "7", "8"])


class TestExperimentSetupConfig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path2 = os.path.join(data_path, "test_setups")
        path3 = os.path.join(data_path, "experiment_setup.txt")
        cls.config1 = SetupConfig(data_path, "standard_setup.txt")
        cls.config2 = SetupConfig(path=path2, fname="setup_internal.txt")
        cls.experiment_config = ExperimentSetupConfig(path3,
                                                      ecohab1=cls.config1,
                                                      ecohab2=cls.config2)

        path1 = os.path.join(data_path, "test_experiment_setups")
        cls.config3 = SetupConfig(path1, "setup1.txt")
        cls.config4 = SetupConfig(path1, "setup2.txt")
        full_path  = os.path.join(path1, "experiment_setup.txt")
        cls.full_exp = ExperimentSetupConfig(full_path,
                                             ecohab_1=cls.config3,
                                             ecohab_2=cls.config4)


    def test_indentity_points(self):
        out = {"ecohab1 cage A": "shared cage 1",
               "ecohab2 cage B": "shared cage 1"}
        self.assertEqual(self.experiment_config.identity_points,
                         out)

    def test_indentity_points_full_exp(self):
        out = {"ecohab_1 cage A": "cage A",
               "ecohab_2 cage A": "cage A",
               "ecohab_1 cage C": "cage C",
               "ecohab_2 cage C": "cage C",}
        self.assertEqual(self.full_exp.identity_points,
                         out)

    def test_all_section_names(self):
        out = sorted(["shared cage 1", "ecohab1 cage B", "ecohab1 cage C",
                      "ecohab1 cage D", "ecohab1 tunnel 1", "ecohab1 tunnel 2",
                      "ecohab1 tunnel 3", "ecohab1 tunnel 4", "ecohab2 cage A",
                      "ecohab2 tunnel 1"])

        self.assertEqual(out, sorted(self.experiment_config.sections()))

    def test_all_section_names_full_exp(self):
        out = sorted(["cage A", "ecohab_1 cage B", "ecohab_2 cage D",
                      "cage C", "ecohab_1 tunnel 1", "ecohab_1 tunnel 2",
                      "ecohab_2 tunnel 1", "ecohab_2 tunnel 2"])

        self.assertEqual(out, sorted(self.full_exp.sections()))

    def test_all_antennas(self):
        self.assertEqual(sorted(self.experiment_config.ALL_ANTENNAS),
                         sorted(['1_ecohab1', '1_ecohab2', '2_ecohab1',
                                 '2_ecohab2', '3_ecohab1', '4_ecohab1',
                                 '5_ecohab1', '6_ecohab1', '7_ecohab1',
                                 '8_ecohab1', '8_ecohab2']))

    def test_all_antennas_full_exp(self):
        self.assertEqual(sorted(self.full_exp.ALL_ANTENNAS),
                         sorted(['1_ecohab_1', '2_ecohab_1',
                                 '3_ecohab_1', '4_ecohab_1',
                                 '5_ecohab_2', '6_ecohab_2',
                                 '7_ecohab_2', '8_ecohab_2']))

    def test_all_cages(self):
        self.assertEqual(sorted(self.experiment_config.cages),
                         sorted(["shared cage 1", "ecohab1 cage B",
                                 "ecohab1 cage C", "ecohab1 cage D",
                                 "ecohab2 cage A"]))

    def test_all_cages_full_exp(self):
        self.assertEqual(sorted(self.full_exp.cages),
                         sorted(["cage A", "cage C",
                                 "ecohab_1 cage B", "ecohab_2 cage D"]))

    def test_all_tunnels(self):
        self.assertEqual(sorted(["ecohab1 tunnel 1", "ecohab1 tunnel 2",
                                 "ecohab1 tunnel 3", "ecohab1 tunnel 4",
                                 "ecohab2 tunnel 1"]),
                         sorted(self.experiment_config.tunnels))

    def test_all_tunnels_full_exp(self):
        self.assertEqual(sorted(["ecohab_1 tunnel 1", "ecohab_1 tunnel 2",
                                 "ecohab_2 tunnel 1", "ecohab_2 tunnel 2"]),
                         sorted(self.full_exp.tunnels))

    def test_get_cages_dict_1(self):
        key =  "shared cage 1"
        out = sorted(["8_ecohab1", "8_ecohab2","1_ecohab1", "2_ecohab2"])
        self.assertEqual(sorted(self.experiment_config.cages_dict[key]),
                         out)

    def test_get_cages_dict_2(self):
        key =  "ecohab1 cage B"
        out = sorted(["2_ecohab1", "3_ecohab1"])
        self.assertEqual(sorted(self.experiment_config.cages_dict[key]),
                         out)

    def test_get_cages_dict_3(self):
        key =  "ecohab1 cage C"
        out = sorted(["4_ecohab1", "5_ecohab1"])
        self.assertEqual(sorted(self.experiment_config.cages_dict[key]),
                         out)

    def test_get_cages_dict_4(self):
        key =  "ecohab1 cage D"
        out = sorted(["6_ecohab1", "7_ecohab1"])
        self.assertEqual(sorted(self.experiment_config.cages_dict[key]),
                         out)

    def test_get_cages_dict_5(self):
        key =  "ecohab2 cage A"
        out = sorted(["1_ecohab2"])
        self.assertEqual(sorted(self.experiment_config.cages_dict[key]),
                         out)

    def test_get_cages_dict_keys(self):
        keys = sorted(["shared cage 1", "ecohab1 cage B", "ecohab1 cage C",
                       "ecohab1 cage D", "ecohab2 cage A"])
        self.assertEqual(sorted(self.experiment_config.cages_dict.keys()),
                         keys)

    def test_get_cages_dict_full_exp_1(self):
        key =  "cage A"
        out = sorted(["1_ecohab_1", "8_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.cages_dict[key]),
                         out)

    def test_get_cages_dict_full_exp_2(self):
        key =  "ecohab_1 cage B"
        out = sorted(["2_ecohab_1", "3_ecohab_1"])
        self.assertEqual(sorted(self.full_exp.cages_dict[key]),
                         out)

    def test_get_cages_dict_full_exp_3(self):
        key =  "ecohab_2 cage D"
        out = sorted(["6_ecohab_2", "7_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.cages_dict[key]),
                         out)

    def test_get_cages_dict_full_exp_4(self):
        key =  "cage C"
        out = sorted(["4_ecohab_1", "5_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.cages_dict[key]),
                         out)

    def test_get_cages_dict_full_exp_keys(self):
        keys =  sorted(["cage A", "ecohab_1 cage B", "cage C",
                       "ecohab_2 cage D"])

        self.assertEqual(sorted(self.full_exp.cages_dict.keys()),
                         keys)

    def test_internal_antennas(self):
        self.assertEqual(self.experiment_config.internal_antennas,
                         ["8_ecohab2"])

    def test_internal_antennas_full_exp(self):
        self.assertEqual(self.full_exp.internal_antennas,
                         [])

    def test_tunnels_dict_keys(self):
        keys = sorted(["ecohab1 tunnel 1", "ecohab1 tunnel 2",
                       "ecohab1 tunnel 3", "ecohab1 tunnel 4",
                       "ecohab2 tunnel 1"])
        act_keys = sorted(self.experiment_config.tunnels_dict.keys())
        self.assertEqual(keys, act_keys)

    def test_tunnels_dict_1(self):
        key = "ecohab1 tunnel 1"
        out = sorted(["1_ecohab1", "2_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.tunnels_dict[key]))

    def test_tunnels_dict_2(self):
        key = "ecohab1 tunnel 2"
        out = sorted(["3_ecohab1", "4_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.tunnels_dict[key]))

    def test_tunnels_dict_3(self):
        key = "ecohab1 tunnel 3"
        out = sorted(["5_ecohab1", "6_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.tunnels_dict[key]))

    def test_tunnels_dict_4(self):
        key = "ecohab1 tunnel 4"
        out = sorted(["7_ecohab1", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.tunnels_dict[key]))

    def test_tunnels_dict_5(self):
        key = "ecohab2 tunnel 1"
        out = sorted(["1_ecohab2", "2_ecohab2"])
        self.assertEqual(out,
                         sorted(self.experiment_config.tunnels_dict[key]))

    def test_tunnels_dict_keys_full_exp(self):
        keys = sorted(["ecohab_1 tunnel 1", "ecohab_1 tunnel 2",
                       "ecohab_2 tunnel 1", "ecohab_2 tunnel 2"])
        act_keys = sorted(self.full_exp.tunnels_dict.keys())
        self.assertEqual(keys, act_keys)

    def test_tunnels_dict_full_exp_1(self):
        key = "ecohab_1 tunnel 1"
        out = sorted(["1_ecohab_1", "2_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.tunnels_dict[key]))

    def test_tunnels_dict_full_exp_2(self):
        key = "ecohab_1 tunnel 2"
        out = sorted(["3_ecohab_1", "4_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.tunnels_dict[key]))

    def test_tunnels_dict_full_exp_3(self):
        key = "ecohab_2 tunnel 2"
        out = sorted(["5_ecohab_2", "6_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.tunnels_dict[key]))

    def test_tunnels_dict_full_exp_4(self):
        key = "ecohab_2 tunnel 1"
        out = sorted(["7_ecohab_2", "8_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.tunnels_dict[key]))

    def test_same_tunnel_keys(self):
        keys = sorted(set(self.experiment_config.ALL_ANTENNAS)
                           - set(self.experiment_config.internal_antennas))
        self.assertEqual(keys,
                         sorted(self.experiment_config.same_tunnel.keys()))

    def test_same_tunnel_1(self):
        key = "1_ecohab1"
        out = sorted(["1_ecohab1", "2_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_2(self):
        key = "2_ecohab1"
        out = sorted(["1_ecohab1", "2_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_3(self):
        key = "3_ecohab1"
        out = sorted(["3_ecohab1", "4_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_4(self):
        key = "4_ecohab1"
        out = sorted(["3_ecohab1", "4_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_5(self):
        key = "5_ecohab1"
        out = sorted(["5_ecohab1", "6_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_6(self):
        key = "6_ecohab1"
        out = sorted(["5_ecohab1", "6_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_7(self):
        key = "7_ecohab1"
        out = sorted(["7_ecohab1", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_8(self):
        key = "8_ecohab1"
        out = sorted(["7_ecohab1", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_9(self):
        key = "1_ecohab2"
        out = sorted(["1_ecohab2", "2_ecohab2"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_10(self):
        key = "2_ecohab2"
        out = sorted(["1_ecohab2", "2_ecohab2"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_tunnel[key]))

    def test_same_tunnel_keys_full_exp(self):
        keys = sorted(set(self.full_exp.ALL_ANTENNAS)
                           - set(self.full_exp.internal_antennas))
        self.assertEqual(keys,
                         sorted(self.full_exp.same_tunnel.keys()))

    def test_same_tunnel_full_exp_1(self):
        key = "1_ecohab_1"
        out = sorted(["1_ecohab_1", "2_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_2(self):
        key = "2_ecohab_1"
        out = sorted(["1_ecohab_1", "2_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_3(self):
        key = "3_ecohab_1"
        out = sorted(["3_ecohab_1", "4_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_4(self):
        key = "4_ecohab_1"
        out = sorted(["3_ecohab_1", "4_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_5(self):
        key = "5_ecohab_2"
        out = sorted(["5_ecohab_2", "6_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_6(self):
        key = "6_ecohab_2"
        out = sorted(["5_ecohab_2", "6_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_7(self):
        key = "7_ecohab_2"
        out = sorted(["7_ecohab_2", "8_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_tunnel_full_exp_8(self):
        key = "8_ecohab_2"
        out = sorted(["7_ecohab_2", "8_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_tunnel[key]))

    def test_same_address_keys(self):
        self.assertEqual(sorted(self.experiment_config.all_antennas),
                         sorted(self.experiment_config.same_address.keys()))

    def test_same_address_11(self):
        key = "1_ecohab1"
        out = sorted(["1_ecohab1", "8_ecohab2", "2_ecohab2", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_81(self):
        key = "8_ecohab1"
        out = sorted(["1_ecohab1", "8_ecohab2", "2_ecohab2", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_22(self):
        key = "2_ecohab2"
        out = sorted(["1_ecohab1", "8_ecohab2", "2_ecohab2", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))
        
    def test_same_address_82(self):
        key = "8_ecohab2"
        out = sorted(["1_ecohab1", "8_ecohab2", "2_ecohab2", "8_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_21(self):
        key = "2_ecohab1"
        out = sorted(["2_ecohab1", "3_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_31(self):
        key = "3_ecohab1"
        out = sorted(["2_ecohab1", "3_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_41(self):
        key = "4_ecohab1"
        out = sorted(["4_ecohab1", "5_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_51(self):
        key = "5_ecohab1"
        out = sorted(["4_ecohab1", "5_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_61(self):
        key = "6_ecohab1"
        out = sorted(["6_ecohab1", "7_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_71(self):
        key = "7_ecohab1"
        out = sorted(["6_ecohab1", "7_ecohab1"])
        self.assertEqual(out,
                         sorted(self.experiment_config.same_address[key]))

    def test_same_address_full_exp_keys(self):
        keys = sorted(self.full_exp.all_antennas)
        self.assertEqual(keys,
                         sorted(self.full_exp.same_address.keys()))

    def test_same_address_full_exp_82(self):
        key = "8_ecohab_2"
        out = sorted(["8_ecohab_2", "1_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_11(self):
        key = "1_ecohab_1"
        out = sorted(["8_ecohab_2", "1_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_21(self):
        key = "2_ecohab_1"
        out = sorted(["2_ecohab_1", "3_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_31(self):
        key = "3_ecohab_1"
        out = sorted(["2_ecohab_1", "3_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_41(self):
        key = "4_ecohab_1"
        out = sorted(["4_ecohab_1", "5_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_52(self):
        key = "5_ecohab_2"
        out = sorted(["5_ecohab_2", "4_ecohab_1"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_62(self):
        key = "6_ecohab_2"
        out = sorted(["6_ecohab_2", "7_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_same_address_full_exp_72(self):
        key = "7_ecohab_2"
        out = sorted(["6_ecohab_2", "7_ecohab_2"])
        self.assertEqual(out,
                         sorted(self.full_exp.same_address[key]))

    def test_opposite_tunnel_dict_full_exp_keys(self):
        self.assertEqual(sorted(self.full_exp.all_antennas),
                         sorted(self.full_exp.opposite_tunnel.keys()))

    def test_opposite_tunnel_dict_keys(self):
        keys = set(self.experiment_config.all_antennas) - set(self.experiment_config.internal_antennas)
        self.assertEqual(sorted(keys),
                         sorted(self.experiment_config.opposite_tunnel.keys()))

    def test_opposite_tunnel_full_exp_1(self):
        key = "1_ecohab_1"
        out = sorted(["5_ecohab_2", "6_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_full_exp_2(self):
        key = "2_ecohab_1"
        out = sorted(["5_ecohab_2", "6_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)


    def test_opposite_tunnel_full_exp_5(self):
        key = "5_ecohab_2"
        out = sorted(["1_ecohab_1", "2_ecohab_1"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_full_exp_6(self):
        key = "6_ecohab_2"
        out = sorted(["1_ecohab_1", "2_ecohab_1"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_full_exp_3(self):
        key = "3_ecohab_1"
        out = sorted(["7_ecohab_2", "8_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_full_exp_4(self):
        key = "4_ecohab_1"
        out = sorted(["7_ecohab_2", "8_ecohab_2"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_full_exp_7(self):
        key = "7_ecohab_2"
        out = sorted(["3_ecohab_1", "4_ecohab_1"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_full_exp_8(self):
        key = "8_ecohab_2"
        out = sorted(["3_ecohab_1", "4_ecohab_1"])
        self.assertEqual(sorted(self.full_exp.opposite_tunnel[key]),
                         out)

    def test_opposite_tunnel_12(self):
        key = "1_ecohab2"
        out = sorted(["5_ecohab1", "6_ecohab1", "3_ecohab1", "4_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_22(self):
        key = "2_ecohab2"
        out = sorted(["5_ecohab1", "6_ecohab1", "3_ecohab1", "4_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_31(self):
        key = "3_ecohab1"
        out = sorted(["1_ecohab2", "2_ecohab2", "7_ecohab1", "8_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_41(self):
        key = "4_ecohab1"
        out = sorted(["1_ecohab2", "2_ecohab2", "7_ecohab1", "8_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_71(self):
        key = "7_ecohab1"
        out = sorted([ "3_ecohab1", "4_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_81(self):
        key = "8_ecohab1"
        out = sorted(["3_ecohab1", "4_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_11(self):
        key = "1_ecohab1"
        out = sorted(["5_ecohab1", "6_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])
    
    def test_opposite_tunnel_21(self):
        key = "2_ecohab1"
        out = sorted(["5_ecohab1", "6_ecohab1"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_51(self):
        key = "5_ecohab1"
        out = sorted(["1_ecohab1", "2_ecohab1", "1_ecohab2", "2_ecohab2"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])

    def test_opposite_tunnel_61(self):
        key = "6_ecohab1"
        out = sorted(["1_ecohab1", "2_ecohab1", "1_ecohab2", "2_ecohab2"])
        self.assertEqual(out, self.experiment_config.opposite_tunnel[key])


    def test_cage_address_dict(self):
        out = {}
        out["1_ecohab2"] = "ecohab2 cage A"
        out["2_ecohab2"] = "shared cage 1"
        out["8_ecohab2"] = "shared cage 1"
        out["1_ecohab1"] = "shared cage 1"
        out["8_ecohab1"] = "shared cage 1"
        out["2_ecohab1"] = "ecohab1 cage B"
        out["3_ecohab1"] = "ecohab1 cage B"
        out["4_ecohab1"] = "ecohab1 cage C"
        out["5_ecohab1"] = "ecohab1 cage C"
        out["6_ecohab1"] = "ecohab1 cage D"
        out["7_ecohab1"] = "ecohab1 cage D"
        self.assertEqual(out, self.experiment_config.address)

    def test_cage_address_full_expdict(self):
        out = {}
        out["1_ecohab_1"] = "cage A"
        out["8_ecohab_2"] = "cage A"
        out["2_ecohab_1"] = "ecohab_1 cage B"
        out["3_ecohab_1"] = "ecohab_1 cage B"
        out["4_ecohab_1"] = "cage C"
        out["5_ecohab_2"] = "cage C"
        out["6_ecohab_2"] = "ecohab_2 cage D"
        out["7_ecohab_2"] = "ecohab_2 cage D"
        self.assertEqual(out,
                         self.full_exp.address)

    def test_address_non_adjacent_full_exp(self):
        out = {
            "1_ecohab_1": "ecohab_1 cage B",
            "2_ecohab_1": "cage A",
            "3_ecohab_1": "cage C",
            "4_ecohab_1": "ecohab_1 cage B",
            "5_ecohab_2": "ecohab_2 cage D",
            "6_ecohab_2": "cage C",
            "7_ecohab_2": "cage A",
            "8_ecohab_2": "ecohab_2 cage D"
        }
        self.assertEqual(out, self.full_exp.address_non_adjacent)

    def test_address_non_adjacent(self):
        out = {
            "1_ecohab2": "shared cage 1",
            "2_ecohab2": "ecohab2 cage A",
            "1_ecohab1": "ecohab1 cage B",
            "2_ecohab1": "shared cage 1",
            "3_ecohab1": "ecohab1 cage C",
            "4_ecohab1": "ecohab1 cage B",
            "5_ecohab1": "ecohab1 cage D",
            "6_ecohab1": "ecohab1 cage C",
            "7_ecohab1": "shared cage 1",
            "8_ecohab1": "ecohab1 cage D",
        }
        
        self.assertEqual(out, self.experiment_config.address_non_adjacent)
        
        
if __name__ == '__main__':
    unittest.main()
