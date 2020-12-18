from __future__ import print_function, division, absolute_import

import unittest
import os
import numpy as np
from pyEcoHAB import trajectories as tr
from pyEcoHAB import utility_functions as uf
from pyEcoHAB import Loader
from pyEcoHAB import Timeline
from pyEcoHAB import data_path


class TestSingleMouseAntennaTransitions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas = ["1", "2", "3", "4", "4", "3", "2", "1", "2", "4"]
        times = [1,    10,  12, 15,  16,  17,  20,  25,  30,   31]
        cls.expected = {"1 2": [9, 5],
                        "2 3": [2],
                        "3 4": [3],
                        "4 4": [1],
                        "4 3": [1],
                        "3 2": [3],
                        "2 1": [5],
                        "2 4": [1]}
        cls.calc = tr.single_mouse_antenna_transitions(antennas, times)

    def test_1(self):
        self.assertEqual(sorted(self.expected.keys()),
                         sorted(self.calc.keys()))


class TestAntennaTransitions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_very_short_3_mice")
        cls.timeline = Timeline(path)
        cls.data = Loader(path)
        cls.expected = {"ALL": {0.0: {}}}
        cls.expected_dark = {"1 dark": {0.0: {}},
                             "1 light": {0.0: {}}}
        for a1 in cls.data.setup_config.all_antennas:
            for a2 in cls.data.setup_config.all_antennas:
                key = "%s %s" % (a1, a2)
                cls.expected["ALL"][0.0][key] = []
                cls.expected_dark["1 dark"][0.0][key] = []
                cls.expected_dark["1 light"][0.0][key] = []
        cls.expected["ALL"][0.0]["4 4"] = [2.348, 0.25, 18.326]
        cls.expected["ALL"][0.0]["4 5"] = [4.38, 3.923, ]
        cls.expected["ALL"][0.0]["5 6"] = [1.963, 19.035, 1.136, 12.458]
        cls.expected["ALL"][0.0]["6 6"] = [3.882, .906, 101.777, 17.995,
                                           11.994]
        cls.expected["ALL"][0.0]["6 5"] = [88.224, 0.734, 0.547, ]
        cls.expected["ALL"][0.0]["5 5"] = [34.989, 67.763, 7.751, .752]
        cls.expected["ALL"][0.0]["5 4"] = [13.382]
        cls.expected["ALL"][0.0]["6 3"] = [16.194]
        cls.expected["ALL"][0.0]["3 3"] = [15.26]
        cls.expected["ALL"][0.0]["6 4"] = [126.495]
        cls.expected["ALL"][0.0]["4 3"] = [0.825]
        cls.expected["ALL"][0.0]["3 2"] = [3.109]
        cls.expected["ALL"][0.0]["2 1"] = [.891]
        cls.expected["ALL"][0.0]["1 1"] = [8.874]
        cls.expected["ALL"][0.0]["1 2"] = [0.651]
        cls.expected_dark["1 dark"][0.0]["4 4"] = [2.348, 0.25, 18.326]
        cls.expected_dark["1 dark"][0.0]["4 5"] = [4.38, 3.923, ]
        cls.expected_dark["1 dark"][0.0]["5 6"] = [1.963, 19.035, 1.136,
                                                   12.458]
        cls.expected_dark["1 dark"][0.0]["6 6"] = [3.882, .906, 101.777,
                                                   17.995, 11.994]
        cls.expected_dark["1 dark"][0.0]["6 5"] = [88.224, 0.734, 0.547]
        cls.expected_dark["1 dark"][0.0]["5 5"] = [34.989, 67.763, 7.751, .752]
        cls.expected_dark["1 dark"][0.0]["5 4"] = [13.382]
        cls.expected_dark["1 dark"][0.0]["6 3"] = [16.194]
        cls.expected_dark["1 dark"][0.0]["3 3"] = [15.26]
        cls.expected_dark["1 dark"][0.0]["6 4"] = [126.495]
        cls.expected_dark["1 dark"][0.0]["4 3"] = [0.825]
        cls.expected_dark["1 dark"][0.0]["3 2"] = [3.109]
        cls.expected_dark["1 dark"][0.0]["2 1"] = [.891]
        cls.expected_dark["1 dark"][0.0]["1 1"] = [8.874]
        cls.expected_dark["1 dark"][0.0]["1 2"] = [0.651]

        cls.calc = tr.get_antenna_transition_durations(cls.data,
                                                       cls.timeline,
                                                       bins="ALL")
        cls.calc_dark = tr.get_antenna_transition_durations(cls.data,
                                                            cls.timeline,
                                                            bins=1800)
        for key in cls.expected["ALL"][0.0]:
            line = cls.expected["ALL"][0.0][key]
            cls.expected["ALL"][0.0][key] = [np.round(a, decimals=3)
                                             for a in line]
            cls.expected_dark["1 dark"][0.0][key] = [np.round(a, decimals=3)
                                                     for a in line]
        for key in cls.calc["ALL"][0.0]:
            line = cls.calc["ALL"][0.0][key]
            cls.calc["ALL"][0.0][key] = [np.round(a, decimals=3)
                                         for a in line]
            cls.calc_dark["1 dark"][0.0][key] = [np.round(a, decimals=3)
                                                 for a in line]

    def test_keys(self):
        self.assertEqual(sorted(self.calc.keys()),
                         sorted(self.expected.keys()))

    def test_empty(self):
        empty_expected = []
        for key in self.expected.keys():
            if not len(self.expected[key]):
                empty_expected.append(key)
        empty_calc = []
        for key in self.calc.keys():
            if not len(self.calc[key]):
                empty_calc.append(key)
        self.assertEqual(sorted(empty_calc), sorted(empty_expected))

    def test_4_4(self):
        self.assertEqual(self.expected["ALL"][0.0]["4 4"],
                         self.calc["ALL"][0.0]["4 4"])

    def test_4_5(self):
        self.assertEqual(self.expected["ALL"][0.0]["4 5"],
                         self.calc["ALL"][0.0]["4 5"])

    def test_5_6(self):
        self.assertEqual(self.expected["ALL"][0.0]["5 6"],
                         self.calc["ALL"][0.0]["5 6"])

    def test_6_6(self):
        self.assertEqual(self.expected["ALL"][0.0]["6 6"],
                         self.calc["ALL"][0.0]["6 6"])

    def test_6_5(self):
        self.assertEqual(self.expected["ALL"][0.0]["6 5"],
                         self.calc["ALL"][0.0]["6 5"])

    def test_5_5(self):
        self.assertEqual(self.expected["ALL"][0.0]["5 5"],
                         self.calc["ALL"][0.0]["5 5"])

    def test_5_4(self):
        self.assertEqual(self.expected["ALL"][0.0]["5 4"],
                         self.calc["ALL"][0.0]["5 4"])

    def test_6_3(self):
        self.assertEqual(self.expected["ALL"][0.0]["6 3"],
                         self.calc["ALL"][0.0]["6 3"])

    def test_3_3(self):
        self.assertEqual(self.expected["ALL"][0.0]["3 3"],
                         self.calc["ALL"][0.0]["3 3"])

    def test_6_4(self):
        self.assertEqual(self.expected["ALL"][0.0]["6 4"],
                         self.calc["ALL"][0.0]["6 4"])

    def test_4_3(self):
        self.assertEqual(self.expected["ALL"][0.0]["4 3"],
                         self.calc["ALL"][0.0]["4 3"])

    def test_3_2(self):
        self.assertEqual(self.expected["ALL"][0.0]["3 2"],
                         self.calc["ALL"][0.0]["3 2"])

    def test_2_1(self):
        self.assertEqual(self.expected["ALL"][0.0]["2 1"],
                         self.calc["ALL"][0.0]["2 1"])

    def test_1_1(self):
        self.assertEqual(self.expected["ALL"][0.0]["1 1"],
                         self.calc["ALL"][0.0]["1 1"])

    def test_1_2(self):
        self.assertEqual(self.expected["ALL"][0.0]["1 2"],
                         self.calc["ALL"][0.0]["1 2"])

    def test_4_4_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["4 4"],
                         self.calc_dark["1 dark"][0.0]["4 4"])

    def test_4_5_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["4 5"],
                         self.calc_dark["1 dark"][0.0]["4 5"])

    def test_5_6_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["5 6"],
                         self.calc_dark["1 dark"][0.0]["5 6"])

    def test_6_6_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["6 6"],
                         self.calc_dark["1 dark"][0.0]["6 6"])

    def test_6_5_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["6 5"],
                         self.calc_dark["1 dark"][0.0]["6 5"])

    def test_5_5_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["5 5"],
                         self.calc_dark["1 dark"][0.0]["5 5"])

    def test_5_4_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["5 4"],
                         self.calc_dark["1 dark"][0.0]["5 4"])

    def test_6_3_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["6 3"],
                         self.calc_dark["1 dark"][0.0]["6 3"])

    def test_3_3_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["3 3"],
                         self.calc_dark["1 dark"][0.0]["3 3"])

    def test_6_4_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["6 4"],
                         self.calc_dark["1 dark"][0.0]["6 4"])

    def test_4_3_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["4 3"],
                         self.calc_dark["1 dark"][0.0]["4 3"])

    def test_3_2_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["3 2"],
                         self.calc_dark["1 dark"][0.0]["3 2"])

    def test_2_1_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["2 1"],
                         self.calc_dark["1 dark"][0.0]["2 1"])

    def test_1_1_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["1 1"],
                         self.calc_dark["1 dark"][0.0]["1 1"])

    def test_1_2_dark(self):
        self.assertEqual(self.expected_dark["1 dark"][0.0]["1 2"],
                         self.calc_dark["1 dark"][0.0]["1 2"])

    def test_1_light(self):
        self.assertEqual(self.expected_dark["1 light"],
                         self.calc_dark["1 light"])


class TestGetAntennaTransitions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_very_short_3_mice")
        cls.data = Loader(path)
        cls.timeline = Timeline(path)

    def test_all_phases(self):
        transitions = tr.get_antenna_transition_durations(self.data,
                                                          self.timeline,
                                                          bins=30*60)

    def test_ALL(self):
        transitions = tr.get_antenna_transition_durations(self.data,
                                                          self.timeline,
                                                          bins="ALL")


class TestTrainsOfSingleAntennaRegistrations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_very_short_3_mice")
        cls.data = Loader(path)
        cls.dur, cls.count = tr.get_registration_trains(cls.data)
        cls.pred_dur = {"ALL": {0: {}}}
        cls.pred_count = {"ALL": {0: {}}}
        for antenna in cls.data.setup_config.all_antennas:
            cls.pred_dur["ALL"][0][antenna] = []
            cls.pred_count["ALL"][0][antenna] = []
        cls.pred_dur["ALL"][0]["4"] = [2.598]
        cls.pred_count["ALL"][0]["4"] = [3]
        cls.pred_dur["ALL"][0]["6"] = [102.683]
        cls.pred_count["ALL"][0]["6"] = [3]
        cls.pred_count["ALL"][0]["5"] = [3]
        cls.pred_dur["ALL"][0]["5"] = [8.503]

    def test_dur_4(self):
        line_pred = [np.round(a, decimals=3)
                     for a in self.pred_dur["ALL"][0]["4"]]
        line_calc = [np.round(a, decimals=3) for a in self.dur["ALL"][0]["4"]]
        self.assertEqual(line_pred, line_calc)

    def test_dur_5(self):
        line_pred = [np.round(a, decimals=3)
                     for a in self.pred_dur["ALL"][0]["5"]]
        line_calc = [np.round(a, decimals=3) for a in self.dur["ALL"][0]["5"]]
        self.assertEqual(line_pred, line_calc)

    def test_dur_6(self):
        line_pred = [np.round(a, decimals=3) for a
                     in self.pred_dur["ALL"][0]["6"]]
        line_calc = [np.round(a, decimals=3) for a in self.dur["ALL"][0]["6"]]
        self.assertEqual(line_pred, line_calc)

    def test_count(self):
        self.assertEqual(self.pred_count["ALL"][0], self.count["ALL"][0])


if __name__ == '__main__':
    unittest.main()
