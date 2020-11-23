from __future__ import print_function, division, absolute_import

import unittest
import os
import numpy as np
from pyEcoHAB import trajectories as tr
from pyEcoHAB import utility_functions as uf
from pyEcoHAB import Loader
from pyEcoHAB import Timeline
from pyEcoHAB import data_path


# class TestSingleMouseAntennaTransitions(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         antennas = ["1", "2", "3", "4", "4", "3", "2", "1", "2", "4"]
#         times =    [1,    10,  12, 15,  16,  17,  20,  25,  30,   31]
#         cls.expected = {"1 2": [9, 5],
#                     "2 3": [2],
#                     "3 4": [3],
#                     "4 4": [1],
#                     "4 3": [1],
#                     "3 2": [3],
#                     "2 1": [5],
#                     "2 4": [1]}
#         cls.calc = tr.single_mouse_antenna_transitions(antennas, times)

#     def test_1(self):
#         self.assertEqual(sorted(self.expected.keys()), sorted(self.calc.keys()))

        
# class TestAntennaTransitions(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         path = os.path.join(data_path, "weird_very_short_3_mice")
#         cls.timeline = Timeline(path)
#         cls.data = Loader(path)
#         cls.expected = {}
#         for a1 in cls.data.setup_config.all_antennas:
#             for a2 in cls.data.setup_config.all_antennas:
#                  key = "%s %s" % (a1, a2)
#                  cls.expected[key] = []
#         cls.expected["4 4"] = [2.348, 0.25, 18.326]
#         cls.expected["4 5"] = [4.38, 3.923, ]
#         cls.expected["5 6"] = [1.963, 19.035, 1.136, 12.458, ]
#         cls.expected["6 6"] = [3.882, .906, 101.777, 17.995, 11.994]
#         cls.expected["6 5"] = [88.224, 0.734, 0.547, ]
#         cls.expected["5 5"] = [34.989, 67.763, 7.751, .752, ]
#         cls.expected["5 4"] = [13.382, ]
#         cls.expected["6 3"] = [16.194]
#         cls.expected["3 3"] = [15.26, ]
#         cls.expected["6 4"] = [126.495, ]
#         cls.expected["4 3"] = [0.825, ]
#         cls.expected["3 2"] = [3.109, ]
#         cls.expected["2 1"] = [.891]
#         cls.expected["1 1"] = [8.874]
#         cls.expected["1 2"] = [0.651]
#         cls.calc =  tr.get_antenna_transitions(cls.data, cls.timeline)
#         for key in cls.expected:
#             line = cls.expected[key]
#             cls.expected[key] = [np.round(a, decimals=3) for a in line]
            
#         for key in cls.calc:
#             line = cls.calc[key]
#             cls.calc[key] = [np.round(a, decimals=3) for a in line]
        

#     def test_keys(self):
#         self.assertEqual(sorted(self.calc.keys()),
#                          sorted(self.expected.keys()))

#     def test_empty(self):
#         empty_expected = []
#         for key in self.expected.keys():
#             if not len(self.expected[key]):
#                 empty_expected.append(key)
#         empty_calc = []
#         for key in self.calc.keys():
#             if not len(self.calc[key]):
#                 empty_calc.append(key)
#         self.assertEqual(sorted(empty_calc), sorted(empty_expected))

#     def test_4_4(self):
#         self.assertEqual(self.expected["4 4"], self.calc["4 4"])
        

#     def test_4_5(self):
#         self.assertEqual(self.expected["4 5"], self.calc["4 5"])

#     def test_5_6(self):
#         self.assertEqual(self.expected["5 6"], self.calc["5 6"])

#     def test_6_6(self):
#         self.assertEqual(self.expected["6 6"], self.calc["6 6"])

#     def test_6_5(self):
#         self.assertEqual(self.expected["6 5"], self.calc["6 5"])

#     def test_5_5(self):
#         self.assertEqual(self.expected["5 5"], self.calc["5 5"])

#     def test_5_4(self):
#         self.assertEqual(self.expected["5 4"], self.calc["5 4"])

#     def test_6_3(self):
#         self.assertEqual(self.expected["6 3"], self.calc["6 3"])

#     def test_3_3(self):
#         self.assertEqual(self.expected["3 3"], self.calc["3 3"])

#     def test_6_4(self):
#         self.assertEqual(self.expected["6 4"], self.calc["6 4"])

#     def test_4_3(self):
#         self.assertEqual(self.expected["4 3"], self.calc["4 3"])

#     def test_3_2(self):
#         self.assertEqual(self.expected["3 2"], self.calc["3 2"])

#     def test_2_1(self):
#         self.assertEqual(self.expected["2 1"], self.calc["2 1"])

#     def test_1_1(self):
#         self.assertEqual(self.expected["1 1"], self.calc["1 1"])

#     def test_1_2(self):
#         self.assertEqual(self.expected["1 2"], self.calc["1 2"])

class TestGetAntennaTransitions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_very_short_3_mice")
        cls.data = Loader(path)
        cls.timeline = Timeline(path)

    def test_all_phases(self):
        transitions = tr.get_antenna_transitions(self.data,
                                                 self.timeline,
                                                 what_phases="")

    def test_ALL(self):
        transitions = tr.get_antenna_transitions(self.data,
                                                 self.timeline,
                                                 what_phases="ALL")


class TestTrainsOfSingleAntennaRegistrations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_very_short_3_mice")
        cls.data = Loader(path)
        cls.dur, cls.count = tr.get_registration_trains(cls.data)
        cls.pred_dur = {}
        cls.pred_count = {}
        for antenna in cls.data.setup_config.all_antennas:
            cls.pred_dur[antenna] = []
            cls.pred_count[antenna] = []
        cls.pred_dur["4"] = [2.598]
        cls.pred_count["4"] = [3]
        cls.pred_dur["6"] = [102.683]
        cls.pred_count["6"] = [3]
        cls.pred_count["5"] = [3]
        cls.pred_dur["5"] = [8.503]

    def test_dur_4(self):
        line_pred = [np.round(a, decimals=3) for a in self.pred_dur["4"]]
        line_calc = [np.round(a, decimals=3) for a in self.dur["4"]]
        self.assertEqual(line_pred, line_calc)

    def test_dur_5(self):
        line_pred = [np.round(a, decimals=3) for a in self.pred_dur["5"]]
        line_calc = [np.round(a, decimals=3) for a in self.dur["5"]]
        self.assertEqual(line_pred, line_calc)

    def test_dur_6(self):
        line_pred = [np.round(a, decimals=3) for a in self.pred_dur["6"]]
        line_calc = [np.round(a, decimals=3) for a in self.dur["6"]]
        self.assertEqual(line_pred, line_calc)

    def test_count(self):
        self.assertEqual(self.pred_count, self.count)



        
if __name__ == '__main__':
    unittest.main()
