from __future__ import print_function, division, absolute_import
import os
import unittest
import numpy as np

from pyEcoHAB import incohort_sociability as ics
from pyEcoHAB import utility_functions as utils
from pyEcoHAB import data_path
from pyEcoHAB import Loader
from pyEcoHAB import ExperimentConfigFile

try:
    basestring
except NameError:
    basestring = str

class TestPrepareMouseIntervals(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.out1 = ics.prepare_mice_intervals(data, 1)
        cls.out2 = ics.prepare_mice_intervals(data, 2)
        cls.out3 = ics.prepare_mice_intervals(data, 3)
        cls.out4 = ics.prepare_mice_intervals(data, 4)
        
    def test_check_mouse1(self):
        out = {
            'mouse1': [
                [2, 14],
                [3, 20],       
            ],
            'mouse2': [
                [0, 13],
                [3, 18]]
            }
        self.assertEqual(self.out1, out)
    
    def test_check_mouse2(self):
        out = {
            'mouse1': [
                [10, 21],
                [12, 28],       
            ],
            'mouse2': [
                [5],
                [6]
            ]
            }
        self.assertEqual(self.out2, out)

    def test_check_mouse3(self):
        out = {
            'mouse1': [
                [8, 31],
                [9, 35],       
            ],
            'mouse2': [
                [8],
                [9]
            ]
            }
        self.assertEqual(self.out3, out)

    def test_check_mouse4(self):
        out = {
            'mouse1': [
                [5, 40],
                [6, 45],       
            ],
            'mouse2': [
                [10, 22],
                [12, 50]]
            }
        self.assertEqual(self.out4, out)
        
class TestCheckInterval(unittest.TestCase):
    def setUp(self):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        mouse3 = [[1, 2, 3.1],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse4 = [[1, 2, 2.5],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]

        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            'mouse3': mouse3,
            'mouse4': mouse4,
            }
        self.out1 = ics.prepare_mice_intervals(data, 1)
        self.out2 = ics.prepare_mice_intervals(data, 2)
        self.out3 = ics.prepare_mice_intervals(data, 3)
        self.out4 = ics.prepare_mice_intervals(data, 4)

    def test_address_mouse1_mouse2_remove_True(self):
        im1 = self.out1['mouse1'][:]
        im2 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertTrue(out)

    def test_address_mouse1_mouse2_im1(self):
        im1 = self.out1['mouse1'][:]
        im2 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[14], [20]])

    def test_address_mouse1_mouse2_im2(self):
        im2 = self.out1['mouse1'][:]
        im1 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[0, 13], [2, 18]])

    def test_address_mouse3_mouse2_remove_False(self):
        im1 = self.out1['mouse3'][:]
        im2 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertFalse(out)

    def test_address_mouse3_mouse2_im1(self):
        im1 = self.out1['mouse3'][:]
        im2 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[3, 14], [3.1, 20]])

    def test_address_mouse3_mouse2_im2(self):
        im2 = self.out1['mouse3'][:]
        im1 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[0, 13], [2, 18]])

    def test_address_mouse4_mouse2_remove_False(self):
        im1 = self.out1['mouse4'][:]
        im2 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertTrue(out)

    def test_address_mouse4_mouse2_im1(self):
        im1 = self.out1['mouse4'][:]
        im2 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[14], [20]])

    def test_address_mouse4_mouse2_im2(self):
        im2 = self.out1['mouse4'][:]
        im1 = self.out1['mouse2'][:]
        out = ics.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[0, 2.5, 13], [2, 3, 18]])

class TestRemoveOverlappingIntervals(unittest.TestCase):
    def setUp(self):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        mouse3 = [[1, 2, 3.1],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse4 = [[1, 2, 2.5],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]

        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            'mouse3': mouse3,
            'mouse4': mouse4,
            }
        self.out1 = ics.prepare_mice_intervals(data, 1)
        self.out2 = ics.prepare_mice_intervals(data, 2)
        self.out3 = ics.prepare_mice_intervals(data, 3)
        self.out4 = ics.prepare_mice_intervals(data, 4)

    def test_mouse1_mouse2_im1(self):
        im1 = self.out1["mouse1"][:]
        im2 = self.out1["mouse2"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[18], [20]])

    def test_mouse1_mouse2_im2(self):
        im2 = self.out1["mouse1"][:]
        im1 = self.out1["mouse2"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[0, 13], [2, 14]])

    def test_mouse3_mouse1_im1(self):
        im1 = self.out1["mouse3"][:]
        im2 = self.out1["mouse1"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[3], [3.1]])

    def test_mouse3_mouse1_im2(self):
        im2 = self.out1["mouse3"][:]
        im1 = self.out1["mouse1"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[], []])

    def test_mouse2_mouse3_im1(self):
        im1 = self.out1["mouse2"][:]
        im2 = self.out1["mouse3"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[0, 13], [2, 14]])

    def test_mouse2_mouse3_im2(self):
        im2 = self.out1["mouse2"][:]
        im1 = self.out1["mouse3"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[3, 18], [3.1, 20]])

    def test_mouse2_mouse4_im1(self):
        im1 = self.out1["mouse2"][:]
        im2 = self.out1["mouse4"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[0, 2.5, 13], [2, 3, 14]])

    def test_mouse2_mouse4_im2(self):
        im2 = self.out1["mouse2"][:]
        im1 = self.out1["mouse4"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[18], [20]])

    def test_mouse1_mouse2_im1_cage_4(self):
        im1 = self.out4["mouse1"][:]
        im2 = self.out4["mouse2"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[5], [6]])

    def test_mouse1_mouse2_im2_cage_4(self):
        im2 = self.out4["mouse1"][:]
        im1 = self.out4["mouse2"][:]
        ics.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[10, 22, 45], [12, 40, 50]])


class TestMouseAlone(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        mouse3 = [[1, 2, 3.1],
                  [4, 5, 6],
                  [3, 7, 10],
                  [2, 11, 15],
                  [1, 16, 25],
                  [2, 27, 35],
                  [3, 38, 45],
                  [4, 50, 52],
                  ]
        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            'mouse3': mouse3,
            }
        cls.out1 = ics.mouse_alone(data, 1)
        cls.out2 = ics.mouse_alone(data, 2)
        cls.out3 = ics.mouse_alone(data, 3)
        cls.out4 = ics.mouse_alone(data, 4)

    def test_address_1_mouse1(self):
        self.assertEqual(self.out1["mouse1"], 0)

    def test_address_1_mouse2(self):
        self.assertEqual(self.out1["mouse2"], 3)

    def test_address_1_mouse3(self):
        self.assertEqual(self.out1["mouse3"], 5.1)

    def test_address_2_mouse1(self):
        self.assertEqual(self.out2["mouse1"], 7)

    def test_address_2_mouse2(self):
        self.assertEqual(self.out2["mouse2"], 1)

    def test_address_2_mouse3(self):
        self.assertEqual(self.out2["mouse3"], 10)

    def test_address_3_mouse1(self):
        self.assertEqual(self.out3["mouse1"], 4)

    def test_address_3_mouse2(self):
        self.assertEqual(self.out3["mouse2"], 0)

    def test_address_3_mouse3(self):
        self.assertEqual(self.out3["mouse3"], 9)

    def test_address_4_mouse1(self):
        self.assertEqual(self.out4["mouse1"], 0)

    def test_address_4_mouse2(self):
        self.assertEqual(self.out4["mouse2"], 25)

    def test_address_4_mouse3(self):
        self.assertEqual(self.out4["mouse3"], 2)


class TestMiceOverlap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }

    def test_mouse1_mouse2_address_1_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 1)
        ints2 = utils.get_intervals(self.data["mouse2"], 1)
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 1)
        ints2 = utils.get_intervals(self.data["mouse2"], 1)
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 1 + 18 - 14)

    def test_mouse1_mouse2_address_2_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 0)

    def test_mouse1_mouse2_address_3_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 1)

    def test_mouse1_mouse2_address_4_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 5)


class TestTimeTogether(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.duration = 100

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 1)
        ints2 = utils.get_intervals(self.data["mouse2"], 1)
        out1 = ics.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 5/self.duration)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = ics.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 0)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = ics.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 1/self.duration)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = ics.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 5/self.duration)



class TestExpectedTimeTogether(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [[1, 2, 3],
                  [4, 5, 6],
                  [3, 8, 9],
                  [2, 10, 12],
                  [1, 14, 20],
                  [2, 21, 28],
                  [3, 31, 35],
                  [4, 40, 45],
                  ]
        mouse2 = [[1, 0, 3],
                  [2, 5, 6],
                  [3, 8, 9],
                  [4, 10, 12],
                  [1, 13, 18],
                  [4, 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.duration = 100
        cls.duration2 = cls.duration**2

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 1)
        ints2 = utils.get_intervals(self.data["mouse2"], 1)
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 56/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 9/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 5/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 6*30/self.duration2)
        self.assertTrue(res)


class TestExpectedTimeTogether(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [["B", 2, 3],
                  ["A", 5, 6],
                  ["D", 8, 9],
                  ["C", 10, 12],
                  ["B", 14, 20],
                  ["C", 21, 28],
                  ["D", 31, 35],
                  ["A", 40, 45],
                  ]
        mouse2 = [["B", 0, 3],
                  ["C", 5, 6],
                  ["D", 8, 9],
                  ["A", 10, 12],
                  ["B", 13, 18],
                  ["A", 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.duration = 100
        cls.out1, cls.out2 = ics.mice_together(cls.data, "mouse1", "mouse2",
                                              ["A", "B", "C", "D"], cls.duration)


    def test_mouse1_mouse2_exp(self):
        dur2 = self.duration**2
        res = np.isclose(self.out2, 250/dur2)
        self.assertTrue(res)

    def test_mouse1_mouse2_measured(self):
        res = np.isclose(self.out1, 11/self.duration)
        self.assertTrue(res)



class TestPrepareFnamesAndTotals(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        cls.data = Loader(path)
        cls.config = ExperimentConfigFile(path)

        cls.all_phases, cls.all_total_time,\
            cls.all_data, cls.all_keys = ics.prepare_fnames_and_totals(cls.data,
                                                                      cls.config, "",
                                                                      "ALL",
                                                                      ["mouse_1"])
        cls.dark_phases, cls.dark_total_time,\
            cls.dark_data, cls.dark_keys = ics.prepare_fnames_and_totals(cls.data,
                                                                      cls.config, "",
                                                                      "DARK",
                                                                      ["mouse_1"])
        cls.light_phases, cls.light_total_time,\
            cls.light_data, cls.light_keys = ics.prepare_fnames_and_totals(cls.data,
                                                                      cls.config, "",
                                                                      "LIGHT",
                                                                      ["mouse_1"])

        cls.phases_100s_bins, cls.total_time_100s_bins,\
            cls.data_100s_bins, cls.keys_100s_bins = ics.prepare_fnames_and_totals(cls.data,
                                                                                  cls.config, "",
                                                                                  100,
                                                                                  ["mouse_1"])
        cls.phases_900s_bins, cls.total_time_900s_bins,\
            cls.data_900s_bins, cls.keys_900s_bins = ics.prepare_fnames_and_totals(cls.data,
                                                                                  cls.config, "",
                                                                                  900,
                                                                                  ["mouse_1"])

    def test_all_phases(self):
        self.assertEqual(self.all_phases, ["ALL"])

    def test_all_time(self):
        time_dict = {"ALL": {0: 3600}}
        self.assertEqual(self.all_total_time, time_dict)

    def test_all_data(self):
        data = {"mouse_1": utils.prepare_data(self.data, ["mouse_1"])}
        all_data = {"ALL": {0: data}}

    def test_all_keys(self):
        keys = [["ALL"], ["0"]]
        self.assertTrue(keys, self.all_keys)
    
    def test_dark_phases(self):
        self.assertEqual(self.dark_phases, ["DARK"])

    def test_dark_time(self):
        time_dict = {"DARK": {0: 1800.0}}
        self.assertEqual(self.dark_total_time, time_dict)

    def test_dark_data(self):
        data = {"mouse_1": utils.prepare_data(self.data, ["mouse_1"])}
        dark_data = {"DARK": {0: data}}

    def test_dark_keys(self):
        keys = [["DARK"], ["0"]]
        self.assertTrue(keys, self.dark_keys)

    def test_light_phases(self):
        self.assertEqual(self.light_phases, ["LIGHT"])

    def test_light_time(self):
        time_dict = {"LIGHT": {0: 1800.0}}
        self.assertEqual(self.light_total_time, time_dict)

    def test_light_data(self):
        data = {"mouse_1": utils.prepare_data(self.data, ["mouse_1"])}
        light_data = {"LIGHT": {0: data}}

    def test_light_keys(self):
        keys = [["LIGHT"], ["0"]]
        self.assertTrue(keys, self.light_keys)

    def test_bins_keys(self):
        keys = [["1 dark"], [i*100/3600. for i in range(18)]]
        self.assertTrue(keys, self.keys_100s_bins)

    def test_bins_data_1st_bin(self):
        self.assertEqual(self.data_100s_bins["1 dark"][0]["mouse_1"], [])

    def test_bins_data_2nd_bin(self):
        self.assertEqual(self.data_100s_bins["1 dark"][100.]["mouse_1"], [])

    def test_bins_data_3rd_bin_len(self):
        data = self.data_100s_bins["1 dark"][200.]["mouse_1"]
        self.assertEqual(len(data), 4)

    def test_bins_data_3rd_bin_last_value(self):
        data = self.data_100s_bins["1 dark"][200.]["mouse_1"][3]
        self.assertEqual(data[-1], 1286701500)

    def test_bins_data_8th_bin_last_value(self):
        data = self.data_100s_bins["1 dark"][600.]["mouse_1"][0]
        self.assertEqual(data[-1], 1286701700)

    def test_bins_data_8th_bin_last_value(self):
        data = self.data_100s_bins["1 dark"][600.]["mouse_1"][0]
        self.assertEqual(data[1], 1286701800)

    def test_bins_data_8th_address(self):
        data = self.data_100s_bins["1 dark"][600.]["mouse_1"][0]
        self.assertEqual(data[0], "B")

class TestSinglePhaseResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.duration = 43200
        path = os.path.join(data_path, "weird_3_mice")
        cls.config = ExperimentConfigFile(path)
        data = Loader(path)
        cls.phases, cls.total_time,\
            cls.data, cls.keys = ics.prepare_fnames_and_totals(data,
                                                              cls.config, "",
                                                              cls.duration,
                                                              ["mouse_1",
                                                               "mouse_2"])
        DD = 984.282*671.526
        AA = 326.757*254.848
        D = 84.312
        A = 154.154
        cls.correct_res = {
            "mouse_1": {
                "mouse_1": 0,
                "mouse_2": (A+D)/cls.duration,
                },
            "mouse_2": {
                "mouse_1": 0,
                "mouse_2": 0,
            }
        }
        cls.correct_exp_res = {
            "mouse_1": {
                "mouse_1": 0,
                "mouse_2": (AA+DD)/cls.duration**2,
            },
            "mouse_2": {
                "mouse_1": 0,
                "mouse_2": 0,
            }
        }

        cls.out, cls.exp = ics.single_phase_results(cls.data["1 dark"][0],
                                                   ["mouse_1", "mouse_2"],
                                                   ["A", "D"],
                                                   cls.duration)
        cls.out_A, cls.exp_A = ics.single_phase_results(cls.data["1 dark"][0],
                                                       ["mouse_1", "mouse_2"],
                                                       ["A"],
                                                       cls.duration)
        cls.correct_res_A= {
            "mouse_1": {
                "mouse_1": 0,
                "mouse_2": A/cls.duration,
            },
            "mouse_2":{
                "mouse_1": 0,
                "mouse_2": 0,
            }
        }
        cls.correct_exp_res_A= {
            "mouse_1": {
                "mouse_1": 0,
                "mouse_2": AA/cls.duration/cls.duration,
            },
            "mouse_2":{
                "mouse_1": 0,
                "mouse_2": 0,
            }
        }


    def test_keys(self):
        self.assertEqual(sorted(self.out.keys()), ["mouse_1", "mouse_2"])

    def test_exp_A_1_1(self):
        self.assertEqual(self.correct_exp_res_A["mouse_1"]["mouse_1"],
                         self.exp_A["mouse_1"]["mouse_1"])

    def test_exp_A_1_2(self):
        self.assertTrue(np.isclose(self.correct_exp_res_A["mouse_1"]["mouse_2"],
                                   self.exp_A["mouse_1"]["mouse_2"]))

    def test_exp_A_2_1(self):
        self.assertEqual(self.correct_exp_res_A["mouse_2"]["mouse_1"],
        self.exp_A["mouse_2"]["mouse_1"])

    def test_exp_A_2_2(self):
        self.assertEqual(self.correct_exp_res_A["mouse_2"]["mouse_2"],
        self.exp_A["mouse_2"]["mouse_2"])

    def test_A_1_1(self):
        self.assertEqual(self.correct_res_A["mouse_1"]["mouse_1"],
                         self.out_A["mouse_1"]["mouse_1"])

    def test_A_1_2(self):
        self.assertTrue(np.isclose(self.correct_res_A["mouse_1"]["mouse_2"],
                                   self.out_A["mouse_1"]["mouse_2"]))

    def test_A_2_1(self):
        self.assertEqual(self.correct_res_A["mouse_2"]["mouse_1"],
        self.out_A["mouse_2"]["mouse_1"])

    def test_A_2_2(self):
        self.assertEqual(self.correct_res_A["mouse_2"]["mouse_2"],
        self.out_A["mouse_2"]["mouse_2"])

    def test_exp_1_1(self):
        self.assertEqual(self.correct_exp_res["mouse_1"]["mouse_1"],
                         self.exp["mouse_1"]["mouse_1"])

    def test_exp_1_2(self):
        self.assertTrue(np.isclose(self.correct_exp_res["mouse_1"]["mouse_2"],
                                   self.exp["mouse_1"]["mouse_2"]))

    def test_exp_2_1(self):
        self.assertEqual(self.correct_exp_res["mouse_2"]["mouse_1"],
        self.exp_A["mouse_2"]["mouse_1"])

    def test_exp_A_2_2(self):
        self.assertEqual(self.correct_exp_res_A["mouse_2"]["mouse_2"],
        self.exp["mouse_2"]["mouse_2"])

    def test_1_1(self):
        self.assertEqual(self.correct_res["mouse_1"]["mouse_1"],
                         self.out["mouse_1"]["mouse_1"])

    def test_1_2(self):
        self.assertTrue(np.isclose(self.correct_res["mouse_1"]["mouse_2"],
                                   self.out["mouse_1"]["mouse_2"]))

    def test_2_1(self):
        self.assertEqual(self.correct_res["mouse_2"]["mouse_1"],
        self.out["mouse_2"]["mouse_1"])

    def test_2_2(self):
        self.assertEqual(self.correct_res["mouse_2"]["mouse_2"],
        self.out["mouse_2"]["mouse_2"])


if __name__ == '__main__':
    unittest.main()
