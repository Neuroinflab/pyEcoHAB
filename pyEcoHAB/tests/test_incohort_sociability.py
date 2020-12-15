from __future__ import print_function, division, absolute_import
import os
import unittest
import numpy as np

from pyEcoHAB import incohort_sociability as ics
from pyEcoHAB import utility_functions as utils
from pyEcoHAB import data_path, sample_data
from pyEcoHAB import Loader
from pyEcoHAB import Timeline


try:
    basestring
except NameError:
    basestring = str


class TestPrepareMouseIntervals(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.out1 = ics.prepare_mice_intervals(data, "cage B")
        cls.out2 = ics.prepare_mice_intervals(data, "cage C")
        cls.out3 = ics.prepare_mice_intervals(data, "cage D")
        cls.out4 = ics.prepare_mice_intervals(data, "cage A")

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
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        mouse3 = [["cage B", 2, 3.1],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse4 = [["cage B", 2, 2.5],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]

        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            'mouse3': mouse3,
            'mouse4': mouse4,
            }
        self.out1 = ics.prepare_mice_intervals(data, "cage B")
        self.out2 = ics.prepare_mice_intervals(data, "cage C")
        self.out3 = ics.prepare_mice_intervals(data, "cage D")
        self.out4 = ics.prepare_mice_intervals(data, "cage A")

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
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        mouse3 = [["cage B", 2, 3.1],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse4 = [["cage B", 2, 2.5],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]

        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            'mouse3': mouse3,
            'mouse4': mouse4,
            }
        self.out1 = ics.prepare_mice_intervals(data, "cage B")
        self.out2 = ics.prepare_mice_intervals(data, "cage C")
        self.out3 = ics.prepare_mice_intervals(data, "cage D")
        self.out4 = ics.prepare_mice_intervals(data, "cage A")

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
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        mouse3 = [["cage B", 2, 3.1],
                  ["cage A", 5, 6],
                  ["cage D", 7, 10],
                  ["cage C", 11, 15],
                  ["cage B", 16, 25],
                  ["cage C", 27, 35],
                  ["cage D", 38, 45],
                  ["cage A", 50, 52],
                  ]
        data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            'mouse3': mouse3,
            }
        cls.out1 = ics.mouse_alone(data, "cage B")
        cls.out2 = ics.mouse_alone(data, "cage C")
        cls.out3 = ics.mouse_alone(data, "cage D")
        cls.out4 = ics.mouse_alone(data, "cage A")

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
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }

    def test_mouse1_mouse2_address_1_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage B")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage B")
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage B")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage B")
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 1 + 18 - 14)

    def test_mouse1_mouse2_address_2_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage C")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage C")
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage C")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage C")
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 0)

    def test_mouse1_mouse2_address_3_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage D")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage D")
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage D")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage D")
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 1)

    def test_mouse1_mouse2_address_4_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage A")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage A")
        out1 = ics.mice_overlap(ints1, ints2)
        out2 = ics.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage A")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage A")
        out1 = ics.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 5)


class TestTimeTogether(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.duration = 100

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage B")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage B")
        out1 = ics.time_fraction_together_one_cage(ints1, ints2,
                                                   self.duration)
        self.assertEqual(out1, 5/self.duration)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage C")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage C")
        out1 = ics.time_fraction_together_one_cage(ints1, ints2,
                                                   self.duration)
        self.assertEqual(out1, 0)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage D")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage D")
        out1 = ics.time_fraction_together_one_cage(ints1, ints2,
                                                   self.duration)
        self.assertEqual(out1, 1/self.duration)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage A")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage A")
        out1 = ics.time_fraction_together_one_cage(ints1, ints2,
                                                   self.duration)
        self.assertEqual(out1, 5/self.duration)


class TestExpectedTimeTogether(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mouse1 = [["cage B", 2, 3],
                  ["cage A", 5, 6],
                  ["cage D", 8, 9],
                  ["cage C", 10, 12],
                  ["cage B", 14, 20],
                  ["cage C", 21, 28],
                  ["cage D", 31, 35],
                  ["cage A", 40, 45],
                  ]
        mouse2 = [["cage B", 0, 3],
                  ["cage C", 5, 6],
                  ["cage D", 8, 9],
                  ["cage A", 10, 12],
                  ["cage B", 13, 18],
                  ["cage A", 22, 50],
                  ]
        cls.data = {
            'mouse1': mouse1,
            'mouse2': mouse2,
            }
        cls.duration = 100
        cls.duration2 = cls.duration**2

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage B")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage B")
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2,
                                                            self.duration)
        res = np.isclose(out1, 56/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage C")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage C")
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2,
                                                            self.duration)
        res = np.isclose(out1, 9/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage D")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage D")
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2,
                                                            self.duration)
        res = np.isclose(out1, 5/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], "cage A")
        ints2 = utils.get_intervals(self.data["mouse2"], "cage A")
        out1 = ics.expected_time_fraction_together_one_cage(ints1, ints2,
                                                            self.duration)
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
                                               ["A", "B", "C", "D"],
                                               cls.duration)

    def test_mouse1_mouse2_exp(self):
        dur2 = self.duration**2
        res = np.isclose(self.out2, 250/dur2)
        self.assertTrue(res)

    def test_mouse1_mouse2_measured(self):
        res = np.isclose(self.out1, 11/self.duration)
        self.assertTrue(res)


class TestSinglePhaseResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.duration = 3*3600
        path = os.path.join(data_path, "weird_3_mice")
        cls.config = Timeline(path)
        data = Loader(path)
        cls.phases, cls.total_time,\
            cls.data, cls.keys = utils.prepare_binned_data(data,
                                                           cls.config,
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
                                                    ["cage A", "cage D"],
                                                    cls.duration)
        cls.out_A, cls.exp_A = ics.single_phase_results(cls.data["1 dark"][0],
                                                        ["mouse_1", "mouse_2"],
                                                        ["cage A"],
                                                        cls.duration)
        cls.correct_res_A = {
            "mouse_1": {
                "mouse_1": 0,
                "mouse_2": A/cls.duration,
            },
            "mouse_2": {
                "mouse_1": 0,
                "mouse_2": 0,
            }
        }
        cls.correct_exp_res_A = {
            "mouse_1": {
                "mouse_1": 0,
                "mouse_2": AA/cls.duration/cls.duration,
            },
            "mouse_2": {
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


class TestGetIncohortSociability(unittest.TestCase):
    def test_run(cls):
        data = Loader(sample_data)
        config = Timeline(sample_data)
        ics.get_incohort_sociability(data, config, 3600)
        ics.get_incohort_sociability(data, config, 24*3600)


if __name__ == '__main__':
    unittest.main()
