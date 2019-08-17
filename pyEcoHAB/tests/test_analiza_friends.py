from __future__ import print_function, division, absolute_import
import os
from pyEcoHAB import analiza_friends as af
from pyEcoHAB import utility_functions as utils
import unittest
import numpy as np

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
        cls.out1 = af.prepare_mice_intervals(data, 1)
        cls.out2 = af.prepare_mice_intervals(data, 2)
        cls.out3 = af.prepare_mice_intervals(data, 3)
        cls.out4 = af.prepare_mice_intervals(data, 4)
        
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
        self.out1 = af.prepare_mice_intervals(data, 1)
        self.out2 = af.prepare_mice_intervals(data, 2)
        self.out3 = af.prepare_mice_intervals(data, 3)
        self.out4 = af.prepare_mice_intervals(data, 4)

    def test_address_mouse1_mouse2_remove_True(self):
        im1 = self.out1['mouse1'][:]
        im2 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertTrue(out)

    def test_address_mouse1_mouse2_im1(self):
        im1 = self.out1['mouse1'][:]
        im2 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[14], [20]])

    def test_address_mouse1_mouse2_im2(self):
        im2 = self.out1['mouse1'][:]
        im1 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[0, 13], [2, 18]])

    def test_address_mouse3_mouse2_remove_False(self):
        im1 = self.out1['mouse3'][:]
        im2 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertFalse(out)

    def test_address_mouse3_mouse2_im1(self):
        im1 = self.out1['mouse3'][:]
        im2 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[3, 14], [3.1, 20]])

    def test_address_mouse3_mouse2_im2(self):
        im2 = self.out1['mouse3'][:]
        im1 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[0, 13], [2, 18]])

    def test_address_mouse4_mouse2_remove_False(self):
        im1 = self.out1['mouse4'][:]
        im2 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertTrue(out)

    def test_address_mouse4_mouse2_im1(self):
        im1 = self.out1['mouse4'][:]
        im2 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
        self.assertEqual(im1, [[14], [20]])

    def test_address_mouse4_mouse2_im2(self):
        im2 = self.out1['mouse4'][:]
        im1 = self.out1['mouse2'][:]
        out = af.check_interval(im1, im2, 0, 0)
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
        self.out1 = af.prepare_mice_intervals(data, 1)
        self.out2 = af.prepare_mice_intervals(data, 2)
        self.out3 = af.prepare_mice_intervals(data, 3)
        self.out4 = af.prepare_mice_intervals(data, 4)

    def test_mouse1_mouse2_im1(self):
        im1 = self.out1["mouse1"][:]
        im2 = self.out1["mouse2"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[18], [20]])

    def test_mouse1_mouse2_im2(self):
        im2 = self.out1["mouse1"][:]
        im1 = self.out1["mouse2"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[0, 13], [2, 14]])

    def test_mouse3_mouse1_im1(self):
        im1 = self.out1["mouse3"][:]
        im2 = self.out1["mouse1"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[3], [3.1]])

    def test_mouse3_mouse1_im2(self):
        im2 = self.out1["mouse3"][:]
        im1 = self.out1["mouse1"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[], []])

    def test_mouse2_mouse3_im1(self):
        im1 = self.out1["mouse2"][:]
        im2 = self.out1["mouse3"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[0, 13], [2, 14]])

    def test_mouse2_mouse3_im2(self):
        im2 = self.out1["mouse2"][:]
        im1 = self.out1["mouse3"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[3, 18], [3.1, 20]])

    def test_mouse2_mouse4_im1(self):
        im1 = self.out1["mouse2"][:]
        im2 = self.out1["mouse4"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[0, 2.5, 13], [2, 3, 14]])

    def test_mouse2_mouse4_im2(self):
        im2 = self.out1["mouse2"][:]
        im1 = self.out1["mouse4"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[18], [20]])

    def test_mouse1_mouse2_im1_cage_4(self):
        im1 = self.out4["mouse1"][:]
        im2 = self.out4["mouse2"][:]
        af.remove_overlapping_intervals(im1, im2)
        self.assertEqual(im1, [[5], [6]])

    def test_mouse1_mouse2_im2_cage_4(self):
        im2 = self.out4["mouse1"][:]
        im1 = self.out4["mouse2"][:]
        af.remove_overlapping_intervals(im1, im2)
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
        cls.out1 = af.mouse_alone(data, 1)
        cls.out2 = af.mouse_alone(data, 2)
        cls.out3 = af.mouse_alone(data, 3)
        cls.out4 = af.mouse_alone(data, 4)

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
        out1 = af.mice_overlap(ints1, ints2)
        out2 = af.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_1(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 1)
        ints2 = utils.get_intervals(self.data["mouse2"], 1)
        out1 = af.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 1 + 18 - 14)

    def test_mouse1_mouse2_address_2_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = af.mice_overlap(ints1, ints2)
        out2 = af.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = af.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 0)

    def test_mouse1_mouse2_address_3_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = af.mice_overlap(ints1, ints2)
        out2 = af.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = af.mice_overlap(ints1, ints2)
        self.assertEqual(out1, 1)

    def test_mouse1_mouse2_address_4_symmetry(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = af.mice_overlap(ints1, ints2)
        out2 = af.mice_overlap(ints2, ints1)
        self.assertEqual(out1, out2)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = af.mice_overlap(ints1, ints2)
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
        out1 = af.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 5/self.duration)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = af.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 0)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = af.time_fraction_together_one_cage(ints1, ints2, self. duration)
        self.assertEqual(out1, 1/self.duration)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = af.time_fraction_together_one_cage(ints1, ints2, self. duration)
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
        out1 = af.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 56/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_2(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 2)
        ints2 = utils.get_intervals(self.data["mouse2"], 2)
        out1 = af.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 9/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_3(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 3)
        ints2 = utils.get_intervals(self.data["mouse2"], 3)
        out1 = af.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 5/self.duration2)
        self.assertTrue(res)

    def test_mouse1_mouse2_address_4(self):
        ints1 = utils.get_intervals(self.data["mouse1"], 4)
        ints2 = utils.get_intervals(self.data["mouse2"], 4)
        out1 = af.expected_time_fraction_together_one_cage(ints1, ints2, self. duration)
        res = np.isclose(out1, 6*30/self.duration2)
        self.assertTrue(res)


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
        cls.out1, cls.out2 = af.mice_together(cls.data, "mouse1", "mouse2", cls.duration)

    def test_mouse1_mouse2_exp(self):
        dur2 = self.duration**2
        res = np.isclose(self.out2, 250/dur2)
        print(self.out2)
        self.assertTrue(res)

    def test_mouse1_mouse2_measured(self):
        res = np.isclose(self.out1, 11/self.duration)
        self.assertTrue(res)


if __name__ == '__main__':
    unittest.main()
