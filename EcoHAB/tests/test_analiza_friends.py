from __future__ import print_function, division, absolute_import
import os
from EcoHAB import analiza_friends as af
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

if __name__ == '__main__':
    unittest.main()
