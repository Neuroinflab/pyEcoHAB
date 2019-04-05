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
    pass


if __name__ == '__main__':
    unittest.main()
