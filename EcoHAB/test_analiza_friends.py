from __future__ import print_function, division, absolute_import
import os
import unittest
import numpy as np

import analiza_friends as af

try:
  basestring
except NameError:
  basestring = str

class testAnalizaFriends(unittest.TestCase):

    def test_interval_overlap_1(self):
        """Incorrect interval
        """
        inte_1 = [34, 45]
        inte_2 = [34, 23]
        self.assertTrue(af.interval_overlap(inte_1, inte_2) == 0)
        
    def test_interval_overlap_2(self):
        """2nd interval shorter than the first one, same start
        """
        inte_1 = [34, 45]
        inte_2 = [34, 43]
        self.assertTrue(af.interval_overlap(inte_1, inte_2) == 9)
        
    def test_interval_overlap_3(self):
        """2nd interval beginning after the 1st interval has finished
        """
        inte_1 = [34, 45]
        inte_2 = [46, 50]
        self.assertTrue(af.interval_overlap(inte_1, inte_2) == 0)
        
    def test_interval_overlap_4(self):
        """1st interval beginning after the 2nd interval has finished
        """
        inte_1 = [46, 50]
        inte_2 = [34, 45]
        self.assertTrue(af.interval_overlap(inte_1, inte_2) == 0)

    def test_interval_overlap_5(self):
        """Overlapping
        """
        inte_1 = [46, 50]
        inte_2 = [34, 48]
        self.assertTrue(af.interval_overlap(inte_1, inte_2) == 2)
        
    def test_interval_overlap_6(self):
        """Overlapping
        """
        inte_1 = [34, 48]
        inte_2 = [46, 50]
        self.assertTrue(af.interval_overlap(inte_1, inte_2) == 2)    

        
if __name__ == '__main__':
    unittest.main()
