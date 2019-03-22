from __future__ import print_function, division, absolute_import
import os
import EcoHAB
import unittest
import numpy as np
import mouse_speed as ms
try:
    basestring
except NameError:
    basestring = str


class test_mouse_speed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(EcoHAB.sample_data_path, 'test_af')
        cls.data = EcoHAB.EcoHabData(path=path, how_many_appearances=1)

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


class TestCheck2ndMouse(unittest.TestCase):
  
        
if __name__ == '__main__':
    unittest.main()
