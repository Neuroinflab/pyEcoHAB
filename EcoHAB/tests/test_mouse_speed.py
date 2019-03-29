from __future__ import print_function, division, absolute_import
from EcoHAB import mouse_speed as ms
import unittest

try:
    basestring
except NameError:
    basestring = str

class TestCheck2ndMouse(unittest.TestCase):
    def test_following(self):
        antenna1 = 1
        antenna2 = 2
        t1 = 15
        threshold = 3
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        out = ms.check_2nd_mouse(antenna1, antenna2, t1, threshold,
                                 antennas2, times2)
        self.assertEqual(out, 1)

    def test_taking_over(self):
        antenna1 = 1
        antenna2 = 2
        t1 = 15
        threshold = 3
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 18, 19.5, 22, 25]
        out = ms.check_2nd_mouse(antenna1, antenna2, t1, threshold,
                                 antennas2, times2)
        self.assertEqual(out, 0)

    def test_not_following(self):
        antenna1 = 1
        antenna2 = 2
        t1 = 15
        threshold = 3
        antennas2 = [8, 1, 2, 3, 4, 5, 6]
        times2 = [16, 19, 19.5, 22, 25]
        out = ms.check_2nd_mouse(antenna1, antenna2, t1, threshold,
                                 antennas2, times2)
        self.assertEqual(out, 0)
        
if __name__ == '__main__':
    unittest.main()
