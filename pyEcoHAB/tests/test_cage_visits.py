from __future__ import print_function, division, absolute_import
from pyEcoHAB import cage_visits as cv
import numpy as np
import unittest

class TestGetVisits(unittest.TestCase):
    def test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, [3, 2, 30, 10])

    def test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, [3, 2, 30, 10])

    def test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, [3, 2, 30, 10, 5])

    def test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, [3, 2, 30, 10, 5])

    def test_one_longer_interval(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 110]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, [])


class TestVisitsDurationsPhase(unittest.TestCase):
    def test_one_bin_test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertTrue(v, [3, 2, 30, 10])

    def test_one_bin_test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertTrue(v, [3, 2, 30, 10])

    def test_one_bin_test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertTrue(v, [3, 2, 30, 10, 35])

    def test_one_bin_test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertTrue(v, [3, 2, 30, 10, 35])

    def test_more_bin_test_intervals_only_in_the_bin(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v = cv.get_visits_in_bins(intervals, t_start, t_stop, 10)
        self.assertTrue(v, [[4], [3, 2], [], [], [30], [], [], [], [10], []])

    def test_more_bin_test_intervals_starting_before_the_bin(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v = cv.get_visits_in_bins(intervals, t_start, t_stop, 10)
        self.assertTrue(v, [[11], [3, 2], [], [], [30], [], [], [], [10], []])

if __name__ == '__main__':
    unittest.main()
