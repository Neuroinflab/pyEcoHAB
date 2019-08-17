from __future__ import print_function, division, absolute_import
from pyEcoHAB import cage_visits as cv
import numpy as np
import unittest

class TestVisitsAndDurations(unittest.TestCase):
    def test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v, d = cv.visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 4)
        self.assertEqual(d, 45)

    def test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, d = cv.visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 4)
        self.assertEqual(d, 45)

    def test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, d = cv.visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 5)
        self.assertEqual(d, 80)

    def test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, d = cv.visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 5)
        self.assertEqual(d, 80)

    def test_one_loonger_interval(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 110]]
        v, d = cv.visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 0)
        self.assertEqual(d, 0)


class TestVisitsDurationsPhase(unittest.TestCase):
    def test_one_bin_test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v, d = cv.visits_and_durations_bins(intervals, t_start, t_stop, 90)
        out_visits = np.array_equal(v, np.array([4]))
        out_durations = np.array_equal(d, np.array([45]))
        self.assertTrue(out_visits)
        self.assertTrue(out_durations)

    def test_one_bin_test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, d = cv.visits_and_durations_bins(intervals, t_start, t_stop, 90)
        out_visits = np.array_equal(v, np.array([4]))
        out_durations = np.array_equal(d, np.array([45]))
        self.assertTrue(out_visits)
        self.assertTrue(out_durations)

    def test_one_bin_test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, d = cv.visits_and_durations_bins(intervals, t_start, t_stop, 90)
        out_visits = np.array_equal(v, np.array([5]))
        out_durations = np.array_equal(d, np.array([80]))
        self.assertTrue(out_visits)
        self.assertTrue(out_durations)

    def test_one_bin_test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, d = cv.visits_and_durations_bins(intervals, t_start, t_stop, 90)
        out_visits = np.array_equal(v, np.array([5]))
        out_durations = np.array_equal(d, np.array([80]))
        self.assertTrue(out_visits)
        self.assertTrue(out_durations)


    def test_more_bin_test_intervals_only_in_the_bin(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v, d = cv.visits_and_durations_bins(intervals, t_start, t_stop, 10)
        out_visits = np.array_equal(v, np.array([1, 2, 0, 0, 1,
                                                 0, 0, 0, 1, 0]))
        out_durations = np.array_equal(d, np.array([4, 5, 0, 0, 30,
                                                    0, 0, 0, 10, 0]))
        self.assertTrue(out_visits)
        self.assertTrue(out_durations)

    def test_more_bin_test_intervals_starting_before_the_bin(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, d = cv.visits_and_durations_bins(intervals, t_start, t_stop, 10)
        out_visits = np.array_equal(v, np.array([1, 2, 0, 0,
                                                 1, 0, 0, 0, 1, 0]))
        out_durations = np.array_equal(d, np.array([10, 5, 0, 0, 30,
                                                    0, 0, 0, 10, 0]))
        self.assertTrue(out_visits)
        self.assertTrue(out_durations)

if __name__ == '__main__':
    unittest.main()
