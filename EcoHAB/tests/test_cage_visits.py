from __future__ import print_function, division, absolute_import
from EcoHAB import cage_visits as cv
import unittest

class TestGetChamberVisitsAndDurations(unittest.TestCase):
    def test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90], [110, 130]]
        v, d = cv.get_chamber_visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 4)
        self.assertEqual(d, 45)

    def test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, d = cv.get_chamber_visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 4)
        self.assertEqual(d, 46)

    def test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, d = cv.get_chamber_visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 5)
        self.assertEqual(d, 50)

    def test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, d = cv.get_chamber_visits_and_durations(intervals, t_start, t_stop)
        self.assertEqual(v, 5)
        self.assertEqual(d, 51)

if __name__ == '__main__':
    unittest.main()
