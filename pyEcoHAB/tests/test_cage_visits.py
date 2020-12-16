from __future__ import print_function, division, absolute_import
import unittest
import numpy as np
from pyEcoHAB import cage_visits as cv


class TestGetVisits(unittest.TestCase):
    def test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90],
                     [110, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, ([3, 2, 30, 10], False))

    def test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, ([1, 3, 2, 30, 10], True))

    def test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, ([3, 2, 30, 10, 5], False))

    def test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, ([1, 3, 2, 30, 10, 5], True))

    def test_one_longer_interval(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 110]]
        v = cv.get_visits(intervals, t_start, t_stop)
        self.assertEqual(v, ([90], True))


class TestVisitsDurationsPhase(unittest.TestCase):
    def test_one_bin_test_intervals_only_in_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90],
                     [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(v, [[3, 2, 30, 10]])

    def test_one_bin_test_intervals_only_in_the_bin_b(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90],
                     [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(a, [False])

    def test_one_bin_test_intervals_starting_before_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(v, [[1, 3, 2, 30, 10]])

    def test_one_bin_test_intervals_starting_before_the_bin_b(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(a, [True])

    def test_one_bin_test_intervals_ending_after_the_bin(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(v, [[3, 2, 30, 10, 5]])

    def test_one_bin_test_intervals_ending_after_the_bin_b(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 4], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(a, [False])

    def test_one_bin_test_intervals(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(v, [[1, 3, 2, 30, 10, 5]])

    def test_one_bin_test_intervals_b(self):
        t_start = 10
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [95, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 90)
        self.assertEqual(a, [True])

    def test_more_bin_test_intervals_only_in_the_bin(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70], [80, 90],
                     [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 10)
        self.assertEqual(v, [[4], [3, 2], [], [], [10], [10], [10],
                             [], [10], []])

    def test_more_bin_test_intervals_only_in_the_bin_b(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 5], [12, 15], [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 10)
        self.assertEqual(a, [False, False, False, False, False, True,
                             True, False, False, False])

    def test_more_bin_test_intervals_starting_before_the_bin(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 10)
        self.assertEqual(v,
                         [[9], [1, 3, 2], [], [],
                          [10], [10], [10], [], [10], []])

    def test_more_bin_test_intervals_starting_before_the_bin_b(self):
        t_start = 0
        t_stop = 100
        intervals = [[1, 11], [12, 15],
                     [18, 20], [40, 70],
                     [80, 90], [110, 130]]
        v, a = cv.get_visits_in_bins(intervals, t_start, t_stop, 10)
        self.assertEqual(a, [False, True, False, False, False, True,
                             True, False, False, False])


class TestCalcVisitsPerMouse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        t_start = 0
        t_stop = 100
        ints = [[1, 11], [12, 15],
                [18, 20], [40, 70],
                [80, 90], [110, 130]]
        binsize = 10
        cls.vis, cls.dur, cls.all_v = cv.calc_visit_per_mouse(ints,
                                                              t_start,
                                                              t_stop,
                                                              binsize)

    def test_visits(self):
        out = [1, 2, 0, 0, 1, 0, 0, 0, 1, 0]
        self.assertEqual(out, self.vis)

    def test_durations(self):
        durations = [9, 6, 0, 0, 10, 10, 10, 0, 10, 0]
        self.assertEqual(durations, self.dur)

    def test_all_visits(self):
        all_vis = [[9], [1, 3, 2], [], [],
                   [10], [10], [10], [], [10], []]
        self.assertEqual(all_vis, self.all_v)


class TestCalculateVisitsDurations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = {
            "mouse_1": [
                ["A", 1, 11],
                ["A", 12, 15],
                ["A", 18, 20],
                ["A", 40, 70],
                ["A", 80, 90],
                ["A", 110, 130]
            ],
            "mouse_2": [
                ["A", 1, 11],
                ["B", 12, 15],
                ["A", 18, 20],
                ["B", 40, 70],
                ["A", 80, 90],
                ["B", 110, 130]
            ]
        }
        t_start = 0
        t_stop = 100
        binsize = 10
        address = "A"
        cls.mice = ["mouse_1", "mouse_2"]
        cls.vis, cls.dur,\
            cls.all_v = cv.calculate_visits_and_durations(data,
                                                          cls.mice,
                                                          address,
                                                          t_start,
                                                          t_stop,
                                                          binsize)
        cls.visB, cls.durB,\
            cls.all_vB = cv.calculate_visits_and_durations(data,
                                                           cls.mice,
                                                           "B",
                                                           t_start,
                                                           t_stop,
                                                           binsize)

    def test_keys_1(self):
        self.assertEqual(sorted(self.vis.keys()), sorted(self.mice))

    def test_keys_2(self):
        self.assertEqual(sorted(self.dur.keys()), sorted(self.mice))

    def test_keys_3(self):
        self.assertEqual(sorted(self.all_v.keys()), sorted(self.mice))

    def test_vis_mouse1(self):
        out = [1, 2, 0, 0, 1, 0, 0, 0, 1, 0]
        self.assertEqual(out, self.vis["mouse_1"])

    def test_durations_mouse1(self):
        durations = [9, 6, 0, 0, 10, 10, 10, 0, 10, 0]
        self.assertEqual(durations, self.dur["mouse_1"])

    def test_all_visits_mouse1(self):
        all_vis = [[9], [1, 3, 2], [], [],
                   [10], [10], [10], [], [10], []]
        self.assertEqual(all_vis, self.all_v["mouse_1"])

    def test_vis_mouse2(self):
        out = [1, 1, 0, 0, 0, 0, 0, 0, 1, 0]
        self.assertEqual(out, self.vis["mouse_2"])

    def test_durations_mouse2(self):
        durations = [9, 3, 0, 0, 0, 0, 0, 0, 10, 0]
        self.assertEqual(durations, self.dur["mouse_2"])

    def test_all_visits_mouse2(self):
        all_vis = [[9], [1, 2], [], [],
                   [], [], [], [], [10], []]
        self.assertEqual(all_vis, self.all_v["mouse_2"])

    def test_vis_B_mouse1(self):
        out = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(out, self.visB["mouse_1"])

    def test_durations_B_mouse1(self):
        durations = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(durations, self.durB["mouse_1"])

    def test_all_B_visits_mouse1(self):
        all_vis = [[], [], [], [], [], [], [], [], [], []]
        self.assertEqual(all_vis, self.all_vB["mouse_1"])

    def test_vis_B_mouse2(self):
        out = [0, 1, 0, 0, 1, 0, 0, 0, 0, 0]
        self.assertEqual(out, self.visB["mouse_2"])

    def test_durations_B_mouse2(self):
        durations = [0, 3, 0, 0, 10, 10, 10, 0, 0, 0]
        self.assertEqual(durations, self.durB["mouse_2"])

    def test_all_visits_B_mouse2(self):
        all_vis = [[], [3], [], [],
                   [10], [10], [10], [], [], []]
        self.assertEqual(all_vis, self.all_vB["mouse_2"])


if __name__ == '__main__':
    unittest.main()
