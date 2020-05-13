from __future__ import print_function, division, absolute_import
import random
import unittest
import os
from pyEcoHAB import following as fol
from pyEcoHAB import utility_functions as uf
from pyEcoHAB import Loader
from pyEcoHAB import ExperimentConfigFile
from pyEcoHAB import sample_data

try:
    basestring
except NameError:
    basestring = str


class TestFollowing2ndMouseInPipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas1 = [1, 2]
        times1 = [15, 16.5]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        dir1 = uf.extract_directions(times1, antennas1, 3)
        dir2 = uf.extract_directions(times2, antennas2, 6)
        res = fol.following_single_pair(dir1, dir2)
        cls.out1, cls.time_together1, cls.intervals1= res

        antennas1 = [1, 2, 3, 4, 5]
        times1 = [15, 16.5, 19, 20, 21]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        dir1 = uf.extract_directions(times1, antennas1, 6)
        dir2 = uf.extract_directions(times2, antennas2, 6)
        res = fol.following_single_pair(dir2, dir1)
        cls.out2, cls.time_together2, cls.intervals2 = res

        antennas1 = [1, 2,   3, 4,  5,  6,  7,   8,  1, 2]
        times1 = [15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34 ]
        antennas2 = [8, 1,   2,    3,  4,  5,  6,   7, 8,   1,   2]
        times2 =   [10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35]
        dir1 = uf.extract_directions(times1, antennas1, 3)
        dir2 = uf.extract_directions(times2, antennas2, 3)
        res = fol.following_single_pair(dir1, dir2)
        cls.out3, cls.time_together3, cls.intervals3 = res

    def test_following_11(self):
        self.assertEqual(self.out1, 1)

    def test_following_12(self):
        self.assertEqual(self.time_together1, 0.5)

    def test_following_13(self):
        self.assertEqual(self.intervals1, [4])

    def test_not_following_more_1(self):
        self.assertEqual(self.out2, 0)

    def test_not_following_more_2(self):
        self.assertEqual(self.time_together2, 0)

    def test_not_following_more_3(self):
        self.assertEqual(self.intervals2, [])

    def test_following_31(self):
        self.assertEqual(self.out3, 3)

    def test_following_32(self):
        self.assertEqual(self.time_together3, 4)

    def test_following_33(self):
        self.assertEqual(self.intervals3, [4, 6, 3])


class TestFollowingMatrices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ta = {"mouse1": [[15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34 ],
                         [1, 2,   3,    4,  5,  6,  7,   8,  1, 2]],
              "mouse2": [[10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35],
                         [8, 1,   2,    3,  4,  5,  6,   7, 8,   1,   2],],
              "mouse3": [[10, 16, 17, 18, 22, 25, 26, 27],
                         [1,  2,   3, 4,  4,  3,  2, 1],]
        }
        mice_list = ["mouse1", "mouse2", "mouse3"]
        directions_dict = {}
        last = {
            "mouse1": 3,
            "mouse2": 3,
            "mouse3": 8
        }
        for mouse in mice_list:
            directions_dict[mouse] = uf.extract_directions(ta[mouse][0],
                                                           ta[mouse][1],
                                                           last[mouse])
        out = fol.following_matrices(directions_dict,
                                    mice_list,
                                    0, 1000)
        cls.following = out[0]
        cls.time_together = out[1]
        cls.interval_details = out[2]

    def test_following_diag_1(self):
        self.assertEqual(self.following["mouse1"]["mouse1"], 0)

    def test_following_diag_2(self):
        self.assertEqual(self.following["mouse2"]["mouse2"], 0)

    def test_following_diag_3(self):
        self.assertEqual(self.following["mouse3"]["mouse3"], 0)

    def test_following_0_1(self):
        self.assertEqual(self.following["mouse1"]["mouse2"], 3)

    def test_following_0_2(self):
        self.assertEqual(self.following["mouse1"]["mouse3"], 0)

    def test_following_1_0(self):
        self.assertEqual(self.following["mouse2"]["mouse1"], 0)

    def test_following_1_2(self):
        self.assertEqual(self.following["mouse2"]["mouse3"], 0)

    def test_following_2_0(self):
        self.assertEqual(self.following["mouse3"]["mouse1"], 1)

    def test_following_2_1(self):
        self.assertEqual(self.following["mouse3"]["mouse2"] , 0)

    def test_time_together_diag_0(self):
        self.assertEqual(self.time_together["mouse1"]["mouse1"], 0)

    def test_time_together_diag_1(self):
        self.assertEqual(self.time_together["mouse2"]["mouse2"], 0)

    def test_time_together_diag_2(self):
        self.assertEqual(self.time_together["mouse3"]["mouse3"], 0)

    def test_time_together_0_1(self):
        self.assertEqual(self.time_together["mouse1"]["mouse2"], 0.004)

    def test_time_together_0_2(self):
        self.assertEqual(self.time_together["mouse1"]["mouse3"], 0)

    def test_time_together_1_0(self):
        self.assertEqual(self.time_together["mouse2"]["mouse1"], 0)

    def test_time_together_1_2(self):
        self.assertEqual(self.time_together["mouse2"]["mouse3"], 0)

    def test_time_together_2_0(self):
        self.assertEqual(self.time_together["mouse3"]["mouse1"], 0.001)

    def test_time_together_2_1(self):
        self.assertEqual(self.time_together["mouse3"]["mouse2"], 0)

    def test_interval_details_empty(self):
        self.assertEqual(self.time_together["mouse3"]["mouse2"], 0)


class TestInsertInterval(unittest.TestCase):
    def insert_first_1(self):
        t_starts = []
        t_ends = []
        out = fol.insert_interval(0, 1, t_starts, t_ends, 10)
        self.assertEqual(out, 1)

    def insert_first_2(self):
        t_starts = []
        t_ends = []
        out = fol.insert_interval(0, 1, t_starts, t_ends, 10)
        self.assertEqual(t_starts, [0])

    def insert_first_3(self):
        t_starts = []
        t_ends = []
        out = fol.insert_interval(0, 1, t_starts, t_ends, 10)
        self.assertEqual(t_ends, [1])

    def no_insert_first(self):
        t_starts = []
        t_ends = []
        out = fol.insert_interval(0, 1, t_starts, t_ends, .5)
        self.assertEqual(out, 0)

    def test_t_start_in_starts(self):
        t_starts = [1, 5, 10]
        t_ends = [2, 8, 11]
        out = fol.insert_interval(1, 1, t_starts, t_ends, 20)
        self.assertEqual(out, 0)

    def test_t_ends_in_ends(self):
        t_starts = [1, 5, 10]
        t_ends = [2, 8, 11]
        out = fol.insert_interval(0, 2, t_starts, t_ends, 20)
        self.assertEqual(out, 0)

    def test_insert_at_the_beginning_1(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(0, 2, t_starts, t_ends, 20)
        self.assertEqual(out, 1)

    def test_insert_at_the_beginning_2(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(0, 2, t_starts, t_ends, 20)
        self.assertEqual(t_starts, [0, 3, 5, 10])

    def test_insert_at_the_beginning_3(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(0, 2, t_starts, t_ends, 20)
        self.assertEqual(t_ends, [2, 4, 8, 11])

    def test_insert_at_the_end_1(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(12, 2, t_starts, t_ends, 20)
        self.assertEqual(out, 1)

    def test_insert_at_the_end_2(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(12, 2, t_starts, t_ends, 20)
        self.assertEqual(t_starts, [3, 5, 10, 12])

    def test_insert_at_the_end_3(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(12, 2, t_starts, t_ends, 20)
        self.assertEqual(t_ends, [4, 8, 11, 14])


    def test_does_not_fit_1(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(4.5, 2, t_starts, t_ends, 20)
        self.assertEqual(out, 0)

    def test_does_not_fit_2(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(8, 2, t_starts, t_ends, 20)
        self.assertEqual(out, 0)

    def test_in_the_middle_1(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(8.5, 1, t_starts, t_ends, 20)
        self.assertEqual(out, 1)

    def test_in_the_middle_2(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(8.5, 1, t_starts, t_ends, 20)
        self.assertEqual(t_starts, [3, 5, 8.5, 10])

    def test_in_the_middle_3(self):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 11]
        out = fol.insert_interval(8.5, 1, t_starts, t_ends, 20)
        self.assertEqual(t_ends, [4, 8, 9.5, 11])


class TestIntervalGeneration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        t_starts = [3, 5, 10]
        t_ends = [4, 8, 12]
        cls.intervals = set(uf.get_interval_durations_2_lists(t_starts,
                                                              t_ends))
        duration = 40
        random.seed(1)
        cls.out1 = fol.generate_intervals(t_starts, t_ends, duration)
        random.seed(100)
        cls.out2 = fol.generate_intervals(t_starts, t_ends, duration)

    def test_length_1(self):
        self.assertEqual(len(self.out1[0]), 3)

    def test_length_2(self):
        self.assertEqual(len(self.out1[1]), 3)

    def test_length_3(self):
        self.assertEqual(len(self.out2[0]), 3)

    def test_length_4(self):
        self.assertEqual(len(self.out2[1]), 3)

    def test_intervals_1(self):
        intervals = set(uf.get_interval_durations_2_lists(self.out1[0],
                                                          self.out1[1]))
        self.assertEqual(intervals, self.intervals)

    def test_intervals_2(self):
        intervals = set(uf.get_interval_durations_2_lists(self.out2[0],
                                                          self.out2[1]))
        self.assertEqual(intervals, self.intervals)

    def test_different_1(self):
        self.assertFalse(self.out1[0] == self.out2[0])

    def test_different_2(self):
        self.assertFalse(self.out1[1] == self.out2[1])

    def test_different_2(self):
        ints1 = uf.get_interval_durations_2_lists(self.out1[0],
                                                  self.out1[1])
        ints2 = uf.get_interval_durations_2_lists(self.out2[0],
                                                  self.out2[1])
        self.assertFalse(ints1 == ints2)


class TestExecution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = Loader(sample_data)
        cls.config = ExperimentConfigFile(sample_data)

    def test_phases(self):
        fol.get_dynamic_interactions(self.data, self.config, 1,
                                     save_distributions=True,
                                     save_figures=False, return_median=False,
                                     delimiter=";",
                                     save_times_following=False)

    def test_ALL(self):
        fol.get_dynamic_interactions(self.data, self.config, 1, binsize="ALL")

    def test_short_bin(self):
        fol.get_dynamic_interactions(self.data, self.config, 1, binsize=3600)

    def test_long_bin(self):
        fol.get_dynamic_interactions(self.data, self.config, 1, binsize=24*3600)


if __name__ == '__main__':
    unittest.main()
