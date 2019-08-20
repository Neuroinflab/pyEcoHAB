#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
import pyEcoHAB.utility_functions as uf
import unittest
import numpy as np


class TestFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lista = ["Dark", "DARK", "dark", "Light", "LIGHT", "light", "ALL",
        "Dark 6", "DARK 5", "dark 4", "Light 3", "LIGHT 2", "light 1", "ALL"]

    def test_light(self):
        out = uf.filter_light(self.lista)
        self.assertEqual(out, ["Light", "LIGHT", "light",
                               "Light 3", "LIGHT 2", "light 1"])

    def test_dark(self):
        out = uf.filter_dark(self.lista)
        self.assertEqual(out, ["Dark", "DARK", "dark",
                               "Dark 6", "DARK 5", "dark 4"])

    def test_all(self):
        out = uf.filter_dark_light(self.lista)
        self.assertEqual(out, ["Dark", "DARK", "dark",
                               "Light", "LIGHT", "light",
                               "Dark 6", "DARK 5", "dark 4",
                               "Light 3", "LIGHT 2", "light 1"])


class TestGetMoreStates(unittest.TestCase):
    
    def test_2nd_tier_len_states(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 211.214,
                 211.98,
                 217.953,
                 218.61,
                 223.347,
                 223.769,
                 225.192,
                 225.942,
                 228.772,
                 228.972]
        idx = 6
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), 3)
    def test_2nd_tier_len_states_times(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 211.214,
                 211.98,
                 217.953,
                 218.61,
                 223.347,
                 223.769,
                 225.192,
                 225.942,
                 228.772,
                 228.972]
        idx = 6
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), len(readouts))

    def test_2nd_tier_equal_len_midx(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 211.214,
                 211.98,
                 217.953,
                 218.61,
                 223.347,
                 223.769,
                 225.192,
                 225.942,
                 228.772,
                 228.972]
        idx = 6
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(midx + 1, idx + len(states))

    def test_catch_threshold_len_states(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 211.214,
                 211.98,
                 217.953,
                 218.61,
                 223.347,
                 223.769,
                 225.192,
                 225.942,
                 228.772,
                 228.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), 4 + idx)
    def test_catch_threshold_equal_len_states_times(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 211.214,
                 211.98,
                 217.953,
                 218.61,
                 223.347,
                 223.769,
                 225.192,
                 225.942,
                 228.772,
                 228.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), len(readouts))

    def test_catch_threshold_equal_len_midx(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 211.214,
                 211.98,
                 217.953,
                 218.61,
                 223.347,
                 223.769,
                 225.192,
                 225.942,
                 228.772,
                 228.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), midx + idx)
    
    def test_catch_3rd_len_states(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 111.214,
                 111.98,
                 117.953,
                 118.61,
                 123.347,
                 123.769,
                 125.192,
                 125.942,
                 128.772,
                 128.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), 7)

    def test_catch_3rd_equal_len_states_times(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 111.214,
                 111.98,
                 117.953,
                 118.61,
                 123.347,
                 123.769,
                 125.192,
                 125.942,
                 128.772,
                 128.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), len(readouts))

    def test_catch_3rd_equal_len_midx(self):
        antennas = [5, 6, 6, 5, 5, 6, 7, 6, 5, 5, 6, 6, 6, 7]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 111.214,
                 111.98,
                 117.953,
                 118.61,
                 123.347,
                 123.769,
                 125.192,
                 125.942,
                 128.772,
                 128.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states) - 1, midx + idx)

    def test_catch_end(self):
        antennas = [5, 6, 6, 5, 5, 6, 5, 6, 5, 5, 6, 6, 6, 5]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 111.214,
                 111.98,
                 117.953,
                 118.61,
                 123.347,
                 123.769,
                 125.192,
                 125.942,
                 128.772,
                 128.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), midx + idx)

    def test_catch_equal_len_antennas_states(self):
        antennas = [5, 6, 6, 5, 5, 6, 5, 6, 5, 5, 6, 6, 6, 5]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 111.214,
                 111.98,
                 117.953,
                 118.61,
                 123.347,
                 123.769,
                 125.192,
                 125.942,
                 128.772,
                 128.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), len(antennas))

    def test_catch_equal_len_antennas_readouts(self):
        antennas = [5, 6, 6, 5, 5, 6, 5, 6, 5, 5, 6, 6, 6, 5]
        times = [101.976,
                 103.148,
                 109.37,
                 109.761,
                 111.214,
                 111.98,
                 117.953,
                 118.61,
                 123.347,
                 123.769,
                 125.192,
                 125.942,
                 128.772,
                 128.972]
        idx = 0
        states, readouts, midx = uf.get_more_states(antennas, times, idx, 50, 3)
        self.assertEqual(len(states), len(readouts))


class TestInPipe(unittest.TestCase):
    def test_pipe_1_2_clockwise(self):
        self.assertTrue(uf.in_tube(1, 2))
    def test_pipe_2_1_counterclockwise(self):
        self.assertTrue(uf.in_tube(2, 1))
    def test_pipe_2_3_clockwise(self):
        self.assertFalse(uf.in_tube(2, 3))
    def test_pipe_3_2_counterclockwise(self):
        self.assertFalse(uf.in_tube(3, 2))
    def test_pipe_7_8_clockwise(self):
        self.assertTrue(uf.in_tube(7, 8))
    def test_pipe_8_7_counterclockwise(self):
        self.assertTrue(uf.in_tube(8, 7))
    def test_pipe_1_8_clockwise(self):
        self.assertFalse(uf.in_tube(1, 8))
    def test_pipe_8_1_counterclockwise(self):
        self.assertFalse(uf.in_tube(8, 1))

class TestInChambers(unittest.TestCase):
    def test_pipe_1_2_clockwise(self):
        self.assertFalse(uf.in_chamber(1, 2))
    def test_pipe_2_1_counterclockwise(self):
        self.assertFalse(uf.in_chamber(2, 1))
    def test_pipe_2_3_clockwise(self):
        self.assertTrue(uf.in_chamber(2, 3))
    def test_pipe_3_2_counterclockwise(self):
        self.assertTrue(uf.in_chamber(3, 2))
    def test_pipe_7_8_clockwise(self):
        self.assertFalse(uf.in_chamber(7, 8))
    def test_pipe_8_7_counterclockwise(self):
        self.assertFalse(uf.in_chamber(8, 7))
    def test_pipe_1_8_clockwise(self):
        self.assertTrue(uf.in_chamber(1, 8))
    def test_pipe_8_1_counterclockwise(self):
        self.assertTrue(uf.in_chamber(8, 1))
    
class TestChangeState(unittest.TestCase):
    def test_check_no_change(self):
        self.assertEqual(len(uf.change_state([1, 1, 1, 1,])), 0)
    def test_check_change_1(self):
        self.assertFalse(len(uf.change_state([1, 1, 1, 2])) == 0)
    def test_check_change_2(self):
        self.assertTrue(uf.change_state([1, 1, 1, 2]), np.array([2], dtype=int))
        
class TestGetIdxPre(unittest.TestCase):
    def test_false(self):
        self.assertEqual(uf.get_idx_pre(2, [3, 4, 5]), None)
    def test_correct(self):
        self.assertEqual(uf.get_idx_pre(2, [-1, 0, 1, 2, 3]), 2)
    def test_empty(self):
        self.assertEqual(uf.get_idx_pre(0, []), None)

class TestGetIdxPost(unittest.TestCase):
    def test_false(self):
        self.assertEqual(uf.get_idx_post(2, [-1, 0, 1]), None)
    def test_correct(self):
        self.assertEqual(uf.get_idx_post(2, [-1, 0, 1, 2, 3]), 4)
    def test_empty(self):
        self.assertEqual(uf.get_idx_post(0, []), None)

class TestGetIdxBetween(unittest.TestCase):
    def test_false(self):
        self.assertEqual(len(uf.get_idx_between(2, 3, [-1, 0, 1])), 0)
        
    def test_correct1(self):
        out = uf.get_idx_between(2, 3, [-1, 0, 1, 2, 3])
        res = np.array([3, 4], dtype=int)
        self.assertEqual(len(out), len(res))


    def test_correct_loop(self):
        out = uf.get_idx_between(2, 3, [-1, 0, 1, 2, 3])
        res = np.array([3, 4], dtype=int)
        for i, x in enumerate(out):
            self.assertEqual(x, res[i])

    def test_empty(self):
        self.assertEqual(len(uf.get_idx_between(0, 2, [])), 0)


class TestMouseGoingForward(unittest.TestCase):
    def test_going_forward(self):
        self.assertTrue(uf.mouse_going_forward([1, 2, 3]))

    def test_not_going_forward(self):
        self.assertFalse(uf.mouse_going_forward([1, 2, 1, 8]))

    def test_going_forward_other_pipe(self):
        self.assertTrue(uf.mouse_going_forward([7, 8, 1]))

    def test_not_going_forward_other_pipe(self):
        self.assertFalse(uf.mouse_going_forward([8, 7, 8, 1]))

    def test_going_forward_other_from_data(self):
        self.assertTrue(uf.mouse_going_forward([3, 4, 4, 3, 3, 4, 4, 5]))

    def test_not_going_forward_from_data(self):
        self.assertFalse(uf.mouse_going_forward([5, 6, 6, 5, 4]))
        
class TestMouseBacking(unittest.TestCase):
    def test_going_forward(self):
        self.assertFalse(uf.mouse_backing_off([1, 2, 3]))

    def test_not_going_forward(self):
        self.assertTrue(uf.mouse_backing_off([1, 2, 1, 8]))

    def test_going_forward_other_pipe(self):
        self.assertFalse(uf.mouse_backing_off([7, 8, 1]))

    def test_not_going_forward_other_pipe(self):
        self.assertTrue(uf.mouse_backing_off([8, 7, 8, 1]))

    def test_going_forward_other_from_data(self):
        self.assertFalse(uf.mouse_backing_off([3, 4, 4, 3, 3, 4, 4, 5]))

    def test_not_going_forward_from_data(self):
        self.assertTrue(uf.mouse_backing_off([5, 6, 6, 5, 4]))
        
    def test_two_same_reading(self):
        self.assertFalse(uf.mouse_backing_off([1, 1]))
        
    def test_two_different_readings(self):
        self.assertFalse(uf.mouse_backing_off([1, 2]))

    def test_only_two_antennas_True(self):
        self.assertTrue(uf.mouse_backing_off([1, 2, 1]))

    def test_only_two_antennas_False(self):
        self.assertFalse(uf.mouse_backing_off([1, 2, 1, 2]))

    def test_one_reading(self):
        self.assertFalse(uf.mouse_backing_off([1]))
        

class TestSkipAntennas(unittest.TestCase):
    def test_no_skipped_antennas_short(self):
        out = uf.skipped_antennas([1, 2])
        self.assertFalse(out)
    def test_skipped_antennas_short(self):
        out = uf.skipped_antennas([4, 2])
        self.assertTrue(out)
    def test_no_skipped_antennas_long(self):
        out = uf.skipped_antennas([1, 2, 2, 3, 4, 4, 5])
        self.assertFalse(out)
    def test_skipped_antennas_long(self):
        out = uf.skipped_antennas([1, 2, 2, 3, 4, 4, 6])
        self.assertTrue(out)
    def test_no_skipped_antennas_long_with_8(self):
        out = uf.skipped_antennas([8, 1, 2, 2, 3, 4, 4, 5])
        self.assertFalse(out)
    def test_skipped_antennas_long_with_8(self):
        out = uf.skipped_antennas([1, 2, 2, 3, 4, 4, 6])
        self.assertTrue(out)

class TestChangeOneToSeven(unittest.TestCase):
    def test_change_seven(self):
        out = uf.change_seven_to_one([1, 0, -1, 7])
        self.assertEqual(out[3], -1)
    def test_change_minus_seven(self):
        out = uf.change_seven_to_one([1, 0, -1, -7])
        self.assertEqual(out[3], 1)
    def test_no_change(self):
        values = [1, 2, 8, 3, 5, 6]
        out = uf.change_seven_to_one(values)
        for i, val in enumerate(values):
            self.assertEqual(val, out[i])


class TestMouseGoingClockwise(unittest.TestCase):
    def test_correct(self):
        positions = [3, 4, 4, 3, 3, 4, 4, 5]
        self.assertTrue(uf.mouse_going_clockwise(positions))
        
    def test_incorrect(self):
        positions = [4, 3, 2, 1]
        self.assertFalse(uf.mouse_going_clockwise(positions))

    def test_incorrect2(self):
        positions = [4, 4, 4, 4]
        self.assertFalse(uf.mouse_going_clockwise(positions))

    def test_correct_with_8(self):
        positions = [7, 8, 1, 2, 3, 4, 4, 5]
        self.assertTrue(uf.mouse_going_clockwise(positions))

class TestMouseGoingCounterClockwise(unittest.TestCase):
    def test_incorrect(self):
        positions = [3, 4, 4, 3, 3, 4, 4, 5]
        self.assertFalse(uf.mouse_going_counterclockwise(positions))
        
    def test_correct(self):
        positions = [4, 3, 2, 1]
        self.assertTrue(uf.mouse_going_counterclockwise(positions))

    def test_incorrect2(self):
        positions = [4, 4, 4, 4]
        self.assertFalse(uf.mouse_going_counterclockwise(positions))

    def test_incorrect_with_8(self):
        positions = [7, 8, 1, 2, 3, 4, 4, 5]
        self.assertFalse(uf.mouse_going_counterclockwise(positions))

    def test_correct_with_8(self):
        positions = [1, 8, 7]
        self.assertTrue(uf.mouse_going_counterclockwise(positions))


class TestGetTimestamp(unittest.TestCase):

    def test_up(self):
        self.assertEqual(uf.get_timestamp(42.2, 45.76, 0.1), 36)

    def test_down(self):
        self.assertEqual(uf.get_timestamp(42.2, 45.71, 0.1), 35) 


class TestGetAnennas(unittest.TestCase):
    def test(self):
        out = uf.get_antennas([0, 1, 4], [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertEqual(out, [1, 2, 5])


class TestGetKeyForFrequency(unittest.TestCase):
    def test_78(self):
        out = uf.get_key_for_frequencies(7, 8)
        self.assertEqual(out, 15)

    def test_87(self):
        out = uf.get_key_for_frequencies(8, 7)
        self.assertEqual(out, 15)

    def test_chamber_18(self):
        out = uf.get_key_for_frequencies(8, 1)
        self.assertEqual(out, None)

    def test_chamber_23(self):
        out = uf.get_key_for_frequencies(2, 3)
        self.assertEqual(out, None)

    def test_stupid_(self):
        out = uf.get_key_for_frequencies(2, 4)
        self.assertEqual(out, None)


class TestIntervalOverlap(unittest.TestCase):

    def test_interval_overlap_1(self):
        """Incorrect interval
        """
        inte_1 = [34, 45]
        inte_2 = [34, 23]
        self.assertTrue(uf.interval_overlap(inte_1, inte_2) == 0)

    def test_interval_overlap_2(self):
        """2nd interval shorter than the first one, same start
        """
        inte_1 = [34, 45]
        inte_2 = [34, 43]
        self.assertTrue(uf.interval_overlap(inte_1, inte_2) == 9)

    def test_interval_overlap_3(self):
        """2nd interval beginning after the 1st interval has finished
        """
        inte_1 = [34, 45]
        inte_2 = [46, 50]
        self.assertTrue(uf.interval_overlap(inte_1, inte_2) == 0)

    def test_interval_overlap_4(self):
        """1st interval beginning after the 2nd interval has finished
        """
        inte_1 = [46, 50]
        inte_2 = [34, 45]
        self.assertTrue(uf.interval_overlap(inte_1, inte_2) == 0)

    def test_interval_overlap_5(self):
        """Overlapping
        """
        inte_1 = [46, 50]
        inte_2 = [34, 48]
        self.assertTrue(uf.interval_overlap(inte_1, inte_2) == 2)

    def test_interval_overlap_6(self):
        """Overlapping
        """
        inte_1 = [34, 48]
        inte_2 = [46, 50]
        self.assertTrue(uf.interval_overlap(inte_1, inte_2) == 2)

class TestGetStatesForEhs(unittest.TestCase):
    def test_chambers(self):
        out = uf.get_states_for_ehs([1, 5], [2, 3], "mouse 1", 2)
        res = [(1, "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_chambers_8(self):
        out = uf.get_states_for_ehs([1, 5], [1, 8], "mouse 1", 2)
        res = [(4, "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_too_fast(self):
        out = uf.get_states_for_ehs([1, 2], [2, 3], "mouse 1", 2)
        self.assertEqual(out, [])

    def test_pipe(self):
        out = uf.get_states_for_ehs([1, 5], [2, 1], "mouse 1", 2)
        self.assertEqual(out, [])

    def test_same_chamber(self):
        out = uf.get_states_for_ehs([1, 5], [1, 1], "mouse 1", 2)
        res = [(4, "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_skipped_antenna(self):
        out = uf.get_states_for_ehs([1, 5], [1, 3], "mouse 1", 2)
        res = [(1, "mouse 1", 1, 5, 4, False)]
        self.assertEqual(out, res)

    def test_skipped_antenna_2(self):
        out = uf.get_states_for_ehs([1, 5], [7, 1], "mouse 1", 2)
        res = [(4, "mouse 1", 1, 5, 4, False)]
        self.assertEqual(out, res)

    def test_opposite_pipe_1(self):
        antenna = 2
        out = uf.get_states_for_ehs([1, 5], [antenna, antenna + 3],
                                    "mouse 1", 2)
        self.assertEqual(out, [])

    def test_not_opposite_pipe_1(self):
        out = uf.get_states_for_ehs([1, 5], [2, 7],
                                    "mouse 1", 2)
        self.assertEqual(out, [(4, "mouse 1", 1, 5, 4, False)])

    def test_longer(self):
        out = uf.get_states_for_ehs([2, 7, 23, 45, 55, 61],
                                    [1, 2, 3, 4, 5, 6],
                                    "mouse 1", 2)
        res = [(1, "mouse 1", 7, 23, 16, True),
               (2, "mouse 1", 45, 55, 10, True)]
        self.assertEqual(out, res)


class TestGetIntervals(unittest.TestCase):

    def test_empty(self):
        self.assertEqual(uf.get_intervals([[1, 2, 3]], 4), [])

    def test_not_empty(self):
        self.assertEqual(uf.get_intervals([[1, 2, 3]], 1), [[2, 3]])

    def test_longer(self):
        test_input = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [1, 10, 12],
        ]
        test_output = [
            [2, 3],
            [10, 12],
        ]
        self.assertEqual(uf.get_intervals(test_input, 1), test_output)

class TestIntervals2Lists(unittest.TestCase):

    def test_empty(self):
        self.assertEqual(uf.intervals2lists([[1, 2, 3]], 4), [[], []])

    def test_not_empty(self):
        self.assertEqual(uf.intervals2lists([[1, 2, 3]], 1), [[2], [3]])

    def test_longer(self):
        test_input = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [1, 10, 12],
        ]
        test_output = [
            [2, 10],
            [3, 12],
        ]
        self.assertEqual(uf.intervals2lists(test_input, 1), test_output)


class TestGetIndices(unittest.TestCase):
    def test_empty(self):
        starts = [1, 5, 10]
        ends = [4, 7, 15]
        out = uf.get_indices(16, 20, starts, ends)
        self.assertEqual(out, [])

    def test_slice_overlap(self):
        starts = [1, 5, 10]
        ends = [4, 7, 15]
        out = uf.get_indices(3, 12, starts, ends)
        self.assertEqual(out, [0, 1, 2])

    def test_slice(self):
        starts = [1, 5, 10]
        ends = [4, 7, 15]
        out = uf.get_indices(4, 10, starts, ends)
        self.assertEqual(out, [0, 1, 2])

class TestGetDurations(unittest.TestCase):
    def test_1(self):
        a = [1, 5, 8]
        b = [2, 7, 10]
        self.assertEqual(uf.get_duration(a, b), 5)

    def test_1(self):
        a = [1, 5, 8]
        b = [2, 7, 10]
        self.assertEqual(uf.get_duration(b, a), 5)


class TestIntervalGetDurations(unittest.TestCase):
    def test_1(self):
        a = [[1, 2],  [5, 7], [8, 10]]
        self.assertEqual(uf.get_interval_durations(a), [1, 2, 2])

class TestCalculateTotalDuration(unittest.TestCase):
    def test_1(self):
        a = [[1, 2],  [5, 7], [8, 10]]
        self.assertEqual(uf.calculate_total_duration(a), 5)

class TestGetMice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lista = ["Zdzisio",
                     "Zbysio",
                     "Henio",
                     "Gienio"]

    def test_check_none(self):
        self.assertEqual(uf.get_mice(self.lista, None),
                         self.lista)

    def test_bogus(self):
        self.assertEqual(uf.get_mice(self.lista, 5),
                         self.lista)

    def test_string(self):
        self.assertEqual(uf.get_mice(self.lista, 'string'),
                         self.lista)

    def test_string_2(self):
        self.assertEqual(uf.get_mice(self.lista, 'Gienio'),
                         ["Zdzisio", "Zbysio", "Henio"])

    def test_string_2(self):
        self.assertEqual(uf.get_mice(self.lista, ['Genio']),
                         self.lista)

    def test_lista_2(self):
        self.assertEqual(uf.get_mice(self.lista, ['Gienio',
                                                  'Zdzisio']),
                         ["Zbysio", "Henio"])
        

class TestAddInfo(unittest.TestCase):
    def test_None(self):
        self.assertEqual(uf.add_info_mice_filename(None), '')

    def test_one(self):
        self.assertEqual(uf.add_info_mice_filename("Zdzisio"),
                         'remove_Zdzisio')
    def test_more(self):
        lista = ["Zdzisio",
                 "Zbysio",
                 "Henio",
                 "Gienio"]
        self.assertEqual(uf.add_info_mice_filename(lista),
                         'remove_Zdzisio_Zbysio_Henio_Gienio')

class TestListOfpairs(unittest.TestCase):
    def test_list(self):
        lista = ["Zdzisio",
                 "Zbysio",
                 "Henio",
                 "Gienio"]
        out = uf.list_of_pairs(lista)
        self.assertEqual(out, ['Zdzisio|Zbysio', 'Zdzisio|Henio',
                               'Zdzisio|Gienio', 'Zbysio|Henio',
                               'Zbysio|Gienio', 'Henio|Gienio'])

    def test_list_2(self):
        lista = ["Zdzisio",
                 "Zbysio",
                 "Henio",
                 "Gienio",
                 "Rysio"]
        out = uf.list_of_pairs(lista)
        self.assertEqual(out, ['Zdzisio|Zbysio', 'Zdzisio|Henio',
                               'Zdzisio|Gienio', 'Zdzisio|Rysio',
                               'Zbysio|Henio', 'Zbysio|Gienio',
                               'Zbysio|Rysio', 'Henio|Gienio',
                               'Henio|Rysio', 'Gienio|Rysio'])

class TestMakeTableOfPairs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lista = ["Zdzisio",
                     "Zbysio",
                     "Henio",
                     "Gienio"]
        cls.phases = ["0", "1", "2"]
        cls.data = np.zeros((len(cls.phases),
                             len(cls.lista),
                             len(cls.lista)))
        cls.data[0, 0, 1] = 1
        cls.data[0, 0, 2] = 2
        cls.data[0, 0, 3] = 3
        cls.data[0, 1, 2] = 4
        cls.data[0, 1, 3] = 5
        cls.data[0, 2, 3] = 6
        cls.data[1, 0, 1] = 7
        cls.data[1, 0, 2] = 8
        cls.data[1, 0, 3] = 9
        cls.data[1, 1, 2] = 10
        cls.data[1, 1, 3] = 11
        cls.data[1, 2, 3] = 12
        cls.data[2, 0, 1] = 1
        cls.data[2, 0, 2] = 2
        cls.data[2, 0, 3] = 3
        cls.data[2, 1, 2] = 4
        cls.data[2, 1, 3] = 5
        cls.data[2, 2, 3] = 6
        cls.out_data, cls.out = uf.make_table_of_pairs(cls.data,
                                                       cls.phases,
                                                       cls.lista)
    def test_shape(self):
        out_lista = uf.list_of_pairs(self.lista)
        self.assertEqual(self.out_data.shape, (len(out_lista), len(self.phases)))
    def test_lista(self):
        out_lista = uf.list_of_pairs(self.lista)
        self.assertEqual(out_lista, self.out)

    def test_1st_phase(self):
        out = [1, 2, 3, 4, 5, 6]
        first_column = self.out_data[:, 0].tolist()
        self.assertEqual(out, first_column)

    def test_2nd_phase(self):
        out = [7, 8, 9, 10, 11, 12]
        first_column = self.out_data[:, 1].tolist()
        self.assertEqual(out, first_column)

    def test_3rd_phase(self):
        out = [1, 2, 3, 4, 5, 6]
        first_column = self.out_data[:, 2].tolist()
        self.assertEqual(out, first_column)


class TestGetLength(unittest.TestCase):
    def test_1(self):
        out = uf.get_length(0, 14, 7)
        self.assertEqual(out, 2)

    def test_2(self):
        out = uf.get_length(0, 13, 7)
        self.assertEqual(out, 2)

    def test_3(self):
        out = uf.get_length(0, 15, 7)
        self.assertEqual(out, 3)


class TestParseFilename(unittest.TestCase):
    def test_normal(self):
        fname = "20190403_120000.txt"
        hour, date, datenext = uf.parse_fname(fname)
        self.assertEqual(hour, "120000")
        self.assertEqual(date, "20190403")
        self.assertEqual(datenext, "20190404")

    def test_weird(self):
        fname = "20190403_120000_0001.txt"
        hour, date, datenext = uf.parse_fname(fname)
        self.assertEqual(hour, "120000")
        self.assertEqual(date, "20190403")
        self.assertEqual(datenext, "20190404")

    def test_throw_exception(self):
        fname = "20190403_120000_0001_kk.txt"
        self.assertRaises(ValueError, uf.parse_fname, fname=fname)


class TestPrintHumanTime(unittest.TestCase):
    def test_date(self):
        tt = 1554247067
        self.assertTrue('Wed Apr  3 00:00:00 2019', uf.print_human_time(tt))


class TestTimeToSec(unittest.TestCase):
    def test_sec(self):
        string = "20190709 20:05:13.333"
        self.assertTrue(1562695513.333, uf.time_to_sec(string))

    def test_sec_2(self):
        string = "20190709 20:05:13"
        self.assertTrue(1562695513., uf.time_to_sec(string))

    def test_sec_raise_1(self):
        string = "2019070920:05:13"
        self.assertRaises(ValueError, uf.time_to_sec, tt=string)

    def test_sec_raise_2(self):
        string = "20190709.20:05:13"
        self.assertRaises(ValueError, uf.time_to_sec, tt=string)


if __name__ == '__main__':
    unittest.main()
