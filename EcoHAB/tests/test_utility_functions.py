#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
from EcoHAB import utility_functions as uf
import unittest
import numpy as np

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

if __name__ == '__main__':
    unittest.main()
