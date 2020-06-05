#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
import os
import unittest
import random
import time
import calendar
from collections import OrderedDict
import numpy as np

import pyEcoHAB.utility_functions as uf
from pyEcoHAB import data_path
from pyEcoHAB import Loader
from pyEcoHAB import ExperimentConfigFile

class TestFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lista = ["Dark", "DARK", "dark",
                     "Light", "LIGHT", "light", "ALL",
                     "Dark 6", "DARK 5", "dark 4",
                     "Light 3", "LIGHT 2", "light 1", "ALL"]

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

    
class TestChangeState(unittest.TestCase):
    def test_check_no_change(self):
        self.assertEqual(len(uf.change_state([1, 1, 1, 1,])), 0)
    def test_check_change_1(self):
        self.assertFalse(len(uf.change_state([1, 1, 1, 2])) == 0)
    def test_check_change_2(self):
        self.assertEqual(uf.change_state([1, 1, 1, 2]).tolist(),
                         np.array([2], dtype=int).tolist())
        
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
        res = np.array([3], dtype=int)
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
        self.assertEqual(out, "78")

    def test_87(self):
        out = uf.get_key_for_frequencies(8, 7)
        self.assertEqual(out, "87")

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
        out = uf.get_animal_position([1, 5], [2, 3], "mouse 1", 2)
        res = [(1, "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_chambers_8(self):
        out = uf.get_animal_position([1, 5], [1, 8], "mouse 1", 2)
        res = [(4, "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_too_fast(self):
        out = uf.get_animal_position([1, 2], [2, 3], "mouse 1", 2)
        self.assertEqual(out, [])

    def test_pipe(self):
        out = uf.get_animal_position([1, 5], [2, 1], "mouse 1", 2)
        self.assertEqual(out, [])

    def test_same_chamber(self):
        out = uf.get_animal_position([1, 5], [1, 1], "mouse 1", 2)
        res = [(4, "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_skipped_antenna(self):
        out = uf.get_animal_position([1, 5], [1, 3], "mouse 1", 2)
        res = [(1, "mouse 1", 1, 5, 4, False)]
        self.assertEqual(out, res)

    def test_skipped_antenna_2(self):
        out = uf.get_animal_position([1, 5], [7, 1], "mouse 1", 2)
        res = [(4, "mouse 1", 1, 5, 4, False)]
        self.assertEqual(out, res)

    def test_opposite_pipe_1(self):
        antenna = 2
        out = uf.get_animal_position([1, 5], [antenna, antenna + 3],
                                    "mouse 1", 2)
        self.assertEqual(out, [])

    def test_not_opposite_pipe_1(self):
        out = uf.get_animal_position([1, 5], [2, 7],
                                    "mouse 1", 2)
        self.assertEqual(out, [(4, "mouse 1", 1, 5, 4, False)])

    def test_longer(self):
        out = uf.get_animal_position([2, 7, 23, 45, 55, 61],
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
        self.assertEqual(out, [0, 1])

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


class TestIntervalGetDurations2_lists(unittest.TestCase):
    def test_1(self):
        a = [1, 5, 8]
        b = [2, 7, 10]
        self.assertEqual(uf.get_interval_durations_2_lists(a, b), [1, 2, 2])


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


class TestAllPairs(unittest.TestCase):
    def test_list_actual(self):
        lista = ["Zdzisio",
                 "Zbysio",
                 "Henio"]
        out = uf.all_pairs(lista)
        expected = ["Zdzisio|Zbysio", "Zdzisio|Henio", "Zbysio|Zdzisio",
                    "Zbysio|Henio", "Henio|Zdzisio", "Henio|Zbysio"]
        self.assertEqual(out, expected)

    def test_list_reverse(self):
        lista = ["Zdzisio",
                 "Zbysio",
                 "Henio",
                 ]
        expected = ["Zbysio|Zdzisio", "Henio|Zdzisio", "Zdzisio|Zbysio",
                    "Henio|Zbysio", "Zdzisio|Henio", "Zbysio|Henio"]
        out = uf.all_pairs(lista, reverse_order="True")
        self.assertEqual(out, expected)


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


class TestMakeTableOfAllPairs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lista = ["Zdzisio",
                     "Zbysio",
                     "Henio"]
        cls.phases = ["0"]
        cls.data = np.zeros((len(cls.phases),
                             len(cls.lista),
                             len(cls.lista)))
        cls.data[0, 0, 1] = 1
        cls.data[0, 0, 2] = 2
        cls.data[0, 1, 2] = 4
        cls.out_data, cls.out = uf.make_table_of_all_pairs(cls.data,
                                                           cls.phases,
                                                           cls.lista,
                                                           reverse_order=False)
        cls.out_data_rev, cls.out_rev = uf.make_table_of_all_pairs(cls.data,
                                                                   cls.phases,
                                                                   cls.lista,
                                                                   reverse_order=True)
    def test_shape(self):
        out_lista = uf.all_pairs(self.lista)
        self.assertEqual(self.out_data.shape,
                         (len(out_lista), len(self.phases)))

    def test_lista(self):
        out_lista = uf.all_pairs(self.lista)
        self.assertEqual(out_lista, self.out)

    def test_lista_rev(self):
        out_lista = uf.all_pairs(self.lista,
                                 reverse_order=True)
        self.assertEqual(out_lista, self.out_rev)

    def test_1st_phase(self):
        out = [1., 2., 0., 4., 0., 0.]
        first_column = self.out_data[:, 0].tolist()
        self.assertEqual(out, first_column)

    def test_1st_phase_rev(self):
        out = [0., 0., 1., 0., 2., 4.]
        first_column = self.out_data_rev[:, 0].tolist()
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



class TestGETEHSData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lista = ["Dark", "DARK", "dark", "Light", "LIGHT", "light", "ALL",
                     "Dark 6", "DARK 5", "dark 4", "Light 3", "LIGHT 2",
                     "light 1", "ALL"]

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
        res = np.array([3], dtype=int)
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
        self.assertEqual(out, "78")

    def test_87(self):
        out = uf.get_key_for_frequencies(8, 7)
        self.assertEqual(out, "87")

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
        out = uf.get_animal_position([1, 5], [2, 3], "mouse 1", 2)
        res = [("B", "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_chambers_8(self):
        out = uf.get_animal_position([1, 5], [1, 8], "mouse 1", 2)
        res = [("A", "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_too_fast(self):
        out = uf.get_animal_position([1, 2], [2, 3], "mouse 1", 2)
        self.assertEqual(out, [])

    def test_pipe(self):
        out = uf.get_animal_position([1, 5], [2, 1], "mouse 1", 2)
        self.assertEqual(out, [])

    def test_same_chamber(self):
        out = uf.get_animal_position([1, 5], [1, 1], "mouse 1", 2)
        res = [("A", "mouse 1", 1, 5, 4, True)]
        self.assertEqual(out, res)

    def test_skipped_antenna(self):
        out = uf.get_animal_position([1, 5], [1, 3], "mouse 1", 2)
        res = [("B", "mouse 1", 1, 5, 4, False)]
        self.assertEqual(out, res)

    def test_skipped_antenna_2(self):
        out = uf.get_animal_position([1, 5], [7, 1], "mouse 1", 2)
        res = [("A", "mouse 1", 1, 5, 4, False)]
        self.assertEqual(out, res)

    def test_opposite_pipe_1(self):
        antenna = 2
        out = uf.get_animal_position([1, 5], [antenna, antenna + 3],
                                    "mouse 1", 2)
        self.assertEqual(out, [])

    def test_not_opposite_pipe_1(self):
        out = uf.get_animal_position([1, 5], [2, 7],
                                    "mouse 1", 2)
        self.assertEqual(out, [("A", "mouse 1", 1, 5, 4, False)])

    def test_longer(self):
        out = uf.get_animal_position([2, 7, 23, 45, 55, 61],
                                    [1, 2, 3, 4, 5, 6],
                                    "mouse 1", 2)
        res = [("B", "mouse 1", 7, 23, 16, True),
               ("C", "mouse 1", 45, 55, 10, True)]
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
        self.assertEqual(out, [0, 1])

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


class TestIntervalGetDurations2_lists(unittest.TestCase):
    def test_1(self):
        a = [1, 5, 8]
        b = [2, 7, 10]
        self.assertEqual(uf.get_interval_durations_2_lists(a, b), [1, 2, 2])


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


class TestAllPairs(unittest.TestCase):
    def test_list(self):
        lista = ["Zdzisio",
                 "Zbysio",
                 "Henio",
                 "Gienio"]
        out = uf.all_pairs(lista)
        self.assertEqual(len(out), 12)


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



class TestGetEHSData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        cls.data = Loader(path)
        cls.t1 = calendar.timegm(uf.to_struck("10.10.201011:04:30"))
        cls.t2 = calendar.timegm(uf.to_struck("10.10.201011:06:20"))
        cls.m_1_a, cls.s1, cls.e1 = uf.get_ehs_data_with_margin(cls.data,
                                                                "mouse_1",
                                                                cls.t1, cls.t2,
                                                                margin=100)

    def test_get_ehs_address(self):
        self.assertEqual(["C", "C", "D", "C", "D", "C", "B", "C", "D",
                          "C", "D", "C", "D"], self.m_1_a)
       
    def test_get_ehs_starttimes(self):
        out = sorted(self.s1)
        starttimes = [1286708667.302, 1286708669.9, 1286708676.243,
                      1286708680.721, 1286708749.62, 1286708768.349,
                      1286708783.809, 1286708800.057, 1286708804.718,
                      1286708817.259, 1286708825.762, 1286708838.969,
                      1286708858.91]
        
        self.assertEqual(out, starttimes)

    def test_get_ehs_endtimes(self):
        out = sorted(self.e1)
        endtimes = [1286708669.65, 1286708674.28, 1286708680.125,
                    1286708748.484, 1286708767.615, 1286708781.731,
                    1286708799.069, 1286708803.98, 1286708816.712,
                    1286708825.01, 1286708838.22, 1286708858.004,
                    1286708960.687]
        self.assertEqual(out, endtimes)


class TestPrepareData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        cls.data = Loader(path)
        cls.t1 = calendar.timegm(uf.to_struck("10.10.201011:04:30"))
        cls.t2 = calendar.timegm(uf.to_struck("10.10.201011:06:20"))
        data1 = uf.prepare_data(cls.data, "mouse_1",
                               [cls.t1, cls.t2])
        cls.out = data1["mouse_1"]
        data2 = uf.prepare_data(cls.data, "mouse_1")
        cls.out_all = data2["mouse_1"]
        cls.t11 = calendar.timegm(uf.to_struck("10.10.201011:09:21"))
        cls.t12 = calendar.timegm(uf.to_struck("10.10.201011:55:55"))

        data3 = uf.prepare_data(cls.data, "mouse_1",
                                [cls.t11, cls.t12])
        cls.out_different = data3["mouse_1"]

        cls.no_mouse = uf.prepare_data(cls.data, "mouse_2",
                               [cls.t1, cls.t2])["mouse_2"]

        path2 = os.path.join(data_path, "weird_3_mice")
        data_longer = Loader(path2)
        
        cls.data_longer = uf.prepare_data(data_longer, ["mouse_1",
                                                        "mouse_2"],
                                          [cls.t11, cls.t12])
        

    def test_get_ehs_address(self):
        self.assertEqual(["C", "D", "C", "D", "C"], [x[0] for x in self.out])

    def test_get_ehs_starttimes(self):
        starttimes = [1286708670, 1286708676.243, 1286708680.721,
                      1286708749.62, 1286708768.349]
        self.assertTrue([x[1] for x in self.out] == starttimes)

    def test_get_ehs_endtimes(self):
        endtimes = [1286708674.28, 1286708680.125, 1286708748.484,
                    1286708767.615, 1286708780]
        self.assertTrue([x[2] for x in self.out] == endtimes)

    def test_get_ehs_starttimes_3(self):
        all_starttimes = set([self.out[i][1] for i in range(len(self.out))])
        self.assertTrue(len(all_starttimes) > 1)

    def test_get_ehs_endtimes_3(self):
        all_endtimes = set([self.out[i][2] for i in range(len(self.out))])
        self.assertTrue(len(all_endtimes) > 1)

    def test_get_all(self):
        self.assertEqual(len(self.out_all),
                         len(self.data.get_visit_addresses("mouse_1")))

    def test_get2_starttimes(self):
        self.assertTrue(self.out_different[0][1] > self.t11)

    def test_get2_endtimes(self):
        self.assertTrue(self.out_different[-1][-1] < self.t12)
    
    def test_no_mouse(self):
        self.assertEqual(self.no_mouse, [])

    def test_three_mice(self):
        mice_list = ["mouse_1", "mouse_2"]
        self.assertEqual(sorted(list(self.data_longer.keys())),
                         sorted(mice_list))

    def test_mouse_1_all_datasets(self):
        self.assertTrue(sorted(list(self.data_longer["mouse_1"])),
                        sorted(self.out_different))
                        
    def test_mouse_2_all_datasets_start(self):
        data_mouse2 = [x[1]< self.t11 for x in self.data_longer["mouse_2"]]
        self.assertEqual(sum(data_mouse2), 0)

    def test_mouse_2_all_datasets_end(self):
        data_mouse2 = [x[2] > self.t12 for x in self.data_longer["mouse_2"]]
        self.assertEqual(sum(data_mouse2), 0)


class TestGetAnimalPositions(unittest.TestCase):
    def test_threshold(self):
        out = uf.get_animal_position([2, 3], [2, 2], "mouse_1", 2)
        self.assertEqual(out, [])

    def test_pipe_1(self):
        out = uf.get_animal_position([2, 3], [1, 2], "mouse_1", 1)
        self.assertEqual(out, [])

    def test_pipe_2(self):
        out = uf.get_animal_position([2, 3], [3, 4], "mouse_1", 1)
        self.assertEqual(out, [])

    def test_pipe_3(self):
        out = uf.get_animal_position([2, 3], [5, 6], "mouse_1", 1)
        self.assertEqual(out, [])
        
    def test_pipe_4(self):
        out = uf.get_animal_position([2, 3], [7, 8], "mouse_1", 1)
        self.assertEqual(out, [])

    def test_chamber_A1(self):
        out = uf.get_animal_position([2, 6], [1, 8], "mouse_1", 2)
        self.assertEqual([("A", "mouse_1", 2, 6, 4, True)], out)
        
    def test_chamber_A2(self):
        out = uf.get_animal_position([2, 6], [8, 1], "mouse_1", 2)
        self.assertEqual([("A", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_B1(self):
        out = uf.get_animal_position([2, 6], [2, 3], "mouse_1", 2)
        self.assertEqual([("B", "mouse_1", 2, 6, 4, True)], out)
        
    def test_chamber_B2(self):
        out = uf.get_animal_position([2, 6], [3, 2], "mouse_1", 2)
        self.assertEqual([("B", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_C1(self):
        out = uf.get_animal_position([2, 6], [4, 5], "mouse_1", 2)
        self.assertEqual([("C", "mouse_1", 2, 6, 4, True)], out)
        
    def test_chamber_C2(self):
        out = uf.get_animal_position([2, 6], [5, 4], "mouse_1", 2)
        self.assertEqual([("C", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_D1(self):
        out = uf.get_animal_position([2, 6], [6, 7], "mouse_1", 2)
        self.assertEqual([("D", "mouse_1", 2, 6, 4, True)], out)
        
    def test_chamber_D2(self):
        out = uf.get_animal_position([2, 6], [7, 6], "mouse_1", 2)
        self.assertEqual([("D", "mouse_1", 2, 6, 4, True)], out)


class TestDictToArray2D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dictio = {"A":{"C": 1, "D": 2}, "B":{"C":3, "D": 4}}
        cls.out = np.array([[1, 2], [3, 4]])
        cls.new_arr = uf.dict_to_array_2D(dictio, ["A", "B"], ["C", "D"])

    def test00(self):
        self.assertEqual(self.out[0, 0], self.new_arr[0, 0])

    def test01(self):
        self.assertEqual(self.out[0, 1], self.new_arr[0, 1])

    def test10(self):
        self.assertEqual(self.out[1, 0], self.new_arr[1, 0])

    def test11(self):
        self.assertEqual(self.out[1, 1], self.new_arr[1, 1])


class TestDictToArray3D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dictio1 = {"A":{"C": 1, "D": 2}, "B":{"C":3, "D": 4}}
        dictio2 = {"A":{"C": 5, "D": 6}, "B":{"C":7, "D": 8}}
        dictio = {"AA": dictio1, "BB":dictio2}
        cls.out = np.array([[[1, 2],
                             [3, 4]],
                            [[5, 6],
                             [7, 8]] ])
        cls.new_arr = uf.dict_to_array_3D(dictio, ["AA", "BB"], ["A", "B"], ["C", "D"])

    def test000(self):
        self.assertEqual(self.out[0, 0, 0], self.new_arr[0, 0, 0])

    def test001(self):
        self.assertEqual(self.out[0, 0, 1], self.new_arr[0, 0, 1])

    def test010(self):
        self.assertEqual(self.out[0, 1, 0], self.new_arr[0, 1, 0])

    def test100(self):
        self.assertEqual(self.out[1, 0, 0], self.new_arr[1, 0, 0])

    def test110(self):
        self.assertEqual(self.out[1, 1, 0], self.new_arr[1, 1, 0])

    def test101(self):
        self.assertEqual(self.out[1, 0, 1], self.new_arr[1, 0, 1])

    def test111(self):
        self.assertEqual(self.out[1, 1, 1], self.new_arr[1, 1, 1])

    def test011(self):
        self.assertEqual(self.out[0, 1, 1], self.new_arr[0, 1, 1])


class TestCalcExcess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dictio11 = {"A":{"C": 1, "D": 2}, "B":{"C":3, "D": 4}}
        dictio12 = {"A":{"C": 5, "D": 4}, "B":{"C":2, "D": 1}}
        dictio13 = {"A":{"C": 8, "D": 9}, "B":{"C":10, "D": 11}}
        dictio1 = {"AA": dictio11, "BB":dictio12}
        dictio2 = {"AA": dictio12, "BB":dictio13}
        cls.out = uf.calc_excess(dictio1, dictio2)
        cls.new_arr = {
            "AA":{
                "A": {"C": -4, "D": -2},
                "B": {"C": 1, "D": 3}
            },
            "BB": {
                "A": {"C": -3, "D": -5},
                "B": {"C": -8, "D": -10}
            }
            }

    def test000(self):
        self.assertEqual(self.out["AA"]["A"]["C"], self.new_arr["AA"]["A"]["C"])

    def test001(self):
        self.assertEqual(self.out["AA"]["A"]["D"], self.new_arr["AA"]["A"]["D"])

    def test010(self):
        self.assertEqual(self.out["AA"]["B"]["C"], self.new_arr["AA"]["B"]["C"])

    def test100(self):
        self.assertEqual(self.out["AA"]["B"]["D"], self.new_arr["AA"]["B"]["D"])

    def test110(self):
        self.assertEqual(self.out["BB"]["A"]["C"], self.new_arr["BB"]["A"]["C"])

    def test101(self):
        self.assertEqual(self.out["BB"]["A"]["D"], self.new_arr["BB"]["A"]["D"])

    def test111(self):
        self.assertEqual(self.out["BB"]["B"]["C"], self.new_arr["BB"]["B"]["C"])

    def test011(self):
        self.assertEqual(self.out["BB"]["B"]["D"], self.new_arr["BB"]["B"]["D"])

class TestPrepareBinnedData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        cls.data = Loader(path)
        cls.config = ExperimentConfigFile(path)

        cls.all_phases, cls.all_total_time,\
            cls.all_data, cls.all_keys = uf.prepare_binned_data(cls.data,
                                                                cls.config,
                                                               "ALL",
                                                                ["mouse_1"])
        cls.dark_phases, cls.dark_total_time,\
            cls.dark_data, cls.dark_keys = uf.prepare_binned_data(cls.data,
                                                                  cls.config,
                                                                  "DARK",
                                                                  ["mouse_1"])
        cls.light_phases, cls.light_total_time,\
            cls.light_data, cls.light_keys = uf.prepare_binned_data(cls.data,
                                                                    cls.config,
                                                                    "LIGHT",
                                                                    ["mouse_1"])

        cls.phases_100s_bins, cls.total_time_100s_bins,\
            cls.data_100s_bins, cls.keys_100s_bins = uf.prepare_binned_data(cls.data, cls.config,100, ["mouse_1"])
        cls.phases_900s_bins, cls.total_time_900s_bins,\
            cls.data_900s_bins, cls.keys_900s_bins = uf.prepare_binned_data(cls.data, cls.config, 900, ["mouse_1"])

        path = os.path.join(data_path, "weird_short_3_mice")
        cls.data2 = Loader(path)
        cls.config2 = ExperimentConfigFile(path)
        cls.phases_24h_bins, cls.total_time_24h_bins,\
        cls.data_24h_bins, cls.keys_24h_bins = uf.prepare_binned_data(cls.data2, cls.config2, 24*3600, ["mouse_1"])

    def test_all_phases(self):
        self.assertEqual(self.all_phases, ["ALL"])

    def test_all_time(self):
        time_dict = {"ALL": {0: 3600}}
        self.assertEqual(self.all_total_time, time_dict)

    def test_all_data(self):
        data = {"mouse_1": uf.prepare_data(self.data, ["mouse_1"])}
        all_data = {"ALL": {0: data}}

    def test_all_keys(self):
        keys = [["ALL"], [0]]
        self.assertEqual(keys, self.all_keys)
    
    def test_dark_phases(self):
        self.assertEqual(self.dark_phases, ["DARK"])

    def test_dark_time(self):
        time_dict = {"DARK": {0: 1800.0}}
        self.assertEqual(self.dark_total_time, time_dict)

    def test_dark_data(self):
        data = {"mouse_1": uf.prepare_data(self.data, ["mouse_1"])}
        dark_data = {"DARK": {0: data}}

    def test_dark_keys(self):
        keys = [["DARK"], [0]]
        self.assertEqual(keys, self.dark_keys)

    def test_light_phases(self):
        self.assertEqual(self.light_phases, ["LIGHT"])

    def test_light_time(self):
        time_dict = {"LIGHT": {0: 1800.0}}
        self.assertEqual(self.light_total_time, time_dict)

    def test_light_data(self):
        data = {"mouse_1": uf.prepare_data(self.data, ["mouse_1"])}
        light_data = {"LIGHT": {0: data}}

    def test_light_keys(self):
        keys = [["LIGHT"], [0]]
        self.assertEqual(keys, self.light_keys)

    def test_bins_keys(self):
        keys = [["1 dark", "1 light"], [i*100 for i in range(12*3600//100)]]
        self.assertEqual(keys, self.keys_100s_bins)

    def test_bins_data_1st_bin(self):
        self.assertEqual(self.data_100s_bins["1 dark"][0]["mouse_1"], [])

    def test_bins_data_2nd_bin(self):
        self.assertEqual(self.data_100s_bins["1 dark"][100.]["mouse_1"], [])

    def test_bins_data_3rd_bin_len(self):
        data = self.data_100s_bins["1 dark"][200.]["mouse_1"]
        self.assertEqual(len(data), 4)

    def test_bins_data_3rd_bin_last_value(self):
        data = self.data_100s_bins["1 dark"][200.]["mouse_1"][3]
        self.assertEqual(data[-1], 1286708700)

    def test_bins_data_8th_bin_last_value(self):
        data = self.data_100s_bins["1 dark"][600.]["mouse_1"][0]
        self.assertEqual(data[-1], 1286708900)

    def test_bins_data_8th_bin_last_value(self):
        data = self.data_100s_bins["1 dark"][600.]["mouse_1"][0]
        self.assertEqual(data[1], 1286709000)

    def test_bins_data_8th_address(self):
        data = self.data_100s_bins["1 dark"][600.]["mouse_1"][0]
        self.assertEqual(data[0], "B")

    def test_bins_24_h_phases(self):
        self.assertEqual(self.phases_24h_bins, ['1_x_24.00h', '2_x_24.00h'])

    def test_total_time_24_h_phases(self):
        self.assertEqual(OrderedDict([('1_x', OrderedDict([(0.0, 86400.0)])),
                                      ('2_x', OrderedDict([(0.0, 86400.0)]))]),
                         self.total_time_24h_bins)

    def test_keys_24_h(self):
        self.assertEqual(self.keys_24h_bins, [['1_x', '2_x'], [0.0]])

    def test_data_1_bin(self):
        self.assertEqual(self.data_24h_bins["1_x"][0.0]["mouse_1"][-1],
                         ('B', 1286711990.827, 1286711995.578))

    def test_data_2_bin(self):
        self.assertEqual(self.data_24h_bins["2_x"][0.0]["mouse_1"],
                         [('C', 1286798036.218, 1286798396.218)])


class TestExtractDirections(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas1 = [1, 2,   3, 4,  5,  6,  7,   8,  1, 2]
        times1 = [15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34 ]
        antennas2 = [8, 1,   2,    1,  2,  3,  4,   5, 6,   7,   8]
        times2 =   [10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35]
        cls.res1 = uf.extract_directions(times1, antennas1, 3)
        cls.res1_1 = uf.extract_directions(times1, antennas1, 1)
        cls.res2 = uf.extract_directions(times2, antennas2, 1)
        cls.res2_1 = uf.extract_directions(times2, antennas2, 7)

    def test_1(self):
        direction_dict = {key:[[], []] for key in uf.KEYS}
        direction_dict["12"][0].append(15)
        direction_dict["12"][1].append(16.5)
        direction_dict["34"][0].append(19)
        direction_dict["34"][1].append(20)
        direction_dict["56"][0].append(21)
        direction_dict["56"][1].append(22)
        direction_dict["78"][0].append(24)
        direction_dict["78"][1].append(25)
        direction_dict["12"][0].append(29)
        direction_dict["12"][1].append(34)

        self.assertEqual(direction_dict, self.res1)

    def test_2(self):
        direction_dict = {key:[[], []] for key in uf.KEYS}
        direction_dict["12"][0].append(19.5)
        direction_dict["12"][1].append(22)
        direction_dict["34"][0].append(25)
        direction_dict["34"][1].append(26)
        direction_dict["56"][0].append(27)
        direction_dict["56"][1].append(28)
        direction_dict["78"][0].append(31)
        direction_dict["78"][1].append(35)

        self.assertEqual(direction_dict, self.res2)

    def test_1_1(self):
        direction_dict = {key:[[], []] for key in uf.KEYS}
        direction_dict["12"][0].append(15)
        direction_dict["12"][1].append(16.5)
        direction_dict["34"][0].append(19)
        direction_dict["34"][1].append(20)
        direction_dict["56"][0].append(21)
        direction_dict["56"][1].append(22)
        direction_dict["78"][0].append(24)
        direction_dict["78"][1].append(25)

        self.assertEqual(direction_dict, self.res1_1)

    def test_2_1(self):
        direction_dict = {key:[[], []] for key in uf.KEYS}
        direction_dict["12"][0].append(19.5)
        direction_dict["12"][1].append(22)
        direction_dict["34"][0].append(25)
        direction_dict["34"][1].append(26)
        direction_dict["56"][0].append(27)
        direction_dict["56"][1].append(28)

        self.assertEqual(direction_dict, self.res2_1)


class TestPrepareRegistrations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        data = Loader(path)
        config = ExperimentConfigFile(path)
        times = config.gettime("1 dark")
        t_start = times[-1] - 3600/3*2
        t_end = times[-1] - 3600/3
        cls.out = uf.prepare_registrations(data, ["mouse_1"],
                                            t_start, t_end)
    def test_1(self):
        self.assertEqual(["mouse_1"], list(self.out.keys()))

    def test_2(self):
        self.assertEqual(len(self.out["mouse_1"]["56"][0]), 3)

    def test_3(self):
        # check if last antenna is working
        self.assertEqual(len(self.out["mouse_1"]["65"][0]), 2)


class TestPrepareBinnedRegistrations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        cls.data = Loader(path)
        cls.config = ExperimentConfigFile(path)
        cls.out_all = uf.prepare_binned_registrations(cls.data,
                                                      cls.config,
                                                      "All",
                                                      ["mouse_1"])
        cls.out_6h = uf.prepare_binned_registrations(cls.data,
                                                     cls.config,
                                                     6*3600,
                                                     ["mouse_1"])
        cls.out_24h = uf.prepare_binned_registrations(cls.data,
                                                      cls.config,
                                                      24*3600,
                                                      ["mouse_1"])

    def test_all_phases(self):
        self.assertEqual(self.out_all[0], ["ALL"])

    def test_all_total_time(self):
        out = OrderedDict()
        out["ALL"] = {0: self.config.gettime("ALL")}
        self.assertEqual(self.out_all[1], out)

    def test_all_data(self):
        out = OrderedDict()
        time = self.config.gettime("ALL")
        out["ALL"] = {0: uf.prepare_registrations(self.data,
                                                  ["mouse_1"],
                                                  *time)}
        self.assertEqual(self.out_all[2], out)

    def test_all_keys(self):
        self.assertEqual(self.out_all[3], [["ALL"], [0.0]])

    def test_6h_phases(self):
        self.assertEqual(self.out_6h[0], ["1_dark_6.00h",
                                          "1_light_6.00h",
                                          "2_dark_6.00h"])

    def test_6h_data(self):
        phases = ["1 dark", "1 light", "2 dark"]
        out = OrderedDict()
        for phase in phases:
            times = self.config.gettime(phase)
            out[phase] = OrderedDict()
            out[phase][0.0] = uf.prepare_registrations(self.data,
                                                     ["mouse_1"],
                                                      times[0],
                                                      times[0]+
                                                      6*3600)
            key = 6*3600.
            out[phase][key] = uf.prepare_registrations(self.data,
                                                     ["mouse_1"],
                                                      times[0]+
                                                      6*3600,
                                                      times[-1])
        self.assertEqual(self.out_6h[2], out)

    def test_6h_total(self):
        out = OrderedDict()
        phases = ["1 dark", "1 light", "2 dark"]
        for phase in phases:
            out[phase] = OrderedDict()
            times = self.config.gettime(phase)
            out[phase][0.0] = (times[0], times[0] + 6*3600)
            out[phase][6*3600.0] = (times[0] + 6*3600, times[-1])
        self.assertEqual(out, self.out_6h[1])

    def test_6h_keys(self):
        data_keys = [["1 dark", "1 light", "2 dark"],
                     [0.0, 6*3600.0]]

    def test_24h_phases(self):
        self.assertEqual(self.out_24h[0], ['1_x_24.00h',
                                           '2_x_24.00h'])

    def test_24h_data(self):
        phases = ["1_x", "2_x"]
        times = self.config.gettime("ALL")
        out = OrderedDict()
        for phase in phases:
            out[phase] = OrderedDict()

        out["1_x"][0.0] = uf.prepare_registrations(self.data,
                                                   ["mouse_1"],
                                                   times[0],
                                                   times[0]+
                                                   24*3600)

        out["2_x"][0.0] = uf.prepare_registrations(self.data,
                                                   ["mouse_1"],
                                                   times[0]+
                                                   24*3600,
                                                   times[0]+
                                                   2*24*3600)
        self.assertEqual(self.out_24h[2], out)

    def test_24h_total(self):
        out = OrderedDict()
        times = self.config.gettime("ALL")
        phases = ["1_x", "2_x"]
        for phase in phases:
            out[phase] = OrderedDict()
        out["1_x"][0.0] = (times[0], times[0] + 24*3600)
        out["2_x"][0.0] = (times[0] + 24*3600,
                           times[0] + 2*24*3600)
        self.assertEqual(out, self.out_24h[1])

    def test_24h_keys(self):
        data_keys = [["1_x", "2_x"],
                     [0.0]]


if __name__ == '__main__':
    unittest.main()
