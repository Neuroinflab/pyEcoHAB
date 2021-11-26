# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
# encoding: utf-8
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
from pyEcoHAB import Timeline
from pyEcoHAB import SetupConfig


SAME_PIPE = {"1": ["1", "2"],
             "2": ["1", "2"],
             "3": ["3", "4"],
             "4": ["3", "4"],
             "5": ["5", "6"],
             "6": ["5", "6"],
             "7": ["7", "8"],
             "8": ["7", "8"]}

SAME_ADDRESS = {
    "1": ["1", "8"],
    "2": ["2", "3"],
    "3": ["2", "3"],
    "4": ["4", "5"],
    "5": ["4", "5"],
    "6": ["6", "7"],
    "7": ["6", "7"],
    "8": ["1", "8"],
}

OPPOSITE_PIPE = {"1": ["5", "6"],
                 "2": ["5", "6"],
                 "3": ["7", "8"],
                 "4": ["7", "8"],
                 "5": ["1", "2"],
                 "6": ["1", "2"],
                 "7": ["3", "4"],
                 "8": ["3", "4"]}

ADDRESS = {
    "1": "cage A",  # "4"
    "2": "cage B",  # 1,
    "3": "cage B",  # 1,
    "4": "cage C",  # 2,
    "5": "cage C",  # 2,
    "6": "cage D",  # "3",
    "7": "cage D",  # "3",
    "8": "cage A",  # "4"
}

ADDRESS_NON_ADJACENT = {
    "1": "cage B",  # 1,
    "2": "cage A",  # "4",
    "3": "cage C",  # 2,
    "4": "cage B",  # 1,
    "5": "cage D",  # "3",
    "6": "cage C",  # 2,
    "7": "cage A",  # "4",
    "8": "cage D",  # "3"
}
# Surrounding: difference between antennas only 2 or "6" -- skipped one antenna
SURROUNDING = {
    ("1", "3"): "cage B",  # 1,
    ("1", "7"): "cage A",  # "4",
    ("2", "4"): "cage B",  # 1,
    ("2", "8"): "cage A",  # "4",
    ("3", "5"): "cage C",  # 2,
    ("4", "6"): "cage C",  # 2,
    ("5", "7"): "cage D",  # "3",
    ("6", "8"): "cage D",  # "3"
}


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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
        self.assertEqual(len(states), len(readouts))


class TestChangeState(unittest.TestCase):
    def test_check_no_change(self):
        self.assertEqual(len(uf.change_state([1, 1, 1, 1])), 0)

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


class TestGetTimestamp(unittest.TestCase):
    def test_up(self):
        self.assertEqual(uf.get_timestamp(42.2, 45.76, 0.1), 36)

    def test_down(self):
        self.assertEqual(uf.get_timestamp(42.2, 45.71, 0.1), 35)


class TestGetAnennas(unittest.TestCase):
    def test(self):
        out = uf.get_antennas([0, 1, 4], [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertEqual(out, [1, 2, 5])


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
        out = uf.all_mouse_pairs(lista)
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
        out = uf.all_mouse_pairs(lista, reverse="True")
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
        self.assertEqual(self.out_data.shape, (len(out_lista),
                                               len(self.phases)))

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
        cls.out_data, cls.out = uf.make_table_of_all_mouse_pairs(cls.data,
                                                                 cls.phases,
                                                                 cls.lista,
                                                                 reverse=False)
        cls.o_d_rev, cls.o_rev = uf.make_table_of_all_mouse_pairs(cls.data,
                                                                  cls.phases,
                                                                  cls.lista,
                                                                  reverse=True)

    def test_shape(self):
        out_lista = uf.all_mouse_pairs(self.lista)
        self.assertEqual(self.out_data.shape,
                         (len(out_lista), len(self.phases)))

    def test_lista(self):
        out_lista = uf.all_mouse_pairs(self.lista)
        self.assertEqual(out_lista, self.out)

    def test_lista_rev(self):
        out_lista = uf.all_mouse_pairs(self.lista,
                                       reverse=True)
        self.assertEqual(out_lista, self.o_rev)

    def test_1st_phase(self):
        out = [1., 2., 0., 4., 0., 0.]
        first_column = self.out_data[:, 0].tolist()
        self.assertEqual(out, first_column)

    def test_1st_phase_rev(self):
        out = [0., 0., 1., 0., 2., 4.]
        first_column = self.o_d_rev[:, 0].tolist()
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
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
        states, readouts, midx = uf.get_more_states(antennas, times, idx,
                                                    50, 3)
        self.assertEqual(len(states), len(readouts))


class TestChangeState(unittest.TestCase):
    def test_check_no_change(self):
        self.assertEqual(len(uf.change_state([1, 1, 1, 1])), 0)

    def test_check_change_1(self):
        self.assertFalse(len(uf.change_state([1, 1, 1, 2])) == 0)

    def test_check_change_2(self):
        self.assertTrue(uf.change_state([1, 1, 1, 2]), np.array([2],
                                                                dtype=int))


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


class TestGetTimestamp(unittest.TestCase):

    def test_up(self):
        self.assertEqual(uf.get_timestamp(42.2, 45.76, 0.1), 36)

    def test_down(self):
        self.assertEqual(uf.get_timestamp(42.2, 45.71, 0.1), 35)


class TestGetAnennas(unittest.TestCase):
    def test(self):
        out = uf.get_antennas([0, 1, 4], [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertEqual(out, [1, 2, 5])


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
        out = uf.all_mouse_pairs(lista)
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
        self.assertEqual(self.out_data.shape, (len(out_lista),
                                               len(self.phases)))

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
        cls.t1 = 1286701470+7200
        cls.t2 = 1286701580+7200
        cls.m_1_a, cls.s1, cls.e1 = uf.get_ecohab_data_with_margin(cls.data,
                                                                   "mouse_1",
                                                                   cls.t1,
                                                                   cls.t2,
                                                                   margin=100)

    def test_get_ecohab_data_address(self):
        self.assertEqual(["cage C", "cage C", "cage D", "cage C", "cage D",
                          "cage C", "cage B", "cage C", "cage D",
                          "cage C", "cage D", "cage C", "cage D"],
                         self.m_1_a)

    def test_get_ecohab_data_starttimes(self):
        out = sorted(self.s1)
        starttimes = [1286708667.302, 1286708669.9, 1286708676.243,
                      1286708680.721, 1286708749.62, 1286708768.349,
                      1286708783.809, 1286708800.057, 1286708804.718,
                      1286708817.259, 1286708825.762, 1286708838.969,
                      1286708858.91]
        self.assertEqual(out, starttimes)

    def test_get_ecohab_data_endtimes(self):
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

    def test_get_ecohab_data_address(self):
        self.assertEqual(["cage C", "cage D", "cage C", "cage D", "cage C"],
                         [x[0] for x in self.out])

    def test_get_ecohab_data_starttimes(self):
        starttimes = [1286708670, 1286708676.243, 1286708680.721,
                      1286708749.62, 1286708768.349]
        self.assertTrue([x[1] for x in self.out] == starttimes)

    def test_get_ecohab_data_endtimes(self):
        endtimes = [1286708674.28, 1286708680.125, 1286708748.484,
                    1286708767.615, 1286708780]
        self.assertTrue([x[2] for x in self.out] == endtimes)

    def test_get_ecohab_data_starttimes_3(self):
        all_starttimes = set([self.out[i][1] for i in range(len(self.out))])
        self.assertTrue(len(all_starttimes) > 1)

    def test_get_ecohab_data_endtimes_3(self):
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
        data_mouse2 = [x[1] < self.t11 for x in self.data_longer["mouse_2"]]
        self.assertEqual(sum(data_mouse2), 0)

    def test_mouse_2_all_datasets_end(self):
        data_mouse2 = [x[2] > self.t12 for x in self.data_longer["mouse_2"]]
        self.assertEqual(sum(data_mouse2), 0)


class TestGetAnimalPositions(unittest.TestCase):
    def test_threshold(self):
        out = uf.get_animal_position([2, 3], ["2", "2"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual(out, [])

    def test_pipe_1(self):
        out = uf.get_animal_position([2, 3], ["1", "2"], "mouse_1", 1,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual(out, [])

    def test_pipe_2(self):
        out = uf.get_animal_position([2, 3], ["3", "4"], "mouse_1", 1,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual(out, [])

    def test_pipe_3(self):
        out = uf.get_animal_position([2, 3], ["5", "6"], "mouse_1", 1,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual(out, [])

    def test_pipe_4(self):
        out = uf.get_animal_position([2, 3], ["7", "8"], "mouse_1", 1,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual(out, [])

    def test_chamber_A1(self):
        out = uf.get_animal_position([2, 6], ["1", "8"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage A", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_A2(self):
        out = uf.get_animal_position([2, 6], ["8", "1"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage A", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_B1(self):
        out = uf.get_animal_position([2, 6], ["2", "3"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage B", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_B2(self):
        out = uf.get_animal_position([2, 6], ["3", "2"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage B", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_C1(self):
        out = uf.get_animal_position([2, 6], ["4", "5"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage C", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_C2(self):
        out = uf.get_animal_position([2, 6], ["5", "4"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage C", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_D1(self):
        out = uf.get_animal_position([2, 6], ["6", "7"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage D", "mouse_1", 2, 6, 4, True)], out)

    def test_chamber_D2(self):
        out = uf.get_animal_position([2, 6], ["7", "6"], "mouse_1", 2,
                                     same_pipe=SAME_PIPE,
                                     same_address=SAME_ADDRESS,
                                     opposite_pipe=OPPOSITE_PIPE,
                                     address=ADDRESS, surrounding=SURROUNDING,
                                     address_not_adjacent=ADDRESS_NON_ADJACENT,
                                     internal_antennas=[])
        self.assertEqual([("cage D", "mouse_1", 2, 6, 4, True)], out)

    def test_chambers_with_internal_antennas(self):
        positions = ["2", "8", "2", "1", "1", "1", "2", "8",
                     "8", "8", "8", "8", "2", "1"]
        times = [49.819, 54.198, 54.902, 65.126, 77.904,
                 92.234, 93.380, 114.173, 114.623,
                 115.214, 116.625, 121.984, 122.934, 123.571]
        positions2 = ["2", "2", "1", "1", "1", "2",  "2", "1"]
        times2 = [49.819,  54.902, 65.126, 77.904, 92.234, 93.380,
                  122.934, 123.571]
        SAME_PIPE = {
            "1": ["1", "2"],
            "2": ["1", "2"]
        }

        SAME_ADDRESS = {
            "1": ["1"],
            "2": ["2", "8"],
            "8": ["2", "8"]
        }
        OPPOSITE_PIPE = {}

        ADDRESS = {
            "1": "cage A",  # "4"
            "2": "cage B",  # 1,
            "8": "cage B",  # 1,
        }

        ADDRESS_NON_ADJACENT = {
            "1": "cage B",  # 1,
            "2": "cage A",  # "4",
        }
        SURROUNDING = {}
        o1 = uf.get_animal_position(times, positions, "mouse_1", 2,
                                    same_pipe=SAME_PIPE,
                                    same_address=SAME_ADDRESS,
                                    opposite_pipe=OPPOSITE_PIPE,
                                    address=ADDRESS, surrounding=SURROUNDING,
                                    address_not_adjacent=ADDRESS_NON_ADJACENT,
                                    internal_antennas=["8"])
        o2 = uf.get_animal_position(times2, positions2, "mouse_1", 2,
                                    same_pipe=SAME_PIPE,
                                    same_address=SAME_ADDRESS,
                                    opposite_pipe=OPPOSITE_PIPE,
                                    address=ADDRESS, surrounding=SURROUNDING,
                                    address_not_adjacent=ADDRESS_NON_ADJACENT,
                                    internal_antennas=[])
        self.assertEqual(o1, o2)


class TestDictToArray2D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dictio = {"A": {"C": 1, "D": 2}, "B": {"C": 3, "D": 4}}
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
        dictio1 = {"A": {"C": 1, "D": 2}, "B": {"C": 3, "D": 4}}
        dictio2 = {"A": {"C": 5, "D": 6}, "B": {"C": 7, "D": 8}}
        dictio = {"AA": dictio1, "BB": dictio2}
        cls.out = np.array([[[1, 2],
                             [3, 4]],
                            [[5, 6],
                             [7, 8]]])
        cls.new_arr = uf.dict_to_array_3D(dictio, ["AA", "BB"],
                                          ["A", "B"], ["C", "D"])

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
        dictio11 = {"A": {"C": 1, "D": 2}, "B": {"C": 3, "D": 4}}
        dictio12 = {"A": {"C": 5, "D": 4}, "B": {"C": 2, "D": 1}}
        dictio13 = {"A": {"C": 8, "D": 9}, "B": {"C": 10, "D": 11}}
        dictio1 = {"AA": dictio11, "BB": dictio12}
        dictio2 = {"AA": dictio12, "BB": dictio13}
        cls.out = uf.calc_excess(dictio1, dictio2)
        cls.new_arr = {
            "AA": {
                "A": {"C": -4, "D": -2},
                "B": {"C": 1, "D": 3}
            },
            "BB": {
                "A": {"C": -3, "D": -5},
                "B": {"C": -8, "D": -10}
            }
        }

    def test000(self):
        self.assertEqual(self.out["AA"]["A"]["C"],
                         self.new_arr["AA"]["A"]["C"])

    def test001(self):
        self.assertEqual(self.out["AA"]["A"]["D"],
                         self.new_arr["AA"]["A"]["D"])

    def test010(self):
        self.assertEqual(self.out["AA"]["B"]["C"],
                         self.new_arr["AA"]["B"]["C"])

    def test100(self):
        self.assertEqual(self.out["AA"]["B"]["D"],
                         self.new_arr["AA"]["B"]["D"])

    def test110(self):
        self.assertEqual(self.out["BB"]["A"]["C"],
                         self.new_arr["BB"]["A"]["C"])

    def test101(self):
        self.assertEqual(self.out["BB"]["A"]["D"],
                         self.new_arr["BB"]["A"]["D"])

    def test111(self):
        self.assertEqual(self.out["BB"]["B"]["C"],
                         self.new_arr["BB"]["B"]["C"])

    def test011(self):
        self.assertEqual(self.out["BB"]["B"]["D"],
                         self.new_arr["BB"]["B"]["D"])


class TestPrepareBinnedData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        cls.data = Loader(path)
        cls.config = Timeline(path)

        cls.all_phases, cls.all_total_time,\
            cls.all_data, cls.all_keys = uf.prepare_binned_data(cls.data,
                                                                cls.config,
                                                                "ALL",
                                                                ["mouse_1"])
        cls.whole_phases, cls.whole_total_time,\
            cls.whole_data, cls.whole_keys = uf.prepare_binned_data(cls.data,
                                                                    cls.config,
                                                                    "whole phase",
                                                                    ["mouse_1"])

        cls.dark_phases, cls.dark_total_time,\
            cls.dark_data, cls.dark_keys = uf.prepare_binned_data(cls.data,
                                                                  cls.config,
                                                                  "DARK",
                                                                  ["mouse_1"])
        cls.light_phases, cls.light_total_time,\
            cls.light_data,\
            cls.light_keys = uf.prepare_binned_data(cls.data,
                                                    cls.config,
                                                    "LIGHT",
                                                    ["mouse_1"])

        cls.phases_100s_bins, cls.total_time_100s_bins,\
            cls.data_100s_bins,\
            cls.keys_100s_bins = uf.prepare_binned_data(cls.data, cls.config,
                                                        100, ["mouse_1"])
        cls.phases_900s_bins, cls.total_time_900s_bins,\
            cls.data_900s_bins,\
            cls.keys_900s_bins = uf.prepare_binned_data(cls.data, cls.config,
                                                        900, ["mouse_1"])

        path = os.path.join(data_path, "weird_short_3_mice")
        cls.data2 = Loader(path)
        cls.config2 = Timeline(path)
        cls.phases_24h_bins, cls.total_time_24h_bins,\
            cls.data24h_b, cls.keys24h_b = uf.prepare_binned_data(cls.data2,
                                                                  cls.config2,
                                                                  24*3600,
                                                                  ["mouse_1"])

    def test_all_phases(self):
        self.assertEqual(self.all_phases, ["ALL"])

    def test_all_time(self):
        time_dict = {"ALL": {0: 3600}}
        self.assertEqual(self.all_total_time, time_dict)

    def test_all_data(self):
        data =  uf.prepare_data(self.data, ["mouse_1"])
        all_data = {"ALL": {0: data}}
        self.assertEqual(self.all_data, all_data)

    def test_all_keys(self):
        keys = [["ALL"], {"ALL": [0]}]
        self.assertEqual(keys, self.all_keys)

    def test_whole_phases(self):
        self.assertEqual(self.whole_phases, ["1_dark", "1_light"])

    def test_whole_time(self):
        time_dict = {"1 dark": {0: 1800}, "1 light": {0: 1800}}
        self.assertEqual(self.whole_total_time, time_dict)

    def test_whole_keys(self):
        keys = [["1 dark", "1 light"],
                {"1 dark": [0], "1 light": [0]}]
        self.assertEqual(keys, self.whole_keys)
    
    def test_dark_phases(self):
        self.assertEqual(self.dark_phases, ["DARK"])

    def test_dark_time(self):
        time_dict = {"DARK": {0: 1800.0}}
        self.assertEqual(self.dark_total_time, time_dict)

    def test_dark_keys(self):
        keys = [["DARK"], {"DARK": [0]}]
        self.assertEqual(keys, self.dark_keys)

    def test_light_phases(self):
        self.assertEqual(self.light_phases, ["LIGHT"])

    def test_light_time(self):
        time_dict = {"LIGHT": {0: 1800.0}}
        self.assertEqual(self.light_total_time, time_dict)

    def test_light_keys(self):
        keys = [["LIGHT"], {"LIGHT": [0]}]
        self.assertEqual(keys, self.light_keys)

    def test_bins_keys(self):
        keys = [["1 dark", "1 light"],
                {"1 dark": [i*100. for i in range(1800//100)],
                 "1 light": [i*100. for i in range(1800//100)]}]
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
        self.assertEqual(data[0], "cage B")

    def test_bins_24_h_phases(self):
        self.assertEqual(self.phases_24h_bins, ['1_x_24.00h', '2_x_24.00h'])

    def test_total_time_24_h_phases(self):
        self.assertEqual(OrderedDict([('1_x', OrderedDict([(0.0, 86400.0)])),
                                      ('2_x', OrderedDict([(0.0, 86400.0)]))]),
                         self.total_time_24h_bins)

    def test_keys_24_h(self):
        self.assertEqual(self.keys24h_b,
                         [['1_x', '2_x'], {'1_x': [0.0], '2_x': [0.0]}])

    def test_data_1_bin(self):
        self.assertEqual(self.data24h_b["1_x"][0.0]["mouse_1"][-1],
                         ('cage B', 1286711990.827, 1286711995.578))

    def test_data_2_bin(self):
        self.assertEqual(self.data24h_b["2_x"][0.0]["mouse_1"],
                         [('cage C', 1286798036.218, 1286798396.218)])


class TestExtractDirections(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas1 = [1, 2, 3, 4, 5, 6, 7, 8, 1, 2]
        times1 = [15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34]
        antennas2 = [8, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8]
        times2 = [10, 16, 19, 19.5, 22, 25, 26, 27, 28, 31, 35]
        config = SetupConfig()
        cls.keys = config.directions
        cls.res1 = uf.extract_directions(times1, antennas1, 3,
                                         config.directions)
        cls.res1_1 = uf.extract_directions(times1, antennas1, 1,
                                           config.directions)
        cls.res2 = uf.extract_directions(times2, antennas2, 1,
                                         config.directions)
        cls.res2_1 = uf.extract_directions(times2, antennas2, 7,
                                           config.directions)

    def test_1(self):
        direction_dict = {key: [[], []] for key in self.keys}
        direction_dict["1 2"][0].append(15)
        direction_dict["1 2"][1].append(16.5)
        direction_dict["3 4"][0].append(19)
        direction_dict["3 4"][1].append(20)
        direction_dict["5 6"][0].append(21)
        direction_dict["5 6"][1].append(22)
        direction_dict["7 8"][0].append(24)
        direction_dict["7 8"][1].append(25)
        direction_dict["1 2"][0].append(29)
        direction_dict["1 2"][1].append(34)
        self.assertEqual(direction_dict, self.res1)

    def test_2(self):
        direction_dict = {key: [[], []] for key in self.keys}
        direction_dict["1 2"][0].append(19.5)
        direction_dict["1 2"][1].append(22)
        direction_dict["3 4"][0].append(25)
        direction_dict["3 4"][1].append(26)
        direction_dict["5 6"][0].append(27)
        direction_dict["5 6"][1].append(28)
        direction_dict["7 8"][0].append(31)
        direction_dict["7 8"][1].append(35)
        self.assertEqual(direction_dict, self.res2)

    def test_1_1(self):
        direction_dict = {key: [[], []] for key in self.keys}
        direction_dict["1 2"][0].append(15)
        direction_dict["1 2"][1].append(16.5)
        direction_dict["3 4"][0].append(19)
        direction_dict["3 4"][1].append(20)
        direction_dict["5 6"][0].append(21)
        direction_dict["5 6"][1].append(22)
        direction_dict["7 8"][0].append(24)
        direction_dict["7 8"][1].append(25)
        self.assertEqual(direction_dict, self.res1_1)

    def test_2_1(self):
        direction_dict = {key: [[], []] for key in self.keys}
        direction_dict["1 2"][0].append(19.5)
        direction_dict["1 2"][1].append(22)
        direction_dict["3 4"][0].append(25)
        direction_dict["3 4"][1].append(26)
        direction_dict["5 6"][0].append(27)
        direction_dict["5 6"][1].append(28)
        self.assertEqual(direction_dict, self.res2_1)


class TestPrepareRegistrations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        data = Loader(path)
        config = Timeline(path)
        times = config.get_time_from_epoch("1 dark")
        t_start = times[-1] - 3600/3*2
        t_end = times[-1] - 3600/3
        cls.out = uf.prepare_registrations(data, ["mouse_1"],
                                           t_start, t_end)

    def test_1(self):
        self.assertEqual(["mouse_1"], list(self.out.keys()))

    def test_2(self):
        self.assertEqual(len(self.out["mouse_1"]["5 6"][0]), 3)

    def test_3(self):
        #  check if last antenna is working
        self.assertEqual(len(self.out["mouse_1"]["6 5"][0]), 2)


class TestPrepareBinnedRegistrations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        cls.data = Loader(path)
        cls.config = Timeline(path)
        cls.out_all = uf.get_registrations_bins(cls.data,
                                                cls.config,
                                                "All",
                                                ["mouse_1"])
        cls.out_6h = uf.get_registrations_bins(cls.data,
                                               cls.config,
                                               6*3600,
                                               ["mouse_1"])
        cls.out_24h = uf.get_registrations_bins(cls.data,
                                                cls.config,
                                                24*3600,
                                                ["mouse_1"])

    def test_all_phases(self):
        self.assertEqual(self.out_all[0], ["ALL"])

    def test_all_total_time(self):
        out = OrderedDict()
        out["ALL"] = {0: self.config.get_time_from_epoch("ALL")}
        self.assertEqual(self.out_all[1], out)

    def test_all_data(self):
        out = OrderedDict()
        time = self.config.get_time_from_epoch("ALL")
        out["ALL"] = {0: uf.prepare_registrations(self.data,
                                                  ["mouse_1"],
                                                  *time)}
        self.assertEqual(self.out_all[2], out)

    def test_all_keys(self):
        self.assertEqual(self.out_all[3], [["ALL"], {"ALL": [0.0]}])

    def test_6h_phases(self):
        self.assertEqual(self.out_6h[0], ["1_dark_6.00h",
                                          "1_light_6.00h",
                                          "2_dark_6.00h"])

    def test_6h_data(self):
        phases = ["1 dark", "1 light", "2 dark"]
        out = OrderedDict()
        for phase in phases:
            times = self.config.get_time_from_epoch(phase)
            out[phase] = OrderedDict()
            out[phase][0.0] = uf.prepare_registrations(self.data,
                                                       ["mouse_1"],
                                                       times[0],
                                                       times[0] +
                                                       6*3600)
            key = 6*3600.
            out[phase][key] = uf.prepare_registrations(self.data,
                                                       ["mouse_1"],
                                                       times[0] +
                                                       6*3600,
                                                       times[-1])
        self.assertEqual(self.out_6h[2], out)

    def test_6h_total(self):
        out = OrderedDict()
        phases = ["1 dark", "1 light", "2 dark"]
        for phase in phases:
            out[phase] = OrderedDict()
            times = self.config.get_time_from_epoch(phase)
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
        times = self.config.get_time_from_epoch("ALL")
        out = OrderedDict()
        for phase in phases:
            out[phase] = OrderedDict()

        out["1_x"][0.0] = uf.prepare_registrations(self.data,
                                                   ["mouse_1"],
                                                   times[0],
                                                   times[0] +
                                                   24*3600)

        out["2_x"][0.0] = uf.prepare_registrations(self.data,
                                                   ["mouse_1"],
                                                   times[0] +
                                                   24*3600,
                                                   times[0] +
                                                   2*24*3600)
        self.assertEqual(self.out_24h[2], out)

    def test_24h_total(self):
        out = OrderedDict()
        times = self.config.get_time_from_epoch("ALL")
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

class TestMath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.phase = ['EMPTY 1 dark']
        cls.binsizes = [1800, 3600, 1.5 * 3600, 7200, 14400, 43200]
        cls.bins = 12 * 3600 / cls.binsizes[4]
        cls.bin_labels = {"EMPTY 1 dark": []}
        for i in np.arange(0, cls.bins, 1):
            cls.bin_labels["EMPTY 1 dark"].append(cls.binsizes[4] * i);
        cls.mice = ['0065-0136651817', '0065-0136653169', '0065-0136655780']


        cls.excess = OrderedDict()
        cls.test_excess = OrderedDict()
        cls.test_sum = OrderedDict()
        cls.activity = OrderedDict()
        cls.test_activity = OrderedDict()
        cls.test = OrderedDict()
        cls.test_error = OrderedDict()


        for key1 in cls.phase:
            cls.excess[key1] = OrderedDict()
            cls.test_excess[key1] = OrderedDict()
            cls.test_sum[key1] = OrderedDict()
            cls.activity[key1] = OrderedDict()
            cls.test_activity[key1] = OrderedDict()
            cls.test[key1] = OrderedDict()
            cls.test_error[key1] = OrderedDict()
            for key2 in cls.bin_labels:
                cls.excess[key1][key2] = OrderedDict()
                cls.test_excess[key1][key2] = OrderedDict()
                cls.test_sum[key1][key2] = OrderedDict()
                cls.activity[key1][key2] = OrderedDict()
                cls.test_activity[key1][key2] = OrderedDict()
                cls.test[key1][key2] = OrderedDict()
                cls.test_error[key1][key2] = OrderedDict()
                for key3 in cls.mice:
                    cls.excess[key1][key2][key3] = OrderedDict()
                    cls.test_excess[key1][key2][key3] = OrderedDict()
                    cls.activity[key1][key2][key3] = OrderedDict()
                    cls.test_activity[key1][key2][key3] = 4/5
                    cls.test_sum[key1][key2][key3] = 4
                    cls.test[key1][key2][key3] = 2
                    cls.test_error[key1][key2][key3] = 0
                    for key4 in cls.mice:
                        cls.excess[key1][key2][key3][key4] = 0.00
                        cls.test_excess[key1][key2][key3][key4] = 2
                        if key3 == key4:
                            cls.test_excess[key1][key2][key3][key4] = 0.00
                        elif key3 < key4:
                            cls.excess[key1][key2][key3][key4] = 2
                    for key5 in range(5):
                        cls.activity[key1][key2][key3][key5] = 'activity'




    def test_diagonal_reflection_of_matrix(self):

        reflected_excess_time = uf.diagonal_reflection(self.excess[self.phase[0]],
                                                       self.mice,
                                                       self.bin_labels)
        self.assertEqual(reflected_excess_time,
                         self.test_excess[self.phase[0]],
                         "False, diagonal reflection test failed")
        return (reflected_excess_time)


    def test_sum_per_mouse(self):
        sum_time = OrderedDict()
        sum_time[self.phase[0]] = self.test_diagonal_reflection_of_matrix()
        sum_time[self.phase[0]] = uf.sum_per_mouse(sum_time[self.phase[0]],
                                                            self.mice,
                                                            self.bin_labels[self.phase[0]],
                                                            "leader", True)
        self.assertEqual(sum_time, self.test_sum, "False, sum test with phase failed")

    def test_divide_sum_activity(self):
        div_result = uf.divide_sum_activity(self.test_sum[self.phase[0]],
                                            self.activity[self.phase[0]],
                                            self.mice, self.bin_labels)
        self.assertEqual(div_result, self.test_activity[self.phase[0]],
                         "False, division test failed")

    def test_mean(self):
        mean_result = uf.mean(self.test_sum_per_mouse(),
                              len(self.mice)-1,
                              self.mice,
                              self.bin_labels)
        self.assertEqual(mean_result, self.test[self.phase[0]],
                         "False, mean test failed")
        return(mean_result)

    def test_standard_error(self):
        a = self.test_diagonal_reflection_of_matrix()
        error_result = uf.standard_error(a,
                                         self.test_mean(),
                                         self.mice, self.bin_labels)
        self.assertEqual(error_result, self.test_error[self.phase[0]],
                         "False, standard error test failed")
        return(error_result)

if __name__ == '__main__':
    unittest.main()
