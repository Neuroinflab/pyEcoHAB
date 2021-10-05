# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import unittest
import numpy as np
from pyEcoHAB import dominance_in_2_cages as dom
from pyEcoHAB import SetupConfig, data_path


class TestGetStates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dt = 0.1
        cls.t_start = 600.0
        cls.t_end = 800.0
        cls.config1 = SetupConfig(path=data_path,
                                  fname="setup_short_2.txt")
        cls.config2 = SetupConfig(path=data_path,
                                  fname="setup_short_3.txt")
        cls.antennas = ["3", "4", "4", "3", "3", "4", "4", "3", "3", "4",
                        "4", "3", "3", "4", "4", "3", "3"]
        cls.times = [641.083,  # 3  # 0
                     642.135,  # 4  # 1
                     675.134,  # 4  # 2
                     675.869,  # 3  # 3
                     681.127,  # 3  # 4
                     681.734,  # 4  # 5
                     692.744,  # 4  # 6
                     693.207,  # 3  # 7
                     701.82,   # 3  # 8
                     702.603,  # 4  # 9
                     703.499,  # 4  # 10
                     703.961,  # 3  # 11
                     723.136,  # 3  # 12
                     725.633,  # 4  # 13
                     734.133,  # 4  # 14
                     734.945,  # 3  # 15
                     783.411]  # 3  # 16
        cls.out_1 = dom.get_states_mouse(cls.antennas,
                                         cls.times,
                                         cls.t_start,
                                         cls.t_end,
                                         cls.config1,
                                         cls.dt)
        cls.out_2 = dom.get_states_mouse(cls.antennas,
                                         cls.times,
                                         cls.t_start,
                                         cls.t_end,
                                         cls.config2,
                                         cls.dt)

    def test_same_length(self):
        self.assertEqual(len(self.out_1), len(self.out_2))

    def test_different_results_for_different_home_antenna(self):
        for i, x in enumerate(self.out_1):
            if x != 1:
                self.assertNotEqual(x, self.out_2[i])

    def test_different_results_for_home_antenna_1(self):
        self.assertFalse(np.all(self.out_1 == self.out_1[0]))

    def test_different_results_for_home_antenna_2(self):
        self.assertFalse(np.all(self.out_2 == self.out_2[0]))

    def test_same_results_pipe(self):
        for i, x in enumerate(self.out_1):
            if x == 1:
                self.assertEqual(x, self.out_2[i])

    def test_different_home_for_home_antenna_1(self):
        self.assertTrue(np.any(self.out_1 == 0))

    def test_different_home_for_home_antenna_2(self):
        self.assertTrue(np.any(self.out_2 == 0))

    def test_threshold_1(self):
        timestamp_1 = int(round((self.times[9] - self.t_start)/self.dt))
        timestamp_2 = int(round((self.times[9] - self.t_start)/self.dt))
        self.assertTrue(np.all(self.out_1[timestamp_1:timestamp_2] == 1))

    def test_threshold_1(self):
        timestamp_1 = int(round((self.times[9] - self.t_start)/self.dt))
        timestamp_2 = int(round((self.times[9] - self.t_start)/self.dt))
        self.assertTrue(np.all(self.out_2[timestamp_1:timestamp_2] == 1))


class TestFindStimulusCageMice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dt = 0.05
        length = int(10/cls.dt)
        out_1 = np.ones((length))
        out_2 = np.ones((length))
        out_3 = np.ones((length))
        out_1[45:100] = 3
        out_2[80:108] = 3
        cls.data = {
            'mouse 1': out_1,
            'mouse 2': out_2,
            'mouse 3': out_3
        }

    def test_1(self):
        t_start = 4.2
        t_stop = 6.05
        beginning = 2.1
        out = dom.find_stimulus_cage_mice(self.data,
                                          t_start, t_stop,
                                          beginning, self.dt)
        self.assertEqual(out, ['mouse 1'])

    def test_2(self):
        t_start = 7.1
        t_stop = 8
        beginning = 2.1
        out = dom.find_stimulus_cage_mice(self.data,
                                          t_start, t_stop,
                                          beginning, self.dt)
        self.assertEqual(out, ['mouse 2'])

    def test_3(self):
        t_start = 6.2
        t_stop = 8
        beginning = 2.1
        out = dom.find_stimulus_cage_mice(self.data,
                                          t_start, t_stop,
                                          beginning, self.dt)
        self.assertTrue('mouse 1' in out)

    def test_4(self):
        t_start = 6.2
        t_stop = 8
        beginning = 2.1
        out = dom.find_stimulus_cage_mice(self.data,
                                          t_start, t_stop,
                                          beginning, self.dt)
        self.assertTrue('mouse 2' in out)

    def test_5(self):
        t_start = 6.2
        t_stop = 8
        beginning = 2.1
        out = dom.find_stimulus_cage_mice(self.data,
                                          t_start, t_stop,
                                          beginning, self.dt)
        self.assertEqual(len(out), 2)


class TestCheckMouse1NotValid(unittest.TestCase):
    def test_home_antenna(self):
        out = dom.check_mouse1_not_valid("4", "4", "4")
        self.assertTrue(out)

    def test_pipe(self):
        out = dom.check_mouse1_not_valid("4", "3", "4")
        self.assertTrue(out)

    def test_True(self):
        out = dom.check_mouse1_not_valid("4", "4", "3")
        self.assertFalse(out)


class TestCheckMouse2NotValid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.antennas = ["3",  "4",   "4",   "3",   "3"]
        cls.times = [11., 12., 13., 14., 21.]

    def test_no_pre(self):
        self.assertTrue(dom.check_mouse2_not_valid(10, 20,
                                                   ["3", "4"],
                                                   [11., 12.],
                                                   "3"))

    def test_no_between(self):
        self.assertTrue(dom.check_mouse2_not_valid(15, 20,
                                                   self.antennas,
                                                   self.times,
                                                   "3"))

    def test_not_home_antenna(self):
        self.assertTrue(dom.check_mouse2_not_valid(11.5, 20,
                                                   self.antennas,
                                                   self.times,
                                                   "4"))

    def test_return_False(self):
        self.assertFalse(dom.check_mouse2_not_valid(11.5, 20,
                                                    self.antennas,
                                                    self.times,
                                                    "3"))


class TestCountAttempts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.antennas = ["3", "4", "4", "4", "3", "3", "3", "4", "3", "3",
                        "3", "3", "4"]
        cls.times = [5., 12., 13., 14., 21., 22., 24., 25.,
                     26., 28., 35., 41., 44.]
        cls.config = SetupConfig(path=data_path,
                                 fname="setup_short_2.txt")

    def test_check_1_true(self):
        out = dom.count_attempts(12.5, 13.5, self.times, self.antennas, "4",
                                 self.config)
        self.assertEqual(out, 1)

    def test_check_1_false(self):
        out = dom.count_attempts(11.5, 12.5, self.times, self.antennas, "3",
                                 self.config)
        self.assertEqual(out, 0)

    def test_check_2_false(self):
        out = dom.count_attempts(11.5, 14.5, self.times, self.antennas, "3",
                                 self.config)
        self.assertEqual(out, 0)

    def test_check_2_true(self):
        out = dom.count_attempts(11.5, 13.5, self.times, self.antennas, "4",
                                 self.config)
        self.assertEqual(out, 1)

    def test_check_more(self):
        out = dom.count_attempts(11.5, 41.5, self.times, self.antennas, "4",
                                 self.config)
        self.assertEqual(out, 1)

    def test_check_more_opposite_antenna(self):
        out = dom.count_attempts(12.5, 41.5, self.times, self.antennas, "3",
                                 self.config)
        self.assertEqual(out, 3)


class TestMouseDefending(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.timestamp_1 = 10
        cls.timestamp_2 = 20
        cls.timestamp_3 = 40
        cls.home_antenna = 4
        cls.antennas1 = ["3", "3", "3", "3"]
        cls.times1 = [4, 12.5, 21.5, 45]
        cls.antennas2 = ["3", "4", "4", "4", "3", "3", "3", "4", "3", "3",
                         "3", "3", "4"]
        cls.times2 = [5., 12., 13., 14., 21., 22., 24., 25., 26., 28.,
                      35., 41., 44.]
        cls.config = SetupConfig(path=data_path,
                                 fname="setup_short_2.txt")

    def test_mouse_1_in_sugar_cage(self):
        out = dom.check_mouse1_defending(self.antennas1, self.times1,
                                         self.antennas2, self.times2, "4",
                                         self.config)
        self.assertEqual(out, 1)

    def test_mouse_2_more_attempst(self):
        antennas1 = ["4", "4", "4", "4"]
        out = dom.check_mouse1_defending(antennas1, self.times1,
                                         self.antennas2, self.times2, "3",
                                         self.config)
        self.assertEqual(out, 3)


if __name__ == '__main__':
    unittest.main()
