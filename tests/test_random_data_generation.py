# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
import unittest
import numpy as np

from pyEcoHAB.utils import random_data_generation as rdg
from pyEcoHAB.utils import general as utils
from pyEcoHAB import data_path, Loader, Timeline


class TestPseudoLoader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        cls.dataset = Loader(path)
        cls.data = cls.dataset.registrations.data

    def test_generation(self):
        PL = rdg.PseudoLoader(self.data, self.dataset.setup_config)
        equal = set(PL.registrations.data == self.data)
        self.assertEqual(set([True]), equal)


class TestGetShifts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mouse_list = ["Zdzisio", "Zbysio", "Dyzio"]
        cls.shifts_dict = rdg.get_shifts(cls.mouse_list)

    def test_keys(self):
        self.assertEqual(sorted(self.mouse_list),
                         sorted(self.shifts_dict.keys()))

    def different_vals(self):
        self.assertEqual(len(set(self.shifts_dict.values())), 3)

    def test_range_1(self):
        equal = set(np.array(list(self.shifts_dict.values())) < 1800)
        self.assertTrue(equal, set([True]))

    def test_range_2(self):
        equal = set(np.array(list(self.shifts_dict.values())) > -1800)
        self.assertTrue(equal, set([True]))


class TestRandomlyShiftData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        cls.dataset = Loader(path)
        cls.data = cls.dataset.registrations.data
        np.random.seed(1)
        cls.shifts_dict = rdg.get_shifts(cls.dataset.mice)
        np.random.seed(1)
        cls.shifted = rdg.randomly_shift_data(cls.data)

    def test_shifting_mouse3(self):
        data_mouse3 = self.data[32]
        shifted_mouse3 = self.shifted[32]
        self.assertTrue(np.isclose(data_mouse3["Time"]
                                   + self.shifts_dict["mouse_3"],
                                   shifted_mouse3["Time"]))

    def test_shifting_mouse2(self):
        idcs = [10, 11]
        equal = []
        for idx in idcs:
            equal.append(np.isclose(self.data[idx]["Time"]
                                    + self.shifts_dict["mouse_2"],
                                    self.shifted[idx]["Time"]))
        self.assertEqual(set([True]), set(equal))

    def test_shifting_mouse1(self):
        idcs = [10, 11, 33]
        equal = []
        for idx, line in enumerate(self.data):
            if idx in idcs:
                continue
            equal.append(np.isclose(line["Time"]
                                    + self.shifts_dict["mouse_1"],
                                    self.shifted[idx]["Time"]))
        self.assertEqual(set([True]), set(equal))


class TestGenerateSurrogateData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        timeline = Timeline(path)
        cls.dataset = Loader(path)
        cls.data = cls.dataset.registrations.data
        cls.N = 2
        np.random.seed(1)
        func = utils.prepare_tube_data
        cls.surrogate = rdg.generate_surrogate_data(cls.dataset,
                                                    timeline,
                                                    "whole_phases",
                                                    cls.dataset.mice,
                                                    cls.N, func)

    def test_length(self):
        self.assertEqual(self.N, len(self.surrogate))

    def test_keys_1(self):
        self.assertEqual(sorted(self.surrogate[0].keys()),
                         sorted(self.surrogate[1].keys()))

    def test_keys_11(self):
        self.assertEqual(sorted(self.surrogate[0].keys()),
                         sorted(["1 dark", "1 light", "2 dark"]))

    def test_assert_different_values(self):
        self.assertFalse(self.surrogate[0]["1 dark"][0] ==
                         self.surrogate[1]["1 dark"][0])


class TestReshapeSurrogateData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        timeline = Timeline(path)
        cls.dataset = Loader(path)
        cls.data = cls.dataset.registrations.data
        cls.N = 2
        np.random.seed(1)
        func = utils.prepare_tube_data
        cls.surrogate = rdg.generate_surrogate_data(cls.dataset,
                                                    timeline,
                                                    "whole_phases",
                                                    cls.dataset.mice,
                                                    cls.N, func)
        cls.reshaped = rdg.reshape_surrogate_data(cls.surrogate)

    def test_reshaped_1(self):
        self.assertEqual(self.reshaped["1 dark"][0][0],
                         self.surrogate[0]["1 dark"][0])

    def test_reshaped_2(self):
        self.assertEqual(self.reshaped["1 dark"][0][1],
                         self.surrogate[1]["1 dark"][0])

    def test_reshaped_3(self):
        self.assertEqual(self.reshaped["1 light"][0][0],
                         self.surrogate[0]["1 light"][0])

    def test_reshaped_4(self):
        self.assertEqual(self.reshaped["1 light"][0][1],
                         self.surrogate[1]["1 light"][0])

    def test_reshaped_5(self):
        self.assertEqual(self.reshaped["2 dark"][0][0],
                         self.surrogate[0]["2 dark"][0])

    def test_reshaped_6(self):
        self.assertEqual(self.reshaped["2 dark"][0][1],
                         self.surrogate[1]["2 dark"][0])


if __name__ == '__main__':
    unittest.main()
