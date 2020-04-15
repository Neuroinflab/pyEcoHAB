#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
import os
import unittest
import random
import numpy as np

import pyEcoHAB.utils.for_loading as uf
from pyEcoHAB import data_path

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
        self.assertEqual('Wed Apr  3 01:17:47 2019', uf.print_human_time(tt))


class TestTimeToSec(unittest.TestCase):
    def test_sec(self):
        string = "20190709 20:05:13.333"
        self.assertEqual(1562695513.333, uf.time_to_sec(string))

    def test_sec_2(self):
        string = "20190709 20:05:13"
        self.assertEqual(1562695513., uf.time_to_sec(string))

    def test_sec_raise_1(self):
        string = "2019070920:05:13"
        self.assertRaises(ValueError, uf.time_to_sec, tt=string)

    def test_sec_raise_2(self):
        string = "20190709.20:05:13"
        self.assertRaises(ValueError, uf.time_to_sec, tt=string)


class TestReformatDateTime(unittest.TestCase):
    def test_1(self):
        date = "2018.07.27"
        time = "16:00:36.926"
        out = "20180727 16:00:36.926"
        self.assertEqual(uf.reformat_date_time(date, time), out)


class TestProcessProcessLineMore(unittest.TestCase):
    def test_1(self):
        line = ["6406",	"2018.07.27", "16:00:29.329",
                "7", "797", "0065-0161980646"]
        out = ["6406", "20180727 16:00:29.329",
               "7", "797", "0065-0161980646"]
        self.assertEqual(uf.process_line_more_elements(line),
                                                       out)

    def test_2(self):
        line = ["6406",	"2018.07.27", "16:00:29.329",
                "7", "0065-0161980646"]
        out = ["6406", "20180727 16:00:29.329",
               "7", "0065-0161980646"]
        self.assertEqual(uf.process_line_more_elements(line),
                                                       out)


class TestProcessLine5(unittest.TestCase):
    def test1(self):
        line = ["6406",	"16:00:29.329",
                "7", "0065-0161980646"]
        date = "20180727"
        out = ["6406", "20180727 16:00:29.329",
               "7", "0065-0161980646"]
        self.assertEqual(uf.process_line_5_elements(line, date), out)


class TestReadInSingleFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        cls.out = uf.read_single_file(path, "20101010_110000.txt")

    def test_1(self):
        self.assertEqual(101, len(self.out))

    def test_2(self):
        lines_len = set([len(line) for line in self.out])
        self.assertEqual(lines_len, set([5]))

    def test_all_mice(self):
        mice = set([line[-1] for line in self.out])
        self.assertEqual(mice, set(["mouse_1"]))

    def test_last_line(self):
        last_line= ["15894", "20101010 11:59:56.218", "4", "307", "mouse_1"]
        self.assertEqual(last_line, self.out[-1])


class TestRemoveAntennas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_2_mice")
        out = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = {}
        cls.data["Antenna"] = [int(o[2]) for o in out]
        cls.data["Tag"] = [o[-1] for o in out]
        cls.data_without_1 = uf.remove_antennas(cls.data, 1)
        cls.data_without_9 = uf.remove_antennas(cls.data, 9)
        cls.data_without_1_9 = uf.remove_antennas(cls.data, [1, 9])
        cls.data_without_A = uf.remove_antennas(cls.data, "A")

    def test_remove_1_1(self):
        self.assertEqual(len(self.data["Antenna"]) - 3,
                         len(self.data_without_1["Antenna"]))

    def test_remove_1_2(self):
        self.assertEqual(len(self.data_without_1["Antenna"]),
                         len(self.data_without_1["Tag"]))

    def test_remove_1_mouse2(self):
        self.assertEqual(self.data_without_1["Tag"].count("mouse_2"),
                         self.data["Tag"].count("mouse_2"))

    def test_remove_1_mouse2_2(self):
        self.assertEqual(self.data_without_1["Tag"].index("mouse_2"),
                         self.data["Tag"].index("mouse_2"))

    def test_remove_1_mouse_3(self):
        self.assertEqual(self.data_without_1["Tag"].count("mouse_3"),
                         self.data["Tag"].count("mouse_3"))

    def test_remove_1_mouse_3_2(self):
        self.assertEqual(self.data_without_1["Tag"].index("mouse_3"),
                         self.data["Tag"].index("mouse_3") - 2)

    def test_remove_9_1(self):
        self.assertEqual(len(self.data["Antenna"]),
                         len(self.data_without_9["Antenna"]))

    def test_remove_9_2(self):
        self.assertEqual(len(self.data_without_9["Antenna"]),
                         len(self.data_without_9["Tag"]))

    def test_remove_9_mouse2(self):
        self.assertEqual(self.data_without_9["Tag"].count("mouse_2"),
                         self.data["Tag"].count("mouse_2"))

    def test_remove_9_mouse2_2(self):
        self.assertEqual(self.data_without_9["Tag"].index("mouse_2"),
                         self.data["Tag"].index("mouse_2"))

    def test_remove_9_mouse_3(self):
        self.assertEqual(self.data_without_9["Tag"].count("mouse_3"),
                         self.data["Tag"].count("mouse_3"))

    def test_remove_9_mouse_3_2(self):
        self.assertEqual(self.data_without_9["Tag"].index("mouse_3"),
                         self.data["Tag"].index("mouse_3"))

    def test_remove_1_9_1(self):
        self.assertEqual(len(self.data["Antenna"]) - 3,
                         len(self.data_without_1_9["Antenna"]))

    def test_remove_1_9_2(self):
        self.assertEqual(len(self.data_without_1_9["Antenna"]),
                         len(self.data_without_1_9["Tag"]))

    def test_remove_1_9_mouse2(self):
        self.assertEqual(self.data_without_1_9["Tag"].count("mouse_2"),
                         self.data["Tag"].count("mouse_2"))

    def test_remove_1_9_mouse2_2(self):
        self.assertEqual(self.data_without_1_9["Tag"].index("mouse_2"),
                         self.data["Tag"].index("mouse_2"))

    def test_remove_1_9_mouse_3(self):
        self.assertEqual(self.data_without_1_9["Tag"].count("mouse_3"),
                         self.data["Tag"].count("mouse_3"))

    def test_remove_1_9_mouse_3_2(self):
        self.assertEqual(self.data_without_1_9["Tag"].index("mouse_3"),
                         self.data["Tag"].index("mouse_3") - 2)

    def test_remove_A_1(self):
        self.assertEqual(len(self.data["Antenna"]),
                         len(self.data_without_A["Antenna"]))

    def test_remove_A_2(self):
        self.assertEqual(len(self.data_without_A["Antenna"]),
                         len(self.data_without_A["Tag"]))

    def test_remove_A_mouse2(self):
        self.assertEqual(self.data_without_A["Tag"].count("mouse_2"),
                         self.data["Tag"].count("mouse_2"))

    def test_remove_A_mouse2_2(self):
        self.assertEqual(self.data_without_A["Tag"].index("mouse_2"),
                         self.data["Tag"].index("mouse_2"))

    def test_remove_A_mouse_3(self):
        self.assertEqual(self.data_without_A["Tag"].count("mouse_3"),
                         self.data["Tag"].count("mouse_3"))

    def test_remove_A_mouse_3_2(self):
        self.assertEqual(self.data_without_A["Tag"].index("mouse_3"),
                         self.data["Tag"].index("mouse_3"))

if __name__ == '__main__':
    unittest.main()
