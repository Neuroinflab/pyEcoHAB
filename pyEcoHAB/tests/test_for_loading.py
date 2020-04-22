#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
import os
import unittest
import random
import numpy as np

import pyEcoHAB.utils.for_loading as uf
from pyEcoHAB import data_path

STANDARD_ANTENNAS = {'1': 1, '2': 2,
                         '3': 3, '4': 4,
                         '5': 5, '6': 6,
                         '7': 7, '8': 8}

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

class TestRemoveGhostTags(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_2_mice")
        cls.data = uf.read_single_file(path, "20101010_110000.txt")

    def test_1(self):
        out = uf.remove_ghost_tags(self.data, 0, 0)
        self.assertEqual(self.data, out)

    def test_1_app(self):
        out = uf.remove_ghost_tags(self.data, 2, 0)
        tags = [o[4] for o in out]
        self.assertTrue("mouse_3" not in tags)

    def test_1_app_2(self):
        out = uf.remove_ghost_tags(self.data, 2, 0)
        line = ["15151", "20101010 11:11:06.748", "3", "204", "mouse_1"]
        self.assertTrue(line, out[32])

    def test_1_app_3(self):
        out = uf.remove_ghost_tags(self.data, 2, 0)
        self.assertEqual(len(out), len(self.data) - 1)

    def test_1_app_4(self):
        out = uf.remove_ghost_tags(self.data, 2, 0)
        tags = [o[4] for o in out]
        self.assertEqual(set(tags), set(["mouse_1", "mouse_2"]))

    def test_2_app(self):
        out = uf.remove_ghost_tags(self.data, 3, 0)
        tags = [o[4] for o in out]
        self.assertEqual(set(tags), set(["mouse_1"]))

    def test_2_app_2(self):
        out = uf.remove_ghost_tags(self.data, 2, 0)
        line = ["15151", "20101010 11:11:06.748", "3", "204", "mouse_1"]
        self.assertTrue(line, out[29])

    def test_2_app_3(self):
        out = uf.remove_ghost_tags(self.data, 2, 0)
        line = ["15045", "20101010 1:06:23.809", "3",	"1383",	"mouse_1"]
        self.assertTrue(line, out[10])

    def test_none_app(self):
        out = uf.remove_ghost_tags(self.data, 200, 0)
        self.assertEqual(out, [])

    def test_removing_tags(self):
         out = uf.remove_ghost_tags(self.data, 0, 0, "mouse_1")
         lines = [
             ["15040", "20101010 11:06:08.349", "5", "307", "mouse_2"],
             ["15043", "20101010 11:06:21.731", "4", "259", "mouse_2"],
             ["15123", "20101010 11:09:39.065", "2", "256", "mouse_3"],
         ]
         self.assertEqual(out, lines)

    def test_removing_tags_2(self):
         out = uf.remove_ghost_tags(self.data, 0, 0, ["mouse_1"])
         lines = [
             ["15040", "20101010 11:06:08.349", "5", "307", "mouse_2"],
             ["15043", "20101010 11:06:21.731", "4", "259", "mouse_2"],
             ["15123", "20101010 11:09:39.065", "2", "256", "mouse_3"],
         ]
         self.assertEqual(out, lines)

    def test_days_1(self):
        out = uf.remove_ghost_tags(self.data, 0, 1, ["mouse_1"])
        self.assertEqual(out, [])

    def test_days_2(self):
        out = uf.remove_ghost_tags(self.data, 0, 1)
        self.assertEqual(len(out), len(self.data) - 3)

    def test_days_3(self):
        out = uf.remove_ghost_tags(self.data, 0, 1)
        line = ["15894", "20101011 11:59:56.218", "4", "307", "mouse_1"]
        self.assertEqual(out[-1], line)

class TestRemoveAntenna(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_2_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(raw_data, STANDARD_ANTENNAS)

    def test_no_antenna(self):
        data = uf.remove_one_antenna(self.data, None)
        self.assertTrue((data==self.data).all())

    def test_single_antenna_1(self):
        data = uf.remove_one_antenna(self.data, 1)
        self.assertEqual(data.shape[0], self.data.shape[0]-3)

    def test_single_antenna_2(self):
        data = uf.remove_one_antenna(self.data, 1)
        self.assertFalse(1 in self.data["Antennas"])

    def test_single_antenna_2(self):
        data = uf.remove_one_antenna(self.data, 1)
        mice = set(data["Tag"])
        all_mice = set(self.data["Tag"])
        self.assertEqual(mice, all_mice)

    def test_nonexistent_antenna(self):
        data = uf.remove_one_antenna(self.data, 9)
        self.assertTrue((data==self.data).all())

    def test_not_antenna(self):
        data = uf.remove_one_antenna(self.data, "gugu")
        self.assertTrue((data==self.data).all())


class TestRemoveAntennas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_2_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(raw_data, STANDARD_ANTENNAS)

    def test_single_antenna_1(self):
        data = uf.remove_antennas(self.data, 1)

    def test_remove_antenna_list_1(self):
        data = uf.remove_antennas(self.data, [1, 9])
        self.assertEqual(data.shape[0], self.data.shape[0]-3)

    def test_remove_antenna_list_1(self):
        data = uf.remove_antennas(self.data, [1, None])
        self.assertEqual(data.shape[0], self.data.shape[0]-3)

    def test_remove_antenna_list_2_1(self):
        data = uf.remove_antennas(self.data, [1, 2])
        self.assertEqual(self.data.shape[0] - 6, data.shape[0])

    def test_remove_antenna_list_3(self):
        data = uf.remove_antennas(self.data, [1, 2])
        mice = set(data["Tag"])
        self.assertEqual(mice,
                         set(["mouse_1", "mouse_2"]))


class TestTransformRaw(unittest.TestCase):
    def test_date_1(self):
        row = [1, "gugu", "2", 222, "AAA"]
        self.assertRaises(ValueError, uf.transform_raw, row,
                          STANDARD_ANTENNAS)

    def test_correct(self):
        row = [1, "20101010 11:00:49.020", '5', '102', 'mouse_3']
        time = uf.time_to_sec(row[1])
        out = uf.transform_raw(row, STANDARD_ANTENNAS)
        out2 = (1, time, 5, 102, "mouse_3")
        self.assertEqual(out, out2)


class TestTransformAllData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_2_mice")
        cls.raw_data = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(cls.raw_data, STANDARD_ANTENNAS)

    def test_len(self):
        self.assertEqual(len(self.raw_data), self.data.shape[0])

    def test_width(self):
        width = set([len(data) for data in self.raw_data])
        self.assertEqual(width, set([len(self.data.dtype)]))

    def test_first_line(self):
        line = uf.transform_raw(self.raw_data[0], STANDARD_ANTENNAS)
        self.assertEqual(line, self.data[0].tolist())



if __name__ == '__main__':
    unittest.main()

