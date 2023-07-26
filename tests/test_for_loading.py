# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
#encoding: utf-8
from __future__ import print_function, division, absolute_import
import os
import glob
import unittest
import numpy as np
import pyEcoHAB.utils.for_loading as uf
import pyEcoHAB.utils.general as ut
from pyEcoHAB import data_path
from pyEcoHAB.SetupConfig import SetupConfig


SAME_PIPE = {
    "1": ["1", "2"],
    "2": ["1", "2"],
    "3": ["3", "4"],
    "4": ["3", "4"],
    "5": ["5", "6"],
    "6": ["5", "6"],
    "7": ["7", "8"],
    "8": ["7", "8"]
}
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

OPPOSITE_PIPE = {
    "1": ["5", "6"],
    "2": ["5", "6"],
    "3": ["7", "8"],
    "4": ["7", "8"],
    "5": ["1", "2"],
    "6": ["1", "2"],
    "7": ["3", "4"],
    "8": ["3", "4"]
}

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

NON_ADJACENT = {
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


class TestParseFilename(unittest.TestCase):
    def test_normal1(self):
        fname = "20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(hour, "120000")

    def test_normal2(self):
        fname = "20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(date, "20190403")

    def test_normal3(self):
        fname = "20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(datenext, "20190404")
        
    def test_normal4(self):
        fname = "20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(setup, "")

        
    def test_weird(self):
        fname = "20190403_120000_0001.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(hour, "120000")
        self.assertEqual(date, "20190403")
        self.assertEqual(datenext, "20190404")

    def test_throw_exception(self):
        fname = "20190403_120000_0001_kk.txt"
        self.assertRaises(ValueError, uf.parse_fname, fname=fname)

    def test_new_format1(self):
        fname = "COM8_20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(hour, "120000")

    def test_new_format2(self):
        fname = "COM8_20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(date, "20190403")

    def test_new_format3(self):
        fname = "COM8_20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(datenext, "20190404")

    def test_new_format4(self):
        fname = "COM8_20190403_120000.txt"
        hour, date, datenext, setup = uf.parse_fname(fname)
        self.assertEqual(setup, "COM8")

        


class TestPrintHumanTime(unittest.TestCase):
    def test_date(self):
        tt = 1554247067
        self.assertEqual('Tue Apr  2 23:17:47 2019', uf.print_human_time(tt))


class TestTimeToSec(unittest.TestCase):
    def test_sec(self):
        string = "20190709 20:05:13.333"
        self.assertEqual(1562702713.333, uf.time_to_sec(string))

    def test_sec_2(self):
        string = "20190709 20:05:13"
        self.assertEqual(1562702713., uf.time_to_sec(string))

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
        cls.out, setup = uf.read_single_file(path, "20101010_110000.txt")

    def test_1(self):
        self.assertEqual(101, len(self.out))

    def test_2(self):
        lines_len = set([len(line) for line in self.out])
        self.assertEqual(lines_len, set([5]))

    def test_all_mice(self):
        mice = set([line[-1] for line in self.out])
        self.assertEqual(mice, set(["mouse_1"]))

    def test_last_line(self):
        last_line = ["15894", "20101010 11:59:56.218", "4",
                     "307", "mouse_1"]
        self.assertEqual(last_line, self.out[-1])


class TestRemoveGhostTags(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        cls.data, setup = uf.read_single_file(path, "20101010_110000.txt")

    def test_removing_tags_1(self):
        out = uf.remove_ghost_tags(self.data, "mouse_1")
        tags = set()
        for line in out:
            tags.add(line[4])
        self.assertEqual(list(tags), ["mouse_1"])

    def test_removing_tags_1(self):
        out = uf.remove_ghost_tags(self.data, ["mouse_1", "mouse_2"])
        tags = set()
        for line in out:
            tags.add(line[4])
        self.assertEqual(sorted(tags), ["mouse_1", "mouse_2"])

    def test_removing_tags_default(self):
        out = uf.remove_ghost_tags(self.data, "ALL")
        tags = set()
        for line in out:
            tags.add(line[4])
        self.assertEqual(sorted(tags), ["mouse_1", "mouse_2", "mouse_3"])


class TestRemoveAntenna(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(raw_data)

    def test_no_antenna(self):
        data = uf.remove_one_antenna(self.data, None)
        self.assertTrue((data == self.data).all())

    def test_single_antenna_1(self):
        data = uf.remove_one_antenna(self.data, "1")
        self.assertEqual(data.shape[0], self.data.shape[0] - 3)

    def test_single_antenna_2(self):
        data = uf.remove_one_antenna(self.data, "1")
        self.assertFalse(1 in self.data["Antennas"])

    def test_single_antenna_2(self):
        data = uf.remove_one_antenna(self.data, "1")
        mice = set(data["Tag"])
        all_mice = set(self.data["Tag"])
        self.assertEqual(mice, all_mice)

    def test_nonexistent_antenna(self):
        data = uf.remove_one_antenna(self.data, "9")
        self.assertTrue((data == self.data).all())

    def test_not_antenna(self):
        data = uf.remove_one_antenna(self.data, "gugu")
        self.assertTrue((data == self.data).all())


class TestRemoveAntennas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(raw_data)

    def test_single_antenna_1(self):
        data = uf.remove_antennas(self.data, "1")

    def test_remove_antenna_list_1(self):
        data = uf.remove_antennas(self.data, ["1", "9"])
        self.assertEqual(data.shape[0], self.data.shape[0]-3)

    def test_remove_antenna_list_1(self):
        data = uf.remove_antennas(self.data, ["1", None])
        self.assertEqual(data.shape[0], self.data.shape[0]-3)

    def test_remove_antenna_list_2_1(self):
        data = uf.remove_antennas(self.data, ["1", "2"])
        self.assertEqual(self.data.shape[0] - 6, data.shape[0])

    def test_remove_antenna_list_3(self):
        data = uf.remove_antennas(self.data, ["1", "2"])
        mice = set(data["Tag"])
        self.assertEqual(mice,
                         set(["mouse_1", "mouse_2"]))


class TestTransformRaw(unittest.TestCase):
    def test_date_1(self):
        row = [1, "gugu", 2, 222, "AAA"]
        self.assertRaises(ValueError, uf.transform_raw, row)

    def test_correct(self):
        row = [1, "20101010 11:00:49.020", 5, '102', 'mouse_3']
        time = uf.time_to_sec(row[1])
        out = uf.transform_raw(row)
        out2 = (1, time,  5, 102, "mouse_3")
        self.assertEqual(out, out2)


class TestTransformAllData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        cls.raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(cls.raw_data)

    def test_len(self):
        self.assertEqual(len(self.raw_data), self.data.shape[0])

    def test_width(self):
        width = set([len(data) for data in self.raw_data])
        self.assertEqual(width, set([len(self.data.dtype)]))

    def test_first_line(self):
        line = uf.transform_raw(self.raw_data[0])
        self.assertEqual(line, self.data[0].tolist())


class TestAntennaMismatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        config = SetupConfig()
        cls.mismatch1 = uf.antenna_mismatch(data, config)

    def test_1(self):
        self.assertEqual(2, self.mismatch1["3 6"])

    def test_2(self):
        self.assertEqual(1, self.mismatch1["1 3"])

    def test_3(self):
        self.assertEqual(1, self.mismatch1["4 6"])

    def test_all_others(self):
        out = sum(list(self.mismatch1.values()))
        self.assertEqual(4, out)

    def test_empty(self):
        self.assertRaises(Exception, uf.antenna_mismatch, [])


class TestSkippedAntennas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        config = SetupConfig()
        cls.mismatch1 = uf.skipped_registrations(data, config)

    def test_1(self):
        self.assertEqual(2, self.mismatch1["skipped two"])

    def test_2(self):
        self.assertEqual(2, self.mismatch1["skipped one"])

    def test_all_others(self):
        out = sum(list(self.mismatch1.values()))
        self.assertEqual(4, out)

    def test_empty(self):
        self.assertRaises(Exception, uf.skipped_registrations, [])


class TestCheckAntennaPresence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.config = SetupConfig()
        cls.presences = uf.check_antenna_presence(data, cls.config, 24*3600)
        cls.end = data["Time"][-1]
        cls.begs = []
        for key in cls.presences.keys():
            if len(cls.presences[key]):
                cls.begs.append(cls.presences[key][0][0])
        cls.begs_data = []
        for antenna in ["1", "2", "3", "4", "5", "6", "7", "8"]:
            if antenna == "4":
                continue
            idx = np.where(data["Antenna"] == antenna)[0][-1]
            cls.begs_data.append(np.round(data["Time"][idx]))

    def test_1(self):
        self.assertEqual(len(self.presences["4"]), 0)

    def test_end(self):
        ends = []
        for key in self.presences.keys():
            if len(self.presences[key]):
                ends.append(self.presences[key][0][-1])
        self.assertEqual(set(ends), set([self.end]))

    def test_beg_1(self):
        self.assertEqual(set(self.begs), set(self.begs_data))

    def test_beg_2(self):
        self.assertEqual(len(set(self.begs)), len(self.begs))

    def test_empty(self):
        self.assertRaises(Exception, uf.check_antenna_presence, [],
                          self.config, 2)


class TestTotalMismatches(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = SetupConfig()
        cls.mismatched_pairs = config.mismatched_pairs
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data1, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data1 = uf.from_raw_data(raw_data1)
        cls.mismatch1 = uf.antenna_mismatch(cls.data1, config)
        path = os.path.join(data_path, "weird_short")
        raw_data2, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data2 = uf.from_raw_data(raw_data2)
        cls.mismatch2 = uf.antenna_mismatch(cls.data2, config)

    def test_no_mismatches(self):
        res = uf.total_mismatches(self.mismatch2)
        correct = {}
        antennas = sorted(set(self.data2["Antenna"]))
        for antenna in antennas:
            correct[antenna] = 0
        self.assertEqual(correct, res)

    def test_mismatches(self):
        res = uf.total_mismatches(self.mismatch1)
        correct = {}
        antennas = sorted(set(self.data1["Antenna"]))
        for antenna in antennas:
            correct[antenna] = 0
        mice = sorted(set(self.data1["Tag"]))
        all_antennas = self.data1["Antenna"]
        for mouse in mice:
            idxs = np.where(self.data1["Tag"] == mouse)[0]
            for i, idx in enumerate(idxs[:-1]):
                idx_next = idxs[i+1]
                a1, a2 = all_antennas[idx], all_antennas[idx_next]
                key = "%s %s" % (min(a1, a2),  max(a1, a2))
                if key in self.mismatched_pairs:
                    correct[a1] += 1
                    correct[a2] += 1
        self.assertEqual(correct, res)


class TestRunDiagnostics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = SetupConfig()
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.mismatch1 = uf.antenna_mismatch(data, config)
        cls.presences1 = uf.check_antenna_presence(data, config, 24*3600)
        res_path = os.path.join(path, "Results")
        files = glob.glob(os.path.join(res_path + "/diagnostics/*.csv"))
        for f in files:
            print("rm ", f)
        cls.length = len(data["Antenna"])
        out1 = uf.run_diagnostics(data, 24*3600, res_path, config)
        cls.str11, cls.str12, cls.str13, cls.str14, cls.str15 = out1
        path = os.path.join(data_path, "weird_short")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.mismatch2 = uf.antenna_mismatch(data, config)
        cls.presences2 = uf.check_antenna_presence(data, config, 24*3600)
        res_path = os.path.join(path, "Results")
        files = glob.glob(os.path.join(res_path + "/diagnostics/*.csv"))
        for f in files:
            os.remove(f)

        out2 = uf.run_diagnostics(data, 24*3600, res_path, config)
        cls.str21, cls.str22, cls.str23, cls.str24, cls.str25 = out2
        path1 = os.path.join(data_path, "weird_very_short_3_mice")
        raw_data1, setup = uf.read_single_file(path1, "20101010_110000.txt")
        data1 = uf.from_raw_data(raw_data1)
        config1 = SetupConfig()
        res_path1 = os.path.join(path1, "Results")
        out = uf.run_diagnostics(data1, 24*3600, res_path1, config1)
        cls.incorr_tunnel = out[-1]

    def test_mismatched_tunnels(self):
        expected = u"tunnel, count, percentage of all passings through the tunnel\n1 2, 0, 0.00 per 100\n3 4, 0, 0.00 per 100\n5 6, 2, 29.00 per 100\n7 8, 0, 0.00 per 100\n"
        self.assertEqual(self.incorr_tunnel, expected)

    def test_no_registration_breaks_file(self):
        path = os.path.join(data_path, "weird_short")
        res_path = os.path.join(path, "Results")
        f_path = os.path.join(res_path +
                              "/diagnostics/breaks_in_registrations.csv")
        f = open(f_path)
        read = f.read()
        self.assertEqual(read, self.str22)

    def test_no_mismatch_antennas_file(self):
        path = os.path.join(data_path, "weird_short")
        res_path = os.path.join(path, "Results")
        f_path = os.path.join(res_path +
                              "/diagnostics/antenna_mismatches.csv")
        f = open(f_path)
        read = f.read()
        self.assertEqual(read, self.str21)

    def test_total_no_mismatch_antennas_file(self):
        path = os.path.join(data_path, "weird_short")
        res_path = os.path.join(path, "Results")
        f_path = os.path.join(res_path +
                              "/diagnostics/incorrect_antenna_transitions.csv")
        f = open(f_path)
        read = f.read()
        self.assertEqual(read, self.str23)

    def test_registration_breaks_file(self):
        path = os.path.join(data_path, "weird_short_3_mice")
        res_path = os.path.join(path, "Results")
        f_path = os.path.join(res_path +
                              "/diagnostics/breaks_in_registrations.csv")
        f = open(f_path)
        read = f.read()
        self.assertEqual(read, self.str12)

    def test_mismatch_antennas_file(self):
        path = os.path.join(data_path, "weird_short_3_mice")
        res_path = os.path.join(path, "Results")
        f_path = os.path.join(res_path +
                              "/diagnostics/antenna_mismatches.csv")
        f = open(f_path)
        read = f.read()
        self.assertEqual(read, self.str11)

    def test_total_mismatch_antennas_file(self):
        path = os.path.join(data_path, "weird_short_3_mice")
        res_path = os.path.join(path, "Results")
        f_path = os.path.join(res_path +
                              "/diagnostics/incorrect_antenna_transitions.csv")
        f = open(f_path)
        read = f.read()
        self.assertEqual(read, self.str13)

    def test_no_mismatch_string(self):
        out = "antenna pair,  count, percentage\n"
        for pair in uf.PAIRS:
            out += "%s, %d, %3.2f per 100\n" % (pair, 0, 0.00)
        self.assertEqual(self.str21, out)

    def test_mismatch_string(self):
        out = "antenna pair,  count, percentage\n"
        for pair in uf.PAIRS:
            if pair in ["3 6", "1 3", "4 6"]:
                exact_mis = np.round(100*self.mismatch1[pair]/self.length)
                out += "%s, %d, %3.2f per 100\n" % (pair,
                                                    self.mismatch1[pair],
                                                    exact_mis)
            else:
                out += "%s, %d, %3.2f per 100\n" % (pair, 0, 0.00)
        self.assertEqual(out, self.str11)

    def test_no_skipped_string(self):
        out = u"type, count, percentage\n"
        for pair in ["skipped one", "skipped two", "skipped more"]:
            out += u"%s, %d, %3.2f per 100\n" % (pair, 0, 0.00)
        self.assertEqual(self.str24, out)

    def test_skipped_string(self):
        out = u"type, count, percentage\n"
        for pair in ["skipped one", "skipped two", "skipped more"]:
            exact_mis = np.round(100*2/self.length)
            if pair == "skipped more":
                out += u"%s, %d, %3.2f per 100\n" % (pair,
                                                     0,
                                                     0.0)
            else:
                out += u"%s, %d, %3.2f per 100\n" % (pair,
                                                     2,
                                                     exact_mis)
        self.assertEqual(out, self.str14)


class TestTransformVisits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.data = ut.get_animal_position(data["Time"],
                                          data["Antenna"],
                                          "mouse_1", 2,
                                          same_pipe=SAME_PIPE,
                                          same_address=SAME_ADDRESS,
                                          opposite_pipe=OPPOSITE_PIPE,
                                          address=ADDRESS,
                                          surrounding=SURROUNDING,
                                          address_not_adjacent=NON_ADJACENT,
                                          internal_antennas=[])
        cls.visits = uf.transform_visits(cls.data)

    def test1(self):
        self.assertEqual(len(self.data), len(self.visits))

    def test2(self):
        self.assertEqual(len(self.visits.dtype), len(self.data[0]))

    def test_mice(self):
        out = set([row[1] for row in self.data])
        self.assertEqual(set(self.visits["Tag"]), out)


class TestRenameAntennas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.old_data = uf.from_raw_data(raw_data)
        cls.data = uf.rename_antennas("setup1", cls.old_data)

    def test_1(self):
        out = []
        for i, row in enumerate(self.data["Antenna"]):
            out.append(row == self.old_data["Antenna"][i]+"_setup1")
        self.assertEqual(set(out), set([True]))

    def test_2(self):
        self.assertEqual(self.data.dtype, self.old_data.dtype)

    def test_3(self):
        self.assertTrue(np.all(self.data["Id"] == self.old_data["Id"]))

    def test_4(self):
        self.assertTrue(np.all(self.data["Time"] == self.old_data["Time"]))

    def test_5(self):
        self.assertTrue(np.all(self.data["Duration"] ==
                               self.old_data["Duration"]))

    def test_6(self):
        self.assertTrue(np.all(self.data["Tag"] == self.old_data["Tag"]))


class TestAppendData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data1 = uf.from_raw_data(raw_data)
        cls.data2 = uf.from_raw_data(raw_data)
        cls.data2["Time"] += 15*60
        cls.combined_data = uf.append_data_sources([cls.data1, cls.data2])

    def test1(self):
        self.assertEqual(len(self.data1) + len(self.data2),
                         len(self.combined_data))

    def test2(self):
        self.assertTrue(np.all((self.combined_data["Time"][1:] -
                                self.combined_data["Time"][:-1]) > 0))

    def test_mouse3(self):
        indx = np.where(self.combined_data["Tag"] == "mouse_3")[0]
        line1 = self.combined_data[indx[0]].copy()
        line2 = self.combined_data[indx[1]].copy()
        line1["Time"] += 15*60
        self.assertTrue(np.all(line1 == line2))


class TestTunnelErrors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas = ["1", "2", "1", "2", "3", "4", "5", "6", "7"]
        times = [1,   2,  2.5,   3,  4.5, 5.5, 6.5, 7.5, 10.5]
        durations = [3, 600,  3,    34,  55,  66, 1999, 200, 100]
        cls.pred_out = {"1 2": 1, "3 4": 0, "5 6": 1, "7 8": 0}
        keys = cls.pred_out.keys()
        cls.out, cls.tot = uf.incorrect_tunnel_single_mouse(keys,
                                                            antennas,
                                                            times, durations)
        cls.pred_tot = {"1 2": 3, "3 4": 1, "5 6": 1, "7 8": 0}
        path = os.path.join(data_path, "weird_very_short_3_mice")
        cls.raw_data, setup = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(cls.raw_data)
        config = SetupConfig()
        cls.out_i, cls.out_tot_i = uf.incorrect_tunnel_registrations(cls.data,
                                                                     config)
        cls.pred_out_i = {"1 2": 0, "3 4": 0, "5 6": 2, "7 8": 0}
        cls.pred_tot_i = {"1 2": 2, "3 4": 1, "5 6": 7, "7 8": 0}

    def test_incorrect(self):
        self.assertEqual(self.pred_out, self.out)

    def test_total(self):
        self.assertEqual(self.pred_tot, self.tot)

    def test_data_incorrect(self):
        self.assertEqual(self.pred_out_i, self.out_i)

    def test_data_total(self):
        self.assertEqual(self.pred_tot_i, self.out_tot_i)


if __name__ == '__main__':
    unittest.main()
