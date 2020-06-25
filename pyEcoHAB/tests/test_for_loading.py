#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
import os
import glob
import unittest
import random
import numpy as np

import pyEcoHAB.utils.for_loading as uf
import pyEcoHAB.utility_functions as utils
from pyEcoHAB import data_path
from pyEcoHAB.SetupConfig import SetupConfig

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

class TestRemoveGhostTags(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
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
        line = ["15151", "20101010 11:11:06.748","3", "204", "mouse_1"]
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
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        cls.data = uf.from_raw_data(raw_data)

    def test_no_antenna(self):
        data = uf.remove_one_antenna(self.data, None)
        self.assertTrue((data==self.data).all())

    def test_single_antenna_1(self):
        data = uf.remove_one_antenna(self.data, "1")
        self.assertEqual(data.shape[0], self.data.shape[0]-3)

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
        self.assertTrue((data==self.data).all())

    def test_not_antenna(self):
        data = uf.remove_one_antenna(self.data, "gugu")
        self.assertTrue((data==self.data).all())


class TestRemoveAntennas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
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
        cls.raw_data = uf.read_single_file(path, "20101010_110000.txt")
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
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        config = SetupConfig() 
        cls.mismatch1 = uf.antenna_mismatch(data, config.mismatched_pairs)

    def test_1(self):
        self.assertEqual(2, self.mismatch1["3 6"])

    def test_2(self):
        self.assertEqual(1, self.mismatch1["1 3"])

    def test_3(self):
        self.assertEqual(1, self.mismatch1["4 6"])

    def test_all_others(self):
        out = sum(list(self.mismatch1.values()))
        self.assertEqual(4, out)


class TestCheckAntennaPresence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.presences = uf.check_antenna_presence(data, 24*3600)
        cls.end = data["Time"][-1]
        cls.begs = []
        for key in cls.presences.keys():
            if len(cls.presences[key]):
                cls.begs.append(cls.presences[key][0][0])
        cls.begs_data = []
        for antenna in ["1", "2", "3", "4", "5", "6", "7", "8"]:
            if antenna == "4":
                continue
            idx = np.where(data["Antenna"]==antenna)[0][-1]
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


class TestRunDiagnostics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = SetupConfig()
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.mismatch1 = uf.antenna_mismatch(data, config.mismatched_pairs)
        cls.presences1 = uf.check_antenna_presence(data, 24*3600)
        res_path = os.path.join(path, "Results")
        files = glob.glob(os.path.join(res_path + "/diagnostics/*.csv"))
        for f in files:
            print("rm ", f)
        cls.length = len(data["Antenna"])
        cls.str11, cls.str12 = uf.run_diagnostics(data, 24*3600, res_path,
        config.mismatched_pairs)
        
        path = os.path.join(data_path, "weird_short")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.mismatch2 = uf.antenna_mismatch(data, config.mismatched_pairs)
        cls.presences2 = uf.check_antenna_presence(data, 24*3600)
        res_path = os.path.join(path, "Results")
        files = glob.glob(os.path.join(res_path + "/diagnostics/*.csv"))

        for f in files:
            os.remove(f)

        cls.str21, cls.str22 = uf.run_diagnostics(data, 24*3600, res_path,
                                                  config.mismatched_pairs)


    # def test_no_registration_breaks(self):
    #     all_antennas_registered = "Breaks in registrations on antennas:\n"
    #     for antenna in range(1, 9):
    #         all_antennas_registered += "%d:\n" % antenna
    #     self.assertEqual(all_antennas_registered, self.str22)

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

    def test_no_mismatch_string(self):
        out = "First reading, consecutive reading,  count, percentage\n"
        for pair in uf.PAIRS:
            out+= "%s,\t%d, %3.2f per 100\n"% (pair, 0, 0.00)
        self.assertEqual(self.str21, out)

    def test_mismatch_string(self):
        out = "First reading, consecutive reading,  count, percentage\n"
        for pair in uf.PAIRS:
            if pair in ["3 6", "1 3", "4 6"]:
                exact_mis = np.round(100*self.mismatch1[pair]/self.length)
                out+= "%s,\t%d, %3.2f per 100\n" % (pair, self.mismatch1[pair],
                                                    exact_mis)
            else:
                out+= "%s,\t%d, %3.2f per 100\n" % (pair, 0, 0.00)
        self.assertEqual(out, self.str11)

    # def test_presence_string(self):
    #     out = u'Breaks in registrations on antennas:\n'
    #     for antenna in range(1, 9):
    #         out += u"%d:\n" % int(antenna)
    #         for breaks in self.presences1[str(antenna)]:
    #             out += u"%s %s, %4.2f h\n" % (uf.print_human_time(breaks[0]),
    #                                          uf.print_human_time(breaks[1]),
    #                                          (breaks[1] - breaks[0])/3600)
    #     self.assertEqual(out, self.str12)


class TestTransformVisits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        data = uf.from_raw_data(raw_data)
        cls.data =  utils.get_animal_position(data["Time"],
                                              data["Antenna"],
                                              "mouse_1", 2)
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
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
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
        self.assertTrue(np.all(self.data["Duration"] == self.old_data["Duration"]))

    def test_6(self):
        self.assertTrue(np.all(self.data["Tag"] == self.old_data["Tag"]))


class TestAppendData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.path.join(data_path, "weird_short_3_mice")
        raw_data = uf.read_single_file(path, "20101010_110000.txt")
        cls.data1 = uf.from_raw_data(raw_data)
        cls.data2 = uf.from_raw_data(raw_data)
        cls.data2["Time"] += 15*60
        cls.combined_data = uf.append_data_sources([cls.data1, cls.data2])

    def test1(self):
        self.assertEqual(len(self.data1) + len(self.data2),
                         len(self.combined_data))
    def test2(self):
        self.assertTrue(np.all((self.combined_data["Time"][1:] - self.combined_data["Time"][:-1]) > 0))

    def test_mouse3(self):
        indx = np.where(self.combined_data["Tag"] == "mouse_3")[0]
        line1 = self.combined_data[indx[0]].copy()
        line2 = self.combined_data[indx[1]].copy()
        line1["Time"] += 15*60
        self.assertTrue(np.all(line1==line2))

if __name__ == '__main__':
    unittest.main()

