import unittest
import os
from datetime import datetime, tzinfo, date, time
from configparser import ConfigParser

import pyEcoHAB.utils.temporal as ut
from pyEcoHAB.utils import for_loading as fl
from pyEcoHAB import sample_data


class TestConvertIntsToTime(unittest.TestCase):
    def test_1(self):
        self.assertEqual(ut.convert_int_to_time(1), "01")

    def test_2(self):
        self.assertEqual(ut.convert_int_to_time(0), "00")

    def test_3(self):
        self.assertEqual(ut.convert_int_to_time(31), "31")


class TestFindLightBeginnig(unittest.TestCase):
    def test_1(self):
        dark_beg = "12:00"
        dark_len = 12
        self.assertEqual(ut.find_light_beginning(dark_beg, dark_len),
                         "00:00")

    def test_2(self):
        dark_beg = "12:40"
        dark_len = 12.5
        self.assertEqual(ut.find_light_beginning(dark_beg, dark_len),
                         "01:10")

    def test_3(self):
        dark_beg = "00:00"
        dark_len = 12.5
        self.assertEqual(ut.find_light_beginning(dark_beg, dark_len),
                         "12:30")


class DealWithDates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filenames = sorted(fl.get_filenames(sample_data))
        cls.first, cls.last = ut.find_first_last(cls.filenames)

    def test_first(self):
        self.assertEqual(self.first, "20140616")

    def test_last(self):
        self.assertEqual(self.last, "20140619_110000 UTC")

    def test_last_day_to_datetime(self):
        last = ut.last_day_to_datetime(self.last)
        expected = datetime(year=2014, month=6, day=19, hour=11,
                            tzinfo=tzinfo("UTF"))

    def test_strtime_to_datetime(self):
        measured = ut.strtime_to_datetime("%s12:00 UTC" % self.first)
        expected = datetime(year=2014, month=6, day=16, hour=12,
                            tzinfo=tzinfo("UTF"))

    def test_get_date(self):
        self.assertEqual(ut.get_date(date(day=16, month=6, year=2014)),
                         "16.06.2014")

    def test_get_time(self):
        self.assertEqual(ut.get_time(time(hour=12)), "12:00")

    def generate_entry(self):
        first = datetime(year=2014, month=6, day=19, hour=12,
                         tzinfo=tzinfo("UTF"))
        last = datetime(year=2014, month=6, day=20, hour=00,
                        tzinfo=tzinfo("UTF"))
        out = ut.make_cofig_entry(first, last)
        expected = {"startdate": "2014.06.19",
                    "starttime": "12:00",
                    "enddate": "2014.06.20",
                    "endtime": "00:00"}
        self.assertEqual(out, expected)


class TestGenerateTimeline(unittest.TestCase):
    def test(self):
        config_gen = ut.gen_timeline(sample_data,
                                     dark_beginning="12:00",
                                     first_phase="dark",
                                     dark_length=12,
                                     light_length=12,
                                     phase_name="EMPTY")
        config_exp = {
            "EMPTY 1 dark": {
                "startdate": "16.06.2014",
                "starttime": "12:00",
                "enddate": "17.06.2014",
                "endtime": "00:00",

            },
            "EMPTY 1 light": {
                "startdate": "17.06.2014",
                "starttime": "00:00",
                "enddate": "17.06.2014",
                "endtime": "12:00",

            },
            "EMPTY 2 dark": {
                "startdate": "17.06.2014",
                "starttime": "12:00",
                "enddate": "18.06.2014",
                "endtime": "00:00",

            },
            "EMPTY 2 light": {
                "startdate": "18.06.2014",
                "starttime": "00:00",
                "enddate": "18.06.2014",
                "endtime": "12:00",

            },
            "EMPTY 3 dark": {
                "startdate": "18.06.2014",
                "starttime": "12:00",
                "enddate": "19.06.2014",
                "endtime": "00:00",

            },
            "EMPTY 3 light": {
                "startdate": "19.06.2014",
                "starttime": "00:00",
                "enddate": "19.06.2014",
                "endtime": "12:00",

            },
            "ALL": {
                "startdate": "16.06.2014",
                "starttime": "12:00",
                "enddate": "19.06.2014",
                "endtime": "12:00"
            }
        }
        self.assertEqual(config_exp, config_gen)


if __name__ == '__main__':
    unittest.main()
