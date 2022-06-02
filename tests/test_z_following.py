# SPDX-License-Identifier: LGPL-2.1-or-later
from __future__ import print_function, division, absolute_import
import unittest
import os
from pyEcoHAB import following as fol
from pyEcoHAB.utils import general as uf
from pyEcoHAB import Loader
from pyEcoHAB import Timeline
from pyEcoHAB import sample_data, data_path
from pyEcoHAB import SetupConfig


try:
    basestring
except NameError:
    basestring = str



class TestFollowing2ndMouseInPipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas1 = ["1", "2"]
        times1 = [15, 16.5]
        antennas2 = ["8", "1", "2", "3", "4", "5"]
        times2 = [10, 16, 19, 19.5, 22, 25]
        config = SetupConfig()
        dir1 = uf.extract_directions(times1, antennas1, 3,
                                     config.directions)
        dir2 = uf.extract_directions(times2, antennas2, 6,
                                     config.directions)
        res = fol.following_single_pair(dir1, dir2)
        cls.out1, cls.time_together1, cls.intervals1 = res

        antennas1 = ["1", "2", "3", "4", "5"]
        times1 = [15, 16.5, 19, 20, 21]
        antennas2 = ["8", "1", "2", "3", "4", "5"]
        times2 = [10, 16, 19, 19.5, 22, 25]
        dir1 = uf.extract_directions(times1, antennas1, 6,
                                     config.directions)
        dir2 = uf.extract_directions(times2, antennas2, 6,
                                     config.directions)
        res = fol.following_single_pair(dir2, dir1)
        cls.out2, cls.time_together2, cls.intervals2 = res

        antennas1 = ["1", "2", "3", "4", "5",  "6", "7",
                     "8", "1", "2"]
        times1 = [15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34]
        antennas2 = ["8", "1", "2", "3", "4", "5", "6",
                     "7", "8", "1", "2"]
        times2 = [10, 16, 19, 19.5, 22, 25, 26, 27, 28, 31, 35]
        dir1 = uf.extract_directions(times1, antennas1, 3,
                                     config.directions)
        dir2 = uf.extract_directions(times2, antennas2, 3,
                                     config.directions)
        res = fol.following_single_pair(dir1, dir2)
        cls.out3, cls.time_together3, cls.intervals3 = res

    def test_following_11(self):
        self.assertEqual(self.out1, 1)

    def test_following_12(self):
        self.assertEqual(self.time_together1, 0.5)

    def test_following_13(self):
        self.assertEqual(self.intervals1, [4])

    def test_not_following_more_1(self):
        self.assertEqual(self.out2, 0)

    def test_not_following_more_2(self):
        self.assertEqual(self.time_together2, 0)

    def test_not_following_more_3(self):
        self.assertEqual(self.intervals2, [])

    def test_following_31(self):
        self.assertEqual(self.out3, 3)

    def test_following_32(self):
        self.assertEqual(self.time_together3, 4)

    def test_following_33(self):
        self.assertEqual(self.intervals3, [4, 6, 3])


class TestFollowingMatrices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ta = {
            "mouse1": [[15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34],
                       ["1", "2",  "3",    "4",  "5",  "6",  "7",
                        "8", "1", "2"]],
            "mouse2": [[10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35],
                       ["8", "1",  "2", "3", "4", "5", "6", "7", "8", "1",
                        "2"], ],
            "mouse3": [[10, 16, 17, 18, 22, 25, 26, 27],
                       ["1", "2", "3", "4",  "4",  "3",  "2", "1"], ]
        }
        mice_list = ["mouse1", "mouse2", "mouse3"]
        directions_dict = {}
        last = {
            "mouse1": "3",
            "mouse2": "3",
            "mouse3": "8",
        }
        config = SetupConfig()

        for mouse in mice_list:
            directions_dict[mouse] = uf.extract_directions(ta[mouse][0],
                                                           ta[mouse][1],
                                                           last[mouse],
                                                           config.directions)
        out = fol.following_matrices(directions_dict,
                                     mice_list,
                                     0, 1000)
        cls.following = out[0]
        cls.time_together = out[1]
        cls.interval_details = out[2]

    def test_following_diag_1(self):
        self.assertEqual(self.following["mouse1"]["mouse1"], 0)

    def test_following_diag_2(self):
        self.assertEqual(self.following["mouse2"]["mouse2"], 0)

    def test_following_diag_3(self):
        self.assertEqual(self.following["mouse3"]["mouse3"], 0)

    def test_following_0_1(self):
        self.assertEqual(self.following["mouse1"]["mouse2"], 3)

    def test_following_0_2(self):
        self.assertEqual(self.following["mouse1"]["mouse3"], 0)

    def test_following_1_0(self):
        self.assertEqual(self.following["mouse2"]["mouse1"], 0)

    def test_following_1_2(self):
        self.assertEqual(self.following["mouse2"]["mouse3"], 0)

    def test_following_2_0(self):
        self.assertEqual(self.following["mouse3"]["mouse1"], 1)

    def test_following_2_1(self):
        self.assertEqual(self.following["mouse3"]["mouse2"], 0)

    def test_time_together_diag_0(self):
        self.assertEqual(self.time_together["mouse1"]["mouse1"], 0)

    def test_time_together_diag_1(self):
        self.assertEqual(self.time_together["mouse2"]["mouse2"], 0)

    def test_time_together_diag_2(self):
        self.assertEqual(self.time_together["mouse3"]["mouse3"], 0)

    def test_time_together_0_1(self):
        self.assertEqual(self.time_together["mouse1"]["mouse2"], 0.004)

    def test_time_together_0_2(self):
        self.assertEqual(self.time_together["mouse1"]["mouse3"], 0)

    def test_time_together_1_0(self):
        self.assertEqual(self.time_together["mouse2"]["mouse1"], 0)

    def test_time_together_1_2(self):
        self.assertEqual(self.time_together["mouse2"]["mouse3"], 0)

    def test_time_together_2_0(self):
        self.assertEqual(self.time_together["mouse3"]["mouse1"], 0.001)

    def test_time_together_2_1(self):
        self.assertEqual(self.time_together["mouse3"]["mouse2"], 0)

    def test_interval_details_empty(self):
        self.assertEqual(self.time_together["mouse3"]["mouse2"], 0)



class TestExecution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = Loader(sample_data)
        cls.config = Timeline(sample_data)
        cls.uneven = Timeline(data_path, "uneven_phases.txt")

    # def test_phases(self):
    #     fol.get_dynamic_interactions(self.data, self.config, 1,
    #                                  save_distributions=True,
    #                                  save_figures=True, return_median=False,
    #                                  delimiter=";",
    #                                  save_times=True)
    # def test_phases_2(self):
    #     fol.get_dynamic_interactions(self.data, self.config, 1,
    #                                  save_distributions=False,
    #                                  save_figures=False, return_median=False,
    #                                  delimiter=";",
    #                                  save_times=False)

    # def test_ALL(self):
    #     fol.get_dynamic_interactions(self.data, self.config, 1,
    #                                  binsize="ALL")

    # def test_ALL(self):
    #     fol.get_dynamic_interactions(self.data, self.config, 1,
    #                                  binsize="ALL",
    #                                  res_dir=os.path.join(self.data.path,
    #                                                       "Resu2"),
    #                                  save_distributions=True,
    #                                  save_figures=True,
    #                                  return_median=True, delimiter=";",
    #                                  save_times=True,
    #                                  full_dir_tree=False)

    def test_short_2(self):
        fol.get_dynamic_interactions(self.data, self.config, 1,
                                     binsize=4800,
                                     res_dir=os.path.join(self.data.path,
                                                          "Resu2"),

                                     save_distributions=True,
                                     save_figures=True,
                                     return_median=True, delimiter=";",
                                     save_times=True,
                                     full_dir_tree=False)

    # def test_short_bin(self):
    #     fol.get_dynamic_interactions(self.data, self.config, 1, binsize=3600)

    # def test_long_bin(self):
    #     fol.get_dynamic_interactions(self.data, self.config, 1,
    #                                  save_distributions=True,
    #                                  save_figures=True,
    #                                  return_median=True, delimiter=";",
    #                                  save_times=True,
    #                                  seed=1,
    #                                  binsize=24*3600)

    # def test_whole_phase_uneven(self):
    #     fol.get_dynamic_interactions(self.data, self.uneven, 1,
    #                                  binsize="whole phase",
    #                                  save_distributions=True,
    #                                  save_figures=True,
    #                                  return_median=True, delimiter=";",
    #                                  save_times=True,
    #                                  seed=1)


if __name__ == '__main__':
    unittest.main()
