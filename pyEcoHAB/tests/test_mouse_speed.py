from __future__ import print_function, division, absolute_import
from pyEcoHAB import mouse_speed as ms
import unittest

try:
    basestring
except NameError:
    basestring = str

class TestCheck2ndMouse(unittest.TestCase):
    def test_following(self):
        antenna1 = 1
        antenna2 = 2
        t1 = 15
        threshold = 3
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        out = ms.check_2nd_mouse(antenna1, antenna2, t1, threshold,
                                 antennas2, times2)
        self.assertEqual(out, 1)

    def test_taking_over(self):
        antenna1 = 1
        antenna2 = 2
        t1 = 15
        threshold = 3
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 18, 19.5, 22, 25]
        out = ms.check_2nd_mouse(antenna1, antenna2, t1, threshold,
                                 antennas2, times2)
        self.assertEqual(out, 0)

    def test_not_following(self):
        antenna1 = 1
        antenna2 = 2
        t1 = 15
        threshold = 3
        antennas2 = [8, 1, 2, 3, 4, 5, 6]
        times2 = [16, 19, 19.5, 22, 25]
        out = ms.check_2nd_mouse(antenna1, antenna2, t1, threshold,
                                 antennas2, times2)
        self.assertEqual(out, 0)


class TestFollowing2ndMouseInPipe(unittest.TestCase):
    KEY_DICT = {
        "12": 0,
        "21": 0,
        "34": 0,
        "43": 0,
        "56": 0,
        "65": 0,
        "78": 0,
        "87": 0,
    }
    @classmethod
    def setUpClass(cls):
        antennas1 = [1, 2]
        times1 = [15, 16.5]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        res = ms.following_2_mice_in_pipe(antennas1,
                                          times1,
                                          antennas2,
                                          times2)
        cls.out1, cls.time_together1, cls.intervals1, cls.deltas_t11, cls.deltas_t21, cls.tot_deltas_t11, cls.followings_mouse21 = res

        antennas1 = [1, 2, 3, 4, 5]
        times1 = [15, 16.5, 19, 20, 21]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        res = ms.following_2_mice_in_pipe(antennas2,
                                          times2,
                                          antennas1,
                                          times1)
        cls.out2, cls.time_together2, cls.intervals2, cls.deltas_t12, cls.deltas_t22, cls.tot_deltas_t12, cls.followings_mouse22 = res

        antennas1 = [1, 2,   3, 4,  5,  6,  7,   8,  1, 2]
        times1 = [15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34 ]
        antennas2 = [8, 1,   2,    3,  4,  5,  6,   7, 8,   1,   2]
        times2 =   [10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35]
        res = ms.following_2_mice_in_pipe(antennas1,
                                          times1,
                                          antennas2,
                                          times2)
        cls.out3, cls.time_together3, cls.intervals3, cls.deltas_t13, cls.deltas_t23, cls.tot_deltas_t13, cls.followings_mouse23 = res

    def test_following_11(self):
        self.assertEqual(self.out1, 1)

    def test_following_12(self):
        self.assertEqual(self.time_together1, 0.5)

    def test_following_13(self):
        self.assertEqual(self.intervals1, [4])

    def test_following_14(self):
        deltas_t1_out = self.KEY_DICT.copy()
        deltas_t1_out["12"] = 1.5
        self.assertEqual(self.deltas_t11, deltas_t1_out)

    def test_following_15(self):
        deltas_t2_out = self.KEY_DICT.copy()
        deltas_t2_out["12"] = 3
        self.assertEqual(self.deltas_t21, deltas_t2_out)

    def test_following16(self):
        tot_times =  self.KEY_DICT.copy()
        tot_times["12"] =  1.5
        self.assertEqual(tot_times, self.tot_deltas_t11)

    def test_following17(self):
        followings_mouse2_out = self.KEY_DICT.copy()
        followings_mouse2_out["12"] = 1
        self.assertEqual(self.followings_mouse21,
                         followings_mouse2_out)

    def test_not_following_more_1(self):
        self.assertEqual(self.out2, 0)

    def test_not_following_more_2(self):
        self.assertEqual(self.time_together2, 0)

    def test_not_following_more_3(self):
        self.assertEqual(self.intervals2, [])

    def test_not_following_more_4(self):
        self.assertEqual(self.deltas_t12, self.KEY_DICT)

    def test_not_following_more_5(self):
        self.assertEqual(self.deltas_t22, self.KEY_DICT)

    def test_not_following_more_6(self):
        tot_deltas_t1 = self.KEY_DICT.copy()
        tot_deltas_t1["12"] = 3
        tot_deltas_t1["34"] = 2.5
        self.assertEqual(self.tot_deltas_t12, tot_deltas_t1)

    def test_not_following_more_7(self):
        self.assertEqual(self.followings_mouse22, self.KEY_DICT)

    def test_following_31(self):
        self.assertEqual(self.out3, 3)

    def test_following_32(self):
        self.assertEqual(self.time_together3, 4)

    def test_following_33(self):
        self.assertEqual(self.intervals3, [4, 3, 6])

    def test_following_34(self):
        deltas_t1_out = self.KEY_DICT.copy()
        deltas_t1_out["12"] = 6.5
        deltas_t1_out["34"] = 1
        self.assertEqual(self.deltas_t13, deltas_t1_out)

    def test_following_35(self):
        deltas_t2_out = self.KEY_DICT.copy()
        deltas_t2_out["12"] = 7
        deltas_t2_out["34"] = 2.5
        self.assertEqual(self.deltas_t23, deltas_t2_out)

    def test_following_36(self):
        tot_time = self.KEY_DICT.copy()
        tot_time["12"] = 6.5
        tot_time["34"] = 1
        tot_time["56"] = 1
        tot_time["78"] = 1
        self.assertEqual(self.tot_deltas_t13, tot_time)

    def test_following_37(self):
        followings_mouse2_out = self.KEY_DICT.copy()
        followings_mouse2_out["12"] = 2
        followings_mouse2_out["34"] = 1
        self.assertEqual(self.followings_mouse23,
                         followings_mouse2_out)


class TestCalculateExpected(unittest.TestCase):
    def test_summing(self):
        wm1 = {'1': 2, '2':5}
        fm2 = {'1':.5, '2':.5}
        out = ms.calculate_expected(wm1, fm2)
        self.assertTrue(out, 3)


class TestAreM2Dicts_correct(unittest.TestCase):
    def test_KeyError(self):
        dict1 = {"mouse1": 0}
        dict2 = {"mouse2": 0}
        self.assertRaises(KeyError, ms.are_mouse_dicts_correct,
                          dict1, dict2)

    def test_ValueError(self):
        dict1 = {"mouse1":0}
        dict2 = {"mouse1":3}
        self.assertRaises(ValueError, ms.are_mouse_dicts_correct,
                          dict1, dict2)

    def test_AssertionError1(self):
        dict1 = {"mouse1":1}
        dict2 = {"mouse1":1, "mouse2":0}
        self.assertRaises(AssertionError, ms.are_mouse_dicts_correct,
                          dict1, dict2)

class TestCalculateExpectedMatrices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dict1 = {"mouse1": {"12":12, "21":0,
                                "34":1, "43":11,
                                "56":0.4, "65":34,
                                "78":0, "87":5},
                     "mouse2": {"12":1, "21":8,
                                "34":1, "43":2,
                                "56":2.4, "65":15,
                                "78":2, "87":0}}
        cls.dict2 = {"mouse1": {"12":2, "21":1,
                                "34":0, "43":6,
                                "56":11, "65":0,
                                "78":7, "87":1},
                     "mouse2": {"12":11, "21":21,
                                "34":0, "43":1,
                                "56":3.4, "65":0,
                                "78":2, "87":0}}
        cls.fm2 = {"mouse1": {"12":1, "21":1,
                              "34":0, "43":1,
                              "56":3, "65":0,
                              "78":2, "87":1},
                   "mouse2": {"12":4, "21":6,
                              "34":0, "43":1,
                              "56":1, "65":0,
                              "78":1, "87":0}}
        cls.duration = 1000
        mice_list = ["mouse1", "mouse2"]
        cls.outf, cls.outt = ms.calculate_expected_matrices(cls.dict1,
                                                            cls.dict2,
                                                            cls.dict1,
                                                            cls.fm2,
                                                            cls.duration,
                                                            mice_list)

    def test_throws_exception(self):
        self.assertRaises(AssertionError, ms.calculate_expected_matrices,
                          self.dict1, self.dict2, self.dict1, self.fm2,
                          self.duration, ["mouse1"])

    def test_throws_exception2(self):
        dict2 = {"mouse1": {"12":2, "21":1,
                            "34":0, "43":6,
                            "56":11, "65":0,
                            "78":7, "87":1},
                 "mouse2": {"12":11, "21":21,
                            "34":0, "43":1,
                            "56":3.4, "65":0,
                            "78":2, "87":0}}
        fm2 = {"mouse1": {"12":1, "21":1,
                          "34":3, "43":0,
                          "56":11, "65":0,
                          "78":7, "87":1},
               "mouse2": {"12":11, "21":21,
                          "34":0, "43":1,
                          "56":8, "65":0,
                          "78":2, "87":0}}
        self.assertRaises(ValueError, ms.calculate_expected_matrices,
                          self.dict1, dict2, fm2, self.dict1,
                          self.duration, ["mouse1", "mouse2"])

    def test_correct_1(self):
        self.assertEqual(self.outf[0, 0], 0)

    def test_correct_2(self):
        self.assertEqual(self.outf[1, 1], 0)

    def test_correct_3(self):
        self.assertEqual(self.outt[0, 0], 0)

    def test_correct_4(self):
        self.assertEqual(self.outt[1, 1], 0)

    def test_correct_5(self):
        self.assertEqual(self.outf[0, 1], 59.4/1000)

    def test_correct_6(self):
        self.assertEqual(self.outf[1, 0], 22.2/1000)

    def test_correct_7(self):
        self.assertEqual(self.outt[0, 1], 144.36/1000**2)

    def test_correct_8(self):
        self.assertEqual(self.outt[1, 0], 62.4/1000**2)


class TestFollowingMatrices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ta = {"mouse1": [[15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34 ],
                         [1, 2,   3,    4,  5,  6,  7,   8,  1, 2]],
              "mouse2": [[10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35],
                         [8, 1,   2,    3,  4,  5,  6,   7, 8,   1,   2],],
              "mouse3": [[10, 16, 17, 18, 22, 25, 26, 27],
                         [1,  2,   3, 4,  4,  3,  2, 1],]
        }
        mice_list = ["mouse1", "mouse2", "mouse3"]
        out = ms.following_matrices(ta, mice_list, 1000)
        cls.following = out[0]
        cls.time_together = out[1]
        cls.interval_details = out[2]
        cls.int_m1 = out[3]
        cls.int_m2 = out[4]
        cls.fm2 = out[5]

    def test_following_diag_1(self):
        self.assertEqual(self.following[0, 0], 0)

    def test_following_diag_2(self):
        self.assertEqual(self.following[1, 1], 0)

    def test_following_diag_3(self):
        self.assertEqual(self.following[2, 2], 0)

    def test_following_0_1(self):
        self.assertEqual(self.following[0, 1], 3)

    def test_following_0_2(self):
        self.assertEqual(self.following[0, 2], 0)

    def test_following_1_0(self):
        self.assertEqual(self.following[1, 0], 0)

    def test_following_1_2(self):
        self.assertEqual(self.following[1, 2], 0)

    def test_following_2_0(self):
        self.assertEqual(self.following[2, 0], 1)

    def test_following_2_1(self):
        self.assertEqual(self.following[2, 1], 1)

    def test_time_together_diag_0(self):
        self.assertEqual(self.time_together[0, 0], 0)

    def test_time_together_diag_1(self):
        self.assertEqual(self.time_together[1, 1], 0)

    def test_time_together_diag_2(self):
        self.assertEqual(self.time_together[2, 2], 0)

    def test_time_together_0_1(self):
        self.assertEqual(self.time_together[0, 1], 0.004)

    def test_time_together_0_2(self):
        self.assertEqual(self.time_together[0, 2], 0)

    def test_time_together_1_0(self):
        self.assertEqual(self.time_together[1, 0], 0)

    def test_time_together_1_2(self):
        self.assertEqual(self.time_together[1, 2], 0)

    def test_time_together_2_0(self):
        self.assertEqual(self.time_together[2, 0], 0.001)

    def test_time_together_2_1(self):
        self.assertEqual(self.time_together[2, 1], 0)

    def test_interval_details_empty(self):
        self.assertEqual(self.time_together[2, 1], 0)

if __name__ == '__main__':
    unittest.main()
