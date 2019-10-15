from __future__ import print_function, division, absolute_import
from pyEcoHAB import mouse_speed as ms
import unittest

try:
    basestring
except NameError:
    basestring = str


class TestTimeInTunnel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        antennas = [8, 1, 2, 3, 4, 5]
        times = [10, 16, 19, 19.5, 22, 25]
        cls.out = ms.calculate_time_in_tunnel(times, antennas)

    def test_12(self):
        self.assertEqual(self.out[3], 3)

    def test_34(self):
        self.assertEqual(self.out[7], 2.5)

    def test_56(self):
        self.assertEqual(self.out[11], 0)

    def test_78(self):
        self.assertEqual(self.out[15], 0)

class TestCalculateExpectedPerTunnel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas = [8,  1,  2,    3,  4,  5,  6, 7,  8,   1,  2,  3,  4]
        times =  [ 10, 16, 19, 19.5, 22, 25, 27, 45, 48, 51, 52, 67, 70]
        bins = [10, 20, 30, 40, 50, 60, 70]
        t_stop = 80
        cls.out1 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    3, bins, t_stop)
        cls.out2 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    7, bins, t_stop)
        cls.out3 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    11, bins, t_stop)
        cls.out4 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    15, bins, t_stop)

        bins2 = [10]
        t_stop = 80
        cls.out11 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    3, bins2, t_stop)
        cls.out12 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    7, bins2, t_stop)
        cls.out13 = ms.calculate_expected_per_tunnel(times, antennas,
                                                    11, bins2, t_stop)
        cls.out14 = ms.calculate_expected_per_tunnel(times, antennas,
                                                     15, bins2, t_stop)

    def test_12_through(self):
        self.assertEqual(self.out1[0], [1, 0, 0, 0, 1, 0, 0])

    def test_12_time(self):
        self.assertEqual(self.out1[1], [3, 0, 0, 0, 1, 0, 0])

    def test_34_through(self):
        self.assertEqual(self.out2[0], [1, 0, 0, 0, 0, 1, 0])

    def test_34_time(self):
        self.assertEqual(self.out2[1], [2.5, 0, 0, 0, 0, 3, 0])

    def test_56_through(self):
        self.assertEqual(self.out3[0], [0, 1, 0, 0, 0, 0, 0])

    def test_56_time(self):
        self.assertEqual(self.out3[1], [0, 2, 0, 0, 0, 0, 0])

    def test_78_through(self):
        self.assertEqual(self.out4[0], [0, 0, 0, 1, 0, 0, 0])

    def test_78_time(self):
        self.assertEqual(self.out4[1], [0, 0, 0, 3, 0, 0, 0])

    def test_12_through2(self):
        self.assertEqual(self.out11[0], [2])

    def test_12_time2(self):
        self.assertEqual(self.out11[1], [4])

    def test_34_through2(self):
        self.assertEqual(self.out12[0], [2])

    def test_34_time2(self):
        self.assertEqual(self.out12[1], [5.5])

    def test_56_through2(self):
        self.assertEqual(self.out13[0], [1])

    def test_56_time2(self):
        self.assertEqual(self.out13[1], [2])

    def test_78_through2(self):
        self.assertEqual(self.out14[0], [1])

    def test_78_time2(self):
        self.assertEqual(self.out14[1], [3])


class TestCalculateExpected(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        antennas1 = [8, 1,  2,  3,   4,   5,  6, 7,  8,   1,  2,  3,  4]
        times1 =  [10, 16, 19, 19.5, 22, 25, 27, 45, 48, 51, 52, 67, 70]
        antennas2 = [8, 1,  2,  3,   4,   5,  6, 7,  8,   1,  2,  3,  4]
        times2 =  [11, 17, 20, 21.5, 23, 28, 29, 42, 47, 50, 51, 62, 65]
        readings_antennas = {"mouse1": [times1, antennas1],
                             "mouse2": [times2, antennas2],}
        t_start = 10
        t_stop = 80
        cls.out_f, cls.out_t = ms.expected_matrices(readings_antennas,
                                                    ["mouse1", "mouse2"],
                                                    t_start, t_stop)
    def test_out_f00(self):
        self.assertEqual(self.out_f[0, 0], 0)

    def test_out_f11(self):
        self.assertEqual(self.out_f[1, 1], 0)

    def test_out_t00(self):
        self.assertEqual(self.out_t[0, 0], 0)

    def test_out_t11(self):
        self.assertEqual(self.out_t[1, 1], 0)

    def test_out_f01(self):
        
        self.assertEqual(self.out_f[0, 1],
                         0.09938857762559292)

    def test_out_f10(self):
        self.assertEqual(self.out_f[1, 0],
                         0.09434996119025715)
        
    def test_out_t01(self):
        self.assertEqual(self.out_t[0, 1],
                         0.003786155756201839)

    def test_out_t10(self):
        self.assertEqual(self.out_t[1, 0],
                         0.003784684414262045)


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
        cls.out1, cls.time_together1, cls.intervals1= res

        antennas1 = [1, 2, 3, 4, 5]
        times1 = [15, 16.5, 19, 20, 21]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        res = ms.following_2_mice_in_pipe(antennas2,
                                          times2,
                                          antennas1,
                                          times1)
        cls.out2, cls.time_together2, cls.intervals2 = res

        antennas1 = [1, 2,   3, 4,  5,  6,  7,   8,  1, 2]
        times1 = [15, 16.5, 19, 20, 21, 22, 24, 25, 29, 34 ]
        antennas2 = [8, 1,   2,    3,  4,  5,  6,   7, 8,   1,   2]
        times2 =   [10, 16, 19, 19.5, 22, 25,  26, 27, 28, 31, 35]
        res = ms.following_2_mice_in_pipe(antennas1,
                                          times1,
                                          antennas2,
                                          times2)
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
        self.assertEqual(self.intervals3, [4, 3, 6])


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
        out = ms.following_matrices(ta, mice_list, 0, 1000)
        cls.following = out[0]
        cls.time_together = out[1]
        cls.interval_details = out[2]

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
        self.assertEqual(self.following[2, 1], 0)

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
