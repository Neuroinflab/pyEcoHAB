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
    def test_following_1(self):
        antennas1 = [1, 2]
        times1 = [15, 16.5]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        out, intervals = ms.following_2_mice_in_pipe(antennas1, times1,
                                                     antennas2, times2)
        self.assertEqual(out, 1)
        self.assertEqual(intervals, [3])

    def test_not_following(self):
        antennas1 = [1, 2]
        times1 = [15, 16.5]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        out, intervals = ms.following_2_mice_in_pipe(antennas2, times2,
                                                     antennas1, times1)
        self.assertEqual(out, 0)
        self.assertEqual(intervals, [])

    def test_following_more(self):
        antennas1 = [1, 2, 3, 4, 5]
        times1 = [15, 16.5, 19, 20, 21]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        out, intervals = ms.following_2_mice_in_pipe(antennas1, times1,
                                                     antennas2, times2)
        self.assertEqual(out, 2)
        self.assertEqual(intervals, [3, 2.5])

    def test_not_following_more(self):
        antennas1 = [1, 2, 3, 4, 5]
        times1 = [15, 16.5, 19, 20, 21]
        antennas2 = [8, 1, 2, 3, 4, 5]
        times2 = [10, 16, 19, 19.5, 22, 25]
        out, intervals = ms.following_2_mice_in_pipe(antennas2, times2,
                                                     antennas1, times1)
        self.assertEqual(out, 0)
        self.assertEqual(intervals, [])


class TestCalculateExpectedFollowings(unittest.TestCase):
    def test_summing(self):
        wm1 = {'1': 2, '2':5}
        fm2 = {'1':.5, '2':.5}
        out = ms.calculate_expected_followings(wm1, fm2)
        self.assertTrue(out, 3)


class TestFrequencyMouseInTube(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.period = 16.
        cls.antennas = [1, 1, 2, 3, 4, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8]
        cls.times = [0, 1, 1.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        cls.freq, cls.window = ms.frequency_mouse_in_tube(cls.times,
                                                          cls.antennas,
                                                          cls.period)
    def test_same_keys(self):
        self.assertEqual(self.freq.keys(), self.window.keys())

    def test_all_keys(self):
        self.assertEqual(set(self.freq.keys()), set(['12', '21', '34',
                                                     '43', '56', '65',
                                                     '78', '87']))

    def test_freq_12(self):
        self.assertEqual(self.freq['12'], 2/16)

    def test_freq_21(self):
        self.assertEqual(self.freq['21'], 1/16)

    def test_freq_34(self):
        self.assertEqual(self.freq['34'], 2/16)

    def test_freq_43(self):
        self.assertEqual(self.freq['43'], 1/16)

    def test_freq_56(self):
        self.assertEqual(self.freq['56'], 1/16)

    def test_freq_65(self):
        self.assertEqual(self.freq['65'], 0)

    def test_freq_78(self):
        self.assertEqual(self.freq['78'], 1/16)

    def test_freq_87(self):
        self.assertEqual(self.freq['87'], 0)


    def test_window_12(self):
        self.assertEqual(self.window['12'], 1.5)

    def test_window_21(self):
        self.assertEqual(self.window['21'], 1)

    def test_window_34(self):
        self.assertEqual(self.window['34'], 2)

    def test_window_43(self):
        self.assertEqual(self.window['43'], 1)

    def test_window_56(self):
        self.assertEqual(self.window['56'], 1)

    def test_window_65(self):
        self.assertEqual(self.window['65'], 0)

    def test_window_78(self):
        self.assertEqual(self.window['78'], 1)

    def test_window_87(self):
        self.assertEqual(self.window['87'], 0)

if __name__ == '__main__':
    unittest.main()
