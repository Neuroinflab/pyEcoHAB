# SPDX-License-Identifier: LGPL-2.1-or-later
#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import unittest
from pyEcoHAB import tube_dominance as tubed
from pyEcoHAB.utils import general as uf
from pyEcoHAB import SetupConfig, Loader, Timeline
from pyEcoHAB import data_path, sample_data


class TestSingleDirection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = SetupConfig(path=data_path, fname="standard_setup.txt")

    def test_mouse2_pushing_mouse1(self):
        m1_antennas = ["5", "6", "5", "5", "5", "5", "5", "6", "6", "6"]
        m1_times = [705.074,  # 5
                    708.091,  # 6
                    710.577,  # 5
                    744.813,  # 5
                    746.72,  # 5
                    758.851,  # 5
                    802.624,  # 5
                    808.095,  # 6
                    809.252,  # 6
                    809.564]  # 6

        m2_antennas = ["6", "6", "6", "6", "6", "5", "5", "5", "6",
                       "5", "4"]
        m2_times = [751.348,  # 6
                    753.771,  # 6
                    755.115,  # 6
                    755.865,  # 6
                    764.666,  # 6
                    768.856,  # 5
                    769.184,  # 5
                    794.072,  # 5
                    796.386,  # 6
                    801.342,  # 5
                    807.86]  # 4

        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           7, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        4, self.config)
        out = tubed.tube_dominance_single_direction(m1_through["5 6"],
                                                    m2_backing["6 6"])
        self.assertEqual(out[0], 0)

    def test_mouse_simple_not_pushing_mouse2(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["5", "4", "4"]
        m2_times = [936.827, 939.892, 943.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_direction(m1_through["3 4"],
                                                    m2_backing["4 4"])
        self.assertEqual(out[0], 0)

    def test_mouse_simple_pushing_mouse2_1(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["5", "4", "4"]
        m2_times = [937.827, 937.892, 938.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_direction(m1_through["3 4"],
                                                    m2_backing["4 4"])
        self.assertEqual(out[0], 1)

    def test_mouse_simple_pushing_mouse2_2(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["5", "4", "4"]
        m2_times = [937.827, 938.292, 938.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_direction(m1_through["3 4"],
                                                    m2_backing["4 4"])
        self.assertEqual(out[0], 1)

    def test_mouse_simple_mice_together_in_pipe2(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["4", "3", "2"]
        m2_times = [936.827, 939.892, 943.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_direction(m1_through["3 4"],
                                                    m2_backing["4 4"])
        self.assertEqual(out[0], 0)

    def test_mouse1_mouse_2_different_directions_1(self):
        m1_antennas = ["1", "2", "2", "1", "8"]
        m1_times = [59.462, 60.447, 64.418, 64.934, 81.723]
        m2_antennas = ["2", "3", "4", "5", "6", "7"]
        m2_times = [57.727, 74.407, 74.922, 77.392, 79.628, 94.855]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           8, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)

        out = tubed.tube_dominance_single_direction(m1_through["1 2"],
                                                    m2_backing["2 2"])
        self.assertEqual(out[0], 0)

    def test_mouse1_mouse_2_different_directions_1(self):
        m1_antennas = ["1", "2", "2", "1", "8"]
        m1_times = [59.462, 60.447, 64.418, 64.934, 81.723]
        m2_antennas = ["2", "3", "4", "5", "6", "7"]
        m2_times = [57.727, 74.407, 74.922, 77.392, 79.628, 94.855]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           8, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)

        out = tubed.tube_dominance_single_direction(m1_through["2 1"],
                                                    m2_backing["1 1"])
        self.assertEqual(out[0], 0)

    def test_mouse_2_does_nothing(self):
        m1_antennas = ["3", "4", "5"]
        m1_times = [7.865, 8.287, 10.523]
        m2_antennas = ["4", "5"]
        m2_times = [5.677, 15.51]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)
        out = tubed.tube_dominance_single_direction(m1_through["3 4"],
                                                    m2_backing["4 4"])
        self.assertEqual(out[0], 0)


class TestTubeDominanceSinglePair(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = SetupConfig(path=data_path, fname="standard_setup.txt")

    def test_lab_pushing_mouse2_pushing_mouse1(self):
        m1_antennas = ["5", "6", "5", "5", "5", "5", "5", "6", "6", "6", "7"]
        m1_times = [705.074,  # 5
                    708.091,  # 6
                    710.577,  # 5
                    744.813,  # 5
                    746.72,  # 5
                    758.851,  # 5
                    802.624,  # 5
                    808.095,  # 6
                    809.252,  # 6
                    809.564,  # 6
                    813.675]  # 7
        m2_antennas = ["6", "6", "6", "6", "6", "5", "5", "5", "6",
                       "5", "4", "4"]
        m2_times = [751.348,  # 6
                    753.771,  # 6
                    755.115,  # 6
                    755.865,  # 6
                    764.666,  # 6
                    768.856,  # 5
                    769.184,  # 5
                    794.072,  # 5
                    796.386,  # 6
                    801.342,  # 5
                    807.86,  # 4
                    814.3]  # 4
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)
        out = tubed.tube_dominance_single_pair(m1_through, m2_backing)
        self.assertEqual(0, out[0])

    def test_mouse_simple_not_pushing_mouse2(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["5", "4", "4"]
        m2_times = [936.827, 939.892, 943.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 0)

    def test_mouse_simple_pushing_mouse2_1(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["5", "4", "4"]
        m2_times = [937.827, 937.892, 938.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 1)

    def test_mouse_simple_pushing_mouse2_2(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["5", "4", "4"]
        m2_times = [937.827, 938.292, 938.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 1)

    def test_mouse_simple_mice_together_in_pipe2(self):
        m1_antennas = ["3", "4", "4", "5"]
        m1_times = [938.187, 939.297, 940.297, 942.267]
        m2_antennas = ["4", "3", "2"]
        m2_times = [936.827, 939.892, 943.486]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        3, self.config)

        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 0)

    def test_mouse1_mouse_2_different_directions_1(self):
        m1_antennas = ["1", "2", "2", "1", "8"]
        m1_times = [59.462, 60.447, 64.418, 64.934, 81.723]
        m2_antennas = ["2", "3", "4", "5", "6", "7"]
        m2_times = [57.727, 74.407, 74.922, 77.392, 79.628, 94.855]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           8, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)

        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 0)

    def test_mouse1_mouse_2_different_directions_1(self):
        m1_antennas = ["1", "2", "2", "1", "8"]
        m1_times = [59.462, 60.447, 64.418, 64.934, 81.723]
        m2_antennas = ["2", "3", "4", "5", "6", "7"]
        m2_times = [57.727, 74.407, 74.922, 77.392, 79.628, 94.855]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           8, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)

        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 0)

    def test_mouse_2_does_nothing(self):
        m1_antennas = ["3", "4", "5"]
        m1_times = [7.865, 8.287, 10.523]
        m2_antennas = ["4", "5"]
        m2_times = [5.677, 15.51]
        m1_through = uf.extract_directions(m1_times, m1_antennas,
                                           6, self.config.directions)
        m2_backing = uf.extract_backing(m2_times, m2_antennas,
                                        7, self.config)
        out = tubed.tube_dominance_single_pair(m1_through,
                                               m2_backing)
        self.assertEqual(out[0], 0)


class TestExecution(unittest.TestCase):
    def test(self):
        data = Loader(sample_data)
        config = Timeline(sample_data)
        tubed.get_tube_dominance(data, config, 1)


if __name__ == '__main__':
    unittest.main()
