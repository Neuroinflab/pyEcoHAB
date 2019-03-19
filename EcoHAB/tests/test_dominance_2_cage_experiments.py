#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
from EcoHAB import dominance_2_cage_experiments as dom
import unittest
import numpy as np

class TestGetStates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        dt = 0.1
        cls.t_start = 600.0
        cls.t_end = 800.0
        cls.home_antenna_1 = 3
        cls.home_antenna_2 = 4

        cls.antennas_1 = [3, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4, 3, 3]
        cls.times_1 = [641.083, #3
                      642.135, #4
                      675.134, #4
                      675.869, #3
                      681.127, #3
                      681.734, #4
                      692.744, #4
                      693.207, #3
                      701.82,  #3
                      702.603, #4
                      703.499, #4
                      703.961, #3
                      723.136, #3
                      725.633, #4
                      734.133, #4
                      734.945, #3
                      783.411]#3
        cls.out_1 = dom.get_states_mouse(cls.antennas_1,
                                         cls.times_1,
                                         cls.t_start,
                                         cls.t_end,
                                         cls.home_antenna_1,
                                         dt)
        cls.out_2 = dom.get_states_mouse(cls.antennas_1,
                                         cls.times_1,
                                         cls.t_start,
                                         cls.t_end,
                                         cls.home_antenna_2,
                                         dt)
        # for i, x in enumerate(cls.out_2):
        #     print(i, x)
    def test_same_length(self):
        self.assertEqual(len(self.out_1),
                         len(self.out_2))
    
    def test_different_results_for_different_home_antenna(self):
        for i, x in enumerate(self.out_1):
            if x != 1:
                self.assertNotEqual(x, self.out_2[i])

    def test_different_results_for_home_antenna_1(self):
        self.assertFalse(np.all(self.out_1 == self.out_1[0]))

    def test_different_results_for_home_antenna_2(self):
        self.assertFalse(np.all(self.out_2 == self.out_2[0]))

    def test_same_results_pipe(self):
        for i, x in enumerate(self.out_1):
            if x == 1:
                self.assertEqual(x, self.out_2[i])
                
    def test_different_home_for_home_antenna_1(self):
        self.assertTrue(np.any(self.out_1 == 0))

    def test_different_home_for_home_antenna_2(self):
        self.assertTrue(np.any(self.out_2 == 0))


if __name__ == '__main__':
    unittest.main()
