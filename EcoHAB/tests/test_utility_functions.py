#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
from EcoHAB import utility_functions as uf
import unittest
import numpy as np

class TestInPipe(unittest.TestCase):
    def test_pipe_1_2_clockwise(self):
        self.assertTrue(uf.in_tube(1, 2))
    def test_pipe_2_1_counterclockwise(self):
        self.assertTrue(uf.in_tube(2, 1))
    def test_pipe_2_3_clockwise(self):
        self.assertFalse(uf.in_tube(2, 3))
    def test_pipe_3_2_counterclockwise(self):
        self.assertFalse(uf.in_tube(3, 2))
    def test_pipe_7_8_clockwise(self):
        self.assertTrue(uf.in_tube(7, 8))
    def test_pipe_8_7_counterclockwise(self):
        self.assertTrue(uf.in_tube(8, 7))
    def test_pipe_1_8_clockwise(self):
        self.assertFalse(uf.in_tube(1, 8))
    def test_pipe_8_1_counterclockwise(self):
        self.assertFalse(uf.in_tube(8, 1))

class TestInChambers(unittest.TestCase):
    def test_pipe_1_2_clockwise(self):
        self.assertFalse(uf.in_chamber(1, 2))
    def test_pipe_2_1_counterclockwise(self):
        self.assertFalse(uf.in_chamber(2, 1))
    def test_pipe_2_3_clockwise(self):
        self.assertTrue(uf.in_chamber(2, 3))
    def test_pipe_3_2_counterclockwise(self):
        self.assertTrue(uf.in_chamber(3, 2))
    def test_pipe_7_8_clockwise(self):
        self.assertFalse(uf.in_chamber(7, 8))
    def test_pipe_8_7_counterclockwise(self):
        self.assertFalse(uf.in_chamber(8, 7))
    def test_pipe_1_8_clockwise(self):
        self.assertTrue(uf.in_chamber(1, 8))
    def test_pipe_8_1_counterclockwise(self):
        self.assertTrue(uf.in_chamber(8, 1))
    
class TestChangeState(unittest.TestCase):
    def test_check_no_change(self):
        self.assertEqual(len(uf.change_state([1, 1, 1, 1,])), 0)
    def test_check_change_1(self):
        self.assertFalse(len(uf.change_state([1, 1, 1, 2])) == 0)
    def test_check_change_2(self):
        self.assertTrue(uf.change_state([1, 1, 1, 2]), np.array([2], dtype=int))
        
class TestGetIdxPre(unittest.TestCase):
    def test_false(self):
        self.assertEqual(uf.get_idx_pre(2, [3, 4, 5]), None)
    def test_correct(self):
        self.assertEqual(uf.get_idx_pre(2, [-1, 0, 1, 2, 3]), 2)
    def test_empty(self):
        self.assertEqual(uf.get_idx_pre(0, []), None)

class TestGetIdxPost(unittest.TestCase):
    def test_false(self):
        self.assertEqual(uf.get_idx_post(2, [-1, 0, 1]), None)
    def test_correct(self):
        self.assertEqual(uf.get_idx_post(2, [-1, 0, 1, 2, 3]), 4)
    def test_empty(self):
        self.assertEqual(uf.get_idx_post(0, []), None)

class TestGetIdxBetween(unittest.TestCase):
    def test_false(self):
        self.assertEqual(len(uf.get_idx_between(2, 3, [-1, 0, 1])), 0)
        
    def test_correct1(self):
        out = uf.get_idx_between(2, 3, [-1, 0, 1, 2, 3])
        res = np.array([3, 4], dtype=int)
        self.assertEqual(len(out), len(res))


    def test_correct_loop(self):
        out = uf.get_idx_between(2, 3, [-1, 0, 1, 2, 3])
        res = np.array([3, 4], dtype=int)
        for i, x in enumerate(out):
            self.assertEqual(x, res[i])

    def test_empty(self):
        self.assertEqual(len(uf.get_idx_between(0, 2, [])), 0)


class TestMouseGoingForward(unittest.TestCase):
    def test_going_forward(self):
        self.assertTrue(uf.mouse_going_forward([1, 2, 3]))

    def test_not_going_forward(self):
        self.assertFalse(uf.mouse_going_forward([1, 2, 1, 8]))

    def test_going_forward_other_pipe(self):
        self.assertTrue(uf.mouse_going_forward([7, 8, 1]))

    def test_not_going_forward_other_pipe(self):
        self.assertFalse(uf.mouse_going_forward([8, 7, 8, 1]))

    def test_going_forward_other_from_data(self):
        self.assertTrue(uf.mouse_going_forward([3, 4, 4, 3, 3, 4, 4, 5]))

    def test_not_going_forward_from_data(self):
        self.assertFalse(uf.mouse_going_forward([5, 6, 6, 5, 4]))
        

if __name__ == '__main__':
    unittest.main()
