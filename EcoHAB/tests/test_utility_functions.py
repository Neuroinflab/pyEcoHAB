#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function, division, absolute_import
from EcoHAB import utility_functions as uf
import unittest

class TestInPipe(unittest.TestCase):
    def test_pipe_1_2_clockwise(self):
        self.assertEqual(uf.in_tube(1, 2), True)
    def test_pipe_2_1_counterclockwise(self):
        self.assertEqual(uf.in_tube(2, 1), True)
    def test_pipe_2_3_clockwise(self):
        self.assertEqual(uf.in_tube(2, 3), False)
    def test_pipe_3_2_counterclockwise(self):
        self.assertEqual(uf.in_tube(3, 2), False)
    def test_pipe_7_8_clockwise(self):
        self.assertEqual(uf.in_tube(7, 8), True)
    def test_pipe_8_7_counterclockwise(self):
        self.assertEqual(uf.in_tube(8, 7), True)
    def test_pipe_1_8_clockwise(self):
        self.assertEqual(uf.in_tube(1, 8), False)
    def test_pipe_8_1_counterclockwise(self):
        self.assertEqual(uf.in_tube(8, 1), False)

class TestInChambers(unittest.TestCase):
    def test_pipe_1_2_clockwise(self):
        self.assertEqual(uf.in_chamber(1, 2), False)
    def test_pipe_2_1_counterclockwise(self):
        self.assertEqual(uf.in_chamber(2, 1), False)
    def test_pipe_2_3_clockwise(self):
        self.assertEqual(uf.in_chamber(2, 3), True)
    def test_pipe_3_2_counterclockwise(self):
        self.assertEqual(uf.in_chamber(3, 2), True)
    def test_pipe_7_8_clockwise(self):
        self.assertEqual(uf.in_chamber(7, 8), False)
    def test_pipe_8_7_counterclockwise(self):
        self.assertEqual(uf.in_chamber(8, 7), False)
    def test_pipe_1_8_clockwise(self):
        self.assertEqual(uf.in_chamber(1, 8), True)
    def test_pipe_8_1_counterclockwise(self):
        self.assertEqual(uf.in_chamber(8, 1), True)
    
        
    

if __name__ == '__main__':
    unittest.main()
